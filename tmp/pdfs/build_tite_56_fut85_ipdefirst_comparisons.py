from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path

import build_tite_56_comparison_tables as base


ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "Presentation 8-10-2026" / "Raw Data"
base.TITE_CSV = RAW_DIR / "TITE_AIDE_phase_I_II-effthr0p2-futcut0p85-mineff3_IDX_1001_to_2000_dose_summary.csv"
base.TASK_MAP_CSV = ROOT / "tmp" / "pdfs" / "tite_all_arrivals_task_map.csv"
base.NON_TITE_CSV = RAW_DIR / "AIDE_phase_I_II_N30_ncycle2_rp0p15x0p85_rate56_28_14d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv"
base.OUT_DIR = ROOT / "Presentation 8-10-2026" / "Table and Plots" / "Phase I-II" / "TITE"


def build_datasets():
    task_map = base.load_rows(base.TASK_MAP_CSV)
    task_by_id = {int(row["Task_ID"]): row for row in task_map}
    if len(task_by_id) != 3552:
        raise ValueError("The full three-rate TITE task map must contain 3,552 tasks.")

    tite_by_key = defaultdict(list)
    for row in base.load_rows(base.TITE_CSV):
        task = task_by_id.get(int(row["Task_ID"]))
        if task is None:
            raise ValueError(f"No recreated task metadata for TITE Task_ID {row['Task_ID']}.")
        if not base.close(base.get_field(row, "Arrival_Rate"), 1 / 56):
            continue
        if task["Model_ID"] != "additive_shared" or not base.close(float(task["Arrival_Rate"]), 1 / 56):
            continue
        alpha_t = float(task["Toxicity_IPDE_Alpha"])
        alpha_e = float(task["Efficacy_IPDE_Alpha"])
        if not base.close(alpha_t, alpha_e) or not any(base.close(alpha_t, alpha) for alpha in base.ALPHAS):
            continue
        key = (int(task["Scenario"]), task["Allocation"], int(task["Utility_Type"]), alpha_t)
        tite_by_key[key].append(row)

    non_tite_by_key = defaultdict(list)
    expected = {
        "CRM_Prior_a": 0.15, "CRM_Prior_b": 0.85,
        "Efficacy_Prior_a": 0.5, "Efficacy_Prior_b": 0.5,
        "Efficacy_Carryover_Prior_a": 0.15, "Efficacy_Carryover_Prior_b": 0.85,
        "Efficacy_Additive_Alpha_Prior_a": 0.15, "Efficacy_Additive_Alpha_Prior_b": 0.85,
        "Efficacy_Threshold": 0.20, "Futility_Cutoff": 0.85, "Min_Eff_N_for_Futility": 0,
    }
    for row in base.load_rows(base.NON_TITE_CSV):
        if row["Enrollment_Scheme"] != "ipde_first" or row["Allocation"] not in {"one_stage", "two_stage"}:
            continue
        if not base.close(base.get_field(row, "Arrival_Rate"), 0.0179):
            continue
        if row["Model_ID"] != "additive_shared" or row["CRM_r_Model"] != "previous_dose":
            continue
        if row["Efficacy_Model"] != "previous_dose_additive" or int(float(row["IPDE_Design"])) != 1:
            continue
        if not all(base.close(base.get_field(row, field), value) for field, value in expected.items()):
            continue
        alpha_t = base.get_field(row, "Toxicity_IPDE_Alpha")
        alpha_e = base.get_field(row, "Efficacy_IPDE_Alpha")
        if not base.close(alpha_t, alpha_e) or not any(base.close(alpha_t, alpha) for alpha in base.ALPHAS):
            continue
        key = (int(row["Scenario"]), row["Allocation"], int(float(row["Utility_Type"])), alpha_t)
        non_tite_by_key[key].append(row)

    expected_keys = {
        (scenario, allocation, utility, alpha)
        for scenario in range(1, 38)
        for allocation, _ in base.ALLOCATIONS
        for utility in (2, 3)
        for alpha in base.ALPHAS
    }
    for name, dataset in (("TITE", tite_by_key), ("non-TITE IPDE-first", non_tite_by_key)):
        if set(dataset) != expected_keys:
            missing = sorted(expected_keys - set(dataset))[:3]
            extra = sorted(set(dataset) - expected_keys)[:3]
            raise ValueError(f"{name} keys do not match requested comparison grid. Missing={missing}; extra={extra}")
        for key, rows in dataset.items():
            if len(rows) != 5:
                raise ValueError(f"{name} {key} has {len(rows)} dose rows; expected 5.")
            rows.sort(key=lambda item: int(item["Dose"]))
    return non_tite_by_key, tite_by_key


def draw_page_chrome(c, alpha, allocation_label, page_num, total_pages):
    c.setFillColor(base.NAVY)
    c.rect(0, base.PAGE_HEIGHT - 50, base.PAGE_WIDTH, 50, stroke=0, fill=1)
    base.draw_text(c, "Phase I/II Operating Characteristics", base.LEFT, base.PAGE_HEIGHT - 31, size=16, color=base.white, font="Helvetica-Bold")
    subtitle = f"IPDE-first non-TITE vs TITE - 56-day accrual - IPDE alpha = {base.fmt(alpha)} - allocation = {allocation_label} - N = 30"
    base.draw_text(c, subtitle, base.LEFT, base.PAGE_HEIGHT - 74, size=10.4, color=base.NAVY, font="Helvetica-Bold")
    for dose, x in enumerate(base.VALUE_X, start=1):
        base.draw_text(c, f"Dose {dose}", x, base.PAGE_HEIGHT - 94, size=8.5, color=base.BLUE, font="Helvetica-Bold")
    base.draw_text(c, "Stop %", base.STOP_X, base.PAGE_HEIGHT - 94, size=8.5, color=base.BLUE, font="Helvetica-Bold")
    base.draw_rule(c, base.PAGE_HEIGHT - 100, color=base.BLUE, width=0.8)
    footer = "Fixed priors: toxicity Beta(0.15, 0.85); efficacy Beta(0.5, 0.5); efficacy additive discount Beta(0.15, 0.85)."
    base.draw_text(c, footer, base.LEFT, 35, size=7.3, color=base.GRAY)
    note = "Non-TITE: IPDE-first enrollment, efficacy threshold = 0.20, futility cutoff = 0.85. TITE: shared previous-dose additive efficacy model, min. efficacy n = 3."
    base.draw_text(c, note, base.LEFT, 22, size=7.3, color=base.GRAY)
    base.draw_right(c, f"Page {page_num} of {total_pages}", base.RIGHT, 22, size=7.3, color=base.GRAY)


def create_pdf(alpha, non_tite, tite, utilities):
    base.OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = base.OUT_DIR / f"phase12_tite_56day_comparison_nontite_ipde_first_fut0p85_{base.alpha_tag(alpha)}_scenarios_1_to_37_tables.pdf"
    c = base.canvas.Canvas(str(path), pagesize=(base.PAGE_WIDTH, base.PAGE_HEIGHT), pageCompression=1)
    total_pages = 38
    page_num = 0
    for allocation, allocation_label in base.ALLOCATIONS:
        for first in range(1, 38, 2):
            page_num += 1
            draw_page_chrome(c, alpha, allocation_label, page_num, total_pages)
            base.draw_scenario(c, base.PAGE_HEIGHT - 126, first, allocation, alpha, non_tite, tite, utilities)
            if first + 1 <= 37:
                base.draw_scenario(c, base.PAGE_HEIGHT - 538, first + 1, allocation, alpha, non_tite, tite, utilities)
            c.showPage()
    c.save()
    return path


def main():
    utilities = base.utilities_from_reference(base.UTILITY_PDF)
    non_tite, tite = build_datasets()
    for alpha in base.ALPHAS:
        print(create_pdf(alpha, non_tite, tite, utilities))


if __name__ == "__main__":
    main()
