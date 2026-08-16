from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import build_shared_vs_dose_specific_alpha_tables as base


ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "Presentation 8-10-2026" / "Raw Data"
CYCLE2_SOURCE = RAW_DIR / (
    "AIDE_phase_I_II_N30_ncycle2_rp0p15x0p85_ep0p5x0p5_cp0p15x0p85_"
    "eth0p2_fut0p85_ap0p15x0p85_0p3x0p7_0p5x0p5_1x1_IDX_0001_to_1000_dose_summary.csv"
)
CYCLE1_ORACLE_SOURCE = RAW_DIR / (
    "AIDE_phase_I_II_modelsadditive_shared_N30_ncycle1_rp0p15x0p85_rate56d_"
    "ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv"
)
OUT_DIR = ROOT / "Presentation 8-10-2026" / "Table and Plots" / "Phase I-II" / "Model"

ROW_H = 9
DOSE_SPECIFIC_SHADE = base.VERY_LIGHT_BLUE
ORACLE_SHADE = base.HexColor("#E9F4E8")
METHODS = [
    (2, "previous_dose_additive", "U2 - shared alpha", "shared"),
    (2, "dose_specific_previous_dose_additive", "U2 - dose-specific alpha", "dose_specific"),
    (2, "cycle1_oracle", "U2 - cycle-1 oracle", "oracle"),
    (3, "previous_dose_additive", "U3 - shared alpha", "shared"),
    (3, "dose_specific_previous_dose_additive", "U3 - dose-specific alpha", "dose_specific"),
    (3, "cycle1_oracle", "U3 - cycle-1 oracle", "oracle"),
]


def select_rows(path: Path, oracle: bool) -> dict[tuple, list[dict[str, str]]]:
    result: dict[tuple, list[dict[str, str]]] = defaultdict(list)
    numeric_requirements = {
        "CRM_Prior_a": 0.15,
        "CRM_Prior_b": 0.85,
        "Efficacy_Prior_a": 0.5,
        "Efficacy_Prior_b": 0.5,
        "Efficacy_Carryover_Prior_a": 0.15,
        "Efficacy_Carryover_Prior_b": 0.85,
        "Efficacy_Additive_Alpha_Prior_a": 0.15,
        "Efficacy_Additive_Alpha_Prior_b": 0.85,
        "Efficacy_Threshold": 0.20,
        "Futility_Cutoff": 0.85,
    }
    valid_models = {"previous_dose_additive"} if oracle else {
        "previous_dose_additive", "dose_specific_previous_dose_additive"
    }
    for row in base.load_rows(path):
        if row["Allocation"] not in {item[0] for item in base.ALLOCATIONS}:
            continue
        if int(base.number(row, "Cycle_Max")) != (1 if oracle else 2):
            continue
        if row["Enrollment_Scheme"] != "continuous" or not base.close(base.number(row, "Arrival_Rate"), 0.0179):
            continue
        if row["CRM_r_Model"] != "previous_dose" or row["Efficacy_Model"] not in valid_models:
            continue
        if row["Model_ID"] != ("additive_shared" if oracle else row["Model_ID"]):
            continue
        if not oracle and row["Model_ID"] not in {"additive_shared", "additive_dose_specific"}:
            continue
        if int(base.number(row, "IPDE_Design")) != 1:
            continue
        if not all(base.close(base.number(row, field), value) for field, value in numeric_requirements.items()):
            continue
        alpha_t = base.number(row, "Toxicity_IPDE_Alpha")
        alpha_e = base.number(row, "Efficacy_IPDE_Alpha")
        if not base.close(alpha_t, alpha_e) or not any(base.close(alpha_t, alpha) for alpha in base.ALPHAS):
            continue
        model = "cycle1_oracle" if oracle else row["Efficacy_Model"]
        key = (int(base.number(row, "Scenario")), row["Allocation"], int(base.number(row, "Utility_Type")), alpha_t, model)
        result[key].append(row)
    return result


def build_data() -> dict[tuple, list[dict[str, str]]]:
    cycle2 = select_rows(CYCLE2_SOURCE, oracle=False)
    oracle = select_rows(CYCLE1_ORACLE_SOURCE, oracle=True)
    result = cycle2 | oracle
    expected = {
        (scenario, allocation, utility, alpha, model)
        for scenario in range(1, 38)
        for allocation, _ in base.ALLOCATIONS
        for utility in (2, 3)
        for alpha in base.ALPHAS
        for model in ("previous_dose_additive", "dose_specific_previous_dose_additive", "cycle1_oracle")
    }
    if set(result) != expected:
        missing = sorted(expected - set(result))[:5]
        extra = sorted(set(result) - expected)[:5]
        raise ValueError(f"The selected model-comparison rows are incomplete. Missing={missing}; extra={extra}")
    for key, rows in result.items():
        if len(rows) != 5:
            raise ValueError(f"{key} has {len(rows)} dose levels; expected five.")
        rows.sort(key=lambda item: int(base.number(item, "Dose")))
    return result


def draw_metric_row(c, label, values, stop, y, kind):
    if kind == "dose_specific":
        c.setFillColor(DOSE_SPECIFIC_SHADE)
        c.rect(base.LEFT, y - 3, base.RIGHT - base.LEFT, ROW_H, stroke=0, fill=1)
    elif kind == "oracle":
        c.setFillColor(ORACLE_SHADE)
        c.rect(base.LEFT, y - 3, base.RIGHT - base.LEFT, ROW_H, stroke=0, fill=1)
    base.draw_text(c, label, base.LABEL_X, y, size=6.35)
    for x, value in zip(base.VALUE_X, values):
        base.draw_right(c, base.fmt(value), x + 42, y, size=6.35)
    base.draw_right(c, base.fmt(stop if stop is not None else float("nan")), base.STOP_X + 48, y, size=6.35)
    return y - ROW_H


def draw_duration_row(c, label, duration, y, kind):
    if kind == "dose_specific":
        c.setFillColor(DOSE_SPECIFIC_SHADE)
        c.rect(base.LEFT, y - 3, base.RIGHT - base.LEFT, ROW_H, stroke=0, fill=1)
    elif kind == "oracle":
        c.setFillColor(ORACLE_SHADE)
        c.rect(base.LEFT, y - 3, base.RIGHT - base.LEFT, ROW_H, stroke=0, fill=1)
    base.draw_text(c, label, base.LABEL_X, y, size=6.35)
    base.draw_right(c, base.fmt(duration), base.STOP_X + 48, y, size=6.35)
    return y - ROW_H


def draw_scenario(c, y, scenario, allocation, alpha, data, truth):
    base_rows = data[(scenario, allocation, 2, alpha, "previous_dose_additive")]
    c.setFillColor(base.LIGHT_BLUE)
    c.rect(base.LEFT, y - 4, base.RIGHT - base.LEFT, 16, stroke=0, fill=1)
    base.draw_text(c, f"Scenario {scenario}", base.LABEL_X, y, size=9.25, color=base.NAVY, font="Helvetica-Bold")
    y -= 19
    truth_rows = [
        ("DLT rate", [base.number(row, "True_DLT_rate") for row in base_rows], False),
        ("Efficacy rate", [base.number(row, "True_Efficacy_rate") for row in base_rows], False),
        ("Utility 2", truth[scenario][2], True),
        ("Utility 3", truth[scenario][3], True),
    ]
    for label, values, utility in truth_rows:
        base.draw_text(c, label, base.LABEL_X, y, size=6.7, color=base.GRAY, font="Helvetica-Bold" if utility else "Helvetica")
        for x, value in zip(base.VALUE_X, values):
            base.draw_right(c, base.fmt(value) if utility else base.fmt_rate(value), x + 42, y, size=6.7, color=base.GRAY)
        y -= ROW_H
    y -= 2
    sections = [
        ("MTD selection (%)", "MTD_Selection_pct", "Early_Stopping_pct"),
        ("OBD selection (%)", "OBD_Selection_pct", "No_OBD_Selection_pct"),
        ("Mean administrations", "Pts_Treated", None),
        ("Mean IPDE administrations", "IPDE_Doses", None),
    ]
    for heading, field, stop_field in sections:
        base.draw_text(c, heading, base.LABEL_X, y, size=6.75, color=base.BLUE, font="Helvetica-Bold")
        y -= ROW_H
        for utility, model, label, kind in METHODS:
            rows = data[(scenario, allocation, utility, alpha, model)]
            stop = base.number(rows[0], stop_field) if stop_field else None
            y = draw_metric_row(c, label, [base.number(row, field) for row in rows], stop, y, kind)
        y -= 1
    base.draw_text(c, "Mean trial duration (days)", base.LABEL_X, y, size=6.75, color=base.BLUE, font="Helvetica-Bold")
    base.draw_right(c, "Days", base.STOP_X + 48, y, size=6.55, color=base.BLUE, font="Helvetica-Bold")
    y -= ROW_H
    for utility, model, label, kind in METHODS:
        rows = data[(scenario, allocation, utility, alpha, model)]
        y = draw_duration_row(c, label, base.number(rows[0], "Duration"), y, kind)
    base.draw_rule(c, y + 4)
    return y - 9


def draw_chrome(c, alpha, allocation_label, page, total):
    c.setFillColor(base.NAVY)
    c.rect(0, base.PAGE_HEIGHT - 50, base.PAGE_WIDTH, 50, stroke=0, fill=1)
    base.draw_text(c, "Phase I/II Operating Characteristics", base.LEFT, base.PAGE_HEIGHT - 31, size=16, color=base.white, font="Helvetica-Bold")
    subtitle = f"Shared vs dose-specific additive-alpha efficacy with cycle-1 oracle - IPDE alpha = {base.fmt(alpha)} - allocation = {allocation_label} - N = 30"
    base.draw_text(c, subtitle, base.LEFT, base.PAGE_HEIGHT - 74, size=9.7, color=base.NAVY, font="Helvetica-Bold")
    for label, x in zip(("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5"), base.VALUE_X):
        base.draw_text(c, label, x, base.PAGE_HEIGHT - 94, size=8.5, color=base.BLUE, font="Helvetica-Bold")
    base.draw_text(c, "No MTD/OBD %", base.STOP_X, base.PAGE_HEIGHT - 94, size=8.5, color=base.BLUE, font="Helvetica-Bold")
    base.draw_rule(c, base.PAGE_HEIGHT - 100, color=base.BLUE, width=0.8)
    footer = "Priors: efficacy beta-binomial Beta(0.5, 0.5); efficacy and toxicity additive alpha Beta(0.15, 0.85)."
    base.draw_text(c, footer, base.LEFT, 35, size=7.1, color=base.GRAY)
    note = "Final column: early stopping % for MTD selection and no-OBD % for OBD selection; blue rows = dose-specific alpha, green rows = cycle-1 oracle."
    base.draw_text(c, note, base.LEFT, 22, size=7.1, color=base.GRAY)
    base.draw_right(c, f"Page {page} of {total}", base.RIGHT, 22, size=7.1, color=base.GRAY)


def make_pdf(alpha, data, truth):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUT_DIR / f"phase12_shared_vs_dose_specific_additive_alpha_with_cycle1_oracle_fut0p85_{base.alpha_tag(alpha)}_scenarios_1_to_37_tables.pdf"
    pdf = base.canvas.Canvas(str(output), pagesize=(base.PAGE_WIDTH, base.PAGE_HEIGHT), pageCompression=1)
    page = 0
    total_pages = 38
    for allocation, allocation_label in base.ALLOCATIONS:
        for scenario in range(1, 38, 2):
            page += 1
            draw_chrome(pdf, alpha, allocation_label, page, total_pages)
            draw_scenario(pdf, base.PAGE_HEIGHT - 126, scenario, allocation, alpha, data, truth)
            if scenario + 1 <= 37:
                draw_scenario(pdf, base.PAGE_HEIGHT - 536, scenario + 1, allocation, alpha, data, truth)
            pdf.showPage()
    pdf.save()
    return output


if __name__ == "__main__":
    comparison_data = build_data()
    utility_truth = base.load_truth()
    for alpha_value in base.ALPHAS:
        print(make_pdf(alpha_value, comparison_data, utility_truth))
