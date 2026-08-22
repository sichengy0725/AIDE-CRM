from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path

from reportlab.lib.pagesizes import landscape, letter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE_BUILDER = ROOT / "tmp" / "pdfs" / "build_tite_vs_nontite_duration_table.py"
TITE_CSV = ROOT / "Presentation 8-23-2026" / "TITE_AIDE_phase_I_II_IDX_1001_to_2000_dose_summary.csv"
NON_TITE_CSV = ROOT / "Presentation 8-23-2026" / (
    "AIDE_phase_I_II_modelsadditive_shared_N30_ncycle2_"
    "rp0p15x0p85_rate56d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_"
    "ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv"
)
OUTPUT = ROOT / "output" / "pdf" / (
    "phase12_tite_vs_nontite_N30_scenarios_1_16_20_24_27_38_tables_"
    "ipde_first_tite_IDX_1001_to_2000.pdf"
)

SCENARIOS = [1, 16, 20, 24, 27, 38]
DOSES = [1, 2, 3, 4, 5]
TITE_TASKS = {
    "two_stage": {1: 1, 16: 2, 20: 3, 24: 4, 27: 5, 38: 6},
    "one_stage": {1: 7, 16: 8, 20: 9, 24: 10, 27: 11, 38: 12},
}
NON_TITE_TASKS = TITE_TASKS
DESIGNS = [("one_stage", "one-stage"), ("two_stage", "two-stage")]


def load_template():
    spec = importlib.util.spec_from_file_location("duration_template", TEMPLATE_BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the 8/17 table-layout builder.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def number(row: dict[str, str], field: str) -> float:
    return float(row[field])


def close(value: str, expected: float) -> bool:
    return math.isclose(float(value), expected, rel_tol=1e-9, abs_tol=1e-9)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def complete_rows(rows: list[dict[str, str]], description: str) -> list[dict[str, str]]:
    rows.sort(key=lambda row: int(row["Dose"]))
    if len(rows) != 5 or [int(row["Dose"]) for row in rows] != DOSES:
        raise RuntimeError(f"{description} is not a complete five-dose row set.")
    return rows


def select_tite(raw: list[dict[str, str]], scenario: int, design_key: str) -> list[dict[str, str]]:
    stage_one_max = 30 if design_key == "one_stage" else 6
    rows = [
        row for row in raw
        if int(row["Task_ID"]) == TITE_TASKS[design_key][scenario]
        and int(row["Scenario"]) == scenario
        and row["Allocation"] == design_key
        and row["Stage2_Allocation"] == "one_stage"
        and row["Model_ID"] == "additive_shared"
        and row["Carryover_Model"] == "additive_shared"
        and row["Efficacy_Model"] == "previous_dose_additive"
        and int(row["Nmax"]) == 30
        and int(row["N_s1"]) == stage_one_max
        and int(row["N_s2"]) == 30
        and int(row["n_eval"]) == 3
        and int(row["Cycle_Max"]) == 2
        and int(row["Utility_Type"]) == 3
        and close(row["Lambda_T"], 1)
        and close(row["Toxicity_IPDE_Alpha"], 0)
        and close(row["Efficacy_IPDE_Alpha"], 0)
        and int(row["ntrial"]) == 1000
        and close(row["Arrival_Rate"], 1 / 56)
    ]
    return complete_rows(rows, f"TITE {design_key}, Scenario {scenario}")


def select_non_tite(raw: list[dict[str, str]], scenario: int, design_key: str) -> list[dict[str, str]]:
    stage_one_max = 30 if design_key == "one_stage" else 6
    rows = [
        row for row in raw
        if int(row["Task_ID"]) == NON_TITE_TASKS[design_key][scenario]
        and int(row["Scenario"]) == scenario
        and row["Allocation"] == design_key
        and int(row["N_s1"]) == stage_one_max
        and row["Enrollment_Scheme"] == "ipde_first"
        and row["Model_ID"] == "additive_shared"
        and int(row["Nmax"]) == 30
        and int(row["N_s2"]) == 30
        and int(row["Utility_Type"]) == 3
        and int(row["Cycle_Max"]) == 2
        and int(row["n_valid"]) == 1000
        and close(row["Toxicity_IPDE_Alpha"], 0)
        and close(row["Efficacy_IPDE_Alpha"], 0)
        and close(row["Efficacy_Threshold"], 0.2)
        and close(row["Futility_Cutoff"], 0.85)
        and int(row["Min_Eff_N_for_Futility"]) == 0
    ]
    return complete_rows(rows, f"IPDE-first non-TITE {design_key}, Scenario {scenario}")


def normalize(rows: list[dict[str, str]], version: str, design_key: str, short_name: str, scenario: int) -> dict:
    return {
        "scenario": scenario,
        "version": version,
        "design": design_key,
        "method": f"{version} | {short_name}",
        "dlt": [number(row, "True_DLT_rate") for row in rows],
        "efficacy": [number(row, "True_Efficacy_rate") for row in rows],
        "selection": [number(row, "OBD_Selection_pct") for row in rows],
        "administrations": [number(row, "Pts_Treated") for row in rows],
        "no_obd": number(rows[0], "No_OBD_Selection_pct"),
        "duration": number(rows[0], "Duration"),
    }


def build_summary() -> list[dict]:
    tite_raw = read_csv(TITE_CSV)
    non_tite_raw = read_csv(NON_TITE_CSV)
    records: list[dict] = []
    for scenario in SCENARIOS:
        scenario_records: list[dict] = []
        for design_key, short_name in DESIGNS:
            scenario_records.extend([
                normalize(select_tite(tite_raw, scenario, design_key), "TITE", design_key, short_name, scenario),
                normalize(select_non_tite(non_tite_raw, scenario, design_key), "IPDE-first non-TITE", design_key, short_name, scenario),
            ])
        truth = (scenario_records[0]["dlt"], scenario_records[0]["efficacy"])
        if any((record["dlt"], record["efficacy"]) != truth for record in scenario_records[1:]):
            raise RuntimeError(f"Dose-level truth differs between source rows for requested Scenario {scenario}.")
        records.extend(scenario_records)
    return records


def record_for(records: list[dict], scenario: int, design: str, version: str) -> dict:
    return next(row for row in records if row["scenario"] == scenario and row["design"] == design and row["version"] == version)


def draw_scenario(pdf, y: float, scenario: int, records: list[dict], template) -> float:
    truth = record_for(records, scenario, "one_stage", "TITE")
    y = template.draw_row(pdf, y, 15, f"Scenario {scenario}", [""] * 6, fill=template.SECTION, bold=True)
    y = template.draw_row(pdf, y, 10.5, "DLT rate", [template.fmt(value, 2) for value in truth["dlt"]] + [""])
    y = template.draw_row(pdf, y, 10.5, "Efficacy rate", [template.fmt(value, 2) for value in truth["efficacy"]] + [""])
    y = template.draw_row(pdf, y, 12, "OBD selection (%)", [""] * 6, fill=template.SECTION, bold=True)
    for design_key, _ in DESIGNS:
        for version, fill, bold in [("TITE", template.TITE_FILL, False), ("IPDE-first non-TITE", template.NON_TITE_FILL, True)]:
            row = record_for(records, scenario, design_key, version)
            values = [template.fmt(value) for value in row["selection"]] + [template.fmt(row["no_obd"])]
            y = template.draw_row(pdf, y, 10.5, row["method"], values, fill=fill, bold=bold)
    y = template.draw_row(pdf, y, 12, "Mean administrations by dose", [""] * 6, fill=template.SECTION, bold=True)
    for design_key, _ in DESIGNS:
        for version, fill, bold in [("TITE", template.TITE_FILL, False), ("IPDE-first non-TITE", template.NON_TITE_FILL, True)]:
            row = record_for(records, scenario, design_key, version)
            y = template.draw_row(pdf, y, 10.5, row["method"], [template.fmt(value) for value in row["administrations"]] + [""], fill=fill, bold=bold)
    y = template.draw_duration_header(pdf, y)
    for design_key, _ in DESIGNS:
        for version, fill, bold in [("TITE", template.TITE_FILL, False), ("IPDE-first non-TITE", template.NON_TITE_FILL, True)]:
            row = record_for(records, scenario, design_key, version)
            y = template.draw_duration_row(pdf, y, row["method"], row["duration"], fill, bold)
    return y


def draw_page(pdf, scenario: int, records: list[dict], template, page_number: int) -> None:
    template.draw_text(pdf, "Phase I/II Operating Characteristics", template.PAGE_WIDTH / 2, 576, size=20, bold=True, align="center")
    template.draw_text(pdf, "TITE versus IPDE-first non-TITE: one-stage and two-stage allocation", template.PAGE_WIDTH / 2, 556, size=10.2, color=template.MUTED, align="center")
    pdf.setStrokeColor(template.NAVY)
    pdf.setLineWidth(0.8)
    pdf.line(template.LEFT, 545, template.RIGHT, 545)
    for index, header in enumerate(["Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "No OBD %"]):
        template.draw_text(pdf, header, template.LEFT + template.LABEL_WIDTH + template.VALUE_WIDTH * (index + 0.5), 531, size=8.8, bold=True, align="center")
    draw_scenario(pdf, 518, scenario, records, template)
    template.draw_text(pdf, "Sources: 8/23 TITE IDX 1001-2000 and IPDE-first non-TITE summaries; all rows represent 1,000 simulations.", template.LEFT, 27, size=6.2, color=template.MUTED)
    template.draw_text(pdf, "TITE and non-TITE rows are matched by scenario truth, allocation, model, IPDE alpha = 0, N = 30, and a 56-day arrival rate.", template.LEFT, 21, size=6.1, color=template.MUTED)
    template.draw_text(pdf, "Duration is mean calendar days to complete the trial.", template.LEFT, 15, size=6.1, color=template.MUTED)
    template.draw_text(pdf, f"Page {page_number} of {len(SCENARIOS)}", template.RIGHT, 15, size=6.2, color=template.MUTED, align="right")


def build() -> Path:
    template = load_template()
    records = build_summary()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    pdf = canvas.Canvas(str(OUTPUT), pagesize=landscape(letter), pageCompression=1)
    for page_number, scenario in enumerate(SCENARIOS, start=1):
        draw_page(pdf, scenario, records, template, page_number)
        pdf.showPage()
    pdf.save()
    return OUTPUT


if __name__ == "__main__":
    print(build())
