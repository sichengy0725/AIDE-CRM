from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path

from reportlab.lib.pagesizes import landscape, letter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE_BUILDER = ROOT / "tmp" / "pdfs" / "build_tite_vs_nontite_duration_table.py"
TITE_CSV = ROOT / "Presentation 8-23-2026" / (
    "TITE_AIDE_phase_I_II-effthr0p2-futcut0p85-mineff3_"
    "IDX_1001_to_2000_dose_summary.csv"
)
NON_TITE_CSV = ROOT / "Presentation 8-23-2026" / (
    "AIDE_phase_I_II_modelsadditive_shared_N30_ncycle2_"
    "rp0p15x0p85_rate56d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_"
    "ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv"
)
OUTPUT = ROOT / "output" / "pdf" / (
    "phase12_tite_vs_nontite_N30_fut0p85_scenarios_"
    "1_16_20_24_27_38_tables_ipde_first_newdata.pdf"
)

SCENARIOS = [1, 16, 20, 24, 27, 38]
DOSES = [1, 2, 3, 4, 5]

# The updated TITE source reports only 37 base scenarios. It has no complete
# truth-vector match for the subsequently added requested Scenario 38.
TITE_SOURCE_SCENARIO = {1: 1, 16: 16, 20: 20, 24: 24, 27: 27, 38: None}
TITE_TASKS = {
    "two_stage": {1: 1, 16: 16, 20: 20, 24: 24, 27: 27, 38: None},
    "one_stage": {1: 38, 16: 53, 20: 57, 24: 61, 27: 64, 38: None},
}
NON_TITE_TASKS = {
    "two_stage": {1: 1, 16: 2, 20: 3, 24: 4, 27: 5, 38: 6},
    "one_stage": {1: 7, 16: 8, 20: 9, 24: 10, 27: 11, 38: 12},
}
DESIGNS = [
    ("one_stage", "one-stage"),
    ("two_stage", "two-stage"),
]


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


def select_tite(raw: list[dict[str, str]], final_scenario: int, design_key: str) -> list[dict[str, str]] | None:
    source_scenario = TITE_SOURCE_SCENARIO[final_scenario]
    target_task = TITE_TASKS[design_key][final_scenario]
    if source_scenario is None or target_task is None:
        return None
    stage_one_max = 30 if design_key == "one_stage" else 6
    rows = [
        row for row in raw
        if int(row["Task_ID"]) == target_task
        and int(row["Scenario"]) == source_scenario
        and row["Allocation"] == design_key
        and int(row["s1_Max"]) == stage_one_max
        and row["Model_ID"] == "additive_shared"
        and row["Efficacy_Model"] == "previous_dose_additive"
        and int(row["Nmax"]) == 30
        and int(row["N_s2"]) == 30
        and int(row["n_eval"]) == 3
        and int(row["Cycle_Max"]) == 2
        and int(row["ntrial"]) == 1000
        and close(row["Arrival_Rate"], 1 / 56)
    ]
    rows.sort(key=lambda row: int(row["Dose"]))
    if len(rows) != 5 or [int(row["Dose"]) for row in rows] != DOSES:
        raise RuntimeError(f"TITE {design_key}, requested Scenario {final_scenario} is not a complete five-dose row set.")
    return rows


def select_non_tite(raw: list[dict[str, str]], final_scenario: int, design_key: str) -> list[dict[str, str]]:
    target_task = NON_TITE_TASKS[design_key][final_scenario]
    stage_one_max = 30 if design_key == "one_stage" else 6
    rows = [
        row for row in raw
        if int(row["Task_ID"]) == target_task
        and int(row["Scenario"]) == final_scenario
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
    rows.sort(key=lambda row: int(row["Dose"]))
    if len(rows) != 5 or [int(row["Dose"]) for row in rows] != DOSES:
        raise RuntimeError(f"IPDE-first non-TITE {design_key}, Scenario {final_scenario} is not a complete five-dose row set.")
    return rows


def normalize(rows: list[dict[str, str]] | None, version: str, design_key: str, short_name: str, final_scenario: int) -> dict:
    if rows is None:
        return {
            "scenario": final_scenario,
            "version": version,
            "design": design_key,
            "method": f"{version} | {short_name}",
            "available": False,
            "dlt": [],
            "efficacy": [],
            "selection": [],
            "administrations": [],
            "no_obd": None,
            "duration": None,
        }
    dlt = [number(row, "True_DLT_rate") for row in rows]
    efficacy = [number(row, "True_Efficacy_rate") for row in rows]
    selection = [number(row, "OBD_Selection_pct") for row in rows]
    return {
        "scenario": final_scenario,
        "version": version,
        "design": design_key,
        "method": f"{version} | {short_name}",
        "available": True,
        "dlt": dlt,
        "efficacy": efficacy,
        "selection": selection,
        "administrations": [number(row, "Pts_Treated") for row in rows],
        "no_obd": 100 - sum(selection) if version == "TITE" else number(rows[0], "No_OBD_Selection_pct"),
        "duration": number(rows[0], "Duration"),
    }


def build_summary() -> list[dict]:
    tite_raw = read_csv(TITE_CSV)
    non_tite_raw = read_csv(NON_TITE_CSV)
    records: list[dict] = []
    for final_scenario in SCENARIOS:
        scenario_records: list[dict] = []
        for design_key, short_name in DESIGNS:
            scenario_records.extend([
                normalize(select_tite(tite_raw, final_scenario, design_key), "TITE", design_key, short_name, final_scenario),
                normalize(select_non_tite(non_tite_raw, final_scenario, design_key), "IPDE-first non-TITE", design_key, short_name, final_scenario),
            ])
        available_records = [record for record in scenario_records if record["available"]]
        truth = (available_records[0]["dlt"], available_records[0]["efficacy"])
        if any((record["dlt"], record["efficacy"]) != truth for record in available_records[1:]):
            raise RuntimeError(f"Dose-level truth differs between source rows for requested Scenario {final_scenario}.")
        records.extend(scenario_records)
    return records


def record_for(records: list[dict], scenario: int, design: str, version: str) -> dict:
    return next(
        row for row in records
        if row["scenario"] == scenario and row["design"] == design and row["version"] == version
    )


def draw_scenario(pdf, y: float, scenario: int, records: list[dict], template) -> float:
    truth = record_for(records, scenario, "one_stage", "IPDE-first non-TITE")
    y = template.draw_row(pdf, y, 15, f"Scenario {scenario}", [""] * 6, fill=template.SECTION, bold=True)
    y = template.draw_row(pdf, y, 10.5, "DLT rate", [template.fmt(value, 2) for value in truth["dlt"]] + [""])
    y = template.draw_row(pdf, y, 10.5, "Efficacy rate", [template.fmt(value, 2) for value in truth["efficacy"]] + [""])
    y = template.draw_row(pdf, y, 12, "OBD selection (%)", [""] * 6, fill=template.SECTION, bold=True)
    for design_key, _ in DESIGNS:
        for version, fill, bold in [("TITE", template.TITE_FILL, False), ("IPDE-first non-TITE", template.NON_TITE_FILL, True)]:
            row = record_for(records, scenario, design_key, version)
            label = row["method"] if row["available"] else f"{row['method']} (N/A)"
            values = [template.fmt(value) for value in row["selection"]] + [template.fmt(row["no_obd"]) ] if row["available"] else ["N/A"] * 6
            y = template.draw_row(pdf, y, 10.5, label, values, fill=fill, bold=bold)
    y = template.draw_row(pdf, y, 12, "Mean administrations by dose", [""] * 6, fill=template.SECTION, bold=True)
    for design_key, _ in DESIGNS:
        for version, fill, bold in [("TITE", template.TITE_FILL, False), ("IPDE-first non-TITE", template.NON_TITE_FILL, True)]:
            row = record_for(records, scenario, design_key, version)
            label = row["method"] if row["available"] else f"{row['method']} (N/A)"
            values = [template.fmt(value) for value in row["administrations"]] + [""] if row["available"] else ["N/A"] * 5 + [""]
            y = template.draw_row(pdf, y, 10.5, label, values, fill=fill, bold=bold)
    y = template.draw_duration_header(pdf, y)
    for design_key, _ in DESIGNS:
        for version, fill, bold in [("TITE", template.TITE_FILL, False), ("IPDE-first non-TITE", template.NON_TITE_FILL, True)]:
            row = record_for(records, scenario, design_key, version)
            label = row["method"] if row["available"] else f"{row['method']} (N/A)"
            y = template.draw_duration_row(pdf, y, label, row["duration"], fill, bold)
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
    template.draw_text(pdf, "Sources: updated 8/23 TITE and IPDE-first non-TITE summaries; all populated rows represent 1,000 simulations.", template.LEFT, 27, size=6.2, color=template.MUTED)
    scenario_note = "The new TITE source has no full Scenario 38 truth-vector match; its two rows are shown as N/A." if scenario == 38 else "Both displayed designs are fully matched across the updated TITE and IPDE-first non-TITE sources."
    template.draw_text(pdf, scenario_note, template.LEFT, 21, size=6.1, color=template.MUTED)
    template.draw_text(pdf, "TITE No OBD % = 100 - sum of dose OBD selections. Duration is mean calendar days to complete the trial.", template.LEFT, 15, size=6.1, color=template.MUTED)
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
