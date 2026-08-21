from pathlib import Path

import numpy as np
import pandas as pd
from reportlab.lib.colors import HexColor, black
from reportlab.lib.pagesizes import landscape, letter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
PRESENTATION = ROOT / "Presentation 8-17-2026"
RAW_DATA = PRESENTATION / "Raw Data"
OUTPUT_DIR = PRESENTATION / "Table and Plots" / "TITE"
OUTPUT_PDF = OUTPUT_DIR / "phase12_tite_vs_nontite_3designs_N30_fut0p85_scenarios_1_16_20_24_27_38_tables.pdf"
TITE_FILE = RAW_DATA / "TITE_AIDE_phase_I_II_IDX_1001_to_2000_dose_summary.csv"
TITE_ONE_STAGE_REPLACEMENT_FILE = RAW_DATA / "TITE_AIDE_phase_I_II_IDX_1001_to_2000_newdesign_dose_summary.csv"
NON_TITE_FILE = RAW_DATA / (
    "AIDE_phase_I_II_modelsadditive_newdesign_shared_N30_ncycle2_"
    "rp0p15x0p85_rate56d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_"
    "ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv"
)
NON_TITE_ONE_STAGE_REPLACEMENT_FILE = RAW_DATA / (
    "AIDE_phase_I_II_modelsadditive_shared_N30_ncycle2_"
    "rp0p15x0p85_rate56d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_"
    "ap0p15x0p85_IDX_0001_to_1000_newdesign_dose_summary.csv"
)

SCENARIOS = [1, 16, 20, 24, 27, 38]
DOSES = [1, 2, 3, 4, 5]
ONE_STAGE_REPLACEMENT_TASKS = {scenario: index for index, scenario in enumerate(SCENARIOS, start=1)}
DESIGNS = [
    {
        "key": "one_stage",
        "short": "one-stage",
        "allocation": "one_stage",
        "stage2": "highest_utility",
        "tite_tasks": {1: 13, 16: 14, 20: 15, 24: 16, 27: 17, 38: 18},
        "non_tite_tasks": {1: 13, 16: 14, 20: 15, 24: 16, 27: 17, 38: 18},
    },
    {
        "key": "two_stage_highest_utility",
        "short": "2-stage highest utility",
        "allocation": "two_stage",
        "stage2": "highest_utility",
        "tite_tasks": {1: 7, 16: 8, 20: 9, 24: 10, 27: 11, 38: 12},
        "non_tite_tasks": {1: 7, 16: 8, 20: 9, 24: 10, 27: 11, 38: 12},
    },
    {
        "key": "two_stage_top2_randomized",
        "short": "2-stage top-2 randomized",
        "allocation": "two_stage",
        "stage2": "top2_randomized",
        "tite_tasks": {1: 1, 16: 2, 20: 3, 24: 4, 27: 5, 38: 6},
        "non_tite_tasks": {1: 1, 16: 2, 20: 3, 24: 4, 27: 5, 38: 6},
    },
]

PAGE_WIDTH, _ = landscape(letter)
LEFT = 36
RIGHT = PAGE_WIDTH - 36
TABLE_WIDTH = RIGHT - LEFT
LABEL_WIDTH = 255
VALUE_WIDTH = (TABLE_WIDTH - LABEL_WIDTH) / 6
NAVY = HexColor("#173F5F")
SECTION = HexColor("#E9ECEF")
TITE_FILL = HexColor("#E5F1FB")
NON_TITE_FILL = HexColor("#E6F3F1")
GRID = HexColor("#C9CDD1")
MUTED = HexColor("#4D5965")
NON_TITE_LABEL = "New non-TITE"
SUBTITLE = "TITE versus new non-TITE: one-stage, two-stage highest utility, and two-stage top-2 randomized"
SOURCE_NOTE = "Sources: updated TITE and new non-TITE one-stage newdesign summaries; existing two-stage summaries."
SIMULATION_NOTE = "Non-TITE rows represent 1,000 simulations; TITE rows use their available 996-1,000 trials."
DURATION_NOTE = "TITE No OBD % = 100 - sum of dose OBD selections. Duration is mean calendar days to complete the trial."


def fmt(value, digits=1):
    if value is None or pd.isna(value):
        return "-"
    text = f"{float(value):.{digits}f}"
    return text.rstrip("0").rstrip(".")


def require_columns(data, columns, source_name):
    missing = set(columns).difference(data.columns)
    if missing:
        raise ValueError(f"{source_name} is missing columns: {sorted(missing)}")


def tite_filter(data):
    return (
        (data["Model_ID"] == "additive_shared")
        & (data["Nmax"] == 30)
        & (data["N_s2"] == 30)
        & (data["n_eval"] == 3)
        & (data["Cycle_Max"] == 2)
        & (data["Utility_Type"] == 3)
        & np.isclose(data["Lambda_T"], 1)
        & np.isclose(data["Arrival_Rate"], 1 / 56)
    )


def non_tite_filter(data):
    return (
        (data["Model_ID"] == "additive_shared")
        & (data["Nmax"] == 30)
        & (data["N_s2"] == 30)
        & (data["Utility_Type"] == 3)
        & (data["Cycle_Max"] == 2)
        & np.isclose(data["Toxicity_IPDE_Alpha"], 0)
        & np.isclose(data["Efficacy_IPDE_Alpha"], 0)
        & np.isclose(data["Efficacy_Threshold"], 0.2)
        & np.isclose(data["Futility_Cutoff"], 0.85)
        & (data["Min_Eff_N_for_Futility"] == 3)
        & (data["n_valid"] == 1000)
    )


def select_rows(data, design, version, task_map=None):
    n_stage1 = 30 if design["key"] == "one_stage" else 6
    if task_map is None:
        task_map = design["tite_tasks"] if version == "TITE" else design["non_tite_tasks"]
    base = (
        data["Scenario"].isin(SCENARIOS)
        & (data["Task_ID"] == data["Scenario"].map(task_map))
        & (data["Allocation"] == design["allocation"])
        & (data["N_s1"] == n_stage1)
    )
    if version == "TITE":
        out = data.loc[base & tite_filter(data)].copy()
    else:
        out = data.loc[base & (data["Stage2_Allocation"] == design["stage2"]) & non_tite_filter(data)].copy()
    if set(out["Scenario"]) != set(SCENARIOS) or not (out.groupby("Scenario").size() == len(DOSES)).all():
        raise ValueError(f"{version} source does not contain complete rows for {design['short']}.")
    return out


def normalise(data, version, design, scenario):
    rows = data.loc[data["Scenario"] == scenario].sort_values("Dose")
    if rows["Dose"].tolist() != DOSES:
        raise ValueError(f"{version}, {design['short']}, scenario {scenario} does not contain doses 1 through 5.")
    no_obd = 100 - rows["OBD_Selection_pct"].sum() if version == "TITE" else rows["No_OBD_Selection_pct"].iloc[0]
    return {
        "Scenario": scenario,
        "Version": version,
        "Design_Key": design["key"],
        "Method": f"{version} | {design['short']}",
        "DLT": rows["True_DLT_rate"].tolist(),
        "Efficacy": rows["True_Efficacy_rate"].tolist(),
        "Selection": rows["OBD_Selection_pct"].tolist(),
        "Administrations": rows["Pts_Treated"].tolist(),
        "No_OBD": no_obd,
        "Duration": rows["Duration"].iloc[0],
    }


def build_summary():
    tite = pd.read_csv(TITE_FILE)
    tite_one_stage = pd.read_csv(TITE_ONE_STAGE_REPLACEMENT_FILE)
    non_tite = pd.read_csv(NON_TITE_FILE)
    non_tite_one_stage = pd.read_csv(NON_TITE_ONE_STAGE_REPLACEMENT_FILE)
    require_columns(tite, {"Task_ID", "Scenario", "Allocation", "N_s1", "N_s2", "Model_ID", "Utility_Type", "Lambda_T", "n_eval", "Cycle_Max", "Arrival_Rate", "Dose", "True_DLT_rate", "True_Efficacy_rate", "OBD_Selection_pct", "Pts_Treated", "Duration", "ntrial"}, TITE_FILE.name)
    require_columns(tite_one_stage, {"Task_ID", "Scenario", "Allocation", "N_s1", "N_s2", "Model_ID", "Utility_Type", "Lambda_T", "n_eval", "Cycle_Max", "Arrival_Rate", "Dose", "True_DLT_rate", "True_Efficacy_rate", "OBD_Selection_pct", "Pts_Treated", "Duration", "ntrial"}, TITE_ONE_STAGE_REPLACEMENT_FILE.name)
    require_columns(non_tite, {"Task_ID", "Scenario", "Allocation", "Stage2_Allocation", "N_s1", "N_s2", "Model_ID", "Utility_Type", "Toxicity_IPDE_Alpha", "Efficacy_IPDE_Alpha", "Efficacy_Threshold", "Futility_Cutoff", "Min_Eff_N_for_Futility", "Cycle_Max", "Dose", "True_DLT_rate", "True_Efficacy_rate", "OBD_Selection_pct", "Pts_Treated", "No_OBD_Selection_pct", "Duration", "n_valid"}, NON_TITE_FILE.name)
    require_columns(non_tite_one_stage, {"Task_ID", "Scenario", "Allocation", "Stage2_Allocation", "N_s1", "N_s2", "Model_ID", "Utility_Type", "Toxicity_IPDE_Alpha", "Efficacy_IPDE_Alpha", "Efficacy_Threshold", "Futility_Cutoff", "Min_Eff_N_for_Futility", "Cycle_Max", "Dose", "True_DLT_rate", "True_Efficacy_rate", "OBD_Selection_pct", "Pts_Treated", "No_OBD_Selection_pct", "Duration", "n_valid"}, NON_TITE_ONE_STAGE_REPLACEMENT_FILE.name)
    summary = []
    for design in DESIGNS:
        if design["key"] == "one_stage":
            selected_tite = select_rows(tite_one_stage, design, "TITE", ONE_STAGE_REPLACEMENT_TASKS)
            selected_non_tite = select_rows(non_tite_one_stage, design, NON_TITE_LABEL, ONE_STAGE_REPLACEMENT_TASKS)
        else:
            selected_tite = select_rows(tite, design, "TITE")
            selected_non_tite = select_rows(non_tite, design, NON_TITE_LABEL)
        for scenario in SCENARIOS:
            summary.append(normalise(selected_tite, "TITE", design, scenario))
            summary.append(normalise(selected_non_tite, NON_TITE_LABEL, design, scenario))
    for scenario in SCENARIOS:
        records = [row for row in summary if row["Scenario"] == scenario]
        if any(row["DLT"] != records[0]["DLT"] or row["Efficacy"] != records[0]["Efficacy"] for row in records[1:]):
            raise ValueError(f"Truth mismatch across TITE and non-TITE rows for scenario {scenario}.")
    return summary


def draw_text(pdf, text, x, y, size=8, bold=False, color=black, align="left"):
    pdf.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    pdf.setFillColor(color)
    if align == "center":
        pdf.drawCentredString(x, y, text)
    elif align == "right":
        pdf.drawRightString(x, y, text)
    else:
        pdf.drawString(x, y, text)


def draw_row(pdf, y_top, height, label, values, fill=None, bold=False):
    if fill is not None:
        pdf.setFillColor(fill)
        pdf.rect(LEFT, y_top - height, TABLE_WIDTH, height, fill=1, stroke=0)
    baseline = y_top - height + (height - 7.4) / 2 + 1.1
    draw_text(pdf, label, LEFT + 5, baseline, size=7.4, bold=bold)
    for index, value in enumerate(values):
        draw_text(pdf, value, LEFT + LABEL_WIDTH + VALUE_WIDTH * (index + 0.5), baseline, size=7.4, bold=bold, align="center")
    pdf.setStrokeColor(GRID)
    pdf.setLineWidth(0.35)
    pdf.line(LEFT, y_top - height, RIGHT, y_top - height)
    return y_top - height


def draw_duration_header(pdf, y_top):
    height = 12
    pdf.setFillColor(SECTION)
    pdf.rect(LEFT, y_top - height, TABLE_WIDTH, height, fill=1, stroke=0)
    baseline = y_top - height + (height - 7.4) / 2 + 1.1
    draw_text(pdf, "Mean trial duration (days)", LEFT + 5, baseline, size=7.4, bold=True)
    draw_text(pdf, "Mean duration (days)", LEFT + LABEL_WIDTH + (TABLE_WIDTH - LABEL_WIDTH) / 2, baseline, size=7.4, bold=True, align="center")
    pdf.setStrokeColor(GRID)
    pdf.setLineWidth(0.35)
    pdf.line(LEFT, y_top - height, RIGHT, y_top - height)
    return y_top - height


def draw_duration_row(pdf, y_top, label, duration, fill, bold):
    height = 10.5
    pdf.setFillColor(fill)
    pdf.rect(LEFT, y_top - height, TABLE_WIDTH, height, fill=1, stroke=0)
    baseline = y_top - height + (height - 7.4) / 2 + 1.1
    draw_text(pdf, label, LEFT + 5, baseline, size=7.4, bold=bold)
    draw_text(pdf, fmt(duration, 1), LEFT + LABEL_WIDTH + (TABLE_WIDTH - LABEL_WIDTH) / 2, baseline, size=7.4, bold=bold, align="center")
    pdf.setStrokeColor(GRID)
    pdf.setLineWidth(0.35)
    pdf.line(LEFT, y_top - height, RIGHT, y_top - height)
    return y_top - height


def scenario_records(summary, scenario, design):
    return {
        row["Version"]: row
        for row in summary
        if row["Scenario"] == scenario and row["Design_Key"] == design["key"]
    }


def draw_scenario(pdf, y_top, scenario, summary):
    first = scenario_records(summary, scenario, DESIGNS[0])[NON_TITE_LABEL]
    y = draw_row(pdf, y_top, 15, f"Scenario {scenario}", [""] * 6, fill=SECTION, bold=True)
    y = draw_row(pdf, y, 10.5, "DLT rate", [fmt(value, 2) for value in first["DLT"]] + [""])
    y = draw_row(pdf, y, 10.5, "Efficacy rate", [fmt(value, 2) for value in first["Efficacy"]] + [""])
    y = draw_row(pdf, y, 12, "OBD selection (%)", [""] * 6, fill=SECTION, bold=True)
    for design in DESIGNS:
        records = scenario_records(summary, scenario, design)
        for version, fill, bold in [("TITE", TITE_FILL, False), (NON_TITE_LABEL, NON_TITE_FILL, True)]:
            row = records[version]
            y = draw_row(pdf, y, 10.5, row["Method"], [fmt(value) for value in row["Selection"]] + [fmt(row["No_OBD"])], fill=fill, bold=bold)
    y = draw_row(pdf, y, 12, "Mean administrations by dose", [""] * 6, fill=SECTION, bold=True)
    for design in DESIGNS:
        records = scenario_records(summary, scenario, design)
        for version, fill, bold in [("TITE", TITE_FILL, False), (NON_TITE_LABEL, NON_TITE_FILL, True)]:
            row = records[version]
            y = draw_row(pdf, y, 10.5, row["Method"], [fmt(value) for value in row["Administrations"]] + [""], fill=fill, bold=bold)
    y = draw_duration_header(pdf, y)
    for design in DESIGNS:
        records = scenario_records(summary, scenario, design)
        for version, fill, bold in [("TITE", TITE_FILL, False), (NON_TITE_LABEL, NON_TITE_FILL, True)]:
            row = records[version]
            y = draw_duration_row(pdf, y, row["Method"], row["Duration"], fill, bold)
    return y


def draw_page(pdf, scenario, summary, page_number):
    draw_text(pdf, "Phase I/II Operating Characteristics", PAGE_WIDTH / 2, 576, size=20, bold=True, align="center")
    draw_text(pdf, SUBTITLE, PAGE_WIDTH / 2, 556, size=10.2, color=MUTED, align="center")
    pdf.setStrokeColor(NAVY)
    pdf.setLineWidth(0.8)
    pdf.line(LEFT, 545, RIGHT, 545)
    for index, header in enumerate(["Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "No OBD %"]):
        draw_text(pdf, header, LEFT + LABEL_WIDTH + VALUE_WIDTH * (index + 0.5), 531, size=8.8, bold=True, align="center")
    draw_scenario(pdf, 518, scenario, summary)
    draw_text(pdf, SOURCE_NOTE, LEFT, 27, size=6.2, color=MUTED)
    draw_text(pdf, SIMULATION_NOTE, LEFT, 21, size=6.1, color=MUTED)
    draw_text(pdf, DURATION_NOTE, LEFT, 15, size=6.1, color=MUTED)
    draw_text(pdf, f"Page {page_number} of {len(SCENARIOS)}", RIGHT, 15, size=6.2, color=MUTED, align="right")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = build_summary()
    pdf = canvas.Canvas(str(OUTPUT_PDF), pagesize=landscape(letter), pageCompression=1)
    for page_number, scenario in enumerate(SCENARIOS, start=1):
        draw_page(pdf, scenario, summary, page_number)
        pdf.showPage()
    pdf.save()
    print(OUTPUT_PDF)


if __name__ == "__main__":
    main()
