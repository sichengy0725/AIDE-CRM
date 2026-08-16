from pathlib import Path

import numpy as np
import pandas as pd
from reportlab.lib.colors import HexColor, black
from reportlab.lib.pagesizes import landscape, letter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
PRESENTATION_DIR = ROOT / "Presentation 8-17-2026"
OUTPUT_DIR = PRESENTATION_DIR / "Table and Plots"
STEM = "phase12_tite_vs_nontite_3designs_N30_fut0p85_scenarios_1_16_20_24_27_38"
OUTPUT_CSV = OUTPUT_DIR / f"{STEM}_comparison.csv"
OUTPUT_PDF = OUTPUT_DIR / f"{STEM}_tables.pdf"
TITE_FILE = PRESENTATION_DIR / "TITE_AIDE_phase_I_II_IDX_1001_to_2000_dose_summary.csv"
NEW_FILE = next(PRESENTATION_DIR.glob("AIDE_phase_I_II_modelsadditive_newdesign_shared_N30_ncycle2_*IDX_0001_to_1000_dose_summary.csv"))

SCENARIOS = [1, 16, 20, 24, 27, 38]
DOSES = [1, 2, 3, 4, 5]
DESIGNS = [
    {
        "key": "one_stage",
        "label": "One-stage (highest utility)",
        "short": "one-stage",
        "tite_tasks": {1: 13, 16: 14, 20: 15, 24: 16, 27: 17, 38: 18},
        "new_tasks": {1: 13, 16: 14, 20: 15, 24: 16, 27: 17, 38: 18},
        "allocation": "one_stage",
        "stage2": "highest_utility",
    },
    {
        "key": "two_stage_highest_utility",
        "label": "Two-stage (highest utility)",
        "short": "2-stage highest utility",
        "tite_tasks": {1: 7, 16: 8, 20: 9, 24: 10, 27: 11, 38: 12},
        "new_tasks": {1: 7, 16: 8, 20: 9, 24: 10, 27: 11, 38: 12},
        "allocation": "two_stage",
        "stage2": "highest_utility",
    },
    {
        "key": "two_stage_top2_randomized",
        "label": "Two-stage (top-2 randomized)",
        "short": "2-stage top-2 randomized",
        "tite_tasks": {1: 1, 16: 2, 20: 3, 24: 4, 27: 5, 38: 6},
        "new_tasks": {1: 1, 16: 2, 20: 3, 24: 4, 27: 5, 38: 6},
        "allocation": "two_stage",
        "stage2": "top2_randomized",
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
NEW_FILL = HexColor("#E6F3F1")
UNAVAILABLE_FILL = HexColor("#F9E6E6")
GRID = HexColor("#C9CDD1")
MUTED = HexColor("#4D5965")


def fmt(value, digits=1):
    if value is None or pd.isna(value):
        return "-"
    text = f"{float(value):.{digits}f}"
    return text.rstrip("0").rstrip(".")


def require_columns(data, columns, source_name):
    missing = set(columns).difference(data.columns)
    if missing:
        raise ValueError(f"{source_name} is missing columns: {sorted(missing)}")


def common_tite_filter(data):
    return (
        (data["Model_ID"] == "additive_shared")
        & (data["Nmax"] == 30)
        & (data["N_s2"] == 30)
        & (data["n_eval"] == 3)
        & (data["Cycle_Max"] == 2)
        & (data["Utility_Type"] == 3)
        & np.isclose(data["Lambda_T"], 1)
        & np.isclose(data["Arrival_Rate"], 1 / 56)
        & (data["ntrial"] == 1000)
    )


def common_new_filter(data):
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


def select_tite(data, design):
    task_by_scenario = design["tite_tasks"]
    out = data.loc[
        data["Scenario"].isin(SCENARIOS)
        & (data["Task_ID"] == data["Scenario"].map(task_by_scenario))
        & (data["Allocation"] == design["allocation"])
        & (data["N_s1"] == (30 if design["key"] == "one_stage" else 6))
        & common_tite_filter(data),
        :,
    ].copy()
    if set(out["Scenario"]) != set(SCENARIOS) or not (out.groupby("Scenario").size() == len(DOSES)).all():
        raise ValueError(f"TITE source does not contain one complete baseline row set for {design['label']}.")
    return out


def select_new(data, design):
    task_by_scenario = design["new_tasks"]
    out = data.loc[
        data["Scenario"].isin(SCENARIOS)
        & (data["Task_ID"] == data["Scenario"].map(task_by_scenario))
        & (data["Allocation"] == design["allocation"])
        & (data["Stage2_Allocation"] == design["stage2"])
        & (data["N_s1"] == (30 if design["key"] == "one_stage" else 6))
        & common_new_filter(data),
        :,
    ].copy()
    if set(out["Scenario"]) != set(SCENARIOS) or not (out.groupby("Scenario").size() == len(DOSES)).all():
        raise ValueError(f"New non-TITE source does not contain one complete baseline row set for {design['label']}.")
    return out


def normalized_row(data, source_type, design, scenario):
    method = "TITE" if source_type == "tite" else "New non-TITE"
    one = data.loc[data["Scenario"] == scenario].sort_values("Dose")
    if one.empty:
        return {
            "Scenario": scenario,
            "Version": method,
            "Design": design["label"],
            "Design_Key": design["key"],
            "Method": f"{method} | {design['short']}",
            "Available": False,
            "Availability_Note": "No matching result is available in the supplied source.",
            "Source": TITE_FILE.name,
            "Trials": np.nan,
            "True_MTD": np.nan,
            "True_OBD": np.nan,
            "No_OBD_Selection_pct": np.nan,
            "Mean_Total_Administrations": np.nan,
            "Mean_Total_Unique_Patients": np.nan,
            "Mean_Duration_days": np.nan,
            **{f"True_DLT_Dose_{dose}": np.nan for dose in DOSES},
            **{f"True_Efficacy_Dose_{dose}": np.nan for dose in DOSES},
            **{f"OBD_Selection_Dose_{dose}_pct": np.nan for dose in DOSES},
            **{f"Mean_Administrations_Dose_{dose}": np.nan for dose in DOSES},
        }
    if one["Dose"].tolist() != DOSES:
        raise ValueError(f"{method}, {design['label']}, scenario {scenario} does not contain doses 1 through 5.")
    no_obd = 100 - one["OBD_Selection_pct"].sum() if source_type == "tite" else one["No_OBD_Selection_pct"].iloc[0]
    row = {
        "Scenario": scenario,
        "Version": method,
        "Design": design["label"],
        "Design_Key": design["key"],
        "Method": f"{method} | {design['short']}",
        "Available": True,
        "Availability_Note": "",
        "Source": TITE_FILE.name if source_type == "tite" else NEW_FILE.name,
        "Trials": one["ntrial"].iloc[0] if source_type == "tite" else one["n_valid"].iloc[0],
        "True_MTD": one["True_MTD"].iloc[0],
        "True_OBD": one["True_OBD"].iloc[0],
        "No_OBD_Selection_pct": no_obd,
        "Mean_Total_Administrations": one["Total_Administrations"].iloc[0],
        "Mean_Total_Unique_Patients": one["Total_Unique_Patients"].iloc[0],
        "Mean_Duration_days": one["Duration"].iloc[0],
    }
    for dose in DOSES:
        dose_row = one.loc[one["Dose"] == dose].iloc[0]
        row.update({
            f"True_DLT_Dose_{dose}": dose_row["True_DLT_rate"],
            f"True_Efficacy_Dose_{dose}": dose_row["True_Efficacy_rate"],
            f"OBD_Selection_Dose_{dose}_pct": dose_row["OBD_Selection_pct"],
            f"Mean_Administrations_Dose_{dose}": dose_row["Pts_Treated"],
        })
    return row


def build_summary():
    tite_data = pd.read_csv(TITE_FILE)
    new_data = pd.read_csv(NEW_FILE)
    require_columns(tite_data, {"Task_ID", "Scenario", "Allocation", "N_s1", "N_s2", "Model_ID", "Utility_Type", "Lambda_T", "n_eval", "Cycle_Max", "Arrival_Rate", "Dose", "True_DLT_rate", "True_Efficacy_rate", "OBD_Selection_pct", "Pts_Treated", "Total_Administrations", "Total_Unique_Patients", "Duration", "ntrial"}, TITE_FILE.name)
    require_columns(new_data, {"Task_ID", "Scenario", "Allocation", "Stage2_Allocation", "N_s1", "N_s2", "Model_ID", "Utility_Type", "Toxicity_IPDE_Alpha", "Efficacy_IPDE_Alpha", "Efficacy_Threshold", "Futility_Cutoff", "Min_Eff_N_for_Futility", "Cycle_Max", "Dose", "True_DLT_rate", "True_Efficacy_rate", "OBD_Selection_pct", "Pts_Treated", "Total_Administrations", "Total_Unique_Patients", "No_OBD_Selection_pct", "Duration", "n_valid"}, NEW_FILE.name)
    rows = []
    for design in DESIGNS:
        tite = select_tite(tite_data, design)
        new = select_new(new_data, design)
        for scenario in SCENARIOS:
            rows.append(normalized_row(tite, "tite", design, scenario))
            rows.append(normalized_row(new, "new", design, scenario))
    summary = pd.DataFrame(rows)
    actual = summary.loc[summary["Available"]]
    for scenario in SCENARIOS:
        truth_columns = [column for column in summary.columns if column.startswith("True_DLT_") or column.startswith("True_Efficacy_")]
        values = actual.loc[actual["Scenario"] == scenario, truth_columns]
        if values.nunique(dropna=False).max() != 1:
            raise ValueError(f"Truth mismatch across design rows for scenario {scenario}.")
    return summary


def draw_text(c, text, x, y, size=8, bold=False, color=black, align="left"):
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.setFillColor(color)
    if align == "center":
        c.drawCentredString(x, y, text)
    elif align == "right":
        c.drawRightString(x, y, text)
    else:
        c.drawString(x, y, text)


def draw_row(c, y_top, height, label, values, fill=None, bold=False):
    if fill:
        c.setFillColor(fill)
        c.rect(LEFT, y_top - height, TABLE_WIDTH, height, fill=1, stroke=0)
    baseline = y_top - height + (height - 7.4) / 2 + 1.1
    draw_text(c, label, LEFT + 5, baseline, size=7.4, bold=bold)
    for index, value in enumerate(values):
        draw_text(c, value, LEFT + LABEL_WIDTH + VALUE_WIDTH * (index + .5), baseline, size=7.4, bold=bold, align="center")
    c.setStrokeColor(GRID)
    c.setLineWidth(.35)
    c.line(LEFT, y_top - height, RIGHT, y_top - height)
    return y_top - height


def values_for(row, prefix, suffix="", digits=1):
    return [fmt(row[f"{prefix}{dose}{suffix}"], digits) for dose in DOSES]


def draw_scenario(c, y_top, scenario, summary):
    y = draw_row(c, y_top, 15, f"Scenario {scenario}", [""] * 6, fill=SECTION, bold=True)
    truth = summary.loc[(summary["Scenario"] == scenario) & (summary["Version"] == "New non-TITE")].iloc[0]
    y = draw_row(c, y, 10.5, "DLT rate", values_for(truth, "True_DLT_Dose_", digits=2) + [""])
    y = draw_row(c, y, 10.5, "Efficacy rate", values_for(truth, "True_Efficacy_Dose_", digits=2) + [""])
    y = draw_row(c, y, 12, "OBD selection (%)", [""] * 6, fill=SECTION, bold=True)
    for design in DESIGNS:
        for version, fill, bold in [("TITE", TITE_FILL, False), ("New non-TITE", NEW_FILL, True)]:
            row = summary.loc[(summary["Scenario"] == scenario) & (summary["Design_Key"] == design["key"]) & (summary["Version"] == version)].iloc[0]
            label = f"{version} | {design['short']}"
            if row["Available"]:
                values = values_for(row, "OBD_Selection_Dose_", "_pct") + [fmt(row["No_OBD_Selection_pct"])]
                y = draw_row(c, y, 10.5, label, values, fill=fill, bold=bold)
            else:
                y = draw_row(c, y, 10.5, f"{label} (rerun)", ["N/A"] * 6, fill=UNAVAILABLE_FILL)
    y = draw_row(c, y, 12, "Mean administrations by dose", [""] * 6, fill=SECTION, bold=True)
    for design in DESIGNS:
        for version, fill, bold in [("TITE", TITE_FILL, False), ("New non-TITE", NEW_FILL, True)]:
            row = summary.loc[(summary["Scenario"] == scenario) & (summary["Design_Key"] == design["key"]) & (summary["Version"] == version)].iloc[0]
            label = f"{version} | {design['short']}"
            if row["Available"]:
                y = draw_row(c, y, 10.5, label, values_for(row, "Mean_Administrations_Dose_") + [""], fill=fill, bold=bold)
            else:
                y = draw_row(c, y, 10.5, f"{label} (rerun)", ["N/A"] * 5 + [""], fill=UNAVAILABLE_FILL)
    return y


def draw_page(c, scenarios, summary, page_number, page_count):
    draw_text(c, "Phase I/II Operating Characteristics", PAGE_WIDTH / 2, 576, size=20, bold=True, align="center")
    draw_text(c, "TITE versus new non-TITE: one-stage, two-stage highest utility, and two-stage top-2 randomized", PAGE_WIDTH / 2, 556, size=10.2, color=MUTED, align="center")
    c.setStrokeColor(NAVY)
    c.setLineWidth(.8)
    c.line(LEFT, 545, RIGHT, 545)
    for index, header in enumerate(["Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "No OBD %"]):
        draw_text(c, header, LEFT + LABEL_WIDTH + VALUE_WIDTH * (index + .5), 531, size=8.8, bold=True, align="center")
    y = 518
    for index, scenario in enumerate(scenarios):
        y = draw_scenario(c, y, scenario, summary)
        if index < len(scenarios) - 1:
            y -= 8
    draw_text(c, "Sources: TITE IDX 1001-2000 and new non-TITE IDX 0001-1000, both under Presentation 8-17. All populated rows represent 1,000 simulations.", LEFT, 27, size=6.2, color=MUTED)
    draw_text(c, "All three designs are present in both sources. TITE No OBD % = 100 - sum of dose OBD selections.", LEFT, 16, size=6.1, color=MUTED)
    draw_text(c, f"Page {page_number} of {page_count}", RIGHT, 16, size=6.2, color=MUTED, align="right")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = build_summary()
    summary.to_csv(OUTPUT_CSV, index=False)
    c = canvas.Canvas(str(OUTPUT_PDF), pagesize=landscape(letter), pageCompression=1)
    groups = [SCENARIOS[index:index + 2] for index in range(0, len(SCENARIOS), 2)]
    for page_number, scenarios in enumerate(groups, start=1):
        draw_page(c, scenarios, summary, page_number, len(groups))
        c.showPage()
    c.save()
    print(OUTPUT_CSV)
    print(OUTPUT_PDF)


if __name__ == "__main__":
    main()
