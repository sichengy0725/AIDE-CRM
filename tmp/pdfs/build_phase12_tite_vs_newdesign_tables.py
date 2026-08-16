from pathlib import Path

import numpy as np
import pandas as pd
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.pagesizes import landscape, letter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
PRESENTATION_DIR = ROOT / "Presentation 8-17-2026"
OUTPUT_DIR = PRESENTATION_DIR / "Table and Plots"
OUTPUT_STEM = "phase12_tite_vs_newdesign_N30_fut0p85_scenarios_1_16_20_24_27_38"
OUTPUT_PDF = OUTPUT_DIR / f"{OUTPUT_STEM}_tables.pdf"
OUTPUT_CSV = OUTPUT_DIR / f"{OUTPUT_STEM}_comparison.csv"

TITE_FILE = PRESENTATION_DIR / "TITE_AIDE_phase_I_II_IDX_1001_to_2000_dose_summary.csv"
NEW_FILE = next(
    PRESENTATION_DIR.glob(
        "AIDE_phase_I_II_modelsadditive_newdesign_shared_N30_ncycle2_*IDX_0001_to_1000_dose_summary.csv"
    )
)

SCENARIOS = [1, 16, 20, 24, 27, 38]
DOSES = [1, 2, 3, 4, 5]
TITE_TASK_BY_SCENARIO = {1: 1, 16: 2, 20: 3, 24: 4, 27: 5, 38: 6}

PAGE_WIDTH, PAGE_HEIGHT = landscape(letter)
LEFT = 36
RIGHT = PAGE_WIDTH - 36
TABLE_WIDTH = RIGHT - LEFT
LABEL_WIDTH = 250
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


def read_tite():
    data = pd.read_csv(TITE_FILE)
    require_columns(
        data,
        {
            "Task_ID", "Scenario", "Allocation", "Model_ID", "Nmax", "N_s1",
            "N_s2", "n_eval", "Cycle_Max", "Arrival_Rate", "Dose",
            "True_DLT_rate", "True_Efficacy_rate", "OBD_Selection_pct",
            "Pts_Treated", "Total_Administrations", "Total_Unique_Patients",
            "Duration", "ntrial",
        },
        TITE_FILE.name,
    )
    out = data.loc[
        # This Presentation 8-17 TITE batch stores the six requested scenarios
        # as Tasks 1--6 in order; remaining tasks are alternate configurations.
        data["Scenario"].isin(SCENARIOS)
        & (data["Task_ID"] == data["Scenario"].map(TITE_TASK_BY_SCENARIO))
        & (data["Allocation"] == "two_stage")
        & (data["Model_ID"] == "additive_shared")
        & (data["Nmax"] == 30)
        & (data["N_s1"] == 6)
        & (data["N_s2"] == 30)
        & (data["n_eval"] == 3)
        & (data["Cycle_Max"] == 2)
        & np.isclose(data["Arrival_Rate"], 1 / 56),
        :,
    ].copy()
    if set(out["Scenario"]) != set(SCENARIOS):
        raise ValueError("The selected TITE source does not contain all six requested scenarios.")
    counts = out.groupby("Scenario").size()
    if not (counts == len(DOSES)).all():
        raise ValueError("Every selected TITE scenario must contain one row per dose.")
    if not (out["ntrial"] == 1000).all():
        raise ValueError("The selected TITE rows do not all have 1,000 trials.")
    return out


def read_new():
    data = pd.read_csv(NEW_FILE)
    require_columns(
        data,
        {
            "Task_ID", "Scenario", "Allocation", "Stage2_Allocation", "Nmax",
            "N_s1", "N_s2", "Utility_Type", "Toxicity_IPDE_Alpha",
            "Efficacy_IPDE_Alpha", "Efficacy_Threshold", "Futility_Cutoff",
            "Min_Eff_N_for_Futility", "Cycle_Max", "Dose", "True_DLT_rate",
            "True_Efficacy_rate", "OBD_Selection_pct", "Pts_Treated",
            "Total_Administrations", "Total_Unique_Patients", "No_OBD_Selection_pct",
            "Duration", "n_valid",
        },
        NEW_FILE.name,
    )
    out = data.loc[
        # The new-design batch contains the requested scenarios as Tasks 1--6
        # in order.  Other tasks are alternate prior/model configurations.
        data["Scenario"].isin(SCENARIOS)
        & (data["Task_ID"].isin(range(1, len(SCENARIOS) + 1)))
        & (data["Allocation"] == "two_stage")
        & (data["Stage2_Allocation"] == "top2_randomized")
        & (data["Nmax"] == 30)
        & (data["N_s1"] == 6)
        & (data["N_s2"] == 30)
        & (data["Utility_Type"] == 3)
        & (data["Cycle_Max"] == 2)
        & np.isclose(data["Toxicity_IPDE_Alpha"], 0)
        & np.isclose(data["Efficacy_IPDE_Alpha"], 0)
        & np.isclose(data["Efficacy_Threshold"], 0.2)
        & np.isclose(data["Futility_Cutoff"], 0.85)
        & (data["Min_Eff_N_for_Futility"] == 3),
        :,
    ].copy()
    if set(out["Scenario"]) != set(SCENARIOS):
        raise ValueError("The selected new-design source does not contain all six requested scenarios.")
    counts = out.groupby("Scenario").size()
    if not (counts == len(DOSES)).all():
        raise ValueError("Every selected new-design scenario must contain one row per dose.")
    if not (out["n_valid"] == 1000).all():
        raise ValueError("The selected new-design rows do not all have 1,000 valid trials.")
    return out


def normalized_rows(tite, new):
    rows = []
    sources = [
        ("TITE design", tite, "tite"),
        ("New non-TITE design", new, "new"),
    ]
    for method, data, source_type in sources:
        for scenario in SCENARIOS:
            one = data.loc[data["Scenario"] == scenario].sort_values("Dose")
            if one.empty:
                rows.append(
                    {
                        "Scenario": scenario,
                        "Method": method,
                        "Available": False,
                        "Source": TITE_FILE.name,
                        "Trials": np.nan,
                        "True_MTD": np.nan,
                        "True_OBD": np.nan,
                        **{f"True_DLT_Dose_{dose}": np.nan for dose in DOSES},
                        **{f"True_Efficacy_Dose_{dose}": np.nan for dose in DOSES},
                        **{f"OBD_Selection_Dose_{dose}_pct": np.nan for dose in DOSES},
                        "No_OBD_Selection_pct": np.nan,
                        **{f"Mean_Administrations_Dose_{dose}": np.nan for dose in DOSES},
                        "Mean_Total_Administrations": np.nan,
                        "Mean_Total_Unique_Patients": np.nan,
                        "Mean_Duration_days": np.nan,
                    }
                )
                continue
            if one["Dose"].tolist() != DOSES:
                raise ValueError(f"{method}, scenario {scenario} does not contain doses 1 through 5.")
            no_obd = (
                100 - one["OBD_Selection_pct"].sum()
                if source_type == "tite"
                else one["No_OBD_Selection_pct"].iloc[0]
            )
            row = {
                "Scenario": scenario,
                "Method": method,
                "Available": True,
                "Source": TITE_FILE.name if source_type == "tite" else NEW_FILE.name,
                "Trials": one["ntrial"].iloc[0] if source_type == "tite" else one["n_valid"].iloc[0],
                "True_MTD": one["True_MTD"].iloc[0],
                "True_OBD": one["True_OBD"].iloc[0] if "True_OBD" in one else np.nan,
                "No_OBD_Selection_pct": no_obd,
                "Mean_Total_Administrations": one["Total_Administrations"].iloc[0],
                "Mean_Total_Unique_Patients": one["Total_Unique_Patients"].iloc[0],
                "Mean_Duration_days": one["Duration"].iloc[0],
            }
            for index, dose in enumerate(DOSES):
                source_row = one.iloc[index]
                row[f"True_DLT_Dose_{dose}"] = source_row["True_DLT_rate"]
                row[f"True_Efficacy_Dose_{dose}"] = source_row["True_Efficacy_rate"]
                row[f"OBD_Selection_Dose_{dose}_pct"] = source_row["OBD_Selection_pct"]
                row[f"Mean_Administrations_Dose_{dose}"] = source_row["Pts_Treated"]
            rows.append(row)
    return pd.DataFrame(rows)


def validate_truth(summary):
    comparable = summary.loc[summary["Available"]].copy()
    for scenario in SCENARIOS[:-1]:
        rows = comparable.loc[comparable["Scenario"] == scenario]
        for prefix in ("True_DLT_Dose_", "True_Efficacy_Dose_"):
            values = rows[[f"{prefix}{dose}" for dose in DOSES]].to_numpy()
            if values.shape != (2, 5) or not np.allclose(values[0], values[1]):
                raise ValueError(f"TITE and new non-TITE truth values differ for scenario {scenario}.")


def draw_text(c, text, x, y, size=8, color=black, bold=False, align="left"):
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    c.setFillColor(color)
    if align == "center":
        c.drawCentredString(x, y, text)
    elif align == "right":
        c.drawRightString(x, y, text)
    else:
        c.drawString(x, y, text)


def draw_row(c, y_top, height, label, values, fill=None, bold=False, label_color=black):
    if fill:
        c.setFillColor(fill)
        c.rect(LEFT, y_top - height, TABLE_WIDTH, height, fill=1, stroke=0)
    baseline = y_top - height + (height - 7.4) / 2 + 1.1
    draw_text(c, label, LEFT + 5, baseline, size=7.4, color=label_color, bold=bold)
    for index, value in enumerate(values):
        x_center = LEFT + LABEL_WIDTH + VALUE_WIDTH * (index + 0.5)
        draw_text(c, value, x_center, baseline, size=7.4, bold=bold, align="center")
    c.setStrokeColor(GRID)
    c.setLineWidth(0.35)
    c.line(LEFT, y_top - height, RIGHT, y_top - height)
    return y_top - height


def values_for(row, prefix, suffix="", digits=1):
    return [fmt(row[f"{prefix}{dose}{suffix}"], digits) for dose in DOSES]


def draw_scenario(c, y_top, scenario, summary):
    y = draw_row(c, y_top, 15, f"Scenario {scenario}", [""] * 6, fill=SECTION, bold=True)
    new = summary.loc[(summary["Scenario"] == scenario) & (summary["Method"] == "New non-TITE design")].iloc[0]
    tite = summary.loc[(summary["Scenario"] == scenario) & (summary["Method"] == "TITE design")].iloc[0]
    y = draw_row(c, y, 10.5, "DLT rate", values_for(new, "True_DLT_Dose_", digits=2) + [""])
    y = draw_row(c, y, 10.5, "Efficacy rate", values_for(new, "True_Efficacy_Dose_", digits=2) + [""])
    y = draw_row(c, y, 12, "OBD selection (%)", [""] * 6, fill=SECTION, bold=True)
    if tite["Available"]:
        tite_values = values_for(tite, "OBD_Selection_Dose_", suffix="_pct") + [fmt(tite["No_OBD_Selection_pct"])]
        y = draw_row(c, y, 10.5, "TITE design", tite_values, fill=TITE_FILL)
    else:
        y = draw_row(c, y, 10.5, "TITE design - unavailable", ["N/A"] * 6, fill=UNAVAILABLE_FILL)
    new_values = values_for(new, "OBD_Selection_Dose_", suffix="_pct") + [fmt(new["No_OBD_Selection_pct"])]
    y = draw_row(c, y, 10.5, "New non-TITE design", new_values, fill=NEW_FILL, bold=True)
    y = draw_row(c, y, 12, "Mean administrations by dose", [""] * 6, fill=SECTION, bold=True)
    if tite["Available"]:
        y = draw_row(c, y, 10.5, "TITE design", values_for(tite, "Mean_Administrations_Dose_", digits=1) + [""], fill=TITE_FILL)
    else:
        y = draw_row(c, y, 10.5, "TITE design - unavailable", ["N/A"] * 5 + [""], fill=UNAVAILABLE_FILL)
    y = draw_row(c, y, 10.5, "New non-TITE design", values_for(new, "Mean_Administrations_Dose_", digits=1) + [""], fill=NEW_FILL, bold=True)
    return y


def draw_page(c, scenarios, summary, page_number, page_count):
    draw_text(c, "Phase I/II Operating Characteristics", PAGE_WIDTH / 2, 576, size=20, bold=True, align="center")
    draw_text(
        c,
        "TITE versus new non-TITE design - N = 30, Cycle_Max = 2, efficacy futility cutoff = 0.85",
        PAGE_WIDTH / 2,
        556,
        size=10.5,
        color=MUTED,
        align="center",
    )
    c.setStrokeColor(NAVY)
    c.setLineWidth(0.8)
    c.line(LEFT, 545, RIGHT, 545)
    headers = ["Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "No OBD %"]
    for index, header in enumerate(headers):
        x_center = LEFT + LABEL_WIDTH + VALUE_WIDTH * (index + 0.5)
        draw_text(c, header, x_center, 531, size=8.8, bold=True, align="center")
    y = 518
    for index, scenario in enumerate(scenarios):
        y = draw_scenario(c, y, scenario, summary)
        if index < len(scenarios) - 1:
            y -= 10
    draw_text(
        c,
        "TITE source: 1,000 simulations (IDX 1001-2000) from Presentation 8-17. New non-TITE source: 1,000 simulations (IDX 0001-1000).",
        LEFT,
        27,
        size=6.4,
        color=MUTED,
    )
    draw_text(
        c,
        "All six scenarios are available in both sources. TITE No OBD % = 100 - sum of dose-specific OBD selections.",
        LEFT,
        16,
        size=6.4,
        color=MUTED,
    )
    draw_text(c, f"Page {page_number} of {page_count}", RIGHT, 16, size=6.4, color=MUTED, align="right")


def write_pdf(summary):
    c = canvas.Canvas(str(OUTPUT_PDF), pagesize=landscape(letter), pageCompression=1)
    page_groups = [SCENARIOS[index:index + 2] for index in range(0, len(SCENARIOS), 2)]
    for page_number, scenarios in enumerate(page_groups, start=1):
        draw_page(c, scenarios, summary, page_number, len(page_groups))
        c.showPage()
    c.save()


def main():
    tite = read_tite()
    new = read_new()
    summary = normalized_rows(tite, new)
    validate_truth(summary)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUTPUT_CSV, index=False)
    write_pdf(summary)
    print(OUTPUT_CSV)
    print(OUTPUT_PDF)


if __name__ == "__main__":
    main()
