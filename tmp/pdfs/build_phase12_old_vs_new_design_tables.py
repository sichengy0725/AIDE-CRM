from pathlib import Path

import numpy as np
import pandas as pd
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.pagesizes import landscape, letter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = ROOT / "Presentation 8-17-2026"
OUTPUT = (
    SOURCE_DIR
    / "Table and Plots"
    / "phase12_old_vs_new_design_N30_fut0p85_scenarios_1_16_20_24_27_38_tables.pdf"
)

OLD_FILE = next(
    SOURCE_DIR.glob(
        "AIDE_phase_I_II_modelsadditive_shared_N30_ncycle2_*IDX_0001_to_1000_dose_summary.csv"
    )
)
NEW_FILE = next(
    SOURCE_DIR.glob(
        "AIDE_phase_I_II_modelsadditive_newdesign_shared_N30_ncycle2_*IDX_0001_to_1000_dose_summary.csv"
    )
)

SCENARIOS = [1, 16, 20, 24, 27, 38]
METHODS = [
    ("one_stage", "highest_utility", "One-stage"),
    ("two_stage", "highest_utility", "Two-stage (highest utility)"),
    ("two_stage", "top2_randomized", "Two-stage (top-2 randomized)"),
]

PAGE_WIDTH, PAGE_HEIGHT = landscape(letter)
LEFT = 36
RIGHT = PAGE_WIDTH - 36
TABLE_WIDTH = RIGHT - LEFT
LABEL_WIDTH = 240
VALUE_WIDTH = (TABLE_WIDTH - LABEL_WIDTH) / 6

NAVY = HexColor("#173F5F")
SECTION = HexColor("#E9ECEF")
OLD_FILL = HexColor("#F6F7F8")
NEW_FILL = HexColor("#E6F3F1")
GRID = HexColor("#C9CDD1")
MUTED = HexColor("#4D5965")


def fmt(value, digits=1):
    if value == "":
        return ""
    text = f"{float(value):.{digits}f}"
    return text.rstrip("0").rstrip(".")


def utility3(toxicity, efficacy):
    return (
        40 * (1 - toxicity) * (1 - efficacy)
        + 100 * (1 - toxicity) * efficacy
        + 60 * toxicity * efficacy
    )


def read_summary(path, design_label):
    data = pd.read_csv(path)
    required = {
        "Scenario",
        "Allocation",
        "Stage2_Allocation",
        "Nmax",
        "Cycle_Max",
        "Utility_Type",
        "Futility_Cutoff",
        "Toxicity_IPDE_Alpha",
        "Efficacy_IPDE_Alpha",
        "Dose",
        "True_DLT_rate",
        "True_Efficacy_rate",
        "OBD_Selection_pct",
        "Pts_Treated",
        "No_OBD_Selection_pct",
        "n_valid",
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"{path.name} is missing columns: {sorted(missing)}")

    selected = data.loc[
        (data["Scenario"].isin(SCENARIOS))
        & (data["Nmax"] == 30)
        & (data["Cycle_Max"] == 2)
        & (data["Utility_Type"] == 3)
        & np.isclose(data["Toxicity_IPDE_Alpha"], 0)
        & np.isclose(data["Efficacy_IPDE_Alpha"], 0)
        & (data["Futility_Cutoff"] == 0.85),
        :,
    ].copy()
    selected["Design"] = design_label

    expected_rows = len(SCENARIOS) * len(METHODS) * 5
    if len(selected) != expected_rows:
        raise ValueError(
            f"{path.name} produced {len(selected)} selected rows; expected {expected_rows}."
        )
    observed = set(zip(selected["Allocation"], selected["Stage2_Allocation"]))
    expected = {(allocation, stage2) for allocation, stage2, _ in METHODS}
    if observed != expected:
        raise ValueError(f"{path.name} has unexpected allocation settings: {observed}")
    counts = selected.groupby(["Scenario", "Allocation", "Stage2_Allocation"]).size()
    if not (counts == 5).all():
        raise ValueError(f"{path.name} does not contain exactly five doses per setting.")
    if not (selected["n_valid"] == 1000).all():
        raise ValueError(f"{path.name} does not contain 1,000 valid trials for every row.")
    return selected


def validate_truth(old, new):
    truth_columns = ["True_DLT_rate", "True_Efficacy_rate"]
    join_columns = ["Scenario", "Allocation", "Stage2_Allocation", "Dose"]
    merged = old[join_columns + truth_columns].merge(
        new[join_columns + truth_columns], on=join_columns, suffixes=("_old", "_new")
    )
    for column in truth_columns:
        if not np.allclose(merged[f"{column}_old"], merged[f"{column}_new"]):
            raise ValueError(f"Truth values differ between old and new summaries for {column}.")


def group_rows(data, scenario, allocation, stage2):
    rows = data.loc[
        (data["Scenario"] == scenario)
        & (data["Allocation"] == allocation)
        & (data["Stage2_Allocation"] == stage2)
    ].sort_values("Dose")
    if rows["Dose"].tolist() != [1, 2, 3, 4, 5]:
        raise ValueError(f"Scenario {scenario} lacks one or more dose rows.")
    return rows


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


def draw_scenario(c, y_top, scenario, old, new):
    y = draw_row(c, y_top, 15, f"Scenario {scenario}", [""] * 6, fill=SECTION, bold=True)

    first = group_rows(old, scenario, "one_stage", "highest_utility")
    dlt = [fmt(value, 2) for value in first["True_DLT_rate"]] + [""]
    efficacy = [fmt(value, 2) for value in first["True_Efficacy_rate"]] + [""]
    utility = [fmt(value, 1) for value in utility3(first["True_DLT_rate"], first["True_Efficacy_rate"])] + [""]
    y = draw_row(c, y, 10.5, "DLT rate", dlt)
    y = draw_row(c, y, 10.5, "Efficacy rate", efficacy)
    y = draw_row(c, y, 10.5, "Utility 3", utility)

    y = draw_row(c, y, 12, "OBD selection (%)", [""] * 6, fill=SECTION, bold=True)
    for allocation, stage2, label in METHODS:
        for design_label, data, fill in (("Old", old, OLD_FILL), ("New", new, NEW_FILL)):
            rows = group_rows(data, scenario, allocation, stage2)
            stop_values = rows["No_OBD_Selection_pct"].to_numpy()
            if not np.allclose(stop_values, stop_values[0]):
                raise ValueError(f"Scenario {scenario} {label} {design_label} has inconsistent stop percentages.")
            values = [fmt(value, 1) for value in rows["OBD_Selection_pct"]] + [fmt(stop_values[0], 1)]
            y = draw_row(c, y, 10.5, f"{label} - {design_label}", values, fill=fill, bold=(design_label == "New"))

    y = draw_row(c, y, 12, "Mean patients treated", [""] * 6, fill=SECTION, bold=True)
    for allocation, stage2, label in METHODS:
        for design_label, data, fill in (("Old", old, OLD_FILL), ("New", new, NEW_FILL)):
            rows = group_rows(data, scenario, allocation, stage2)
            values = [fmt(value, 1) for value in rows["Pts_Treated"]] + [""]
            y = draw_row(c, y, 10.5, f"{label} - {design_label}", values, fill=fill, bold=(design_label == "New"))
    return y


def draw_page(c, scenarios, old, new, page_number, page_count):
    draw_text(c, "Phase I/II Operating Characteristics", PAGE_WIDTH / 2, 576, size=20, bold=True, align="center")
    draw_text(
        c,
        "Old versus new allocation design - N = 30, Cycle_Max = 2, efficacy futility cutoff = 0.85",
        PAGE_WIDTH / 2,
        556,
        size=10.5,
        color=MUTED,
        align="center",
    )
    c.setStrokeColor(NAVY)
    c.setLineWidth(0.8)
    c.line(LEFT, 545, RIGHT, 545)
    headers = ["Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "Stop %"]
    for index, header in enumerate(headers):
        x_center = LEFT + LABEL_WIDTH + VALUE_WIDTH * (index + 0.5)
        draw_text(c, header, x_center, 531, size=8.8, bold=True, align="center")

    y = 518
    for index, scenario in enumerate(scenarios):
        y = draw_scenario(c, y, scenario, old, new)
        if index < len(scenarios) - 1:
            y -= 10

    draw_text(
        c,
        "Source: paired old/new 1,000-trial dose summaries; shared-alpha additive model, continuous enrollment, IPDE alpha = 0.",
        LEFT,
        27,
        size=6.7,
        color=MUTED,
    )
    draw_text(
        c,
        "Utility 3 scores: no toxicity/no efficacy = 40; no toxicity/efficacy = 100; toxicity/no efficacy = 0; toxicity/efficacy = 60. Stop % is no-OBD selection.",
        LEFT,
        16,
        size=6.7,
        color=MUTED,
    )
    draw_text(c, f"Page {page_number} of {page_count}", RIGHT, 16, size=6.7, color=MUTED, align="right")


def main():
    old = read_summary(OLD_FILE, "Old")
    new = read_summary(NEW_FILE, "New")
    validate_truth(old, new)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    c = canvas.Canvas(str(OUTPUT), pagesize=landscape(letter), pageCompression=1)
    page_groups = [SCENARIOS[index : index + 2] for index in range(0, len(SCENARIOS), 2)]
    for page_number, scenarios in enumerate(page_groups, start=1):
        draw_page(c, scenarios, old, new, page_number, len(page_groups))
        c.showPage()
    c.save()
    print(OUTPUT)


if __name__ == "__main__":
    main()
