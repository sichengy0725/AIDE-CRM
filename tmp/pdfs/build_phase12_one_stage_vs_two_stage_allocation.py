"""Build the August 17 Phase I/II allocation-strategy comparison tables."""

from __future__ import annotations

import csv
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
SOURCE = next((ROOT / "Presentation 8-17-2026").glob("*dose_summary.csv"))
OUT_DIR = ROOT / "Presentation 8-17-2026" / "Table and Plots"
OUTPUT = OUT_DIR / (
    "phase12_one_stage_vs_two_stage_allocation_N30_fut0p85_"
    "scenarios_1_16_20_24_27_38_tables.pdf"
)

SCENARIOS = (1, 16, 20, 24, 27, 38)
METHODS = (
    ("AIDE one-stage", "one_stage", "highest_utility"),
    ("AIDE two-stage (highest utility)", "two_stage", "highest_utility"),
    ("AIDE two-stage (top-2 randomized)", "two_stage", "top2_randomized"),
)

PAGE_WIDTH, PAGE_HEIGHT = 18 * inch, 14 * inch
LEFT = 80
LABEL_WIDTH = 330
DOSE_WIDTH = 130
STOP_WIDTH = 156
TABLE_WIDTH = LABEL_WIDTH + 5 * DOSE_WIDTH + STOP_WIDTH


def number(value: str) -> float:
    return float(value)


def display(value: float, decimals: int = 1) -> str:
    text = f"{value:.{decimals}f}"
    return text.rstrip("0").rstrip(".")


def utility3(toxicity: float, efficacy: float) -> float:
    return (
        40 * (1 - toxicity) * (1 - efficacy)
        + 100 * (1 - toxicity) * efficacy
        + 60 * toxicity * efficacy
    )


def read_records() -> tuple[dict[int, dict[str, list[float]]], dict[str, dict[int, dict[str, list[float] | float]]]]:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))

    truth: dict[int, dict[str, list[float]]] = {}
    records: dict[str, dict[int, dict[str, list[float] | float]]] = {
        label: {} for label, _, _ in METHODS
    }

    for scenario in SCENARIOS:
        scenario_rows = [row for row in rows if int(row["Scenario"]) == scenario]
        alpha_zero_rows = [
            row for row in scenario_rows
            if number(row["Toxicity_IPDE_Alpha"]) == 0
            and number(row["Efficacy_IPDE_Alpha"]) == 0
        ]
        first = alpha_zero_rows[0]
        dlt = [number(next(row["True_DLT_rate"] for row in alpha_zero_rows if int(row["Dose"]) == dose)) for dose in range(1, 6)]
        efficacy = [number(next(row["True_Efficacy_rate"] for row in alpha_zero_rows if int(row["Dose"]) == dose)) for dose in range(1, 6)]
        truth[scenario] = {
            "dlt": dlt,
            "efficacy": efficacy,
            "utility3": [utility3(p, e) for p, e in zip(dlt, efficacy)],
        }

        for label, allocation, stage2_allocation in METHODS:
            selected = [
                row for row in alpha_zero_rows
                if row["Allocation"] == allocation
                and row["Stage2_Allocation"] == stage2_allocation
            ]
            selected = sorted(selected, key=lambda row: int(row["Dose"]))
            records[label][scenario] = {
                "selected": [number(row["OBD_Selection_pct"]) for row in selected],
                "treated": [number(row["Pts_Treated"]) for row in selected],
                "stop": number(selected[0]["No_OBD_Selection_pct"]),
                "ntrial": number(first["ntrial_from_files"]),
            }
    return truth, records


def draw_row(pdf: canvas.Canvas, y: float, label: str, values: list[float], stop: float | None, decimals: int = 1) -> float:
    row_height = 13.2
    pdf.setFillColor(colors.HexColor("#202020"))
    pdf.setFont("Helvetica", 9.5)
    pdf.drawString(LEFT + 8, y - 9.6, label)
    x = LEFT + LABEL_WIDTH
    for value in values:
        pdf.drawCentredString(x + DOSE_WIDTH / 2, y - 9.6, display(value, decimals))
        x += DOSE_WIDTH
    if stop is not None:
        pdf.drawCentredString(x + STOP_WIDTH / 2, y - 9.6, display(stop))
    pdf.setStrokeColor(colors.HexColor("#D0D0D0"))
    pdf.setLineWidth(0.45)
    pdf.line(LEFT, y - row_height, LEFT + TABLE_WIDTH, y - row_height)
    return y - row_height


def draw_group_header(pdf: canvas.Canvas, y: float, label: str) -> float:
    height = 14
    pdf.setFillColor(colors.HexColor("#EFEFEF"))
    pdf.rect(LEFT, y - height + 2, TABLE_WIDTH, height, fill=1, stroke=0)
    pdf.setFillColor(colors.HexColor("#202020"))
    pdf.setFont("Helvetica-Bold", 10.5)
    pdf.drawString(LEFT + 8, y - 9, label)
    return y - height


def draw_scenario(
    pdf: canvas.Canvas,
    y: float,
    scenario: int,
    truth: dict[str, list[float]],
    records: dict[str, dict[str, list[float] | float]],
) -> float:
    height = 18
    pdf.setFillColor(colors.HexColor("#EFEFEF"))
    pdf.rect(LEFT, y - height + 2, TABLE_WIDTH, height, fill=1, stroke=0)
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 13)
    pdf.drawString(LEFT + 8, y - 12, f"Scenario {scenario}")
    y -= 19
    y = draw_row(pdf, y, "DLT rate", truth["dlt"], None, decimals=2)
    y = draw_row(pdf, y, "Efficacy rate", truth["efficacy"], None, decimals=2)
    y = draw_row(pdf, y, "Utility 3", truth["utility3"], None)
    y -= 2
    y = draw_group_header(pdf, y, "OBD selection (%)")
    for label, _, _ in METHODS:
        record = records[label]
        y = draw_row(pdf, y, label, record["selected"], record["stop"])
    y = draw_group_header(pdf, y, "Mean patients treated")
    for label, _, _ in METHODS:
        record = records[label]
        y = draw_row(pdf, y, label, record["treated"], None)
    return y - 9


def draw_header(pdf: canvas.Canvas, page_number: int, page_count: int) -> None:
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 28)
    pdf.drawCentredString(PAGE_WIDTH / 2, PAGE_HEIGHT - 56, "Phase I/II Operating Characteristics")
    pdf.setFont("Helvetica", 14)
    pdf.drawCentredString(
        PAGE_WIDTH / 2,
        PAGE_HEIGHT - 78,
        "One-stage versus two-stage allocation strategies - N = 30, Cycle_Max = 2, efficacy futility cutoff = 0.85",
    )
    pdf.setStrokeColor(colors.HexColor("#999999"))
    pdf.setLineWidth(0.8)
    pdf.line(LEFT, PAGE_HEIGHT - 92, PAGE_WIDTH - LEFT, PAGE_HEIGHT - 92)
    x = LEFT + LABEL_WIDTH
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 13)
    for dose in range(1, 6):
        pdf.drawCentredString(x + DOSE_WIDTH / 2, PAGE_HEIGHT - 110, f"Dose {dose}")
        x += DOSE_WIDTH
    pdf.drawCentredString(x + STOP_WIDTH / 2, PAGE_HEIGHT - 110, "Stop %")
    pdf.setFont("Helvetica", 7.1)
    pdf.drawString(
        LEFT,
        38,
        "All methods use the shared-alpha additive model, continuous enrollment, and toxicity/efficacy IPDE alpha = 0.",
    )
    pdf.drawString(
        LEFT,
        25,
        "Utility 3 scores: no toxicity/no efficacy = 40; no toxicity/efficacy = 100; toxicity/no efficacy = 0; toxicity/efficacy = 60. Stop % is no-OBD selection.",
    )
    pdf.drawRightString(PAGE_WIDTH - LEFT, 25, f"Page {page_number} of {page_count}")


def main() -> None:
    truth, records = read_records()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pdf = canvas.Canvas(str(OUTPUT), pagesize=(PAGE_WIDTH, PAGE_HEIGHT), pageCompression=1)
    page_count = len(SCENARIOS) // 2
    for page_number, start in enumerate(range(0, len(SCENARIOS), 2), start=1):
        draw_header(pdf, page_number, page_count)
        y = PAGE_HEIGHT - 130
        for scenario in SCENARIOS[start:start + 2]:
            y = draw_scenario(pdf, y, scenario, truth[scenario], {
                label: records[label][scenario] for label, _, _ in METHODS
            })
        pdf.showPage()
    pdf.save()
    print(OUTPUT)


if __name__ == "__main__":
    main()
