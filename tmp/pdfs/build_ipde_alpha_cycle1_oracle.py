"""Add cycle-1 oracle rows to the August 3 IPDE-alpha comparison tables.

The IPDE-alpha rows are read directly from the existing August 3 comparison
PDFs.  The oracle rows are the matching AIDE one-stage/two-stage records in
the August 3 all-method PDF, whose Cycle_Max is 1.  Cycle-1 has no IPDE
administrations, so its mean-IPDE row is zero by design.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import pdfplumber
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = ROOT / "Presentation 8-03-2026" / "Table and Plots" / "Phase I-II"
OUT_DIR = ROOT / "Presentation 8-10-2026" / "Table and Plots" / "Phase I-II"
ALL_METHODS_SOURCE = SOURCE_DIR / "phase12_all_methods_N30_scenarios_1_to_37_tables.pdf"
ALPHAS = (0.0, 0.3, 0.6, 0.9)
SCENARIOS = set(range(1, 38))

PAGE_WIDTH, PAGE_HEIGHT = 18 * inch, 14 * inch
LEFT = 80
LABEL_WIDTH = 330
DOSE_WIDTH = 130
STOP_WIDTH = 156
TABLE_WIDTH = LABEL_WIDTH + 5 * DOSE_WIDTH + STOP_WIDTH
ROW_HEIGHT = 9.6
ROW_FONT = 7.4
ROW_BASELINE = 7.2


def numbers(text: str) -> list[float]:
    return [float(value) for value in re.findall(r"-?(?:\d+\.?\d*|\.\d+)", text)]


def display(value: float, decimals: int = 1) -> str:
    if not math.isfinite(value):
        return "-"
    text = f"{value:.{decimals}f}"
    return text.rstrip("0").rstrip(".")


def alpha_label(alpha: float) -> str:
    return display(alpha, 1)


def pdf_text(path: Path) -> str:
    with pdfplumber.open(path) as document:
        return "\n".join(page.extract_text() or "" for page in document.pages)


def scenario_blocks(text: str) -> dict[int, list[str]]:
    starts = list(re.finditer(r"^Scenario\s+(\d+)\s*$", text, flags=re.MULTILINE))
    blocks: dict[int, list[str]] = {}
    for index, match in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
        blocks[int(match.group(1))] = [line.strip() for line in text[match.end():end].splitlines() if line.strip()]
    if set(blocks) != SCENARIOS:
        raise ValueError(f"Expected scenarios 1-37; found {sorted(blocks)}.")
    return blocks


def section(lines: list[str], start: str, end: str | None) -> list[str]:
    try:
        first = lines.index(start) + 1
    except ValueError as error:
        raise ValueError(f"Missing section {start!r}.") from error
    last = lines.index(end, first) if end is not None else len(lines)
    return lines[first:last]


def values_after(lines: list[str], label: str, expected: int) -> list[float]:
    for line in lines:
        if line.startswith(label + " "):
            values = numbers(line[len(label):])
            if len(values) != expected:
                raise ValueError(f"Expected {expected} values for {label!r}; got {values}.")
            return values
    raise ValueError(f"Missing row {label!r}.")


def parse_truth(lines: list[str]) -> dict[str, list[float]]:
    return {
        "dlt": values_after(lines, "DLT rate", 5),
        "efficacy": values_after(lines, "Efficacy rate", 5),
        "utility2": values_after(lines, "Utility 2", 5),
        "utility3": values_after(lines, "Utility 3", 5),
    }


def parse_ipde_alpha_source(path: Path) -> tuple[dict[int, dict[str, list[float]]], dict[int, dict[int, dict[float, dict[str, Any]]]]]:
    truth: dict[int, dict[str, list[float]]] = {}
    records: dict[int, dict[int, dict[float, dict[str, Any]]]] = {}
    for scenario, lines in scenario_blocks(pdf_text(path)).items():
        truth[scenario] = parse_truth(lines)
        selected_lines = section(lines, "OBD selection (%)", "Mean patients treated")
        treated_lines = section(lines, "Mean patients treated", "Mean IPDE patients")
        ipde_lines = section(lines, "Mean IPDE patients", None)
        records[scenario] = {2: {}, 3: {}}
        for source_lines, metric, expected_count in (
            (selected_lines, "selected", 6),
            (treated_lines, "treated", 5),
            (ipde_lines, "ipde_patients", 5),
        ):
            for line in source_lines:
                matched = re.fullmatch(r"(U[23]) - IPDE alpha = (0(?:\.\d+)?)\s+(.+)", line)
                if not matched:
                    continue
                utility = int(matched.group(1)[1])
                alpha = float(matched.group(2))
                values = numbers(matched.group(3))
                if len(values) != expected_count:
                    raise ValueError(f"Unexpected {metric} row in {path.name}: {line}")
                record = records[scenario][utility].setdefault(alpha, {})
                record[metric] = values[:5]
                if metric == "selected":
                    record["stop"] = values[5]
        for utility in (2, 3):
            if set(records[scenario][utility]) != set(ALPHAS):
                raise ValueError(f"Incomplete alpha rows for scenario {scenario}, U{utility} in {path.name}.")
            for alpha in ALPHAS:
                if set(records[scenario][utility][alpha]) != {"selected", "treated", "ipde_patients", "stop"}:
                    raise ValueError(f"Incomplete metrics for scenario {scenario}, U{utility}, alpha {alpha}.")
    return truth, records


def parse_cycle1_oracles(path: Path) -> dict[int, dict[str, dict[int, dict[str, Any]]]]:
    oracles: dict[int, dict[str, dict[int, dict[str, Any]]]] = {}
    for scenario, lines in scenario_blocks(pdf_text(path)).items():
        selected_lines = section(lines, "OBD selection (%)", "Mean patients treated")
        treated_lines = section(lines, "Mean patients treated", None)
        oracles[scenario] = {"one_stage": {}, "two_stage": {}}
        for allocation in ("one_stage", "two_stage"):
            for utility in (2, 3):
                label = f"AIDE {allocation} (U{utility})"
                selected = values_after(selected_lines, label, 6)
                treated = values_after(treated_lines, label, 5)
                oracles[scenario][allocation][utility] = {
                    "selected": selected[:5],
                    "stop": selected[5],
                    "treated": treated,
                    "ipde_patients": [0.0] * 5,
                }
    return oracles


def draw_row(pdf: canvas.Canvas, y: float, label: str, values: list[float], stop: float | None, *, shaded: bool = False, decimals: int = 1) -> float:
    if shaded:
        pdf.setFillColor(colors.HexColor("#F6F6F6"))
        pdf.rect(LEFT, y - ROW_HEIGHT + 2, TABLE_WIDTH, ROW_HEIGHT, fill=1, stroke=0)
    pdf.setFillColor(colors.HexColor("#202020"))
    pdf.setFont("Helvetica", ROW_FONT)
    pdf.drawString(LEFT + 8, y - ROW_BASELINE, label)
    x = LEFT + LABEL_WIDTH
    for value in values:
        pdf.drawCentredString(x + DOSE_WIDTH / 2, y - ROW_BASELINE, display(value, decimals))
        x += DOSE_WIDTH
    if stop is not None:
        pdf.drawCentredString(x + STOP_WIDTH / 2, y - ROW_BASELINE, display(stop))
    pdf.setStrokeColor(colors.HexColor("#D0D0D0"))
    pdf.setLineWidth(0.45)
    pdf.line(LEFT, y - ROW_HEIGHT, LEFT + TABLE_WIDTH, y - ROW_HEIGHT)
    return y - ROW_HEIGHT


def draw_group_header(pdf: canvas.Canvas, y: float, label: str) -> float:
    pdf.setFillColor(colors.HexColor("#E7E7E7"))
    pdf.rect(LEFT, y - ROW_HEIGHT + 2, TABLE_WIDTH, ROW_HEIGHT, fill=1, stroke=0)
    pdf.setFillColor(colors.HexColor("#202020"))
    pdf.setFont("Helvetica-Bold", 7.8)
    pdf.drawString(LEFT + 8, y - ROW_BASELINE, label)
    pdf.setStrokeColor(colors.HexColor("#D0D0D0"))
    pdf.setLineWidth(0.45)
    pdf.line(LEFT, y - ROW_HEIGHT, LEFT + TABLE_WIDTH, y - ROW_HEIGHT)
    return y - ROW_HEIGHT


def draw_scenario(pdf: canvas.Canvas, y: float, scenario: int, truth: dict[str, list[float]], alpha_records: dict[int, dict[float, dict[str, Any]]], oracle_records: dict[int, dict[str, Any]]) -> float:
    pdf.setFillColor(colors.HexColor("#EFEFEF"))
    pdf.rect(LEFT, y - 16 + 2, TABLE_WIDTH, 16, fill=1, stroke=0)
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 11)
    pdf.drawString(LEFT + 8, y - 10, f"Scenario {scenario}")
    y -= 18
    y = draw_row(pdf, y, "DLT rate", truth["dlt"], None, decimals=2)
    y = draw_row(pdf, y, "Efficacy rate", truth["efficacy"], None, decimals=2)
    y = draw_row(pdf, y, "Utility 2", truth["utility2"], None)
    y = draw_row(pdf, y, "Utility 3", truth["utility3"], None)
    y -= 1
    for group_label, metric in (
        ("OBD selection (%)", "selected"),
        ("Mean patients treated", "treated"),
        ("Mean IPDE patients", "ipde_patients"),
    ):
        y = draw_group_header(pdf, y, group_label)
        for utility in (2, 3):
            for alpha in ALPHAS:
                record = alpha_records[utility][alpha]
                y = draw_row(
                    pdf,
                    y,
                    f"U{utility} - IPDE alpha = {alpha_label(alpha)}",
                    record[metric],
                    record["stop"] if metric == "selected" else None,
                    shaded=metric == "selected",
                )
            oracle = oracle_records[utility]
            y = draw_row(
                pdf,
                y,
                f"U{utility} - Oracle (cycle = 1)",
                oracle[metric],
                oracle["stop"] if metric == "selected" else None,
                shaded=metric == "selected",
            )
    return y - 3


def draw_header(pdf: canvas.Canvas, allocation: str, page_number: int, page_count: int) -> None:
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 28)
    pdf.drawCentredString(PAGE_WIDTH / 2, PAGE_HEIGHT - 56, "Phase I/II Operating Characteristics")
    pdf.setFont("Helvetica", 12.5)
    pdf.drawCentredString(
        PAGE_WIDTH / 2,
        PAGE_HEIGHT - 78,
        f"IPDE alpha comparison + cycle = 1 oracle - allocation = {allocation} - futility cutoff = 0.95 - N = 30",
    )
    pdf.setStrokeColor(colors.HexColor("#999999"))
    pdf.setLineWidth(0.8)
    pdf.line(LEFT, PAGE_HEIGHT - 92, PAGE_WIDTH - LEFT, PAGE_HEIGHT - 92)
    x = LEFT + LABEL_WIDTH
    pdf.setFont("Helvetica-Bold", 13)
    for dose in range(1, 6):
        pdf.drawCentredString(x + DOSE_WIDTH / 2, PAGE_HEIGHT - 110, f"Dose {dose}")
        x += DOSE_WIDTH
    pdf.drawCentredString(x + STOP_WIDTH / 2, PAGE_HEIGHT - 110, "Stop %")
    pdf.setFont("Helvetica", 7.0)
    pdf.drawString(LEFT, 38, "IPDE-alpha rows use the August 3 fixed-prior results. U2: Lambda_T = 0.3; U3: Lambda_T = 1.0; both IPDE alphas equal the row value.")
    pdf.drawString(LEFT, 25, "Oracle rows are the Cycle_Max = 1 AIDE results from the August 3 all-method table; mean IPDE patients = 0 by design.")
    pdf.drawRightString(PAGE_WIDTH - LEFT, 25, f"Page {page_number} of {page_count}")


def build(allocation: str) -> Path:
    source = SOURCE_DIR / f"phase12_ipde_alpha_comparison_{allocation}_scenarios_1_to_37_tables.pdf"
    truth, alpha_records = parse_ipde_alpha_source(source)
    oracles = parse_cycle1_oracles(ALL_METHODS_SOURCE)
    if set(truth) != SCENARIOS or set(oracles) != SCENARIOS:
        raise ValueError("The source PDFs do not cover scenarios 1 through 37.")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUT_DIR / f"phase12_ipde_alpha_comparison_{allocation}_fut0p95_scenarios_1_to_37_tables.pdf"
    pdf = canvas.Canvas(str(output), pagesize=(PAGE_WIDTH, PAGE_HEIGHT), pageCompression=1)
    page_count = 19
    page = 0
    for start in range(1, 38, 2):
        page += 1
        draw_header(pdf, allocation, page, page_count)
        y = PAGE_HEIGHT - 130
        y = draw_scenario(pdf, y, start, truth[start], alpha_records[start], oracles[start][allocation])
        if start + 1 <= 37:
            draw_scenario(pdf, y, start + 1, truth[start + 1], alpha_records[start + 1], oracles[start + 1][allocation])
        pdf.showPage()
    pdf.save()
    return output


def main() -> None:
    for allocation in ("one_stage", "two_stage"):
        print(build(allocation))


if __name__ == "__main__":
    main()
