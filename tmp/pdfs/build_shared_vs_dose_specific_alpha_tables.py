from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "Presentation 8-10-2026" / "Raw Data" / (
    "AIDE_phase_I_II_N30_ncycle2_rp0p15x0p85_ep0p5x0p5_cp0p15x0p85_"
    "eth0p2_fut0p85_ap0p15x0p85_0p3x0p7_0p5x0p5_1x1_IDX_0001_to_1000_dose_summary.csv"
)
TRUTH_SOURCE = ROOT / "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv"
OUT_DIR = ROOT / "Presentation 8-10-2026" / "Table and Plots" / "Phase I-II" / "Model"

PAGE_WIDTH = 18 * inch
PAGE_HEIGHT = 14 * inch
LEFT = 48
RIGHT = PAGE_WIDTH - 44
LABEL_X = 66
VALUE_X = [545, 665, 785, 905, 1025]
STOP_X = 1145
ROW_H = 11
NAVY = HexColor("#16365C")
BLUE = HexColor("#1F4E79")
LIGHT_BLUE = HexColor("#D9EAF7")
VERY_LIGHT_BLUE = HexColor("#F1F6FB")
GRID = HexColor("#AAB7C4")
GRAY = HexColor("#5B6570")
ALPHAS = [0.0, 0.3, 0.6, 0.9]
ALLOCATIONS = [("two_stage", "two-stage"), ("one_stage", "one-stage")]
METHODS = [
    (2, "previous_dose_additive", "U2 - shared alpha"),
    (2, "dose_specific_previous_dose_additive", "U2 - dose-specific alpha"),
    (3, "previous_dose_additive", "U3 - shared alpha"),
    (3, "dose_specific_previous_dose_additive", "U3 - dose-specific alpha"),
]


def number(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float("nan") if value in ("", "NA", "NaN", None) else float(value)


def close(a: float, b: float, tol: float = 1e-10) -> bool:
    return math.isclose(a, b, rel_tol=tol, abs_tol=tol)


def fmt(value: float) -> str:
    if not math.isfinite(value):
        return "-"
    text = f"{value:.1f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def fmt_rate(value: float) -> str:
    return "-" if not math.isfinite(value) else f"{value:.2f}"


def alpha_tag(alpha: float) -> str:
    return "alpha" + ("0" if alpha == 0 else str(alpha).replace(".", "p"))


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def load_truth() -> dict[int, dict[int, list[float]]]:
    truth: dict[int, dict[int, list[float]]] = {}
    for row in load_rows(TRUTH_SOURCE):
        scenario = int(number(row, "Scenario"))
        truth[scenario] = {
            2: [number(row, f"Utility2_Dose{dose}") for dose in range(1, 6)],
            3: [number(row, f"Utility3_Dose{dose}") for dose in range(1, 6)],
        }
    return truth


def selected_rows() -> dict[tuple, list[dict[str, str]]]:
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
    valid_models = {"previous_dose_additive", "dose_specific_previous_dose_additive"}
    for row in load_rows(SOURCE):
        if row["Allocation"] not in {item[0] for item in ALLOCATIONS}:
            continue
        if row["CRM_r_Model"] != "previous_dose" or row["Efficacy_Model"] not in valid_models:
            continue
        if int(number(row, "IPDE_Design")) != 1:
            continue
        if not all(close(number(row, field), value) for field, value in numeric_requirements.items()):
            continue
        alpha_t = number(row, "Toxicity_IPDE_Alpha")
        alpha_e = number(row, "Efficacy_IPDE_Alpha")
        if not close(alpha_t, alpha_e) or not any(close(alpha_t, alpha) for alpha in ALPHAS):
            continue
        key = (
            int(number(row, "Scenario")), row["Allocation"], int(number(row, "Utility_Type")),
            alpha_t, row["Efficacy_Model"],
        )
        result[key].append(row)
    expected = {
        (scenario, allocation, utility, alpha, model)
        for scenario in range(1, 38)
        for allocation, _ in ALLOCATIONS
        for utility in (2, 3)
        for alpha in ALPHAS
        for model in ("previous_dose_additive", "dose_specific_previous_dose_additive")
    }
    if set(result) != expected:
        raise ValueError("The selected summary rows do not form the required complete comparison grid.")
    for rows in result.values():
        if len(rows) != 5:
            raise ValueError("Each scenario-model comparison row must contain five dose levels.")
        rows.sort(key=lambda item: int(number(item, "Dose")))
    return result


def draw_text(c: canvas.Canvas, text: str, x: float, y: float, size: float = 8, color=black, font: str = "Helvetica") -> None:
    c.setFillColor(color)
    c.setFont(font, size)
    c.drawString(x, y, text)


def draw_right(c: canvas.Canvas, text: str, x: float, y: float, size: float = 8, color=black, font: str = "Helvetica") -> None:
    c.setFillColor(color)
    c.setFont(font, size)
    c.drawRightString(x, y, text)


def draw_rule(c: canvas.Canvas, y: float, color=GRID, width: float = 0.45) -> None:
    c.setStrokeColor(color)
    c.setLineWidth(width)
    c.line(LEFT, y, RIGHT, y)


def rows_for(data: dict[tuple, list[dict[str, str]]], scenario: int, allocation: str, utility: int, alpha: float, model: str) -> list[dict[str, str]]:
    return data[(scenario, allocation, utility, alpha, model)]


def draw_metric_row(c: canvas.Canvas, label: str, values: list[float], stop: float | None, y: float, shaded: bool) -> float:
    if shaded:
        c.setFillColor(VERY_LIGHT_BLUE)
        c.rect(LEFT, y - 3, RIGHT - LEFT, ROW_H, stroke=0, fill=1)
    draw_text(c, label, LABEL_X, y, size=7.15)
    for x, value in zip(VALUE_X, values):
        draw_right(c, fmt(value), x + 42, y, size=7.15)
    draw_right(c, fmt(stop if stop is not None else float("nan")), STOP_X + 48, y, size=7.15)
    return y - ROW_H


def draw_duration_row(c: canvas.Canvas, label: str, duration: float, y: float, shaded: bool) -> float:
    if shaded:
        c.setFillColor(VERY_LIGHT_BLUE)
        c.rect(LEFT, y - 3, RIGHT - LEFT, ROW_H, stroke=0, fill=1)
    draw_text(c, label, LABEL_X, y, size=7.15)
    draw_right(c, fmt(duration), STOP_X + 48, y, size=7.15)
    return y - ROW_H


def draw_scenario(c: canvas.Canvas, y: float, scenario: int, allocation: str, alpha: float, data: dict[tuple, list[dict[str, str]]], truth: dict[int, dict[int, list[float]]]) -> float:
    base = rows_for(data, scenario, allocation, 2, alpha, "previous_dose_additive")
    c.setFillColor(LIGHT_BLUE)
    c.rect(LEFT, y - 4, RIGHT - LEFT, 16, stroke=0, fill=1)
    draw_text(c, f"Scenario {scenario}", LABEL_X, y, size=9.25, color=NAVY, font="Helvetica-Bold")
    y -= 19
    truth_rows = [
        ("DLT rate", [number(row, "True_DLT_rate") for row in base], False),
        ("Efficacy rate", [number(row, "True_Efficacy_rate") for row in base], False),
        ("Utility 2", truth[scenario][2], True),
        ("Utility 3", truth[scenario][3], True),
    ]
    for label, values, utility in truth_rows:
        draw_text(c, label, LABEL_X, y, size=7.5, color=GRAY, font="Helvetica-Bold" if utility else "Helvetica")
        for x, value in zip(VALUE_X, values):
            draw_right(c, fmt(value) if utility else fmt_rate(value), x + 42, y, size=7.5, color=GRAY)
        y -= ROW_H
    y -= 2
    sections = [
        ("MTD selection (%)", "MTD_Selection_pct", "Early_Stopping_pct"),
        ("OBD selection (%)", "OBD_Selection_pct", "No_OBD_Selection_pct"),
        ("Mean administrations", "Pts_Treated", None),
        ("Mean IPDE administrations", "IPDE_Doses", None),
    ]
    for heading, field, stop_field in sections:
        draw_text(c, heading, LABEL_X, y, size=7.65, color=BLUE, font="Helvetica-Bold")
        y -= ROW_H
        for utility, model, label in METHODS:
            rows = rows_for(data, scenario, allocation, utility, alpha, model)
            stop = number(rows[0], stop_field) if stop_field else None
            y = draw_metric_row(c, label, [number(row, field) for row in rows], stop, y, shaded=(model != "previous_dose_additive"))
        y -= 1
    draw_text(c, "Mean trial duration (days)", LABEL_X, y, size=7.65, color=BLUE, font="Helvetica-Bold")
    draw_right(c, "Days", STOP_X + 48, y, size=7.4, color=BLUE, font="Helvetica-Bold")
    y -= ROW_H
    for utility, model, label in METHODS:
        rows = rows_for(data, scenario, allocation, utility, alpha, model)
        y = draw_duration_row(c, label, number(rows[0], "Duration"), y, shaded=(model != "previous_dose_additive"))
    draw_rule(c, y + 4)
    return y - 9


def draw_chrome(c: canvas.Canvas, alpha: float, allocation_label: str, page: int, total: int) -> None:
    c.setFillColor(NAVY)
    c.rect(0, PAGE_HEIGHT - 50, PAGE_WIDTH, 50, stroke=0, fill=1)
    draw_text(c, "Phase I/II Operating Characteristics", LEFT, PAGE_HEIGHT - 31, size=16, color=white, font="Helvetica-Bold")
    subtitle = f"Shared vs dose-specific additive-alpha efficacy - IPDE alpha = {fmt(alpha)} - allocation = {allocation_label} - N = 30"
    draw_text(c, subtitle, LEFT, PAGE_HEIGHT - 74, size=10.2, color=NAVY, font="Helvetica-Bold")
    for label, x in zip(("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5"), VALUE_X):
        draw_text(c, label, x, PAGE_HEIGHT - 94, size=8.5, color=BLUE, font="Helvetica-Bold")
    draw_text(c, "No MTD/OBD %", STOP_X, PAGE_HEIGHT - 94, size=8.5, color=BLUE, font="Helvetica-Bold")
    draw_rule(c, PAGE_HEIGHT - 100, color=BLUE, width=0.8)
    footer = "Priors: efficacy beta-binomial Beta(0.5, 0.5); efficacy additive alpha Beta(0.15, 0.85); toxicity additive alpha Beta(0.15, 0.85)."
    draw_text(c, footer, LEFT, 35, size=7.1, color=GRAY)
    draw_text(c, "Final column: early stopping % for MTD selection and no-OBD % for OBD selection; shaded rows are the dose-specific alpha efficacy model.", LEFT, 22, size=7.1, color=GRAY)
    draw_right(c, f"Page {page} of {total}", RIGHT, 22, size=7.1, color=GRAY)


def make_pdf(alpha: float, data: dict[tuple, list[dict[str, str]]], truth: dict[int, dict[int, list[float]]]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUT_DIR / f"phase12_shared_vs_dose_specific_additive_alpha_fut0p85_{alpha_tag(alpha)}_scenarios_1_to_37_tables.pdf"
    pdf = canvas.Canvas(str(output), pagesize=(PAGE_WIDTH, PAGE_HEIGHT), pageCompression=1)
    page = 0
    total_pages = 38
    for allocation, allocation_label in ALLOCATIONS:
        for scenario in range(1, 38, 2):
            page += 1
            draw_chrome(pdf, alpha, allocation_label, page, total_pages)
            draw_scenario(pdf, PAGE_HEIGHT - 126, scenario, allocation, alpha, data, truth)
            if scenario + 1 <= 37:
                draw_scenario(pdf, PAGE_HEIGHT - 536, scenario + 1, allocation, alpha, data, truth)
            pdf.showPage()
    pdf.save()
    return output


if __name__ == "__main__":
    comparison_data = selected_rows()
    utility_truth = load_truth()
    for alpha_value in ALPHAS:
        print(make_pdf(alpha_value, comparison_data, utility_truth))
