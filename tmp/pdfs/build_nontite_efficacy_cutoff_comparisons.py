from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import pdfplumber
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
SOURCE_FILES = {
    0.75: ROOT / "Presentation 8-10-2026" / "AIDE_phase_I_II_N30_ncycle2_rp0p15x0p85_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p75_ap0p15x0p85_0p3x0p7_0p5x0p5_1x1_IDX_0001_to_1000_dose_summary.csv",
    0.85: ROOT / "Presentation 8-10-2026" / "AIDE_phase_I_II_N30_ncycle2_rp0p15x0p85_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_ap0p15x0p85_0p3x0p7_0p5x0p5_1x1_IDX_0001_to_1000_dose_summary.csv",
    0.95: ROOT / "Presentation 8-03-2026" / "AIDE_phase_I_II_IDX_0001_to_1000_dose_summary.csv",
}
UTILITY_PDF = ROOT / "Presentation 8-03-2026" / "Table and Plots" / "Phase I-II" / "phase12_ipde_alpha_comparison_two_stage_scenarios_1_to_37_tables.pdf"
OUT_DIR = ROOT / "Presentation 8-10-2026" / "Table and Plots" / "Phase I-II"

PAGE_WIDTH = 18 * inch
PAGE_HEIGHT = 14 * inch
LEFT, RIGHT = 48, PAGE_WIDTH - 44
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
CUTOFFS = [0.75, 0.85, 0.95]
ROWS = [(utility, cutoff) for utility in (2, 3) for cutoff in CUTOFFS]


def close(a: float, b: float, tol: float = 1e-10) -> bool:
    return math.isclose(float(a), float(b), abs_tol=tol, rel_tol=tol)


def get_field(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value in ("", "NA", "NaN", None):
        return float("nan")
    return float(value)


def fmt(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "-"
    text = f"{value:.1f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def fmt_rate(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "-"
    return f"{value:.2f}"


def alpha_tag(alpha: float) -> str:
    return "alpha" + ("0" if alpha == 0 else str(alpha).replace(".", "p"))


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def load_utilities(path: Path) -> dict[int, dict[int, list[float]]]:
    utilities: dict[int, dict[int, list[float]]] = {}
    with pdfplumber.open(path) as pdf:
        for page in pdf.pages:
            scenario = None
            for raw in (page.extract_text() or "").splitlines():
                line = raw.strip()
                match = re.fullmatch(r"Scenario\s+(\d+)", line)
                if match:
                    scenario = int(match.group(1))
                    utilities.setdefault(scenario, {})
                    continue
                if scenario is None:
                    continue
                for utility_type in (2, 3):
                    prefix = f"Utility {utility_type} "
                    if line.startswith(prefix):
                        values = [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", line[len(prefix):])]
                        if len(values) != 5:
                            raise ValueError(f"Could not parse Utility {utility_type}: {line}")
                        utilities[scenario][utility_type] = values
    if set(utilities) != set(range(1, 38)) or any(set(values) != {2, 3} for values in utilities.values()):
        raise ValueError("The August 3 utility reference is incomplete.")
    return utilities


def selected_rows(path: Path, cutoff: float) -> dict[tuple, list[dict[str, str]]]:
    result: dict[tuple, list[dict[str, str]]] = defaultdict(list)
    priors = {
        "CRM_Prior_a": 0.15,
        "CRM_Prior_b": 0.85,
        "Efficacy_Prior_a": 0.5,
        "Efficacy_Prior_b": 0.5,
        "Efficacy_Carryover_Prior_a": 0.15,
        "Efficacy_Carryover_Prior_b": 0.85,
        "Efficacy_Additive_Alpha_Prior_a": 0.15,
        "Efficacy_Additive_Alpha_Prior_b": 0.85,
    }
    for row in load_rows(path):
        if row["Allocation"] not in {"one_stage", "two_stage"}:
            continue
        if row["CRM_r_Model"] != "previous_dose" or row["Efficacy_Model"] != "previous_dose_additive":
            continue
        if int(float(row["IPDE_Design"])) != 1 or not all(close(get_field(row, field), value) for field, value in priors.items()):
            continue
        # The two updated August 10 sources expose their efficacy settings; the
        # requested August 3 comparator is the legacy 0.95-cutoff source.
        if "Futility_Cutoff" in row:
            if not close(get_field(row, "Efficacy_Threshold"), 0.20) or not close(get_field(row, "Futility_Cutoff"), cutoff):
                continue
        elif cutoff != 0.95:
            continue
        alpha_t = get_field(row, "Toxicity_IPDE_Alpha")
        alpha_e = get_field(row, "Efficacy_IPDE_Alpha")
        if not close(alpha_t, alpha_e) or not any(close(alpha_t, alpha) for alpha in ALPHAS):
            continue
        key = (int(row["Scenario"]), row["Allocation"], int(float(row["Utility_Type"])), alpha_t)
        result[key].append(row)
    expected = {
        (scenario, allocation, utility, alpha)
        for scenario in range(1, 38)
        for allocation, _ in ALLOCATIONS
        for utility in (2, 3)
        for alpha in ALPHAS
    }
    if set(result) != expected:
        raise ValueError(f"Cutoff {cutoff} source has an incomplete comparison grid.")
    for key, rows in result.items():
        if len(rows) != 5:
            raise ValueError(f"Cutoff {cutoff}, {key} has {len(rows)} dose rows; expected 5.")
        rows.sort(key=lambda item: int(item["Dose"]))
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


def rows_for(data: dict[float, dict[tuple, list[dict[str, str]]]], cutoff: float, scenario: int, allocation: str, utility: int, alpha: float) -> list[dict[str, str]]:
    return data[cutoff][(scenario, allocation, utility, alpha)]


def draw_metric_row(c: canvas.Canvas, label: str, values: list[float], stop: float | None, y: float, shaded: bool) -> float:
    if shaded:
        c.setFillColor(VERY_LIGHT_BLUE)
        c.rect(LEFT, y - 3, RIGHT - LEFT, ROW_H, stroke=0, fill=1)
    draw_text(c, label, LABEL_X, y, size=7.25)
    for x, value in zip(VALUE_X, values):
        draw_right(c, fmt(value), x + 42, y, size=7.25)
    draw_right(c, fmt(stop), STOP_X + 48, y, size=7.25)
    return y - ROW_H


def draw_duration_row(c: canvas.Canvas, label: str, duration: float, y: float, shaded: bool) -> float:
    if shaded:
        c.setFillColor(VERY_LIGHT_BLUE)
        c.rect(LEFT, y - 3, RIGHT - LEFT, ROW_H, stroke=0, fill=1)
    draw_text(c, label, LABEL_X, y, size=7.25)
    draw_right(c, fmt(duration), STOP_X + 48, y, size=7.25)
    return y - ROW_H


def draw_scenario(c: canvas.Canvas, y: float, scenario: int, allocation: str, alpha: float, data: dict[float, dict[tuple, list[dict[str, str]]]], utilities: dict[int, dict[int, list[float]]]) -> float:
    base = rows_for(data, 0.75, scenario, allocation, 2, alpha)
    c.setFillColor(LIGHT_BLUE)
    c.rect(LEFT, y - 4, RIGHT - LEFT, 16, stroke=0, fill=1)
    draw_text(c, f"Scenario {scenario}", LABEL_X, y, size=9.5, color=NAVY, font="Helvetica-Bold")
    y -= 20
    truth = [
        ("DLT rate", [get_field(row, "True_DLT_rate") for row in base]),
        ("Efficacy rate", [get_field(row, "True_Efficacy_rate") for row in base]),
        ("Utility 2", utilities[scenario][2]),
        ("Utility 3", utilities[scenario][3]),
    ]
    for label, values in truth:
        draw_text(c, label, LABEL_X, y, size=7.7, color=GRAY, font="Helvetica-Bold" if label.startswith("Utility") else "Helvetica")
        for x, value in zip(VALUE_X, values):
            value_text = fmt_rate(value) if label in {"DLT rate", "Efficacy rate"} else fmt(value)
            draw_right(c, value_text, x + 42, y, size=7.7, color=GRAY)
        y -= ROW_H
    y -= 2
    for heading, field, has_stop in [
        ("OBD selection (%)", "OBD_Selection_pct", True),
        ("Mean administrations", "Pts_Treated", False),
        ("Mean IPDE administrations", "IPDE_Doses", False),
    ]:
        draw_text(c, heading, LABEL_X, y, size=7.8, color=BLUE, font="Helvetica-Bold")
        y -= ROW_H
        for utility, cutoff in ROWS:
            rows = rows_for(data, cutoff, scenario, allocation, utility, alpha)
            stop = get_field(rows[0], "No_OBD_Selection_pct") if has_stop else None
            y = draw_metric_row(c, f"U{utility} - cutoff = {cutoff:.2f}", [get_field(row, field) for row in rows], stop, y, shaded=(cutoff != 0.95))
        y -= 2
    draw_text(c, "Mean trial duration (days)", LABEL_X, y, size=7.8, color=BLUE, font="Helvetica-Bold")
    draw_right(c, "Days", STOP_X + 48, y, size=7.5, color=BLUE, font="Helvetica-Bold")
    y -= ROW_H
    for utility, cutoff in ROWS:
        rows = rows_for(data, cutoff, scenario, allocation, utility, alpha)
        y = draw_duration_row(c, f"U{utility} - cutoff = {cutoff:.2f}", get_field(rows[0], "Duration"), y, shaded=(cutoff != 0.95))
    draw_rule(c, y + 4)
    return y - 10


def draw_chrome(c: canvas.Canvas, alpha: float, allocation: str, page: int, total: int) -> None:
    c.setFillColor(NAVY)
    c.rect(0, PAGE_HEIGHT - 50, PAGE_WIDTH, 50, stroke=0, fill=1)
    draw_text(c, "Phase I/II Operating Characteristics", LEFT, PAGE_HEIGHT - 31, size=16, color=white, font="Helvetica-Bold")
    draw_text(c, f"Non-TITE efficacy cutoff comparison - IPDE alpha = {fmt(alpha)} - allocation = {allocation} - N = 30", LEFT, PAGE_HEIGHT - 74, size=10.4, color=NAVY, font="Helvetica-Bold")
    for label, x in zip(("Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5"), VALUE_X):
        draw_text(c, label, x, PAGE_HEIGHT - 94, size=8.5, color=BLUE, font="Helvetica-Bold")
    draw_text(c, "Stop %", STOP_X, PAGE_HEIGHT - 94, size=8.5, color=BLUE, font="Helvetica-Bold")
    draw_rule(c, PAGE_HEIGHT - 100, color=BLUE, width=0.8)
    draw_text(c, "Priors: toxicity alpha Beta(0.15, 0.85); efficacy beta-binomial Beta(0.5, 0.5); efficacy alpha Beta(0.15, 0.85).", LEFT, 35, size=7.1, color=GRAY)
    draw_text(c, "Cutoffs 0.75 and 0.85 use August 10 data with efficacy threshold = 0.20; cutoff 0.95 uses the specified August 3 source.", LEFT, 22, size=7.1, color=GRAY)
    draw_right(c, f"Page {page} of {total}", RIGHT, 22, size=7.1, color=GRAY)


def make_pdf(alpha: float, data: dict[float, dict[tuple, list[dict[str, str]]]], utilities: dict[int, dict[int, list[float]]]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"phase12_nontite_efficacy_cutoff_comparison_{alpha_tag(alpha)}_scenarios_1_to_37_tables.pdf"
    pdf = canvas.Canvas(str(path), pagesize=(PAGE_WIDTH, PAGE_HEIGHT), pageCompression=1)
    total_pages, page = 38, 0
    for allocation, allocation_label in ALLOCATIONS:
        for scenario in range(1, 38, 2):
            page += 1
            draw_chrome(pdf, alpha, allocation_label, page, total_pages)
            draw_scenario(pdf, PAGE_HEIGHT - 126, scenario, allocation, alpha, data, utilities)
            if scenario + 1 <= 37:
                draw_scenario(pdf, PAGE_HEIGHT - 536, scenario + 1, allocation, alpha, data, utilities)
            pdf.showPage()
    pdf.save()
    return path


def main() -> None:
    utilities = load_utilities(UTILITY_PDF)
    data = {cutoff: selected_rows(path, cutoff) for cutoff, path in SOURCE_FILES.items()}
    for alpha in ALPHAS:
        print(make_pdf(alpha, data, utilities))


if __name__ == "__main__":
    main()
