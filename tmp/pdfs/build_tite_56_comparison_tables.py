from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import pdfplumber
from reportlab.lib.colors import Color, HexColor, black, white
from reportlab.lib.units import inch
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
TITE_CSV = ROOT / "Presentation 8-10-2026" / "TITE_AIDE_phase_I_II_IDX_1001_to_2000_dose_summary.csv"
TASK_MAP_CSV = ROOT / "tmp" / "pdfs" / "tite_56_task_map.csv"
NON_TITE_CSV = ROOT / "Presentation 8-10-2026" / "AIDE_phase_I_II_N30_ncycle2_rp0p15x0p85_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_ap0p15x0p85_0p3x0p7_0p5x0p5_1x1_IDX_0001_to_1000_dose_summary.csv"
UTILITY_PDF = ROOT / "Presentation 8-03-2026" / "Table and Plots" / "Phase I-II" / "phase12_ipde_alpha_comparison_two_stage_scenarios_1_to_37_tables.pdf"
OUT_DIR = ROOT / "Presentation 8-10-2026" / "Table and Plots" / "Phase I-II"

PAGE_WIDTH = 18 * inch
PAGE_HEIGHT = 14 * inch
LEFT = 48
LABEL_X = 66
VALUE_X = [545, 665, 785, 905, 1025]
STOP_X = 1145
RIGHT = PAGE_WIDTH - 44
NAVY = HexColor("#16365C")
BLUE = HexColor("#1F4E79")
LIGHT_BLUE = HexColor("#D9EAF7")
VERY_LIGHT_BLUE = HexColor("#F1F6FB")
GRID = HexColor("#AAB7C4")
GRAY = HexColor("#5B6570")

ALPHAS = [0.0, 0.3, 0.6, 0.9]
ALLOCATIONS = [("two_stage", "two-stage"), ("one_stage", "one-stage")]
METHODS = [
    (2, "Original non-TITE"),
    (2, "TITE, 56-day accrual"),
    (3, "Original non-TITE"),
    (3, "TITE, 56-day accrual"),
]


def num(value: str) -> float:
    return float(value)


def close(a: float, b: float, tol: float = 1e-10) -> bool:
    return math.isclose(float(a), float(b), abs_tol=tol, rel_tol=tol)


def fmt(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "-"
    text = f"{value:.1f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def alpha_tag(alpha: float) -> str:
    return "alpha" + ("0" if alpha == 0 else str(alpha).replace(".", "p"))


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def utilities_from_reference(path: Path) -> dict[int, dict[int, list[float]]]:
    utilities: dict[int, dict[int, list[float]]] = {}
    with pdfplumber.open(path) as pdf:
        for page in pdf.pages:
            scenario = None
            text = page.extract_text() or ""
            for raw in text.splitlines():
                line = raw.strip()
                scenario_match = re.fullmatch(r"Scenario\s+(\d+)", line)
                if scenario_match:
                    scenario = int(scenario_match.group(1))
                    utilities.setdefault(scenario, {})
                    continue
                if scenario is None:
                    continue
                for utility_type in (2, 3):
                    prefix = f"Utility {utility_type} "
                    if line.startswith(prefix):
                        values = [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", line[len(prefix):])]
                        if len(values) != 5:
                            raise ValueError(f"Could not parse five Utility {utility_type} values from: {line}")
                        utilities[scenario][utility_type] = values
    if set(utilities) != set(range(1, 38)):
        raise ValueError("Utility reference table does not cover scenarios 1-37.")
    for scenario, values in utilities.items():
        if set(values) != {2, 3}:
            raise ValueError(f"Utility reference is incomplete for scenario {scenario}.")
    return utilities


def get_field(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value in ("", "NA", "NaN", None):
        return float("nan")
    return float(value)


def build_datasets() -> tuple[dict[tuple, list[dict[str, str]]], dict[tuple, list[dict[str, str]]]]:
    task_map = load_rows(TASK_MAP_CSV)
    task_by_id = {int(row["Task_ID"]): row for row in task_map}
    if len(task_by_id) != 1184:
        raise ValueError("The recreated 56-day task map must contain 1,184 tasks.")

    tite_by_key: dict[tuple, list[dict[str, str]]] = defaultdict(list)
    for row in load_rows(TITE_CSV):
        task = task_by_id.get(int(row["Task_ID"]))
        if task is None:
            raise ValueError(f"No recreated task metadata for TITE Task_ID {row['Task_ID']}.")
        if not close(get_field(row, "Arrival_Rate"), 1 / 56):
            raise ValueError("The available TITE summary is not a 56-day accrual run.")
        if task["Model_ID"] != "additive_shared":
            continue
        key = (
            int(task["Scenario"]), task["Allocation"], int(task["Utility_Type"]),
            float(task["Toxicity_IPDE_Alpha"]),
        )
        tite_by_key[key].append(row)

    non_tite_by_key: dict[tuple, list[dict[str, str]]] = defaultdict(list)
    for row in load_rows(NON_TITE_CSV):
        if row["Allocation"] not in {"one_stage", "two_stage"}:
            continue
        if row["CRM_r_Model"] != "previous_dose" or row["Efficacy_Model"] != "previous_dose_additive":
            continue
        if int(float(row["IPDE_Design"])) != 1:
            continue
        expected = {
            "CRM_Prior_a": 0.15, "CRM_Prior_b": 0.85,
            "Efficacy_Prior_a": 0.5, "Efficacy_Prior_b": 0.5,
            "Efficacy_Carryover_Prior_a": 0.15, "Efficacy_Carryover_Prior_b": 0.85,
            "Efficacy_Additive_Alpha_Prior_a": 0.15, "Efficacy_Additive_Alpha_Prior_b": 0.85,
            "Efficacy_Threshold": 0.20, "Futility_Cutoff": 0.85, "Min_Eff_N_for_Futility": 0,
        }
        if not all(close(get_field(row, field), value) for field, value in expected.items()):
            continue
        alpha_t = get_field(row, "Toxicity_IPDE_Alpha")
        alpha_e = get_field(row, "Efficacy_IPDE_Alpha")
        if not close(alpha_t, alpha_e) or not any(close(alpha_t, alpha) for alpha in ALPHAS):
            continue
        key = (int(row["Scenario"]), row["Allocation"], int(float(row["Utility_Type"])), alpha_t)
        non_tite_by_key[key].append(row)

    expected_keys = {
        (scenario, allocation, utility, alpha)
        for scenario in range(1, 38)
        for allocation, _ in ALLOCATIONS
        for utility in (2, 3)
        for alpha in ALPHAS
    }
    for name, dataset in (("TITE", tite_by_key), ("non-TITE", non_tite_by_key)):
        if set(dataset) != expected_keys:
            missing = sorted(expected_keys - set(dataset))[:3]
            extra = sorted(set(dataset) - expected_keys)[:3]
            raise ValueError(f"{name} keys do not match requested comparison grid. Missing={missing}; extra={extra}")
        for key, rows in dataset.items():
            if len(rows) != 5:
                raise ValueError(f"{name} {key} has {len(rows)} dose rows; expected 5.")
            rows.sort(key=lambda item: int(item["Dose"]))
    return non_tite_by_key, tite_by_key


def draw_text(c: canvas.Canvas, text: str, x: float, y: float, size: float = 10, color=black, font: str = "Helvetica") -> None:
    c.setFillColor(color)
    c.setFont(font, size)
    c.drawString(x, y, text)


def draw_right(c: canvas.Canvas, text: str, x: float, y: float, size: float = 10, color=black, font: str = "Helvetica") -> None:
    c.setFillColor(color)
    c.setFont(font, size)
    c.drawRightString(x, y, text)


def draw_rule(c: canvas.Canvas, y: float, color=GRID, width: float = 0.45) -> None:
    c.setStrokeColor(color)
    c.setLineWidth(width)
    c.line(LEFT, y, RIGHT, y)


def draw_value_row(c: canvas.Canvas, label: str, values: list[float], stop: float | None, y: float, shaded: bool = False) -> float:
    if shaded:
        c.setFillColor(VERY_LIGHT_BLUE)
        c.rect(LEFT, y - 4, RIGHT - LEFT, 15, stroke=0, fill=1)
    draw_text(c, label, LABEL_X, y, size=8.6)
    for x, value in zip(VALUE_X, values):
        draw_right(c, fmt(value), x + 42, y, size=8.6)
    draw_right(c, fmt(stop), STOP_X + 48, y, size=8.6)
    return y - 15


def draw_duration_row(c: canvas.Canvas, label: str, duration: float, y: float, shaded: bool = False) -> float:
    if shaded:
        c.setFillColor(VERY_LIGHT_BLUE)
        c.rect(LEFT, y - 4, RIGHT - LEFT, 15, stroke=0, fill=1)
    draw_text(c, label, LABEL_X, y, size=8.6)
    draw_right(c, fmt(duration), STOP_X + 48, y, size=8.6)
    return y - 15


def scenario_rows(dataset: dict[tuple, list[dict[str, str]]], scenario: int, allocation: str, utility: int, alpha: float) -> list[dict[str, str]]:
    key = (scenario, allocation, utility, alpha)
    if key not in dataset:
        raise ValueError(f"Missing data for {key}")
    return dataset[key]


def tite_stop(rows: list[dict[str, str]]) -> float:
    return max(0.0, 100.0 - sum(get_field(row, "OBD_Selection_pct") for row in rows))


def draw_scenario(
    c: canvas.Canvas,
    y: float,
    scenario: int,
    allocation: str,
    alpha: float,
    non_tite: dict[tuple, list[dict[str, str]]],
    tite: dict[tuple, list[dict[str, str]]],
    utilities: dict[int, dict[int, list[float]]],
) -> float:
    base = scenario_rows(non_tite, scenario, allocation, 2, alpha)
    dlt = [get_field(row, "True_DLT_rate") for row in base]
    efficacy = [get_field(row, "True_Efficacy_rate") for row in base]

    c.setFillColor(LIGHT_BLUE)
    c.rect(LEFT, y - 5, RIGHT - LEFT, 18, stroke=0, fill=1)
    draw_text(c, f"Scenario {scenario}", LABEL_X, y, size=10.5, color=NAVY, font="Helvetica-Bold")
    y -= 22

    for label, values in (("DLT rate", dlt), ("Efficacy rate", efficacy), ("Utility 2", utilities[scenario][2]), ("Utility 3", utilities[scenario][3])):
        draw_text(c, label, LABEL_X, y, size=8.6, color=GRAY, font="Helvetica-Bold" if label.startswith("Utility") else "Helvetica")
        for x, value in zip(VALUE_X, values):
            draw_right(c, fmt(value), x + 42, y, size=8.6, color=GRAY)
        y -= 15
    y -= 3

    sections = [
        ("OBD selection (%)", "OBD_Selection_pct", True),
        ("Mean administrations", "Pts_Treated", False),
        ("Mean IPDE administrations", "IPDE_Doses", False),
    ]
    for section, field, show_stop in sections:
        draw_text(c, section, LABEL_X, y, size=8.8, color=BLUE, font="Helvetica-Bold")
        y -= 15
        for utility, method in METHODS:
            rows = scenario_rows(non_tite if method == "Original non-TITE" else tite, scenario, allocation, utility, alpha)
            values = [get_field(row, field) for row in rows]
            stop = None
            if show_stop:
                stop = get_field(rows[0], "No_OBD_Selection_pct") if method == "Original non-TITE" else tite_stop(rows)
            y = draw_value_row(c, f"U{utility} - {method}", values, stop, y, shaded=(method != "Original non-TITE"))
        y -= 3

    draw_text(c, "Mean trial duration (days)", LABEL_X, y, size=8.8, color=BLUE, font="Helvetica-Bold")
    draw_right(c, "Days", STOP_X + 48, y, size=8.2, color=BLUE, font="Helvetica-Bold")
    y -= 15
    for utility, method in METHODS:
        rows = scenario_rows(non_tite if method == "Original non-TITE" else tite, scenario, allocation, utility, alpha)
        duration = get_field(rows[0], "Duration")
        y = draw_duration_row(c, f"U{utility} - {method}", duration, y, shaded=(method != "Original non-TITE"))
    draw_rule(c, y + 5)
    return y - 15


def draw_page_chrome(c: canvas.Canvas, alpha: float, allocation_label: str, page_num: int, total_pages: int) -> None:
    c.setFillColor(NAVY)
    c.rect(0, PAGE_HEIGHT - 50, PAGE_WIDTH, 50, stroke=0, fill=1)
    draw_text(c, "Phase I/II Operating Characteristics", LEFT, PAGE_HEIGHT - 31, size=16, color=white, font="Helvetica-Bold")
    subtitle = f"Original non-TITE vs TITE - 56-day accrual - IPDE alpha = {fmt(alpha)} - allocation = {allocation_label} - N = 30"
    draw_text(c, subtitle, LEFT, PAGE_HEIGHT - 74, size=10.4, color=NAVY, font="Helvetica-Bold")
    draw_text(c, "Dose 1", VALUE_X[0], PAGE_HEIGHT - 94, size=8.5, color=BLUE, font="Helvetica-Bold")
    draw_text(c, "Dose 2", VALUE_X[1], PAGE_HEIGHT - 94, size=8.5, color=BLUE, font="Helvetica-Bold")
    draw_text(c, "Dose 3", VALUE_X[2], PAGE_HEIGHT - 94, size=8.5, color=BLUE, font="Helvetica-Bold")
    draw_text(c, "Dose 4", VALUE_X[3], PAGE_HEIGHT - 94, size=8.5, color=BLUE, font="Helvetica-Bold")
    draw_text(c, "Dose 5", VALUE_X[4], PAGE_HEIGHT - 94, size=8.5, color=BLUE, font="Helvetica-Bold")
    draw_text(c, "Stop %", STOP_X, PAGE_HEIGHT - 94, size=8.5, color=BLUE, font="Helvetica-Bold")
    draw_rule(c, PAGE_HEIGHT - 100, color=BLUE, width=0.8)
    footer = "Fixed priors: CRM toxicity Beta(0.15, 0.85); efficacy Beta(0.5, 0.5); efficacy additive discount Beta(0.15, 0.85)."
    draw_text(c, footer, LEFT, 35, size=7.3, color=GRAY)
    draw_text(c, "Non-TITE: efficacy threshold = 0.20, futility cutoff = 0.85. TITE uses the shared previous-dose additive efficacy model; 28- and 14-day results are pending.", LEFT, 22, size=7.3, color=GRAY)
    draw_right(c, f"Page {page_num} of {total_pages}", RIGHT, 22, size=7.3, color=GRAY)


def create_pdf(alpha: float, non_tite: dict[tuple, list[dict[str, str]]], tite: dict[tuple, list[dict[str, str]]], utilities: dict[int, dict[int, list[float]]]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"phase12_tite_56day_comparison_nontite_fut0p85_{alpha_tag(alpha)}_scenarios_1_to_37_tables.pdf"
    c = canvas.Canvas(str(path), pagesize=(PAGE_WIDTH, PAGE_HEIGHT), pageCompression=1)
    total_pages = 38
    page_num = 0
    for allocation, allocation_label in ALLOCATIONS:
        for first in range(1, 38, 2):
            page_num += 1
            draw_page_chrome(c, alpha, allocation_label, page_num, total_pages)
            draw_scenario(c, PAGE_HEIGHT - 126, first, allocation, alpha, non_tite, tite, utilities)
            if first + 1 <= 37:
                draw_scenario(c, PAGE_HEIGHT - 538, first + 1, allocation, alpha, non_tite, tite, utilities)
            c.showPage()
    c.save()
    return path


def main() -> None:
    utilities = utilities_from_reference(UTILITY_PDF)
    non_tite, tite = build_datasets()
    paths = [create_pdf(alpha, non_tite, tite, utilities) for alpha in ALPHAS]
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
