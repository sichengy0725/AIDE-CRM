"""Build August 10 Phase I/II all-method comparison tables by futility cutoff.

The layout mirrors the August 3 N=30 all-method table.  Comparator results
come from the new August 10 Raw Data files.  AIDE rows are restricted to the
requested shared-alpha model/prior configuration and the alpha=0 baseline used
by the August 3 all-method table.
"""

from __future__ import annotations

import csv
import html
import math
import re
from pathlib import Path
from typing import Any

import pdfplumber
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "Presentation 8-10-2026" / "Raw Data"
OUT_DIR = ROOT / "Presentation 8-10-2026" / "Table and Plots" / "Phase I-II"
TRUTH_CSV = ROOT / "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv"

PAGE_WIDTH, PAGE_HEIGHT = 18 * inch, 14 * inch
LEFT = 80
LABEL_WIDTH = 330
DOSE_WIDTH = 130
STOP_WIDTH = 156
TABLE_WIDTH = LABEL_WIDTH + 5 * DOSE_WIDTH + STOP_WIDTH
CUTOFFS = (0.75, 0.85)
SCENARIOS = set(range(1, 38))


def number(value: str | float | int | None) -> float:
    if value in (None, "", "NA", "NaN"):
        raise ValueError(f"Missing numeric value: {value!r}")
    return float(value)


def close(value: str | float, expected: float, tolerance: float = 1e-10) -> bool:
    return math.isclose(float(value), expected, abs_tol=tolerance, rel_tol=tolerance)


def display(value: float, decimals: int = 1) -> str:
    if not math.isfinite(value):
        return "-"
    text = f"{value:.{decimals}f}"
    return text.rstrip("0").rstrip(".")


def read_csv_dict(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def read_truth() -> dict[int, dict[str, list[float]]]:
    truth: dict[int, dict[str, list[float]]] = {}
    for row in read_csv_dict(TRUTH_CSV):
        scenario = int(row["Scenario"])
        truth[scenario] = {
            "dlt": [number(row[f"Tox_Dose{dose}"]) for dose in range(1, 6)],
            "efficacy": [number(row[f"Eff_Dose{dose}"]) for dose in range(1, 6)],
            "utility2": [number(row[f"Utility2_Dose{dose}"]) for dose in range(1, 6)],
            "utility3": [number(row[f"Utility3_Dose{dose}"]) for dose in range(1, 6)],
        }
    if set(truth) != SCENARIOS:
        raise ValueError("The truth file must contain scenarios 1 through 37.")
    return truth


def parse_boin12_csv(path: Path) -> dict[int, dict[str, Any]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.reader(stream))
    records: dict[int, dict[str, Any]] = {}
    current: int | None = None
    for row in rows[1:]:
        label = row[0].strip()
        scenario_match = re.fullmatch(r"Scenario\s+(\d+)", label)
        if scenario_match:
            current = int(scenario_match.group(1))
            records[current] = {}
            continue
        if current is None:
            continue
        if label == "No. Pts treated":
            records[current]["treated"] = [number(value) for value in row[1:6]]
        elif label == "Select %":
            records[current]["selected"] = [number(value) for value in row[1:6]]
            records[current]["stop"] = number(row[6])
    if set(records) != SCENARIOS or any(set(record) != {"treated", "selected", "stop"} for record in records.values()):
        raise ValueError(f"BOIN12 source is incomplete: {path.name}")
    return records


def parse_operating_characteristics_pdf(path: Path, scenario_offset: int) -> dict[int, dict[str, Any]]:
    with pdfplumber.open(path) as document:
        text = "\n".join(page.extract_text() or "" for page in document.pages)
    blocks = list(re.finditer(r"^Scenario\s+(\d+)\s*$", text, flags=re.MULTILINE))
    records: dict[int, dict[str, Any]] = {}
    for index, match in enumerate(blocks):
        block_end = blocks[index + 1].start() if index + 1 < len(blocks) else len(text)
        block = text[match.end():block_end]

        def values_for(label: str, count: int) -> list[float]:
            found = re.search(re.escape(label) + r"\s+([^\n]+)", block)
            if not found:
                raise ValueError(f"Could not find {label} for scenario {match.group(1)} in {path.name}.")
            values = [float(value) for value in re.findall(r"-?(?:\d+\.?\d*|\.\d+)", found.group(1))]
            if len(values) < count:
                raise ValueError(f"Could not parse {count} {label} values for scenario {match.group(1)} in {path.name}.")
            return values[-count:]

        selected = values_for("Select %", 6)
        records[int(match.group(1)) + scenario_offset] = {
            "treated": values_for("# Pts treated", 5),
            "selected": selected[:5],
            "stop": selected[5],
        }
    if not records:
        raise ValueError(f"No scenarios found in {path.name}.")
    return records


def parse_uboin(cutoff: float) -> dict[int, dict[str, Any]]:
    records: dict[int, dict[str, Any]] = {}
    groups = (("1-10", 10, 0), ("11-20", 10, 10), ("21-30", 10, 20), ("31-37", 7, 30))
    for label, expected_count, offset in groups:
        path = RAW_DIR / f"U-BOIN Sce {label} fut {cutoff:.2f}.pdf"
        parsed = parse_operating_characteristics_pdf(path, offset)
        if len(parsed) != expected_count:
            raise ValueError(f"Expected {expected_count} scenarios in {path.name}; found {len(parsed)}.")
        records.update(parsed)
    if set(records) != SCENARIOS:
        raise ValueError(f"U-BOIN source is incomplete at cutoff {cutoff:.2f}.")
    return records


def parse_efftox(path: Path) -> dict[int, dict[str, Any]]:
    content = path.read_text(encoding="utf-8")
    starts = list(re.finditer(r'<tr><td\s+colspan\s*=\s*"8"><b>\s*(\d+)\s*</b></tr>', content, re.IGNORECASE))
    records: dict[int, dict[str, Any]] = {}
    for index, start in enumerate(starts):
        block_end = starts[index + 1].start() if index + 1 < len(starts) else len(content)
        block = content[start.start():block_end]

        def values_after(label: str, terminator: str) -> list[float]:
            found = re.search(label + r"\s*</b>\s*</td>([\s\S]*?)" + terminator, block, re.IGNORECASE)
            if not found:
                return []
            cells = re.findall(r"<td[^>]*>([\s\S]*?)</td>", found.group(1), re.IGNORECASE)
            values: list[float] = []
            for cell in cells:
                value = re.sub(r"\s+", " ", re.sub(r"<[^>]+>|&nbsp;", " ", html.unescape(cell))).strip()
                values.append(float("nan") if value in {"", "-", "--"} else float(value))
            return values

        selected = values_after(r"%\s*selected", r"<tr\b")
        treated = values_after(r"#\s*Patients\s*Treated", r"</tr>")
        if len(selected) != 6 or len(treated) != 6:
            raise ValueError(f"Could not parse five doses and Stop % for EffTox scenario {start.group(1)}.")
        records[int(start.group(1))] = {
            "treated": treated[:5],
            "selected": selected[:5],
            "stop": selected[5],
        }
    if set(records) != SCENARIOS:
        raise ValueError(f"EffTox source is incomplete: {path.name}")
    return records


def select_aide_source(cutoff: float) -> Path:
    matches = list(RAW_DIR.glob(f"AIDE*ep0p5x0p5*eth0p2*fut0p{int(cutoff * 100):02d}*dose_summary.csv"))
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one matching AIDE file for cutoff {cutoff:.2f}; found {matches}.")
    return matches[0]


def extract_aide_records(path: Path, cutoff: float, allocation: str, utility_type: int, lambda_t: float, truth: dict[int, dict[str, list[float]]]) -> dict[int, dict[str, Any]]:
    requested_priors = {
        "CRM_Prior_a": 0.15,
        "CRM_Prior_b": 0.85,
        "Efficacy_Prior_a": 0.5,
        "Efficacy_Prior_b": 0.5,
        "Efficacy_Carryover_Prior_a": 0.15,
        "Efficacy_Carryover_Prior_b": 0.85,
        "Efficacy_Additive_Alpha_Prior_a": 0.15,
        "Efficacy_Additive_Alpha_Prior_b": 0.85,
        "Toxicity_IPDE_Alpha": 0.0,
        "Efficacy_IPDE_Alpha": 0.0,
        "Efficacy_Threshold": 0.20,
        "Futility_Cutoff": cutoff,
    }
    selected: list[dict[str, str]] = []
    for row in read_csv_dict(path):
        if (
            row["Model_ID"] != "additive_shared"
            or row["Carryover_Model"] != "additive_shared"
            or row["CRM_r_Model"] != "previous_dose"
            or row["Efficacy_Model"] != "previous_dose_additive"
            or row["Allocation"] != allocation
            or int(row["Nmax"]) != 30
            or int(row["Cycle_Max"]) != 2
            or int(row["IPDE_Design"]) != 1
            or int(row["Flexible_IPDE"]) != 1
            or row["Enrollment_Scheme"] != "continuous"
            or int(row["Utility_Type"]) != utility_type
            or not close(row["Lambda_T"], lambda_t)
            or not all(close(row[field], value) for field, value in requested_priors.items())
        ):
            continue
        selected.append(row)

    expected_rows = 37 * 5
    if len(selected) != expected_rows:
        raise ValueError(
            f"Expected {expected_rows} AIDE rows for cutoff={cutoff}, allocation={allocation}, U{utility_type}; found {len(selected)}."
        )

    records: dict[int, dict[str, Any]] = {}
    for scenario in range(1, 38):
        rows = sorted((row for row in selected if int(row["Scenario"]) == scenario), key=lambda row: int(row["Dose"]))
        if len(rows) != 5 or [int(row["Dose"]) for row in rows] != [1, 2, 3, 4, 5]:
            raise ValueError(f"AIDE rows are incomplete for scenario {scenario} at cutoff {cutoff:.2f}.")
        if any(int(row["n_valid"]) != int(row["ntrial_from_files"]) or int(row["n_valid"]) <= 0 for row in rows):
            raise ValueError(f"AIDE replicate counts are invalid for scenario {scenario} at cutoff {cutoff:.2f}.")
        for key, truth_key in (("True_DLT_rate", "dlt"), ("True_Efficacy_rate", "efficacy")):
            observed = [number(row[key]) for row in rows]
            if any(not close(value, expected) for value, expected in zip(observed, truth[scenario][truth_key], strict=True)):
                raise ValueError(f"AIDE truth mismatch in scenario {scenario}, field {key}.")
        records[scenario] = {
            "treated": [number(row["Pts_Treated"]) for row in rows],
            "selected": [number(row["OBD_Selection_pct"]) for row in rows],
            "stop": number(rows[0]["No_OBD_Selection_pct"]),
        }
    return records


def draw_row(pdf: canvas.Canvas, y: float, label: str, values: list[float], stop: float | None, *, decimals: int = 1) -> float:
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


def draw_scenario(pdf: canvas.Canvas, y: float, scenario: int, truth: dict[str, list[float]], methods: list[tuple[str, dict[str, Any]]]) -> float:
    height = 18
    pdf.setFillColor(colors.HexColor("#EFEFEF"))
    pdf.rect(LEFT, y - height + 2, TABLE_WIDTH, height, fill=1, stroke=0)
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 13)
    pdf.drawString(LEFT + 8, y - 12, f"Scenario {scenario}")
    y -= 19
    y = draw_row(pdf, y, "DLT rate", truth["dlt"], None, decimals=2)
    y = draw_row(pdf, y, "Efficacy rate", truth["efficacy"], None, decimals=2)
    y = draw_row(pdf, y, "Utility 2", truth["utility2"], None)
    y = draw_row(pdf, y, "Utility 3", truth["utility3"], None)
    y -= 2
    y = draw_group_header(pdf, y, "OBD selection (%)")
    for label, record in methods:
        y = draw_row(pdf, y, label, record["selected"], record["stop"])
    y = draw_group_header(pdf, y, "Mean patients treated")
    for label, record in methods:
        y = draw_row(pdf, y, label, record["treated"], None)
    return y - 9


def draw_header(pdf: canvas.Canvas, cutoff: float, page_number: int, page_count: int) -> None:
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 28)
    pdf.drawCentredString(PAGE_WIDTH / 2, PAGE_HEIGHT - 56, "Phase I/II Operating Characteristics")
    pdf.setFont("Helvetica", 14)
    pdf.drawCentredString(
        PAGE_WIDTH / 2,
        PAGE_HEIGHT - 78,
        f"Method comparison - N = 30, Cycle_Max = 2, efficacy futility cutoff = {cutoff:.2f}",
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
        "AIDE: shared-alpha additive model; toxicity Beta(0.15, 0.85); efficacy beta-binomial Beta(0.5, 0.5); efficacy alpha Beta(0.15, 0.85).",
    )
    pdf.drawString(
        LEFT,
        25,
        "AIDE rows use IPDE alpha = 0 and continuous enrollment. U2: Type = 2, Lambda_T = 0.3; U3: Type = 3, Lambda_T = 1.",
    )
    pdf.drawRightString(PAGE_WIDTH - LEFT, 25, f"Page {page_number} of {page_count}")


def build_pdf(cutoff: float, truth: dict[int, dict[str, list[float]]]) -> Path:
    aide_path = select_aide_source(cutoff)
    sources = {
        "BOIN12": parse_boin12_csv(next(RAW_DIR.glob(f"BOIN12*_fut{cutoff:.2f}.csv"))),
        "U-BOIN": parse_uboin(cutoff),
        "EffTox": parse_efftox(RAW_DIR / f"EffTox N = 30 Fut {cutoff:.2f}.html"),
        "AIDE one-stage (U2)": extract_aide_records(aide_path, cutoff, "one_stage", 2, 0.3, truth),
        "AIDE one-stage (U3)": extract_aide_records(aide_path, cutoff, "one_stage", 3, 1.0, truth),
        "AIDE two-stage (U2)": extract_aide_records(aide_path, cutoff, "two_stage", 2, 0.3, truth),
        "AIDE two-stage (U3)": extract_aide_records(aide_path, cutoff, "two_stage", 3, 1.0, truth),
    }
    for method, records in sources.items():
        if set(records) != SCENARIOS:
            raise ValueError(f"{method} does not contain every scenario at cutoff {cutoff:.2f}.")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tag = str(cutoff).replace(".", "p")
    output = OUT_DIR / f"phase12_all_methods_N30_fut{tag}_scenarios_1_to_37_tables.pdf"
    pdf = canvas.Canvas(str(output), pagesize=(PAGE_WIDTH, PAGE_HEIGHT), pageCompression=1)
    page_count = 19
    method_order = tuple(sources)
    page = 0
    for start in range(1, 38, 2):
        page += 1
        draw_header(pdf, cutoff, page, page_count)
        y = PAGE_HEIGHT - 130
        y = draw_scenario(pdf, y, start, truth[start], [(name, sources[name][start]) for name in method_order])
        if start + 1 <= 37:
            draw_scenario(pdf, y, start + 1, truth[start + 1], [(name, sources[name][start + 1]) for name in method_order])
        pdf.showPage()
    pdf.save()
    return output


def main() -> None:
    truth = read_truth()
    for cutoff in CUTOFFS:
        print(build_pdf(cutoff, truth))


if __name__ == "__main__":
    main()
