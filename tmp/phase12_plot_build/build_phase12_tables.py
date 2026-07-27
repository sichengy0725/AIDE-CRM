"""Create BOIN12-style Phase I/II operating-characteristic comparison tables.

This script consumes the already-validated PDF extraction JSON created during
the 7-27 comparison run, the EffTox HTML summaries, and the AIDE dose summary.
It produces one 10-page landscape PDF for N=30 and one for N=60.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "Presentation 7-27-2026"
OUT_DIR = INPUT_DIR / "Comparison Figures"
DEBUG_PATH = ROOT / "tmp" / "phase12_plot_build" / "pdf_parse_debug.json"
PAGE_SIZE = (18 * inch, 14 * inch)
PAGE_WIDTH, PAGE_HEIGHT = PAGE_SIZE


def number(value: str | float | int | None) -> float:
    if value in (None, "", "NA", "NaN"):
        raise ValueError(f"Missing numeric value: {value!r}")
    return float(value)


def display(value: float, decimals: int = 1) -> str:
    text = f"{value:.{decimals}f}"
    return text.rstrip("0").rstrip(".")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_efftox(path: Path) -> dict[int, dict[str, object]]:
    html = path.read_text(encoding="utf-8")
    starts = list(re.finditer(r'<tr><td\s+colspan\s*=\s*"8"><b>\s*(\d+)\s*</b></tr>', html, re.IGNORECASE))
    records: dict[int, dict[str, object]] = {}

    for index, start in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(html)
        block = html[start.start():end]

        def cells_for(label: str) -> list[str]:
            found = re.search(label, block, re.IGNORECASE)
            if not found:
                return []
            close = block.find("</tr>", found.start())
            next_row = block.find("<tr", found.start() + 2)
            end_candidates = [value for value in (close, next_row) if value >= 0]
            fragment_end = min(end_candidates) if end_candidates else len(block)
            fragment = block[found.start():fragment_end]
            raw_cells = re.findall(r"<td[^>]*>([\s\S]*?)</td>", fragment, re.IGNORECASE)
            return [re.sub(r"\s+", " ", re.sub(r"<[^>]+>|&nbsp;", " ", cell)).strip() for cell in raw_cells]

        selected_cells = cells_for(r"%\s*selected")
        treated_cells = cells_for(r"#\s*Patients\s*Treated")
        selected_values = [number(value) for value in selected_cells[-6:]]
        treated_values = [number(value) for value in treated_cells[-6:-1]]
        scenario = int(start.group(1))
        records[scenario] = {
            "treated": treated_values,
            "selected": selected_values[:5],
            "stop": selected_values[5],
        }

    if set(records) != set(range(1, 38)):
        raise ValueError(f"EffTox source should contain scenarios 1-37; found {sorted(records)}")
    return records


def read_truth() -> dict[int, dict[str, list[float]]]:
    rows = read_csv(INPUT_DIR / "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv")
    truth = {}
    for row in rows:
        scenario = int(row["Scenario"])
        truth[scenario] = {
            "dlt": [number(row[f"Tox_Dose{dose}"]) for dose in range(1, 6)],
            "efficacy": [number(row[f"Eff_Dose{dose}"]) for dose in range(1, 6)],
            "utility": [number(row[f"Utility3_Dose{dose}"]) for dose in range(1, 6)],
        }
    if set(truth) != set(range(1, 38)):
        raise ValueError("Truth summary must contain scenarios 1-37.")
    return truth


def read_aide(nmax: int, allocation: str) -> dict[int, dict[str, object]]:
    rows = read_csv(INPUT_DIR / "AIDE_phase_I_II_IDX_1001_to_2000_dose_summary.csv")
    output: dict[int, dict[str, object]] = {}
    for scenario in range(1, 38):
        selected_rows = [
            row for row in rows
            if int(row["Scenario"]) == scenario
            and int(row["Nmax"]) == nmax
            and row["Allocation"] == allocation
            and int(row["Utility_Type"]) == 3
            and abs(number(row["Lambda_T"]) - 1) < 1e-12
        ]
        if len(selected_rows) != 5:
            raise ValueError(
                f"Expected five AIDE rows for scenario {scenario}, N={nmax}, allocation={allocation}; found {len(selected_rows)}."
            )
        selected_rows.sort(key=lambda row: int(row["Dose"]))
        output[scenario] = {
            "treated": [number(row["Pts_Treated"]) for row in selected_rows],
            "selected": [number(row["OBD_Selection_pct"]) for row in selected_rows],
            "stop": number(selected_rows[0]["No_OBD_Selection_pct"]),
        }
    return output


def draw_header(pdf: canvas.Canvas, nmax: int, page_number: int, page_count: int) -> None:
    left = 80
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 28)
    pdf.drawCentredString(PAGE_WIDTH / 2, PAGE_HEIGHT - 56, "Phase I/II Operating Characteristics")
    pdf.setFont("Helvetica", 14)
    pdf.drawCentredString(PAGE_WIDTH / 2, PAGE_HEIGHT - 78, f"Method comparison - N = {nmax}")
    pdf.setStrokeColor(colors.HexColor("#999999"))
    pdf.setLineWidth(0.8)
    pdf.line(left, PAGE_HEIGHT - 92, PAGE_WIDTH - left, PAGE_HEIGHT - 92)

    label_width = 330
    dose_width = 130
    stop_width = 156
    x = left + label_width
    header_y = PAGE_HEIGHT - 110
    pdf.setFont("Helvetica-Bold", 13)
    for dose in range(1, 6):
        pdf.drawCentredString(x + dose_width / 2, header_y, f"Dose {dose}")
        x += dose_width
    pdf.drawCentredString(x + stop_width / 2, header_y, "Stop %")

    pdf.setFont("Helvetica", 9)
    footer = "AIDE phase I = one-stage; AIDE two-stage = Phase I/II. AIDE rows use Utility_Type = 3 and Lambda_T = 1."
    pdf.drawString(left, 38, footer)
    pdf.drawRightString(PAGE_WIDTH - left, 38, f"Page {page_number} of {page_count}")


def draw_row(
    pdf: canvas.Canvas,
    y: float,
    label: str,
    values: list[float],
    stop: float | None,
    shaded: bool = False,
    decimals: int = 1,
) -> float:
    left = 80
    label_width = 330
    dose_width = 130
    stop_width = 156
    row_height = 13.2
    table_width = label_width + 5 * dose_width + stop_width

    if shaded:
        pdf.setFillColor(colors.HexColor("#F6F6F6"))
        pdf.rect(left, y - row_height + 2, table_width, row_height, fill=1, stroke=0)
    pdf.setFillColor(colors.HexColor("#202020"))
    pdf.setFont("Helvetica", 9.5)
    pdf.drawString(left + 8, y - 9.6, label)

    x = left + label_width
    for value in values:
        pdf.drawCentredString(x + dose_width / 2, y - 9.6, display(value, decimals))
        x += dose_width
    if stop is not None:
        pdf.drawCentredString(x + stop_width / 2, y - 9.6, display(stop))

    pdf.setStrokeColor(colors.HexColor("#D0D0D0"))
    pdf.setLineWidth(0.45)
    pdf.line(left, y - row_height, left + table_width, y - row_height)
    return y - row_height


def draw_scenario(
    pdf: canvas.Canvas,
    y: float,
    scenario: int,
    truth: dict[str, list[float]],
    methods: list[tuple[str, dict[str, object]]],
) -> float:
    left = 80
    label_width = 330
    dose_width = 130
    stop_width = 156
    table_width = label_width + 5 * dose_width + stop_width
    pdf.setFillColor(colors.HexColor("#EFEFEF"))
    pdf.rect(left, y - 16, table_width, 18, fill=1, stroke=0)
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 13)
    pdf.drawString(left + 8, y - 12, f"Scenario {scenario}")
    y -= 19

    y = draw_row(pdf, y, "DLT rate", truth["dlt"], None, decimals=2)
    y = draw_row(pdf, y, "Efficacy rate", truth["efficacy"], None, decimals=2)
    y = draw_row(pdf, y, "Utility", truth["utility"], None)
    y -= 2

    for method_label, record in methods:
        y = draw_row(pdf, y, f"{method_label} - No. Pts treated", record["treated"], None)
        y = draw_row(pdf, y, f"{method_label} - Select %", record["selected"], float(record["stop"]), shaded=True)
    y -= 9
    return y


def write_table_pdf(nmax: int, source_rows: dict[str, dict[int, dict[str, object]]], truth: dict[int, dict[str, list[float]]]) -> None:
    groups = [list(range(start, min(start + 4, 38))) for start in range(1, 38, 4)]
    output = OUT_DIR / f"phase12_all_methods_N{nmax}_scenarios_1_to_37_tables.pdf"
    pdf = canvas.Canvas(str(output), pagesize=PAGE_SIZE)
    for page_index, scenarios in enumerate(groups, start=1):
        draw_header(pdf, nmax, page_index, len(groups))
        y = PAGE_HEIGHT - 130
        for scenario in scenarios:
            methods = [
                ("BOIN12", source_rows["BOIN12"][scenario]),
                ("U-BOIN", source_rows["U-BOIN"][scenario]),
                ("EffTox", source_rows["EffTox"][scenario]),
                ("AIDE phase I", source_rows["AIDE phase I"][scenario]),
                ("AIDE two-stage", source_rows["AIDE two-stage"][scenario]),
            ]
            y = draw_scenario(pdf, y, scenario, truth[scenario], methods)
        pdf.showPage()
    pdf.save()


def write_readme() -> None:
    (OUT_DIR / "README.txt").write_text(
        "Generated 2026-07-26 from the Presentation 7-27-2026 result files.\n\n"
        "Random-prior figures: four six-panel figures (alpha 0, 0.3, 0.6, 0.9) comparing the four Beta priors for Discount CRM only.\n"
        "All-method figures: four six-panel figures (same alphas) comparing BOIN, CRM, Alpha-CRM, IPCRM, and all four Discount CRM prior settings; CFO is excluded.\n"
        "Phase I/II tables: one BOIN12-style operating-characteristic table PDF for N=30 and one for N=60, covering scenarios 1-37. Each table presents truth (DLT, efficacy, Utility3) and treated/selection results for BOIN12, U-BOIN, EffTox, AIDE phase I, and AIDE two-stage.\n\n"
        "Phase I/II table conventions:\n"
        "- Stop % is the source method's no-dose-selection probability.\n"
        "- AIDE phase I uses Allocation=one_stage; AIDE two-stage uses Allocation=two_stage.\n"
        "- AIDE rows use Utility_Type=3 and Lambda_T=1, matching the supplied utility example.\n",
        encoding="utf-8",
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    debug = json.loads(DEBUG_PATH.read_text(encoding="utf-8"))
    truth = read_truth()

    for nmax, suffix in ((30, "30"), (60, "60")):
        boin = {int(record["scenario"]): record for record in debug[f"boin{suffix}"]}
        uboin = {int(record["scenario"]): record for record in debug[f"uboin{suffix}"]}
        efftox = parse_efftox(INPUT_DIR / f"EffTox N = {nmax}.html")
        for source_name, records in (("BOIN12", boin), ("U-BOIN", uboin), ("EffTox", efftox)):
            if set(records) != set(range(1, 38)):
                raise ValueError(f"{source_name} N={nmax} is missing one or more scenarios.")
        write_table_pdf(
            nmax,
            {
                "BOIN12": boin,
                "U-BOIN": uboin,
                "EffTox": efftox,
                "AIDE phase I": read_aide(nmax, "one_stage"),
                "AIDE two-stage": read_aide(nmax, "two_stage"),
            },
            truth,
        )
    write_readme()


if __name__ == "__main__":
    main()
