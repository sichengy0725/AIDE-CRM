from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path.cwd()
PRESENTATION = ROOT / "Presentation 8-17-2026"
TABLES_DIR = PRESENTATION / "Table and Plots"
NEW_DESIGN_CSV = PRESENTATION / (
    "AIDE_phase_I_II_modelsadditive_newdesign_shared_N30_ncycle2_"
    "rp0p15x0p85_rate56d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_"
    "ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv"
)
BOIN12_CSV = PRESENTATION / "BOIN12_v1.4.2.0_Operating Characteristics_2026-08-09 170536.992978_fut0.85.csv"
ALPHA_PDF = TABLES_DIR / "nontite_newdesign_alpha_0_0p3_0p6_0p9_scenarios_1_16_20_24_27_38_tables.pdf"
METHOD_PDF = TABLES_DIR / "nontite_alpha0_newdesign_vs_uboin_boin12_scenarios_1_16_20_24_27_38_tables.pdf"

SCENARIOS = [1, 16, 20, 24, 27, 38]
ALPHAS = [0.0, 0.3, 0.6, 0.9]
DESIGN_ORDER = ["one", "high", "top2"]

NAVY = colors.HexColor("#173F5F")
TEAL = colors.HexColor("#0F766E")
MID_GRAY = colors.HexColor("#E5E7EB")
LIGHT_GRAY = colors.HexColor("#F8FAFC")
NOTE_BG = colors.HexColor("#FFF7ED")
ALPHA_BG = {
    0.0: colors.HexColor("#E2F0D9"),
    0.3: colors.HexColor("#DBEAFE"),
    0.6: colors.HexColor("#FEF3C7"),
    0.9: colors.HexColor("#F3E8FF"),
}


def num(value: str | float | int | None) -> float | None:
    if value in (None, "", "NA"):
        return None
    return float(value)


def fmt(value: float | None, digits: int = 1) -> str:
    if value is None:
        return "-"
    return f"{value:.{digits}f}"


def design_info(row: dict[str, str]) -> tuple[str, str, str]:
    if row["Allocation"] == "one_stage":
        return "one", "One-stage", "Highest utility throughout"
    if row["Stage2_Allocation"] == "highest_utility":
        return "high", "Two-stage highest utility", "Stage II: highest utility"
    if row["Stage2_Allocation"] == "top2_randomized":
        return "top2", "Two-stage top-2 randomized", "Stage II: top-2 randomized"
    raise ValueError(f"Unexpected design row: {row['Allocation']} / {row['Stage2_Allocation']}")


def read_new_design() -> list[dict]:
    with NEW_DESIGN_CSV.open(newline="", encoding="utf-8-sig") as handle:
        raw = list(csv.DictReader(handle))
    grouped: dict[tuple[int, str, float], list[dict[str, str]]] = defaultdict(list)
    for row in raw:
        scenario = int(row["Scenario"])
        tox_alpha = num(row["Toxicity_IPDE_Alpha"])
        eff_alpha = num(row["Efficacy_IPDE_Alpha"])
        if scenario not in SCENARIOS or tox_alpha not in ALPHAS or tox_alpha != eff_alpha:
            continue
        design_key, _, _ = design_info(row)
        grouped[(scenario, design_key, tox_alpha)].append(row)

    if len(grouped) != 72 or any(len(rows) != 5 for rows in grouped.values()):
        raise RuntimeError("Expected a complete six-scenario x three-design x four-alpha grid.")

    records: list[dict] = []
    for (scenario, key, alpha), dose_rows in grouped.items():
        dose_rows = sorted(dose_rows, key=lambda item: int(item["Dose"]))
        row = dose_rows[0]
        _, label, rule = design_info(row)
        records.append(
            {
                "scenario": scenario,
                "key": key,
                "design": label,
                "rule": rule,
                "alpha": alpha,
                "true_obd": int(row["True_OBD"]),
                "tox": [num(item["True_DLT_rate"]) for item in dose_rows],
                "eff": [num(item["True_Efficacy_rate"]) for item in dose_rows],
                "obd": [num(item["OBD_Selection_pct"]) for item in dose_rows],
                "treated": [num(item["Pts_Treated"]) for item in dose_rows],
                "no_obd": num(row["No_OBD_Selection_pct"]),
                "early_stop": num(row["Early_Stopping_pct"]),
                "total_admin": num(row["Total_Administrations"]),
                "unique_patients": num(row["Total_Unique_Patients"]),
                "duration": num(row["Duration"]),
                "ntrial": int(row["ntrial_from_files"]),
            }
        )
    records.sort(key=lambda item: (item["scenario"], DESIGN_ORDER.index(item["key"]), item["alpha"]))
    return records


UBOIN = {
    1: {"treated": [9, 5.4, 4.7, 3.9, 3.7], "selection": [15.6, 17.4, 20.6, 18.7, 13.7], "outcome": 14, "outcome_label": "Stop %"},
    16: {"treated": [6.1, 6.2, 6.4, 6.4, 5.4], "selection": [28.1, 26.8, 21.1, 16.3, 7.7], "outcome": 0, "outcome_label": "Stop %"},
    20: {"treated": [5.9, 10.8, 7, 4, 2.4], "selection": [8.5, 47.5, 18.6, 12.4, 13], "outcome": 0, "outcome_label": "Stop %"},
    24: {"treated": [3.5, 4, 5.1, 9.8, 7.9], "selection": [1.5, 3.2, 8.2, 60.1, 27], "outcome": 0, "outcome_label": "Stop %"},
    27: {"treated": [6.9, 8.1, 8.8, 4.3, 1.8], "selection": [8.7, 20.9, 38.9, 16.5, 14.5], "outcome": 0.5, "outcome_label": "Stop %"},
    38: {"treated": [5.2, 5.5, 5.6, 6.3, 7.8], "selection": [6.2, 8.9, 13, 13.3, 58.6], "outcome": 0, "outcome_label": "Stop %"},
}


def read_boin12() -> dict[int, dict]:
    parsed: dict[int, dict] = {}
    current: int | None = None
    with BOIN12_CSV.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            first = row[0].strip() if row else ""
            if first.startswith("Scenario "):
                current = int(first.split()[1])
                parsed.setdefault(current, {})
            elif current is not None and first == "No. Pts treated":
                parsed[current]["treated"] = [num(value) for value in row[1:6]]
            elif current is not None and first == "Select %":
                parsed[current]["selection"] = [num(value) for value in row[1:6]]
                parsed[current]["outcome"] = num(row[6])
                parsed[current]["outcome_label"] = "No selection %"
    parsed[38] = {
        "treated": [4.2, 4.7, 5.2, 6.1, 9.6],
        "selection": [13.2, 13.8, 16.5, 13.6, 41.5],
        "outcome": 1.4,
        "outcome_label": "No selection %",
    }
    if any(scenario not in parsed for scenario in SCENARIOS):
        raise RuntimeError("BOIN12 comparator rows are missing.")
    return parsed


styles = getSampleStyleSheet()
TITLE = ParagraphStyle("Title", parent=styles["Heading1"], fontName="Helvetica-Bold", fontSize=16.5, leading=19, textColor=NAVY, alignment=TA_CENTER, spaceAfter=1)
SUBTITLE = ParagraphStyle("Subtitle", parent=styles["Normal"], fontName="Helvetica", fontSize=9.8, leading=12, textColor=colors.HexColor("#475569"), alignment=TA_CENTER, spaceAfter=4)
SCENARIO = ParagraphStyle("Scenario", parent=styles["Heading2"], fontName="Helvetica-Bold", fontSize=11.5, leading=14, textColor=colors.HexColor("#111827"), spaceBefore=3, spaceAfter=3)
SECTION = ParagraphStyle("Section", parent=styles["Heading3"], fontName="Helvetica-Bold", fontSize=9, leading=11, textColor=NAVY, spaceBefore=4, spaceAfter=2)
NOTE = ParagraphStyle("Note", parent=styles["Normal"], fontName="Helvetica-Oblique", fontSize=7.6, leading=9, textColor=colors.HexColor("#7C2D12"), alignment=TA_LEFT)
FOOT = ParagraphStyle("Foot", parent=styles["Normal"], fontName="Helvetica", fontSize=6.8, leading=8, textColor=colors.HexColor("#475569"), alignment=TA_LEFT)


def base_table(data: list[list[str]], widths: list[float], header: bool = True, font_size: float = 7.2) -> Table:
    table = Table(data, colWidths=widths, repeatRows=1 if header else 0, hAlign="LEFT")
    style = [
        ("FONTNAME", (0, 0), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), font_size),
        ("LEADING", (0, 0), (-1, -1), font_size + 1.2),
        ("ALIGN", (1, 0), (-1, -1), "RIGHT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.0),
        ("TOPPADDING", (0, 0), (-1, -1), 2.0),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#CBD5E1")),
    ]
    if header:
        style.extend([
            ("BACKGROUND", (0, 0), (-1, 0), TEAL),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("ALIGN", (0, 0), (-1, 0), "CENTER"),
        ])
    table.setStyle(TableStyle(style))
    return table


def alpha_color_rows(table: Table, alphas_by_row: list[float], alpha_column: int = 1) -> None:
    commands = []
    for index, alpha in enumerate(alphas_by_row, start=1):
        commands.append(("BACKGROUND", (alpha_column, index), (alpha_column, index), ALPHA_BG[alpha]))
        if index % 2 == 0:
            commands.append(("BACKGROUND", (0, index), (-1, index), colors.HexColor("#FAFAFA")))
            commands.append(("BACKGROUND", (alpha_column, index), (alpha_column, index), ALPHA_BG[alpha]))
    table.setStyle(TableStyle(commands))


def page_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(colors.HexColor("#475569"))
    canvas.drawString(20 * mm, 11 * mm, "Presentation 8-17-2026 | All populated new-design rows represent 1,000 simulations.")
    canvas.drawRightString(277 * mm, 11 * mm, f"Page {doc.page}")
    canvas.restoreState()


def truth_table(record: dict) -> Table:
    data = [
        ["Dose", "D1", "D2", "D3", "D4", "D5"],
        ["DLT rate", *[fmt(value, 2) for value in record["tox"]]],
        ["Efficacy rate", *[fmt(value, 2) for value in record["eff"]]],
    ]
    table = base_table(data, [38, 48, 48, 48, 48, 48], font_size=7.3)
    table.setStyle(TableStyle([("ALIGN", (0, 1), (0, -1), "LEFT")]))
    return table


def alpha_story(records: list[dict]) -> list:
    story = []
    for position, scenario in enumerate(SCENARIOS):
        scenario_records = [record for record in records if record["scenario"] == scenario]
        truth = scenario_records[0]
        story.extend([
            Paragraph("Non-TITE New-Design Alpha Sensitivity", TITLE),
            Paragraph("Comparison of IPDE alpha = 0, 0.3, 0.6, and 0.9 (toxicity alpha = efficacy alpha)", SUBTITLE),
            Paragraph(f"Scenario {scenario} | True OBD: Dose {truth['true_obd']}", SCENARIO),
            truth_table(truth),
            Spacer(1, 3),
            Paragraph("OBD selection (%)", SECTION),
        ])
        selection_data = [["New-design allocation", "Alpha", "D1", "D2", "D3", "D4", "D5", "No OBD"]]
        for record in scenario_records:
            selection_data.append([
                record["design"], fmt(record["alpha"]), *[fmt(value) for value in record["obd"]], fmt(record["no_obd"])
            ])
        selection_table = base_table(selection_data, [150, 35, 53, 53, 53, 53, 53, 60], font_size=7.0)
        alpha_color_rows(selection_table, [record["alpha"] for record in scenario_records])
        story.extend([selection_table, Spacer(1, 4), Paragraph("Mean patients treated by dose and run-level results", SECTION)])
        allocation_data = [["New-design allocation", "Alpha", "D1", "D2", "D3", "D4", "D5", "Total admin.", "Unique pts", "Duration", "Early stop"]]
        for record in scenario_records:
            allocation_data.append([
                record["design"], fmt(record["alpha"]), *[fmt(value) for value in record["treated"]],
                fmt(record["total_admin"]), fmt(record["unique_patients"]), fmt(record["duration"]), fmt(record["early_stop"]),
            ])
        allocation_table = base_table(allocation_data, [120, 30, 39, 39, 39, 39, 39, 59, 50, 55, 48], font_size=6.55)
        alpha_color_rows(allocation_table, [record["alpha"] for record in scenario_records])
        story.extend([
            allocation_table,
            Spacer(1, 5),
            Paragraph("Source: AIDE_phase_I_II_modelsadditive_newdesign_shared_N30_ncycle2...IDX_0001_to_1000_dose_summary.csv.", FOOT),
        ])
        if position < len(SCENARIOS) - 1:
            story.append(PageBreak())
    return story


def methods_story(records: list[dict], boin12: dict[int, dict]) -> list:
    story = []
    for page_index, scenario in enumerate(SCENARIOS):
        story.extend([
            Paragraph("Alpha = 0: New Design vs U-BOIN vs BOIN12", TITLE),
            Paragraph("New-design results use IPDE alpha = 0; U-BOIN and BOIN12 are shown from Presentation 8-17-2026 sources.", SUBTITLE),
            Paragraph("Outcome labels remain source-specific: New Design = No OBD selection; U-BOIN = Stop %; BOIN12 = No selection %.", NOTE),
            Spacer(1, 3),
        ])
        scenario_records = [record for record in records if record["scenario"] == scenario and record["alpha"] == 0]
        truth = scenario_records[0]
        story.extend([
            Paragraph(f"Scenario {scenario} | True OBD: Dose {truth['true_obd']}", SCENARIO),
            truth_table(truth),
            Spacer(1, 2),
            Paragraph("Selection (%)", SECTION),
        ])
        entries = [
            (record["design"], record["obd"], record["no_obd"], "No OBD") for record in scenario_records
        ] + [
            ("U-BOIN", UBOIN[scenario]["selection"], UBOIN[scenario]["outcome"], "Stop"),
            ("BOIN12", boin12[scenario]["selection"], boin12[scenario]["outcome"], "No selection"),
        ]
        select_data = [["Method", "D1", "D2", "D3", "D4", "D5", "Outcome", "Definition"]]
        for label, selection, outcome, definition in entries:
            select_data.append([label, *[fmt(value) for value in selection], fmt(outcome), definition])
        selection_table = base_table(select_data, [150, 55, 55, 55, 55, 55, 60, 70], font_size=7.0)
        selection_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 1), (-1, 3), colors.HexColor("#E2F0D9")),
            ("BACKGROUND", (0, 4), (-1, 5), colors.HexColor("#DBEAFE")),
        ]))
        story.extend([selection_table, Spacer(1, 3), Paragraph("Mean patients treated by dose", SECTION)])
        allocation_entries = [(record["design"], record["treated"]) for record in scenario_records] + [
            ("U-BOIN", UBOIN[scenario]["treated"]),
            ("BOIN12", boin12[scenario]["treated"]),
        ]
        allocation_data = [["Method", "D1", "D2", "D3", "D4", "D5"]]
        allocation_data.extend([[label, *[fmt(value) for value in treated]] for label, treated in allocation_entries])
        allocation_table = base_table(allocation_data, [170, 68, 68, 68, 68, 68], font_size=7.0)
        allocation_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 1), (-1, 3), colors.HexColor("#E2F0D9")),
            ("BACKGROUND", (0, 4), (-1, 5), colors.HexColor("#DBEAFE")),
        ]))
        story.extend([allocation_table, Spacer(1, 4)])
        story.append(Paragraph("Sources: specified new-design summary; U-BOIN presentation PDFs; BOIN12 operating-characteristics CSV (Scenarios 1, 16, 20, 24, 27) and BOIN12 Sce38 PDF.", FOOT))
        if page_index < len(SCENARIOS) - 1:
            story.append(PageBreak())
    return story


def build_document(path: Path, story: list) -> None:
    doc = SimpleDocTemplate(
        str(path),
        pagesize=landscape(A4),
        leftMargin=14 * mm,
        rightMargin=14 * mm,
        topMargin=8 * mm,
        bottomMargin=14 * mm,
        title=path.stem,
        author="Codex",
    )
    doc.build(story, onFirstPage=page_footer, onLaterPages=page_footer)


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    records = read_new_design()
    boin12 = read_boin12()
    build_document(ALPHA_PDF, alpha_story(records))
    build_document(METHOD_PDF, methods_story(records, boin12))
    print(ALPHA_PDF)
    print(METHOD_PDF)


if __name__ == "__main__":
    main()
