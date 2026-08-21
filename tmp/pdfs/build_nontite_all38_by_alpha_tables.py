"""Build all-scenario New non-TITE, U-BOIN, and BOIN12 comparison tables."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pdfplumber
from reportlab.lib.colors import HexColor, black
from reportlab.lib.pagesizes import landscape, letter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
PRESENTATION = ROOT / "Presentation 8-17-2026"
RAW_DATA = PRESENTATION / "Raw Data"
SOURCE_FILE = RAW_DATA / (
    "AIDE_phase_I_II_Sc1-38_N30_ncycle2_rp0p15x0p85_rate56d_ep0p5x0p5_"
    "cp0p15x0p85_eth0p2_fut0p85_ap0p15x0p85_IDX_0001_to_1000_"
    "newdesign_dose_summary.csv"
)
BOIN12_CSV = RAW_DATA / "BOIN12_v1.4.2.0_Operating Characteristics_2026-08-09 170536.992978_fut0.85.csv"
BOIN12_SCENARIO_38 = RAW_DATA / "BOIN12 Sce38.pdf"
UBOIN_GROUPS = (
    ("U-BOIN Sce 1-10 fut 0.85.pdf", 0, 10),
    ("U-BOIN Sce 11-20 fut 0.85.pdf", 10, 10),
    ("U-BOIN Sce 21-30 fut 0.85.pdf", 20, 10),
    ("U-BOIN Sce 31-37 fut 0.85.pdf", 30, 7),
)
UBOIN_SCENARIO_38 = RAW_DATA / "UBOIN Sce38.pdf"
OUTPUT_DIR = PRESENTATION / "Table and Plots" / "New Design" / "All Scenarios 1-38 By Alpha"

SCENARIOS = tuple(range(1, 39))
ALPHAS = (0.0, 0.3, 0.6, 0.9)
DESIGNS = (
    ("one", "one-stage"),
    ("high", "2-stage highest utility"),
    ("top2", "2-stage top-2 randomized"),
)

PAGE_WIDTH, _ = landscape(letter)
LEFT = 36
RIGHT = PAGE_WIDTH - 36
TABLE_WIDTH = RIGHT - LEFT
LABEL_WIDTH = 255
VALUE_WIDTH = (TABLE_WIDTH - LABEL_WIDTH) / 6
NAVY = HexColor("#173F5F")
SECTION = HexColor("#E9ECEF")
NEW_FILL = HexColor("#E6F3F1")
UBOIN_FILL = HexColor("#E5F1FB")
BOIN_FILL = HexColor("#EEF2FF")
GRID = HexColor("#C9CDD1")
MUTED = HexColor("#4D5965")


def alpha_tag(alpha: float) -> str:
    return "0" if alpha == 0 else f"{alpha:.1f}".replace(".", "p")


def fmt(value: float | None, digits: int = 1) -> str:
    if value is None:
        return "-"
    text = f"{float(value):.{digits}f}"
    return text.rstrip("0").rstrip(".")


def draw_text(
    pdf: canvas.Canvas,
    text: str,
    x: float,
    y: float,
    *,
    size: float = 8,
    bold: bool = False,
    color=black,
    align: str = "left",
) -> None:
    pdf.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    pdf.setFillColor(color)
    if align == "center":
        pdf.drawCentredString(x, y, text)
    elif align == "right":
        pdf.drawRightString(x, y, text)
    else:
        pdf.drawString(x, y, text)


def design_info(allocation: str, stage2_allocation: str) -> tuple[str, str]:
    if allocation == "one_stage":
        return "one", "one-stage"
    if allocation == "two_stage" and stage2_allocation == "highest_utility":
        return "high", "2-stage highest utility"
    if allocation == "two_stage" and stage2_allocation == "top2_randomized":
        return "top2", "2-stage top-2 randomized"
    raise ValueError(f"Unexpected design: {allocation} / {stage2_allocation}")


def read_new_design() -> dict[tuple[int, str, float], dict]:
    data = pd.read_csv(SOURCE_FILE)
    required = {
        "Scenario", "Nmax", "Cycle_Max", "Utility_Type", "Futility_Cutoff", "Allocation",
        "Stage2_Allocation", "Toxicity_IPDE_Alpha", "Efficacy_IPDE_Alpha", "Dose",
        "True_DLT_rate", "True_Efficacy_rate", "OBD_Selection_pct", "Pts_Treated",
        "No_OBD_Selection_pct", "n_valid", "ntrial_from_files",
    }
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError(f"New-design source is missing: {missing}")

    data = data.loc[
        data["Scenario"].isin(SCENARIOS)
        & (data["Nmax"] == 30)
        & (data["Cycle_Max"] == 2)
        & (data["Utility_Type"] == 3)
        & np.isclose(data["Futility_Cutoff"], 0.85)
        & data["Toxicity_IPDE_Alpha"].isin(ALPHAS)
        & np.isclose(data["Toxicity_IPDE_Alpha"], data["Efficacy_IPDE_Alpha"])
        & (data["n_valid"] == 1000)
        & (data["ntrial_from_files"] == 1000)
    ].copy()
    allowed = {
        ("one_stage", "highest_utility"),
        ("two_stage", "highest_utility"),
        ("two_stage", "top2_randomized"),
    }
    data = data.loc[
        data.apply(lambda row: (row["Allocation"], row["Stage2_Allocation"]) in allowed, axis=1)
    ].copy()

    expected = {(scenario, key, alpha) for scenario in SCENARIOS for key, _ in DESIGNS for alpha in ALPHAS}
    records: dict[tuple[int, str, float], dict] = {}
    for (scenario, alpha, allocation, stage2), rows in data.groupby(
        ["Scenario", "Toxicity_IPDE_Alpha", "Allocation", "Stage2_Allocation"], sort=False
    ):
        key, short = design_info(allocation, stage2)
        rows = rows.sort_values("Dose")
        if rows["Dose"].tolist() != [1, 2, 3, 4, 5]:
            raise ValueError(f"Expected five dose rows for scenario {scenario}, alpha {alpha}, {short}.")
        none_values = rows["No_OBD_Selection_pct"].unique()
        if len(none_values) != 1:
            raise ValueError(f"No OBD selection is inconsistent for scenario {scenario}, alpha {alpha}, {short}.")
        records[(int(scenario), key, float(alpha))] = {
            "short": short,
            "tox": rows["True_DLT_rate"].astype(float).tolist(),
            "eff": rows["True_Efficacy_rate"].astype(float).tolist(),
            "selection": rows["OBD_Selection_pct"].astype(float).tolist(),
            "treated": rows["Pts_Treated"].astype(float).tolist(),
            "none": float(none_values[0]),
        }
    if set(records) != expected:
        missing = sorted(expected.difference(records))
        extra = sorted(set(records).difference(expected))
        raise ValueError(f"New-design grid is incomplete. Missing={missing[:5]}, extra={extra[:5]}")

    for scenario in SCENARIOS:
        reference = records[(scenario, "one", 0.0)]
        for key, _ in DESIGNS:
            for alpha in ALPHAS:
                record = records[(scenario, key, alpha)]
                if not np.allclose(record["tox"], reference["tox"]) or not np.allclose(record["eff"], reference["eff"]):
                    raise ValueError(f"Truth rows disagree for scenario {scenario}.")
    return records


def parse_operating_characteristics_pdf(path: Path) -> dict[int, dict]:
    with pdfplumber.open(path) as document:
        text = "\n".join(page.extract_text() or "" for page in document.pages)
    blocks = list(re.finditer(r"^Scenario\s+(\d+)\s*$", text, flags=re.MULTILINE))
    records: dict[int, dict] = {}
    for index, match in enumerate(blocks):
        block_end = blocks[index + 1].start() if index + 1 < len(blocks) else len(text)
        block = text[match.end():block_end]

        def values_for(labels: tuple[str, ...], count: int) -> list[float]:
            for label in labels:
                found = re.search(re.escape(label) + r"\s+([^\n]+)", block)
                if found:
                    values = [float(value) for value in re.findall(r"-?(?:\d+\.?\d*|\.\d+)", found.group(1))]
                    if len(values) >= count:
                        return values[-count:]
            raise ValueError(f"Could not parse {labels} in {path.name}, scenario {match.group(1)}.")

        selection = values_for(("Select %",), 6)
        records[int(match.group(1))] = {
            "treated": values_for(("No. Pts treated", "# Pts treated"), 5),
            "selection": selection[:5],
            "none": selection[5],
        }
    if not records:
        raise ValueError(f"No comparator scenarios found in {path.name}.")
    return records


def read_uboin() -> dict[int, dict]:
    records: dict[int, dict] = {}
    for filename, offset, expected_count in UBOIN_GROUPS:
        parsed = parse_operating_characteristics_pdf(RAW_DATA / filename)
        if len(parsed) != expected_count or set(parsed) != set(range(1, expected_count + 1)):
            raise ValueError(f"Unexpected U-BOIN scenarios in {filename}.")
        records.update({scenario + offset: row for scenario, row in parsed.items()})
    scenario_38 = parse_operating_characteristics_pdf(UBOIN_SCENARIO_38)
    if set(scenario_38) != {1}:
        raise ValueError("U-BOIN scenario 38 source must contain exactly local scenario 1.")
    records[38] = scenario_38[1]
    if set(records) != set(SCENARIOS):
        raise ValueError("U-BOIN sources are incomplete for scenarios 1-38.")
    return records


def read_boin12() -> dict[int, dict]:
    records: dict[int, dict] = {}
    current: int | None = None
    with BOIN12_CSV.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.reader(stream):
            first = row[0].strip() if row else ""
            if first.startswith("Scenario "):
                current = int(first.split()[1])
                records[current] = {}
            elif current is not None and first == "No. Pts treated":
                records[current]["treated"] = [float(value) for value in row[1:6]]
            elif current is not None and first == "Select %":
                records[current]["selection"] = [float(value) for value in row[1:6]]
                records[current]["none"] = float(row[6])
    if set(records) != set(range(1, 38)) or any(set(row) != {"treated", "selection", "none"} for row in records.values()):
        raise ValueError("BOIN12 CSV is incomplete for scenarios 1-37.")
    scenario_38 = parse_operating_characteristics_pdf(BOIN12_SCENARIO_38)
    if set(scenario_38) != {1}:
        raise ValueError("BOIN12 scenario 38 source must contain exactly local scenario 1.")
    records[38] = scenario_38[1]
    return records


def draw_row(
    pdf: canvas.Canvas,
    y_top: float,
    height: float,
    label: str,
    values: list[str],
    *,
    fill=None,
    bold: bool = False,
) -> float:
    if fill is not None:
        pdf.setFillColor(fill)
        pdf.rect(LEFT, y_top - height, TABLE_WIDTH, height, fill=1, stroke=0)
    baseline = y_top - height + (height - 7.4) / 2 + 1.1
    draw_text(pdf, label, LEFT + 5, baseline, size=7.4, bold=bold)
    for index, value in enumerate(values):
        draw_text(pdf, value, LEFT + LABEL_WIDTH + VALUE_WIDTH * (index + 0.5), baseline, size=7.4, bold=bold, align="center")
    pdf.setStrokeColor(GRID)
    pdf.setLineWidth(0.35)
    pdf.line(LEFT, y_top - height, RIGHT, y_top - height)
    return y_top - height


def draw_truth(pdf: canvas.Canvas, y: float, scenario: int, record: dict) -> float:
    y = draw_row(pdf, y, 15, f"Scenario {scenario}", [""] * 6, fill=SECTION, bold=True)
    y = draw_row(pdf, y, 10.5, "DLT rate", [fmt(value, 2) for value in record["tox"]] + [""])
    return draw_row(pdf, y, 10.5, "Efficacy rate", [fmt(value, 2) for value in record["eff"]] + [""])


def draw_page_header(pdf: canvas.Canvas, alpha: float) -> float:
    draw_text(pdf, "Phase I/II Operating Characteristics", PAGE_WIDTH / 2, 576, size=20, bold=True, align="center")
    draw_text(
        pdf,
        f"New non-TITE alpha = {fmt(alpha)} versus U-BOIN and BOIN12 | scenarios 1-38",
        PAGE_WIDTH / 2,
        556,
        size=10.2,
        color=MUTED,
        align="center",
    )
    pdf.setStrokeColor(NAVY)
    pdf.setLineWidth(0.8)
    pdf.line(LEFT, 545, RIGHT, 545)
    for index, header in enumerate(["Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "No OBD %"]):
        draw_text(pdf, header, LEFT + LABEL_WIDTH + VALUE_WIDTH * (index + 0.5), 531, size=8.8, bold=True, align="center")
    return 518


def record_lookup(records: dict[tuple[int, str, float], dict], scenario: int, key: str, alpha: float) -> dict:
    return records[(scenario, key, alpha)]


def draw_method_scenario(
    pdf: canvas.Canvas,
    y: float,
    scenario: int,
    alpha: float,
    new_design: dict[tuple[int, str, float], dict],
    uboin: dict[int, dict],
    boin12: dict[int, dict],
) -> float:
    reference = record_lookup(new_design, scenario, "one", alpha)
    y = draw_truth(pdf, y, scenario, reference)
    entries: list[tuple[str, dict, object, bool]] = []
    for key, _ in DESIGNS:
        row = record_lookup(new_design, scenario, key, alpha)
        entries.append((f"New non-TITE | {row['short']}", row, NEW_FILL, True))
    entries.extend([
        ("U-BOIN", uboin[scenario], UBOIN_FILL, False),
        ("BOIN12", boin12[scenario], BOIN_FILL, False),
    ])
    y = draw_row(pdf, y, 12, "OBD selection (%)", [""] * 6, fill=SECTION, bold=True)
    for label, row, fill, bold in entries:
        y = draw_row(pdf, y, 10.5, label, [fmt(value) for value in row["selection"]] + [fmt(row["none"])], fill=fill, bold=bold)
    y = draw_row(pdf, y, 12, "Mean administrations by dose", [""] * 6, fill=SECTION, bold=True)
    for label, row, fill, bold in entries:
        y = draw_row(pdf, y, 10.5, label, [fmt(value) for value in row["treated"]] + [""], fill=fill, bold=bold)
    return y


def draw_footer(pdf: canvas.Canvas, alpha: float, page_number: int, page_count: int) -> None:
    draw_text(
        pdf,
        f"Sources: Scenarios 1-38 new-design alpha = {fmt(alpha)}, U-BOIN futility 0.85 PDFs, and BOIN12 operating-characteristics CSV/PDF (1,000 simulations each).",
        LEFT,
        27,
        size=6.2,
        color=MUTED,
    )
    draw_text(
        pdf,
        "Outcome column is source-specific: New non-TITE = No OBD; U-BOIN = Stop; BOIN12 = No selection.",
        LEFT,
        16,
        size=6.1,
        color=MUTED,
    )
    draw_text(pdf, f"Page {page_number} of {page_count}", RIGHT, 16, size=6.2, color=MUTED, align="right")


def build_pdf(
    alpha: float,
    new_design: dict[tuple[int, str, float], dict],
    uboin: dict[int, dict],
    boin12: dict[int, dict],
) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUTPUT_DIR / f"nontite_newdesign_alpha_{alpha_tag(alpha)}_scenarios_1_to_38_tables.pdf"
    groups = [SCENARIOS[index:index + 2] for index in range(0, len(SCENARIOS), 2)]
    pdf = canvas.Canvas(str(output), pagesize=landscape(letter), pageCompression=1)
    pdf.setTitle(f"New non-TITE alpha {fmt(alpha)} versus U-BOIN and BOIN12")
    for page_number, scenarios in enumerate(groups, start=1):
        y = draw_page_header(pdf, alpha)
        for index, scenario in enumerate(scenarios):
            y = draw_method_scenario(pdf, y, scenario, alpha, new_design, uboin, boin12)
            if index < len(scenarios) - 1:
                y -= 8
        draw_footer(pdf, alpha, page_number, len(groups))
        pdf.showPage()
    pdf.save()
    return output


def main() -> None:
    new_design = read_new_design()
    uboin = read_uboin()
    boin12 = read_boin12()
    for alpha in ALPHAS:
        print(build_pdf(alpha, new_design, uboin, boin12))


if __name__ == "__main__":
    main()
