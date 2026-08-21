from __future__ import annotations

import csv
import html
import re
from collections import defaultdict
from pathlib import Path

from reportlab.lib.colors import HexColor, black
from reportlab.lib.pagesizes import landscape, letter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
PRESENTATION = ROOT / "Presentation 8-17-2026"
OUTPUT_DIR = PRESENTATION / "Table and Plots" / "New Design"
BASE_NEW_DESIGN_CSV = PRESENTATION / "Raw Data" / (
    "AIDE_phase_I_II_modelsadditive_newdesign_shared_N30_ncycle2_"
    "rp0p15x0p85_rate56d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_"
    "ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv"
)
ONE_STAGE_REPLACEMENT_CSV = PRESENTATION / "Raw Data" / (
    "AIDE_phase_I_II_modelsadditive_shared_N30_ncycle2_"
    "rp0p15x0p85_rate56d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_"
    "ap0p15x0p85_IDX_0001_to_1000_newdesign_dose_summary.csv"
)
ORACLE_ONE_STAGE_REPLACEMENT_CSV = PRESENTATION / "Raw Data" / (
    "AIDE_phase_N30_ncycle1_rp0p15x0p85_rate56d_ep0p5x0p5_"
    "cp0p15x0p85_eth0p2_fut0p85_ap0p15x0p85_IDX_0001_to_1000_"
    "newdesign_continuous_dose_summary.csv"
)
BOIN12_CSV = PRESENTATION / "Raw Data" / "BOIN12_v1.4.2.0_Operating Characteristics_2026-08-09 170536.992978_fut0.85.csv"
EFFTOX_MAIN = PRESENTATION / "Raw Data" / "EffTox N = 30.html"
EFFTOX_38 = PRESENTATION / "Raw Data" / "Efftox N = 30 Sce 38.html"

ALPHA_PDF = OUTPUT_DIR / "nontite_newdesign_alpha_0_0p3_0p6_0p9_scenarios_1_16_20_24_27_38_tables.pdf"
METHOD_PDF = OUTPUT_DIR / "nontite_alpha0_newdesign_vs_uboin_boin12_scenarios_1_16_20_24_27_38_tables.pdf"

SCENARIOS = [1, 16, 20, 24, 27, 38]
DOSES = [1, 2, 3, 4, 5]
ALPHAS = [0.0, 0.3, 0.6, 0.9]
DESIGN_ORDER = ["one", "high", "top2"]

# Retained oracle two-stage rows from the existing table.  The requested
# update replaces only the oracle one-stage rows.
ORACLE_TWO_STAGE_ROWS = {
    1: {
        "high": {"selection": [17.6, 5.7, 0.9, 0.4, 0.0], "none": 75.4, "treated": [11.2, 4.2, 2.3, 0.8, 0.2]},
        "top2": {"selection": [16.6, 5.7, 1.5, 0.5, 0.1], "none": 75.6, "treated": [11.1, 4.2, 2.2, 0.9, 0.2]},
    },
    16: {
        "high": {"selection": [32.5, 22.7, 16.3, 17.8, 10.6], "none": 0.1, "treated": [7.5, 6.9, 5.8, 5.6, 4.2]},
        "top2": {"selection": [31.5, 22.3, 17.1, 16.1, 13.0], "none": 0.0, "treated": [6.6, 6.3, 5.7, 6.1, 5.3]},
    },
    20: {
        "high": {"selection": [13.5, 51.2, 24.4, 6.6, 0.7], "none": 3.6, "treated": [5.4, 11.6, 8.0, 3.5, 1.0]},
        "top2": {"selection": [9.2, 50.0, 25.3, 9.8, 2.0], "none": 3.7, "treated": [5.7, 9.8, 8.0, 4.4, 1.5]},
    },
    24: {
        "high": {"selection": [2.0, 1.1, 6.8, 56.7, 33.4], "none": 0.0, "treated": [3.4, 3.5, 4.4, 10.9, 7.7]},
        "top2": {"selection": [0.3, 0.3, 2.9, 56.3, 40.2], "none": 0.0, "treated": [3.6, 3.8, 4.9, 9.5, 8.4]},
    },
    27: {
        "high": {"selection": [16.4, 25.1, 43.9, 6.8, 0.6], "none": 7.2, "treated": [7.8, 8.2, 9.4, 3.0, 0.7]},
        "top2": {"selection": [12.0, 24.5, 45.2, 10.6, 1.1], "none": 6.6, "treated": [7.6, 8.4, 8.5, 3.7, 0.9]},
    },
    38: {
        "high": {"selection": [19.2, 10.3, 8.6, 7.4, 52.8], "none": 1.7, "treated": [4.9, 4.8, 5.5, 5.6, 9.1]},
        "top2": {"selection": [12.0, 5.6, 7.4, 7.8, 65.8], "none": 1.4, "treated": [4.8, 4.9, 5.4, 6.8, 7.8]},
    },
}

PAGE_WIDTH, _ = landscape(letter)
LEFT = 36
RIGHT = PAGE_WIDTH - 36
TABLE_WIDTH = RIGHT - LEFT
LABEL_WIDTH = 255
VALUE_WIDTH = (TABLE_WIDTH - LABEL_WIDTH) / 6
SUMMARY_LABEL_WIDTH = 240
SUMMARY_VALUE_WIDTH = (TABLE_WIDTH - SUMMARY_LABEL_WIDTH) / 7
NAVY = HexColor("#173F5F")
SECTION = HexColor("#E9ECEF")
NEW_FILL = HexColor("#E6F3F1")
ORACLE_FILL = HexColor("#F6F7F8")
UBOIN_FILL = HexColor("#E5F1FB")
BOIN_FILL = HexColor("#EEF2FF")
EFFTOX_FILL = HexColor("#FFF4D6")
GRID = HexColor("#C9CDD1")
MUTED = HexColor("#4D5965")


def num(value: str | float | int | None) -> float | None:
    if value in (None, "", "NA", "-"):
        return None
    return float(value)


def fmt(value: float | None, digits: int = 1) -> str:
    if value is None:
        return "-"
    text = f"{float(value):.{digits}f}"
    return text.rstrip("0").rstrip(".")


def design_info(row: dict[str, str]) -> tuple[str, str]:
    if row["Allocation"] == "one_stage":
        return "one", "one-stage"
    if row["Stage2_Allocation"] == "highest_utility":
        return "high", "2-stage highest utility"
    if row["Stage2_Allocation"] == "top2_randomized":
        return "top2", "2-stage top-2 randomized"
    raise ValueError(f"Unexpected design row: {row['Allocation']} / {row['Stage2_Allocation']}")


def read_design_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        raw = list(csv.DictReader(handle))
    grouped: dict[tuple[int, str, float], list[dict[str, str]]] = defaultdict(list)
    for row in raw:
        scenario = int(row["Scenario"])
        tox_alpha = num(row["Toxicity_IPDE_Alpha"])
        eff_alpha = num(row["Efficacy_IPDE_Alpha"])
        if scenario not in SCENARIOS or tox_alpha not in ALPHAS or tox_alpha != eff_alpha:
            continue
        design_key, _ = design_info(row)
        grouped[(scenario, design_key, tox_alpha)].append(row)
    if any(len(rows) != 5 for rows in grouped.values()):
        raise RuntimeError(f"{path.name} does not contain five dose rows for every selected setting.")

    records: list[dict] = []
    for (scenario, key, alpha), dose_rows in grouped.items():
        dose_rows.sort(key=lambda item: int(item["Dose"]))
        first = dose_rows[0]
        _, short = design_info(first)
        records.append({
            "scenario": scenario,
            "key": key,
            "short": short,
            "alpha": alpha,
            "true_obd": int(first["True_OBD"]),
            "tox": [num(row["True_DLT_rate"]) for row in dose_rows],
            "eff": [num(row["True_Efficacy_rate"]) for row in dose_rows],
            "selection": [num(row["OBD_Selection_pct"]) for row in dose_rows],
            "treated": [num(row["Pts_Treated"]) for row in dose_rows],
            "none": num(first["No_OBD_Selection_pct"]),
            "unique_patients": num(first["Total_Unique_Patients"]),
        })
    return records


def read_new_design() -> list[dict]:
    records = read_design_csv(BASE_NEW_DESIGN_CSV)
    if len(records) != 72:
        raise RuntimeError("The base New Design CSV does not contain the complete 6 x 3 x 4 alpha grid.")

    replacements = {
        (record["scenario"], record["key"], record["alpha"]): record
        for record in read_design_csv(ONE_STAGE_REPLACEMENT_CSV)
    }
    required = {(scenario, "one", alpha) for scenario in SCENARIOS for alpha in ALPHAS}
    if not required.issubset(replacements):
        raise RuntimeError("The requested one-stage alpha-grid replacement rows are incomplete.")

    return [
        replacements[(record["scenario"], record["key"], record["alpha"])]
        if (record["scenario"], record["key"], record["alpha"]) in required
        else record
        for record in records
    ]


def read_oracle_design() -> list[dict]:
    replacements = read_design_csv(ORACLE_ONE_STAGE_REPLACEMENT_CSV)
    required = {(scenario, "one", 0.0) for scenario in SCENARIOS}
    observed = {(record["scenario"], record["key"], record["alpha"]) for record in replacements}
    if observed != required:
        raise RuntimeError("The requested cycle-1 oracle one-stage replacement rows are incomplete.")

    records = list(replacements)
    short_labels = {"high": "2-stage highest utility", "top2": "2-stage top-2 randomized"}
    for scenario, rows in ORACLE_TWO_STAGE_ROWS.items():
        for key, row in rows.items():
            records.append({
                "scenario": scenario,
                "key": key,
                "short": short_labels[key],
                "alpha": 0.0,
                "true_obd": 0,
                "tox": [],
                "eff": [],
                "selection": row["selection"],
                "treated": row["treated"],
                "none": row["none"],
                "unique_patients": None,
            })
    return records


def clean_html(value: str) -> str:
    value = re.sub(r"<[^>]+>", "", value)
    return html.unescape(value).replace("\xa0", " ").strip()


def extract_efftox_file(path: Path, scenario_map: dict[int, int] | None = None) -> dict[int, dict]:
    raw = path.read_text(encoding="utf-8")
    pattern = re.compile(
        r'<tr><td colspan\s*=\s*"8"><b>\s*(\d+)\s*</b>(?:</td>)?</tr>(.*?)(?=<tr><td colspan\s*=\s*"8"|</table>)',
        re.S,
    )
    result: dict[int, dict] = {}
    for match in pattern.finditer(raw):
        scenario = int(match.group(1))
        scenario = (scenario_map or {}).get(scenario, scenario)
        block = match.group(2)
        selection_match = re.search(r"%\s*selected.*?(?=<tr>|$)", block, re.S | re.I)
        treated_match = re.search(r"#\s*Patients\s*Treated.*?(?=<tr>|$)", block, re.S | re.I)
        truth_match = re.search(r"True\s*pT,\s*pE.*?(?=<tr>|$)", block, re.S | re.I)
        if not selection_match or not treated_match or not truth_match:
            raise ValueError(f"EffTox scenario {scenario} could not be parsed from {path.name}.")
        selection_cells = [clean_html(item) for item in re.findall(r"<td[^>]*>(.*?)</td>", selection_match.group(0), re.S)]
        treated_cells = [clean_html(item) for item in re.findall(r"<td[^>]*>(.*?)</td>", treated_match.group(0), re.S)]
        truth_cells = [clean_html(item) for item in re.findall(r"<td[^>]*>(.*?)</td>", truth_match.group(0), re.S)]
        selection = [num(item) for item in selection_cells[-6:-1]]
        none = num(selection_cells[-1])
        treated = [num(item) for item in treated_cells[-6:-1]]
        truth_pairs = truth_cells[-6:-1]
        tox = [float(item.split(",")[0]) for item in truth_pairs]
        eff = [float(item.split(",")[1]) for item in truth_pairs]
        if len(selection) != 5 or len(treated) != 5:
            raise ValueError(f"EffTox scenario {scenario} does not contain five dose-level values.")
        result[scenario] = {"selection": selection, "treated": treated, "none": none, "tox": tox, "eff": eff}
    return result


def read_efftox() -> dict[int, dict]:
    data = extract_efftox_file(EFFTOX_MAIN)
    data.update(extract_efftox_file(EFFTOX_38, {1: 38}))
    missing = set(SCENARIOS).difference(data)
    if missing:
        raise ValueError(f"EffTox results are missing scenarios: {sorted(missing)}")
    return {scenario: data[scenario] for scenario in SCENARIOS}


UBOIN = {
    1: {"treated": [9, 5.4, 4.7, 3.9, 3.7], "selection": [15.6, 17.4, 20.6, 18.7, 13.7], "none": 14},
    16: {"treated": [6.1, 6.2, 6.4, 6.4, 5.4], "selection": [28.1, 26.8, 21.1, 16.3, 7.7], "none": 0},
    20: {"treated": [5.9, 10.8, 7, 4, 2.4], "selection": [8.5, 47.5, 18.6, 12.4, 13], "none": 0},
    24: {"treated": [3.5, 4, 5.1, 9.8, 7.9], "selection": [1.5, 3.2, 8.2, 60.1, 27], "none": 0},
    27: {"treated": [6.9, 8.1, 8.8, 4.3, 1.8], "selection": [8.7, 20.9, 38.9, 16.5, 14.5], "none": 0.5},
    38: {"treated": [5.2, 5.5, 5.6, 6.3, 7.8], "selection": [6.2, 8.9, 13, 13.3, 58.6], "none": 0},
}


def read_boin12() -> dict[int, dict]:
    parsed: dict[int, dict] = {}
    current: int | None = None
    with BOIN12_CSV.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            first = row[0].strip() if row else ""
            if first.startswith("Scenario "):
                current = int(first.split()[1])
            elif current is not None and first == "No. Pts treated":
                parsed.setdefault(current, {})["treated"] = [num(value) for value in row[1:6]]
            elif current is not None and first == "Select %":
                parsed.setdefault(current, {})["selection"] = [num(value) for value in row[1:6]]
                parsed[current]["none"] = num(row[6])
    parsed[38] = {
        "treated": [4.2, 4.7, 5.2, 6.1, 9.6],
        "selection": [13.2, 13.8, 16.5, 13.6, 41.5],
        "none": 1.4,
    }
    missing = set(SCENARIOS).difference(parsed)
    if missing:
        raise ValueError(f"BOIN12 results are missing scenarios: {sorted(missing)}")
    return {scenario: parsed[scenario] for scenario in SCENARIOS}


def draw_text(pdf, text: str, x: float, y: float, size: float = 8, bold: bool = False, color=black, align: str = "left") -> None:
    pdf.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    pdf.setFillColor(color)
    if align == "center":
        pdf.drawCentredString(x, y, text)
    elif align == "right":
        pdf.drawRightString(x, y, text)
    else:
        pdf.drawString(x, y, text)


def draw_row(pdf, y_top: float, height: float, label: str, values: list[str], fill=None, bold: bool = False) -> float:
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


def draw_truth(pdf, y: float, record: dict) -> float:
    y = draw_row(pdf, y, 15, f"Scenario {record['scenario']}", [""] * 6, fill=SECTION, bold=True)
    y = draw_row(pdf, y, 10.5, "DLT rate", [fmt(value, 2) for value in record["tox"]] + [""])
    return draw_row(pdf, y, 10.5, "Efficacy rate", [fmt(value, 2) for value in record["eff"]] + [""])


def draw_page_header(pdf, subtitle: str) -> float:
    draw_text(pdf, "Phase I/II Operating Characteristics", PAGE_WIDTH / 2, 576, size=20, bold=True, align="center")
    draw_text(pdf, subtitle, PAGE_WIDTH / 2, 556, size=10.2, color=MUTED, align="center")
    pdf.setStrokeColor(NAVY)
    pdf.setLineWidth(0.8)
    pdf.line(LEFT, 545, RIGHT, 545)
    for index, header in enumerate(["Dose 1", "Dose 2", "Dose 3", "Dose 4", "Dose 5", "No OBD %"]):
        draw_text(pdf, header, LEFT + LABEL_WIDTH + VALUE_WIDTH * (index + 0.5), 531, size=8.8, bold=True, align="center")
    return 518


def draw_summary_header(pdf) -> float:
    draw_text(pdf, "Mean Unique Patients by Method", PAGE_WIDTH / 2, 576, size=20, bold=True, align="center")
    draw_text(pdf, "Mean number of unique patients per simulated trial", PAGE_WIDTH / 2, 556, size=10.2, color=MUTED, align="center")
    pdf.setStrokeColor(NAVY)
    pdf.setLineWidth(0.8)
    pdf.line(LEFT, 545, RIGHT, 545)
    headers = ["Scenario 1", "16", "20", "24", "27", "38", "Mean"]
    for index, header in enumerate(headers):
        draw_text(
            pdf,
            header,
            LEFT + SUMMARY_LABEL_WIDTH + SUMMARY_VALUE_WIDTH * (index + 0.5),
            531,
            size=8.8,
            bold=True,
            align="center",
        )
    return 518


def record_lookup(records: list[dict], scenario: int, key: str, alpha: float) -> dict:
    return next(record for record in records if record["scenario"] == scenario and record["key"] == key and record["alpha"] == alpha)


def draw_alpha_scenario(
    pdf, y: float, scenario: int, records: list[dict], oracle_records: list[dict]
) -> float:
    reference = record_lookup(records, scenario, "one", 0.0)
    y = draw_truth(pdf, y, reference)
    y = draw_row(pdf, y, 12, "OBD selection (%)", [""] * 6, fill=SECTION, bold=True)
    for key in DESIGN_ORDER:
        oracle = record_lookup(oracle_records, scenario, key, 0.0)
        y = draw_row(
            pdf,
            y,
            10.5,
            f"Oracle | {oracle['short']} | cycle 1, alpha=0",
            [fmt(value) for value in oracle["selection"]] + [fmt(oracle["none"])],
            fill=ORACLE_FILL,
        )
        for alpha in ALPHAS:
            row = record_lookup(records, scenario, key, alpha)
            y = draw_row(
                pdf, y, 10.5, f"New non-TITE | {row['short']} | alpha={fmt(alpha)}",
                [fmt(value) for value in row["selection"]] + [fmt(row["none"])], fill=NEW_FILL,
            )
    y = draw_row(pdf, y, 12, "Mean administrations by dose", [""] * 6, fill=SECTION, bold=True)
    for key in DESIGN_ORDER:
        oracle = record_lookup(oracle_records, scenario, key, 0.0)
        y = draw_row(
            pdf,
            y,
            10.5,
            f"Oracle | {oracle['short']} | cycle 1, alpha=0",
            [fmt(value) for value in oracle["treated"]] + [""],
            fill=ORACLE_FILL,
        )
        for alpha in ALPHAS:
            row = record_lookup(records, scenario, key, alpha)
            y = draw_row(pdf, y, 10.5, f"New non-TITE | {row['short']} | alpha={fmt(alpha)}", [fmt(value) for value in row["treated"]] + [""], fill=NEW_FILL)
    return y


def draw_method_scenario(pdf, y: float, scenario: int, records: list[dict], efftox: dict[int, dict], boin12: dict[int, dict]) -> float:
    reference = record_lookup(records, scenario, "one", 0.0)
    y = draw_truth(pdf, y, reference)
    entries = []
    for key in DESIGN_ORDER:
        row = record_lookup(records, scenario, key, 0.0)
        entries.append((f"New non-TITE | {row['short']}", row, NEW_FILL, True))
    entries.extend([
        ("U-BOIN", UBOIN[scenario], UBOIN_FILL, False),
        ("BOIN12", boin12[scenario], BOIN_FILL, False),
        ("EffTox", efftox[scenario], EFFTOX_FILL, False),
    ])
    y = draw_row(pdf, y, 12, "OBD selection (%)", [""] * 6, fill=SECTION, bold=True)
    for label, row, fill, bold in entries:
        y = draw_row(pdf, y, 10.5, label, [fmt(value) for value in row["selection"]] + [fmt(row["none"])], fill=fill, bold=bold)
    y = draw_row(pdf, y, 12, "Mean administrations by dose", [""] * 6, fill=SECTION, bold=True)
    for label, row, fill, bold in entries:
        y = draw_row(pdf, y, 10.5, label, [fmt(value) for value in row["treated"]] + [""], fill=fill, bold=bold)
    return y


def draw_summary_row(pdf, y_top: float, height: float, label: str, values: list[float], fill=None, bold: bool = False) -> float:
    if fill is not None:
        pdf.setFillColor(fill)
        pdf.rect(LEFT, y_top - height, TABLE_WIDTH, height, fill=1, stroke=0)
    baseline = y_top - height + (height - 7.8) / 2 + 1.2
    draw_text(pdf, label, LEFT + 5, baseline, size=7.8, bold=bold)
    for index, value in enumerate(values):
        draw_text(
            pdf,
            fmt(value),
            LEFT + SUMMARY_LABEL_WIDTH + SUMMARY_VALUE_WIDTH * (index + 0.5),
            baseline,
            size=7.8,
            bold=bold,
            align="center",
        )
    pdf.setStrokeColor(GRID)
    pdf.setLineWidth(0.35)
    pdf.line(LEFT, y_top - height, RIGHT, y_top - height)
    return y_top - height


def draw_unique_patients_summary(pdf, records: list[dict], efftox: dict[int, dict], boin12: dict[int, dict]) -> None:
    y = draw_summary_header(pdf)
    entries = []
    for key in DESIGN_ORDER:
        values = [record_lookup(records, scenario, key, 0.0)["unique_patients"] for scenario in SCENARIOS]
        label = f"New non-TITE | {record_lookup(records, SCENARIOS[0], key, 0.0)['short']}"
        entries.append((label, values, NEW_FILL, True))
    entries.extend([
        ("U-BOIN", [sum(UBOIN[scenario]["treated"]) for scenario in SCENARIOS], UBOIN_FILL, False),
        ("BOIN12", [sum(boin12[scenario]["treated"]) for scenario in SCENARIOS], BOIN_FILL, False),
        ("EffTox", [sum(efftox[scenario]["treated"]) for scenario in SCENARIOS], EFFTOX_FILL, False),
    ])
    for label, values, fill, bold in entries:
        y = draw_summary_row(pdf, y, 18, label, values + [sum(values) / len(values)], fill=fill, bold=bold)
    draw_footer(
        pdf,
        "New non-TITE values are reported Total_Unique_Patients. U-BOIN, BOIN12, and EffTox have no IPDE recycling, so unique patients equal the sum of dose-level means.",
        "The updated New non-TITE | one-stage alpha = 0 rows use the specified modelsadditive_shared ... newdesign_dose_summary.csv source.",
        len([SCENARIOS[index:index + 2] for index in range(0, len(SCENARIOS), 2)]) + 1,
        len([SCENARIOS[index:index + 2] for index in range(0, len(SCENARIOS), 2)]) + 1,
    )


def draw_footer(pdf, source_note: str, outcome_note: str, page_number: int, page_count: int) -> None:
    draw_text(pdf, source_note, LEFT, 27, size=6.2, color=MUTED)
    draw_text(pdf, outcome_note, LEFT, 16, size=6.1, color=MUTED)
    draw_text(pdf, f"Page {page_number} of {page_count}", RIGHT, 16, size=6.2, color=MUTED, align="right")


def build_alpha_pdf(records: list[dict], oracle_records: list[dict]) -> None:
    pdf = canvas.Canvas(str(ALPHA_PDF), pagesize=landscape(letter), pageCompression=1)
    for page_number, scenario in enumerate(SCENARIOS, start=1):
        y = draw_page_header(pdf, "New non-TITE alpha comparison with cycle-1 oracle (alpha = 0)")
        draw_alpha_scenario(pdf, y, scenario, records, oracle_records)
        draw_footer(
            pdf,
            "Sources: cycle-2 new-design alpha grid; cycle-1 continuous one-stage oracle; existing two-stage oracle rows (1,000 simulations per populated row).",
            "Oracle rows use Cycle_Max = 1 and alpha = 0; No OBD % is the source's no-OBD selection percentage.",
            page_number,
            len(SCENARIOS),
        )
        pdf.showPage()
    pdf.save()


def build_methods_pdf(records: list[dict], efftox: dict[int, dict], boin12: dict[int, dict]) -> None:
    pdf = canvas.Canvas(str(METHOD_PDF), pagesize=landscape(letter), pageCompression=1)
    groups = [SCENARIOS[index:index + 2] for index in range(0, len(SCENARIOS), 2)]
    page_count = len(groups) + 1
    for page_number, scenarios in enumerate(groups, start=1):
        y = draw_page_header(pdf, "New non-TITE alpha = 0 versus U-BOIN, BOIN12, and EffTox | updated one-stage rows")
        for index, scenario in enumerate(scenarios):
            y = draw_method_scenario(pdf, y, scenario, records, efftox, boin12)
            if index < len(scenarios) - 1:
                y -= 8
        draw_footer(
            pdf,
            "Sources: updated New non-TITE one-stage alpha = 0, existing two-stage rows, U-BOIN, BOIN12, and EffTox N = 30 (1,000 simulations each).",
            "No OBD % is source-specific: New non-TITE = No OBD; U-BOIN = Stop; BOIN12 = No selection; EffTox = None selected.",
            page_number,
            page_count,
        )
        pdf.showPage()
    draw_unique_patients_summary(pdf, records, efftox, boin12)
    pdf.showPage()
    pdf.save()


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    records = read_new_design()
    oracle_records = read_oracle_design()
    build_alpha_pdf(records, oracle_records)
    print(ALPHA_PDF)


if __name__ == "__main__":
    main()
