from __future__ import annotations

import csv
import importlib.util
from collections import defaultdict
from pathlib import Path

from reportlab.lib.colors import HexColor
from reportlab.lib.pagesizes import landscape, letter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
SOURCE_BUILDER = ROOT / "tmp" / "pdfs" / "build_nontite_alpha_and_methods_efftox_tables.py"
LATEST_CSV = ROOT / "Presentation 8-17-2026" / "Raw Data" / (
    "AIDE_phase_I_II_modelsadditive_shared_N30_ncycle2_"
    "rp0p15x0p85_rate56d_ep0p5x0p5_cp0p15x0p85_eth0p2_fut0p85_"
    "ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv"
)
OUTPUT = ROOT / "output" / "pdf" / (
    "nontite_alpha0_current_and_previous_nontite_vs_uboin_boin12_"
    "scenarios_1_16_20_24_27_38_tables.pdf"
)
SCENARIOS = [1, 16, 20, 24, 27, 38]

# Gray identifies retained results from the earlier table; teal identifies
# the requested continuous-enrollment rows from the updated non-TITE source.
PREVIOUS_FILL = HexColor("#F6F7F8")
CURRENT_FILL = HexColor("#E6F3F1")


def load_source_builder():
    spec = importlib.util.spec_from_file_location("source_builder", SOURCE_BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load the source table builder.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def as_number(value: str | float | int | None) -> float | None:
    if value in (None, "", "NA", "-"):
        return None
    return float(value)


def read_latest_records() -> list[dict]:
    with LATEST_CSV.open(newline="", encoding="utf-8-sig") as stream:
        raw = list(csv.DictReader(stream))

    grouped: dict[tuple[int, str], list[dict[str, str]]] = defaultdict(list)
    for row in raw:
        if int(row["Scenario"]) not in SCENARIOS:
            continue
        if row["Enrollment_Scheme"] != "continuous":
            continue
        if row["Allocation"] not in {"two_stage", "one_stage"}:
            continue
        if row["Stage2_Allocation"] != "highest_utility":
            continue
        if as_number(row["Toxicity_IPDE_Alpha"]) != 0 or as_number(row["Efficacy_IPDE_Alpha"]) != 0:
            continue
        grouped[(int(row["Scenario"]), row["Allocation"])].append(row)

    required = {(scenario, allocation) for scenario in SCENARIOS for allocation in ("two_stage", "one_stage")}
    if set(grouped) != required or any(len(rows) != 5 for rows in grouped.values()):
        raise RuntimeError("The 8/23 continuous alpha = 0 source does not contain the required 12 five-dose rows.")

    label_by_allocation = {
        "two_stage": "2-stage highest utility",
        "one_stage": "one-stage",
    }
    records: list[dict] = []
    for (scenario, allocation), dose_rows in grouped.items():
        dose_rows.sort(key=lambda row: int(row["Dose"]))
        first = dose_rows[0]
        records.append({
            "scenario": scenario,
            "key": "high" if allocation == "two_stage" else "one",
            "short": label_by_allocation[allocation],
            "true_obd": int(first["True_OBD"]),
            "tox": [as_number(row["True_DLT_rate"]) for row in dose_rows],
            "eff": [as_number(row["True_Efficacy_rate"]) for row in dose_rows],
            "selection": [as_number(row["OBD_Selection_pct"]) for row in dose_rows],
            "treated": [as_number(row["Pts_Treated"]) for row in dose_rows],
            "none": as_number(first["No_OBD_Selection_pct"]),
            "unique_patients": as_number(first["Total_Unique_Patients"]),
        })
    return records


def record_lookup(records: list[dict], scenario: int, key: str) -> dict:
    return next(record for record in records if record["scenario"] == scenario and record["key"] == key)


def entries_for(scenario: int, previous: list[dict], latest: list[dict], source) -> list[tuple[str, dict, object, bool]]:
    previous_one = source.record_lookup(previous, scenario, "one", 0.0)
    previous_high = source.record_lookup(previous, scenario, "high", 0.0)
    latest_high = record_lookup(latest, scenario, "high")
    latest_one = record_lookup(latest, scenario, "one")
    return [
        ("Previous non-TITE | one-stage", previous_one, PREVIOUS_FILL, True),
        ("Non-TITE | one-stage", latest_one, CURRENT_FILL, True),
        ("Previous non-TITE | 2-stage highest utility", previous_high, PREVIOUS_FILL, True),
        ("Non-TITE | 2-stage highest utility", latest_high, CURRENT_FILL, True),
        ("U-BOIN", source.UBOIN[scenario], source.UBOIN_FILL, False),
        ("BOIN12", source.read_boin12()[scenario], source.BOIN_FILL, False),
        ("EffTox", source.read_efftox()[scenario], source.EFFTOX_FILL, False),
    ]


def draw_scenario(pdf, y: float, scenario: int, previous: list[dict], latest: list[dict], source, comparators: dict) -> float:
    reference = record_lookup(latest, scenario, "one")
    y = source.draw_truth(pdf, y, reference)
    entries = [
        ("Previous non-TITE | one-stage", source.record_lookup(previous, scenario, "one", 0.0), PREVIOUS_FILL, True),
        ("Non-TITE | one-stage", record_lookup(latest, scenario, "one"), CURRENT_FILL, True),
        ("Previous non-TITE | 2-stage highest utility", source.record_lookup(previous, scenario, "high", 0.0), PREVIOUS_FILL, True),
        ("Non-TITE | 2-stage highest utility", record_lookup(latest, scenario, "high"), CURRENT_FILL, True),
        ("U-BOIN", source.UBOIN[scenario], source.UBOIN_FILL, False),
        ("BOIN12", comparators["boin12"][scenario], source.BOIN_FILL, False),
        ("EffTox", comparators["efftox"][scenario], source.EFFTOX_FILL, False),
    ]
    y = source.draw_row(pdf, y, 12, "OBD selection (%)", [""] * 6, fill=source.SECTION, bold=True)
    for label, row, fill, bold in entries:
        y = source.draw_row(
            pdf, y, 10.5, label,
            [source.fmt(value) for value in row["selection"]] + [source.fmt(row["none"])],
            fill=fill, bold=bold,
        )
    y = source.draw_row(pdf, y, 12, "Mean administrations by dose", [""] * 6, fill=source.SECTION, bold=True)
    for label, row, fill, bold in entries:
        y = source.draw_row(
            pdf, y, 10.5, label,
            [source.fmt(value) for value in row["treated"]] + [""],
            fill=fill, bold=bold,
        )
    return y


def draw_summary(pdf, previous: list[dict], latest: list[dict], source, comparators: dict, page_number: int, page_count: int) -> None:
    y = source.draw_summary_header(pdf)
    entries = [
        ("Previous non-TITE | one-stage", [source.record_lookup(previous, scenario, "one", 0.0)["unique_patients"] for scenario in SCENARIOS], PREVIOUS_FILL, True),
        ("Non-TITE | one-stage", [record_lookup(latest, scenario, "one")["unique_patients"] for scenario in SCENARIOS], CURRENT_FILL, True),
        ("Previous non-TITE | 2-stage highest utility", [source.record_lookup(previous, scenario, "high", 0.0)["unique_patients"] for scenario in SCENARIOS], PREVIOUS_FILL, True),
        ("Non-TITE | 2-stage highest utility", [record_lookup(latest, scenario, "high")["unique_patients"] for scenario in SCENARIOS], CURRENT_FILL, True),
        ("U-BOIN", [sum(source.UBOIN[scenario]["treated"]) for scenario in SCENARIOS], source.UBOIN_FILL, False),
        ("BOIN12", [sum(comparators["boin12"][scenario]["treated"]) for scenario in SCENARIOS], source.BOIN_FILL, False),
        ("EffTox", [sum(comparators["efftox"][scenario]["treated"]) for scenario in SCENARIOS], source.EFFTOX_FILL, False),
    ]
    for label, values, fill, bold in entries:
        y = source.draw_summary_row(pdf, y, 18, label, values + [sum(values) / len(values)], fill=fill, bold=bold)
    source.draw_footer(
        pdf,
        "Previous non-TITE rows are retained from the original table. Current Non-TITE rows use the requested modelsadditive_shared continuous-enrollment summary (1,000 simulations).",
        "Non-TITE unique patients are Total_Unique_Patients. Comparator patient counts are the sums of dose-level means; No OBD % remains source-specific.",
        page_number,
        page_count,
    )


def build() -> Path:
    source = load_source_builder()
    previous = source.read_new_design()
    latest = read_latest_records()
    comparators = {"boin12": source.read_boin12(), "efftox": source.read_efftox()}
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    pdf = canvas.Canvas(str(OUTPUT), pagesize=landscape(letter), pageCompression=1)
    scenario_pairs = [SCENARIOS[index:index + 2] for index in range(0, len(SCENARIOS), 2)]
    page_count = len(scenario_pairs) + 1
    for page_number, scenarios in enumerate(scenario_pairs, start=1):
        y = source.draw_page_header(
            pdf,
            "Current and previous non-TITE (alpha = 0) versus U-BOIN, BOIN12, and EffTox",
        )
        for index, scenario in enumerate(scenarios):
            y = draw_scenario(pdf, y, scenario, previous, latest, source, comparators)
            if index < len(scenarios) - 1:
                y -= 8
        source.draw_footer(
            pdf,
            "Sources: retained previous non-TITE rows; current Non-TITE rows use AIDE_phase_I_II_modelsadditive_shared_N30_ncycle2...dose_summary; U-BOIN, BOIN12, and EffTox N = 30.",
            "All non-TITE rows use IPDE alpha = 0. Current Non-TITE rows use continuous enrollment; no top-2 randomized rows are shown.",
            page_number,
            page_count,
        )
        pdf.showPage()
    draw_summary(pdf, previous, latest, source, comparators, page_count, page_count)
    pdf.showPage()
    pdf.save()
    return OUTPUT


if __name__ == "__main__":
    print(build())
