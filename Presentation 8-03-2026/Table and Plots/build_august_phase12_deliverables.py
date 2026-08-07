"""Build August Phase I/II prior-sensitivity tables, plots, and the refreshed N=30 comparison table.

The two August summaries intentionally vary different model components:
  * toxicity prior: CRM_Prior_a / CRM_Prior_b;
  * efficacy prior: Efficacy_Additive_Alpha_Prior_a / Efficacy_Additive_Alpha_Prior_b.

The N=30 all-method comparison keeps BOIN12, U-BOIN, and EffTox values from
the July raw-result files and refreshes only the AIDE rows from the August
toxicity-prior run at its baseline setting.
"""

from __future__ import annotations

import html
import re
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import pdfplumber
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
AUGUST_DIR = ROOT / "Presentation 8-03-2026"
JULY_DIR = ROOT / "Presentation 7-27-2026"
JULY_RAW_DIR = JULY_DIR / "Raw Result"
TRUTH_CSV = ROOT / "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv"

PAGE_SIZE = (18 * inch, 14 * inch)
PAGE_WIDTH, PAGE_HEIGHT = PAGE_SIZE
LEFT = 80
LABEL_WIDTH = 330
DOSE_WIDTH = 130
STOP_WIDTH = 156
TABLE_WIDTH = LABEL_WIDTH + 5 * DOSE_WIDTH + STOP_WIDTH
PRIORS = ((0.15, 0.85), (0.30, 0.70), (0.50, 0.50), (1.00, 1.00))
ALPHAS = (0.0, 0.3, 0.6, 0.9)


@dataclass(frozen=True)
class PriorRun:
    key: str
    summary_csv: Path
    output_dir: Path
    prior_a_column: str
    prior_b_column: str
    display_label: str
    table_prefix: str
    one_stage_table_prefix: str
    plot_prefix: str
    README_prior_description: str


TOXICITY_RUN = PriorRun(
    key="toxicity",
    summary_csv=AUGUST_DIR / (
        "AIDE_phase_I_II_modelsprevious_dose_additive_N30_"
        "rp0p15x0p85_0p3x0p7_0p5x0p5_1x1_ep0p5x0p5_"
        "cp0p15x0p85_ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv"
    ),
    output_dir=AUGUST_DIR / "Table and Plots" / "Toxicity Prior",
    prior_a_column="CRM_Prior_a",
    prior_b_column="CRM_Prior_b",
    display_label="Toxicity prior",
    table_prefix="phase12_toxicity_prior",
    one_stage_table_prefix="phase12_one_stage_toxicity_prior",
    plot_prefix="phase12_toxicity_prior",
    README_prior_description="the CRM toxicity prior",
)

EFFICACY_RUN = PriorRun(
    key="efficacy",
    summary_csv=AUGUST_DIR / (
        "AIDE_phase_I_II_modelsprevious_dose_additive_N30_"
        "rp0p15x0p85_ep0p5x0p5_cp0p15x0p85_ap0p15x0p85_"
        "0p3x0p7_0p5x0p5_1x1_IDX_0001_to_1000_dose_summary.csv"
    ),
    output_dir=AUGUST_DIR / "Table and Plots" / "Efficacy Prior",
    prior_a_column="Efficacy_Additive_Alpha_Prior_a",
    prior_b_column="Efficacy_Additive_Alpha_Prior_b",
    display_label="Efficacy additive alpha prior",
    table_prefix="phase12_efficacy_prior",
    one_stage_table_prefix="phase12_one_stage_efficacy_prior",
    plot_prefix="phase12_efficacy_prior",
    README_prior_description="the additive-efficacy-alpha prior",
)

N30_COMPARISON_DIR = AUGUST_DIR / "Table and Plots" / "Phase I-II"
RANDOM_MODEL_DIR = AUGUST_DIR / "Table and Plots" / "Random Model"
ALPHA_MODEL_DIR = AUGUST_DIR / "Table and Plots" / "Alpha Model"
COMBINED_IPDE_ALPHA_CSV = AUGUST_DIR / (
    "AIDE_phase_I_II_modelsprevious_dose_additive_N30_"
    "rp0p15x0p85_0p3x0p7_0p5x0p5_1x1_ep0p5x0p5_"
    "cp0p15x0p85_ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv"
)
ALL_METHODS_CYCLE1_CSV = AUGUST_DIR / (
    "AIDE_phase_I_II_modelsrandom_carryoever_previous_dose_additive_N30_ncycle1_"
    "rp0p15x0p85_ep0p5x0p5_cp0p15x0p85_ap0p15x0p85_IDX_0001_to_1000_"
    "dose_summary.csv"
)

RANDOM_TOXICITY_RUN = PriorRun(
    key="toxicity",
    summary_csv=AUGUST_DIR / (
        "AIDE_phase_I_II_modelsrandom_carryoever_N30_ncycle2_"
        "rp0p15x0p85_0p3x0p7_0p5x0p5_1x1_ep0p5x0p5_"
        "cp0p15x0p85_ap0p15x0p85_IDX_0001_to_1000_dose_summary.csv"
    ),
    output_dir=RANDOM_MODEL_DIR,
    prior_a_column="CRM_Prior_a",
    prior_b_column="CRM_Prior_b",
    display_label="Toxicity prior",
    table_prefix="phase12_toxicity_prior",
    one_stage_table_prefix="phase12_toxicity_prior_one_stage",
    plot_prefix="phase12_toxicity_prior",
    README_prior_description="the CRM toxicity prior",
)

RANDOM_EFFICACY_RUN = PriorRun(
    key="efficacy",
    summary_csv=AUGUST_DIR / (
        "AIDE_phase_I_II_modelsrandom_carryoever_N30_ncycle2_"
        "rp0p15x0p85_ep0p5x0p5_cp0p15x0p85_ap0p15x0p85_"
        "0p3x0p7_0p5x0p5_1x1_IDX_0001_to_1000_dose_summary.csv"
    ),
    output_dir=RANDOM_MODEL_DIR,
    prior_a_column="Efficacy_Additive_Alpha_Prior_a",
    prior_b_column="Efficacy_Additive_Alpha_Prior_b",
    display_label="Efficacy additive alpha prior",
    table_prefix="phase12_efficacy_prior",
    one_stage_table_prefix="phase12_efficacy_prior_one_stage",
    plot_prefix="phase12_efficacy_prior",
    README_prior_description="the additive-efficacy-alpha prior",
)

ALPHA_TOXICITY_RUN = PriorRun(
    key=TOXICITY_RUN.key,
    summary_csv=TOXICITY_RUN.summary_csv,
    output_dir=ALPHA_MODEL_DIR,
    prior_a_column=TOXICITY_RUN.prior_a_column,
    prior_b_column=TOXICITY_RUN.prior_b_column,
    display_label=TOXICITY_RUN.display_label,
    table_prefix="phase12_toxicity_prior",
    one_stage_table_prefix="phase12_toxicity_prior_one_stage",
    plot_prefix=TOXICITY_RUN.plot_prefix,
    README_prior_description=TOXICITY_RUN.README_prior_description,
)

ALPHA_EFFICACY_RUN = PriorRun(
    key=EFFICACY_RUN.key,
    summary_csv=EFFICACY_RUN.summary_csv,
    output_dir=ALPHA_MODEL_DIR,
    prior_a_column=EFFICACY_RUN.prior_a_column,
    prior_b_column=EFFICACY_RUN.prior_b_column,
    display_label=EFFICACY_RUN.display_label,
    table_prefix="phase12_efficacy_prior",
    one_stage_table_prefix="phase12_efficacy_prior_one_stage",
    plot_prefix=EFFICACY_RUN.plot_prefix,
    README_prior_description=EFFICACY_RUN.README_prior_description,
)

DESIGNS = (
    {"allocation": "two_stage", "n_s1": 6, "label": "allocation = two_stage"},
    {"allocation": "one_stage", "n_s1": 30, "label": "allocation = one_stage"},
)
IPDE_ALPHA_UTILITIES = (
    (2, 0.3, "U2"),
    (3, 1.0, "U3"),
)


def display(value: float, decimals: int = 1) -> str:
    """Format an operating-characteristic value without unnecessary zeroes."""
    if pd.isna(value):
        return "-"
    rounded = round(float(value), decimals)
    if rounded == 0:
        rounded = 0.0
    text = f"{rounded:.{decimals}f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def alpha_label(alpha: float) -> str:
    return display(alpha, 1)


def prior_label(run: PriorRun, prior: tuple[float, float]) -> str:
    a, b = prior
    return f"{run.display_label} Beta({display(a, 2)}, {display(b, 2)})"


def close_to(series: pd.Series, value: float) -> pd.Series:
    return (series.astype(float) - value).abs() < 1e-10


def load_truth(path: Path = TRUTH_CSV) -> dict[int, dict[str, list[float]]]:
    """Load the five-dose truth vectors shared by the July and August results."""
    truth_df = pd.read_csv(path)
    scenarios = sorted(truth_df["Scenario"].unique())
    if scenarios != list(range(1, 38)):
        raise ValueError("Truth summary must contain scenarios 1 through 37 exactly once.")
    return {
        int(row["Scenario"]): {
            "dlt": [float(row[f"Tox_Dose{dose}"]) for dose in range(1, 6)],
            "efficacy": [float(row[f"Eff_Dose{dose}"]) for dose in range(1, 6)],
            "utility2": [float(row[f"Utility2_Dose{dose}"]) for dose in range(1, 6)],
            "utility3": [float(row[f"Utility3_Dose{dose}"]) for dose in range(1, 6)],
        }
        for _, row in truth_df.iterrows()
    }


def validate_prior_source(data: pd.DataFrame, run: PriorRun) -> None:
    required = {
        "Allocation", "Nmax", "N_s1", "Utility_Type", "Lambda_T", "Dose",
        "Scenario", "Toxicity_IPDE_Alpha", "Efficacy_IPDE_Alpha", "n_valid",
        "ntrial_from_files", "OBD_Selection_pct", "Pts_Treated", "IPDE_Doses",
        "No_OBD_Selection_pct", run.prior_a_column, run.prior_b_column,
    }
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"{run.summary_csv.name} is missing required columns: {missing}")
    observed_priors = {
        (float(a), float(b))
        for a, b in zip(data[run.prior_a_column], data[run.prior_b_column], strict=True)
    }
    if observed_priors != set(PRIORS):
        raise ValueError(
            f"{run.summary_csv.name} has unexpected {run.display_label.lower()} values: "
            f"{sorted(observed_priors)}"
        )


def extract_prior_records(
    data: pd.DataFrame,
    run: PriorRun,
    alpha: float,
    design: dict[str, Any],
) -> dict[int, dict[tuple[float, float], dict[str, Any]]]:
    """Select the exact N=30 utility-2 records that appear in sensitivity outputs."""
    subset = data[
        close_to(data["Toxicity_IPDE_Alpha"], alpha)
        & close_to(data["Efficacy_IPDE_Alpha"], alpha)
        & (data["Allocation"] == design["allocation"])
        & (data["Nmax"] == 30)
        & (data["N_s1"] == design["n_s1"])
        & (data["Utility_Type"] == 2)
        & close_to(data["Lambda_T"], 0.3)
    ].copy()

    records: dict[int, dict[tuple[float, float], dict[str, Any]]] = {}
    for scenario in range(1, 38):
        records[scenario] = {}
        for prior in PRIORS:
            a, b = prior
            rows = subset[
                (subset["Scenario"] == scenario)
                & close_to(subset[run.prior_a_column], a)
                & close_to(subset[run.prior_b_column], b)
            ].sort_values("Dose")
            if len(rows) != 5 or rows["Dose"].astype(int).tolist() != [1, 2, 3, 4, 5]:
                raise ValueError(
                    f"Expected five dose rows for {run.key}, scenario {scenario}, "
                    f"alpha {alpha}, {design['label']}, Beta({a}, {b}); found {len(rows)}."
                )
            if (
                (rows["n_valid"].astype(float) <= 0).any()
                or not (rows["n_valid"].astype(float) == rows["ntrial_from_files"].astype(float)).all()
            ):
                raise ValueError(
                    f"Invalid replicate count for {run.key}, scenario {scenario}, alpha {alpha}, "
                    f"{design['label']}, Beta({a}, {b})."
                )
            records[scenario][prior] = {
                "selected": rows["OBD_Selection_pct"].astype(float).tolist(),
                "treated": rows["Pts_Treated"].astype(float).tolist(),
                "ipde_patients": rows["IPDE_Doses"].astype(float).tolist(),
                "stop": float(rows["No_OBD_Selection_pct"].iloc[0]),
                "n_valid": int(rows["n_valid"].iloc[0]),
            }
    return records


def scenario_groups(size: int) -> list[list[int]]:
    return [list(range(start, min(start + size, 38))) for start in range(1, 38, size)]


def table_scenario_groups() -> list[list[int]]:
    """Use up to three scenarios per page so the two operating-characteristic blocks remain legible."""
    return [list(range(start, start + 3)) for start in range(1, 34, 3)] + [[34, 35], [36, 37]]


def make_pdf_canvas(output: Path) -> tuple[canvas.Canvas, Path]:
    """Write ReportLab output outside OneDrive, then copy it into the requested folder.

    ReportLab occasionally receives an invalid-argument failure when it opens a
    OneDrive-synchronised PDF directly. Building in tmp first avoids that
    transient filesystem condition without changing the requested output path.
    """
    temporary_dir = ROOT / "tmp" / "pdfs" / "phase12_build"
    temporary_dir.mkdir(parents=True, exist_ok=True)
    temporary_path = temporary_dir / output.name
    return canvas.Canvas(str(temporary_path), pagesize=PAGE_SIZE), temporary_path


def finalize_pdf(pdf: canvas.Canvas, temporary_path: Path, output: Path) -> None:
    pdf.save()
    output.parent.mkdir(parents=True, exist_ok=True)
    last_error: OSError | None = None
    for attempt in range(4):
        try:
            shutil.copyfile(temporary_path, output)
            return
        except OSError as error:
            last_error = error
            time.sleep(0.3 * (attempt + 1))
    raise last_error or OSError(f"Could not copy {temporary_path} to {output}")


def draw_table_header(
    pdf: canvas.Canvas,
    title: str,
    subtitle: str,
    page_number: int,
    page_count: int,
    footer: str,
) -> None:
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 28)
    pdf.drawCentredString(PAGE_WIDTH / 2, PAGE_HEIGHT - 56, title)
    pdf.setFont("Helvetica", 14)
    pdf.drawCentredString(PAGE_WIDTH / 2, PAGE_HEIGHT - 78, subtitle)
    pdf.setStrokeColor(colors.HexColor("#999999"))
    pdf.setLineWidth(0.8)
    pdf.line(LEFT, PAGE_HEIGHT - 92, PAGE_WIDTH - LEFT, PAGE_HEIGHT - 92)

    x = LEFT + LABEL_WIDTH
    header_y = PAGE_HEIGHT - 110
    pdf.setFont("Helvetica-Bold", 13)
    for dose in range(1, 6):
        pdf.drawCentredString(x + DOSE_WIDTH / 2, header_y, f"Dose {dose}")
        x += DOSE_WIDTH
    pdf.drawCentredString(x + STOP_WIDTH / 2, header_y, "Stop %")

    pdf.setFont("Helvetica", 8.3)
    pdf.drawString(LEFT, 38, footer)
    pdf.drawRightString(PAGE_WIDTH - LEFT, 38, f"Page {page_number} of {page_count}")


def draw_table_row(
    pdf: canvas.Canvas,
    y: float,
    label: str,
    values: list[float],
    stop: float | None,
    *,
    shaded: bool = False,
    decimals: int = 1,
    compact: bool = False,
    row_height: float | None = None,
    font_size: float | None = None,
    baseline: float | None = None,
) -> float:
    row_height = row_height if row_height is not None else (9.8 if compact else 13.2)
    font_size = font_size if font_size is not None else (7.4 if compact else 9.5)
    baseline = baseline if baseline is not None else (7.2 if compact else 9.6)
    if shaded:
        pdf.setFillColor(colors.HexColor("#F6F6F6"))
        pdf.rect(LEFT, y - row_height + 2, TABLE_WIDTH, row_height, fill=1, stroke=0)
    pdf.setFillColor(colors.HexColor("#202020"))
    pdf.setFont("Helvetica", font_size)
    pdf.drawString(LEFT + 8, y - baseline, label)
    x = LEFT + LABEL_WIDTH
    for value in values:
        pdf.drawCentredString(x + DOSE_WIDTH / 2, y - baseline, display(value, decimals))
        x += DOSE_WIDTH
    if stop is not None:
        pdf.drawCentredString(x + STOP_WIDTH / 2, y - baseline, display(stop))
    pdf.setStrokeColor(colors.HexColor("#D0D0D0"))
    pdf.setLineWidth(0.45)
    pdf.line(LEFT, y - row_height, LEFT + TABLE_WIDTH, y - row_height)
    return y - row_height


def draw_group_header(
    pdf: canvas.Canvas,
    y: float,
    label: str,
    *,
    compact: bool,
    row_height: float | None = None,
    font_size: float | None = None,
    baseline: float | None = None,
) -> float:
    row_height = row_height if row_height is not None else (9.8 if compact else 13.2)
    font_size = font_size if font_size is not None else (7.6 if compact else 9.6)
    baseline = baseline if baseline is not None else (7.2 if compact else 9.6)
    pdf.setFillColor(colors.HexColor("#E7E7E7"))
    pdf.rect(LEFT, y - row_height + 2, TABLE_WIDTH, row_height, fill=1, stroke=0)
    pdf.setFillColor(colors.HexColor("#202020"))
    pdf.setFont("Helvetica-Bold", font_size)
    pdf.drawString(LEFT + 8, y - baseline, label)
    pdf.setStrokeColor(colors.HexColor("#D0D0D0"))
    pdf.setLineWidth(0.45)
    pdf.line(LEFT, y - row_height, LEFT + TABLE_WIDTH, y - row_height)
    return y - row_height


def draw_prior_table_scenario(
    pdf: canvas.Canvas,
    y: float,
    scenario: int,
    truth: dict[str, list[float]],
    records: dict[tuple[float, float], dict[str, Any]],
    run: PriorRun,
    *,
    compact: bool,
) -> float:
    scenario_height = 14 if compact else 18
    pdf.setFillColor(colors.HexColor("#EFEFEF"))
    pdf.rect(LEFT, y - scenario_height + 2, TABLE_WIDTH, scenario_height, fill=1, stroke=0)
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 10.5 if compact else 13)
    pdf.drawString(LEFT + 8, y - (9 if compact else 12), f"Scenario {scenario}")
    y -= 15 if compact else 19

    y = draw_table_row(pdf, y, "DLT rate", truth["dlt"], None, decimals=2, compact=compact)
    y = draw_table_row(pdf, y, "Efficacy rate", truth["efficacy"], None, decimals=2, compact=compact)
    y = draw_table_row(pdf, y, "Utility 2", truth["utility2"], None, compact=compact)
    y = draw_table_row(pdf, y, "Utility 3", truth["utility3"], None, compact=compact)
    y -= 1 if compact else 2
    y = draw_group_header(pdf, y, "OBD selection (%)", compact=compact)
    for prior in PRIORS:
        record = records[prior]
        y = draw_table_row(
            pdf, y, prior_label(run, prior), record["selected"], record["stop"],
            shaded=True, compact=compact,
        )
    y = draw_group_header(pdf, y, "Mean patients treated", compact=compact)
    for prior in PRIORS:
        y = draw_table_row(
            pdf, y, prior_label(run, prior), records[prior]["treated"], None,
            compact=compact,
        )
    y = draw_group_header(pdf, y, "Mean IPDE patients", compact=compact)
    for prior in PRIORS:
        y = draw_table_row(
            pdf, y, prior_label(run, prior), records[prior]["ipde_patients"], None,
            compact=compact,
        )
    # Preserve room for three scenarios while keeping the footer clear.
    return y - 3


def write_prior_table(
    run: PriorRun,
    alpha: float,
    design: dict[str, Any],
    truth: dict[int, dict[str, list[float]]],
    records: dict[int, dict[tuple[float, float], dict[str, Any]]],
) -> Path:
    prefix = run.one_stage_table_prefix if design["allocation"] == "one_stage" else run.table_prefix
    filename = f"{prefix}_alpha{alpha_label(alpha).replace('.', 'p')}_scenarios_1_to_37_tables.pdf"
    output = run.output_dir / filename
    groups = table_scenario_groups()
    pdf, temporary_path = make_pdf_canvas(output)
    for page_number, scenarios in enumerate(groups, start=1):
        draw_table_header(
            pdf,
            "Phase I/II Operating Characteristics",
            f"{run.display_label} comparison - {design['label']} - IPDE alpha = {alpha_label(alpha)} - N = 30",
            page_number,
            len(groups),
            (
                f"AIDE Phase I/II: {design['label']}, Utility 2 (Lambda_T = 0.3). "
                "Selection % is OBD selection; Stop % is no OBD selection; IPDE patients are mean recycled patients. "
                "Toxicity and efficacy IPDE alphas are both the listed alpha."
            ),
        )
        y = PAGE_HEIGHT - 130
        compact = len(scenarios) == 4
        for scenario in scenarios:
            y = draw_prior_table_scenario(
                pdf, y, scenario, truth[scenario], records[scenario], run, compact=compact,
            )
        pdf.showPage()
    finalize_pdf(pdf, temporary_path, output)
    return output


def extract_prior_both_utility_records(
    data: pd.DataFrame,
    run: PriorRun,
    alpha: float,
    design: dict[str, Any],
) -> dict[int, dict[tuple[float, float], dict[int, dict[str, Any]]]]:
    """Select the N=30 U2 and U3 records for one prior-sensitivity table."""
    subset = data[
        close_to(data["Toxicity_IPDE_Alpha"], alpha)
        & close_to(data["Efficacy_IPDE_Alpha"], alpha)
        & (data["Allocation"] == design["allocation"])
        & (data["Nmax"] == 30)
        & (data["N_s1"] == design["n_s1"])
    ].copy()
    records: dict[int, dict[tuple[float, float], dict[int, dict[str, Any]]]] = {}
    for scenario in range(1, 38):
        records[scenario] = {}
        for prior in PRIORS:
            a, b = prior
            records[scenario][prior] = {}
            for utility_type, lambda_t, _ in IPDE_ALPHA_UTILITIES:
                rows = subset[
                    (subset["Scenario"] == scenario)
                    & close_to(subset[run.prior_a_column], a)
                    & close_to(subset[run.prior_b_column], b)
                    & (subset["Utility_Type"] == utility_type)
                    & close_to(subset["Lambda_T"], lambda_t)
                ].sort_values("Dose")
                if len(rows) != 5 or rows["Dose"].astype(int).tolist() != [1, 2, 3, 4, 5]:
                    raise ValueError(
                        f"Expected five dose rows for {run.key}, scenario {scenario}, "
                        f"utility {utility_type}, alpha {alpha}, {design['label']}, "
                        f"Beta({a}, {b}); found {len(rows)}."
                    )
                if (
                    (rows["n_valid"].astype(float) <= 0).any()
                    or not (rows["n_valid"].astype(float) == rows["ntrial_from_files"].astype(float)).all()
                ):
                    raise ValueError(
                        f"Invalid replicate count for {run.key}, scenario {scenario}, "
                        f"utility {utility_type}, alpha {alpha}, {design['label']}, Beta({a}, {b})."
                    )
                records[scenario][prior][utility_type] = {
                    "selected": rows["OBD_Selection_pct"].astype(float).tolist(),
                    "treated": rows["Pts_Treated"].astype(float).tolist(),
                    "ipde_patients": rows["IPDE_Doses"].astype(float).tolist(),
                    "stop": float(rows["No_OBD_Selection_pct"].iloc[0]),
                }
    return records


def draw_prior_both_utility_table_scenario(
    pdf: canvas.Canvas,
    y: float,
    scenario: int,
    truth: dict[str, list[float]],
    records: dict[tuple[float, float], dict[int, dict[str, Any]]],
    run: PriorRun,
) -> float:
    """Draw a dense but legible prior-sensitivity scenario with U2 and U3 rows."""
    row_height = 10.8
    font_size = 8.2
    baseline = 8.0
    pdf.setFillColor(colors.HexColor("#EFEFEF"))
    pdf.rect(LEFT, y - 16, TABLE_WIDTH, 16, fill=1, stroke=0)
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 12)
    pdf.drawString(LEFT + 8, y - 11, f"Scenario {scenario}")
    y -= 18
    y = draw_table_row(
        pdf, y, "DLT rate", truth["dlt"], None, decimals=2, compact=False,
        row_height=row_height, font_size=font_size, baseline=baseline,
    )
    y = draw_table_row(
        pdf, y, "Efficacy rate", truth["efficacy"], None, decimals=2, compact=False,
        row_height=row_height, font_size=font_size, baseline=baseline,
    )
    y = draw_table_row(
        pdf, y, "Utility 2", truth["utility2"], None, compact=False,
        row_height=row_height, font_size=font_size, baseline=baseline,
    )
    y = draw_table_row(
        pdf, y, "Utility 3", truth["utility3"], None, compact=False,
        row_height=row_height, font_size=font_size, baseline=baseline,
    )
    y -= 1
    for group_label, metric in (
        ("OBD selection (%)", "selected"),
        ("Mean patients treated", "treated"),
        ("Mean IPDE patients", "ipde_patients"),
    ):
        y = draw_group_header(
            pdf, y, group_label, compact=False, row_height=row_height,
            font_size=font_size + 0.2, baseline=baseline,
        )
        for prior in PRIORS:
            for utility_type, _, utility_label in IPDE_ALPHA_UTILITIES:
                record = records[prior][utility_type]
                y = draw_table_row(
                    pdf,
                    y,
                    f"{utility_label} - {prior_label(run, prior)}",
                    record[metric],
                    record["stop"] if metric == "selected" else None,
                    shaded=True,
                    compact=False,
                    row_height=row_height,
                    font_size=font_size,
                    baseline=baseline,
                )
    return y - 4


def write_prior_both_utility_table(
    run: PriorRun,
    alpha: float,
    design: dict[str, Any],
    truth: dict[int, dict[str, list[float]]],
    records: dict[int, dict[tuple[float, float], dict[int, dict[str, Any]]]],
) -> Path:
    """Write a prior-sensitivity PDF with both utility definitions in every scenario."""
    prefix = run.one_stage_table_prefix if design["allocation"] == "one_stage" else run.table_prefix
    filename = f"{prefix}_alpha{alpha_label(alpha).replace('.', 'p')}_scenarios_1_to_37_tables.pdf"
    output = run.output_dir / filename
    groups = combined_ipde_alpha_scenario_groups()
    pdf, temporary_path = make_pdf_canvas(output)
    for page_number, scenarios in enumerate(groups, start=1):
        draw_table_header(
            pdf,
            "Phase I/II Operating Characteristics",
            (
                f"{run.display_label} comparison - {design['label']} - IPDE alpha = "
                f"{alpha_label(alpha)} - both utilities - N = 30"
            ),
            page_number,
            len(groups),
            (
                "U2: Lambda_T = 0.3; U3: Lambda_T = 1.0. Selection % is OBD selection; "
                "Stop % is no OBD selection; IPDE patients are mean recycled patients. "
                "Toxicity and efficacy IPDE alphas are both the listed alpha."
            ),
        )
        y = PAGE_HEIGHT - 130
        for scenario in scenarios:
            y = draw_prior_both_utility_table_scenario(
                pdf, y, scenario, truth[scenario], records[scenario], run,
            )
        pdf.showPage()
    finalize_pdf(pdf, temporary_path, output)
    return output


def write_prior_both_utility_readme(run: PriorRun, table_outputs: list[Path]) -> None:
    """Document the refreshed prior-sensitivity tables without duplicating the Phase I-II alpha tables."""
    tables = "\n".join(f"  - {path.name}" for path in table_outputs)
    (run.output_dir / "README.txt").write_text(
        "Built from the August Phase I/II N=30 summary listed below.\n\n"
        f"Source: {run.summary_csv.name}\n\n"
        "Tables: eight 19-page landscape PDFs covering scenarios 1-37 at IPDE alpha 0, 0.3, 0.6, and 0.9. "
        "For every alpha, one table is supplied per allocation design (allocation = one_stage and allocation = two_stage). "
        "Each table includes both utilities: U2 (Utility_Type=2, Lambda_T=0.3) and U3 "
        "(Utility_Type=3, Lambda_T=1.0). Every scenario shows OBD selection, mean patients treated, "
        "mean IPDE patients, and Stop %.\n\n"
        f"Prior varied: {run.README_prior_description}. The four settings are Beta(0.15,0.85), "
        "Beta(0.30,0.70), Beta(0.50,0.50), and Beta(1,1).\n"
        "Toxicity and efficacy IPDE alpha are both set to the listed alpha in each output.\n\n"
        "Table files:\n"
        f"{tables}\n",
        encoding="utf-8",
    )


def build_prior_both_utility_outputs(run: PriorRun) -> list[Path]:
    """Refresh all named prior-sensitivity PDFs with U2 and U3 in the same table."""
    if not run.summary_csv.exists():
        raise FileNotFoundError(run.summary_csv)
    run.output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(run.summary_csv)
    validate_prior_source(data, run)
    truth = load_truth()
    outputs = [
        write_prior_both_utility_table(
            run,
            alpha,
            design,
            truth,
            extract_prior_both_utility_records(data, run, alpha, design),
        )
        for alpha in ALPHAS
        for design in DESIGNS
    ]
    write_prior_both_utility_readme(run, outputs)
    return outputs


def validate_model_prior_source(
    data: pd.DataFrame,
    run: PriorRun,
    expected_model_id: str,
) -> None:
    """Verify a model-specific N=30 prior run before publishing its tables."""
    validate_prior_source(data, run)
    required = {"Model_ID", "Cycle_Max"}
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"{run.summary_csv.name} is missing required columns: {missing}")
    observed_models = set(data["Model_ID"].astype(str).unique())
    if observed_models != {expected_model_id}:
        raise ValueError(
            f"{run.summary_csv.name} has Model_ID {sorted(observed_models)}; "
            f"expected {expected_model_id}."
        )
    observed_cycles = set(data["Cycle_Max"].astype(int).unique())
    if observed_cycles != {2}:
        raise ValueError(
            f"{run.summary_csv.name} has Cycle_Max {sorted(observed_cycles)}; expected 2."
        )


def write_model_prior_both_utility_table(
    run: PriorRun,
    model_label: str,
    alpha: float,
    design: dict[str, Any],
    truth: dict[int, dict[str, list[float]]],
    records: dict[int, dict[tuple[float, float], dict[int, dict[str, Any]]]],
) -> Path:
    """Write a model-specific prior table with explicit one_stage/two_stage filenames."""
    filename = (
        f"phase12_{run.key}_prior_{design['allocation']}_"
        f"alpha{alpha_label(alpha).replace('.', 'p')}_scenarios_1_to_37_tables.pdf"
    )
    output = run.output_dir / filename
    groups = combined_ipde_alpha_scenario_groups()
    pdf, temporary_path = make_pdf_canvas(output)
    for page_number, scenarios in enumerate(groups, start=1):
        draw_table_header(
            pdf,
            "Phase I/II Operating Characteristics",
            (
                f"{model_label} - {run.display_label} comparison - {design['label']} - "
                f"IPDE alpha = {alpha_label(alpha)} - both utilities - N = 30"
            ),
            page_number,
            len(groups),
            (
                "U2: Lambda_T = 0.3; U3: Lambda_T = 1.0. Selection % is OBD selection; "
                "Stop % is no OBD selection; IPDE patients are mean recycled patients. "
                "Toxicity and efficacy IPDE alphas are both the listed alpha."
            ),
        )
        y = PAGE_HEIGHT - 130
        for scenario in scenarios:
            y = draw_prior_both_utility_table_scenario(
                pdf, y, scenario, truth[scenario], records[scenario], run,
            )
        pdf.showPage()
    finalize_pdf(pdf, temporary_path, output)
    return output


def write_model_prior_readme(
    output_dir: Path,
    model_label: str,
    model_description: str,
    runs: tuple[PriorRun, PriorRun],
    table_outputs: list[Path],
) -> None:
    """Describe the separate toxicity- and efficacy-prior tables for one model family."""
    tables = "\n".join(f"  - {path.name}" for path in table_outputs)
    sources = "\n".join(f"- Presentation 8-03-2026/{run.summary_csv.name}" for run in runs)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "README.txt").write_text(
        f"{model_label} prior-sensitivity tables for the August Phase I/II N=30, Cycle_Max=2 runs.\n\n"
        f"Model: {model_description}\n\n"
        "Sources:\n"
        f"{sources}\n\n"
        "Tables: 16 landscape PDFs covering toxicity-prior and efficacy-prior sensitivity at IPDE alpha 0, 0.3, "
        "0.6, and 0.9. Filenames explicitly identify one_stage or two_stage. Each table contains both utilities: "
        "U2 (Utility_Type=2, Lambda_T=0.3) and U3 (Utility_Type=3, Lambda_T=1.0), with scenarios 1-37 "
        "shown two per page.\n\n"
        "The four prior settings are Beta(0.15,0.85), Beta(0.30,0.70), Beta(0.50,0.50), and Beta(1,1). "
        "Every table reports OBD selection, mean patients treated, mean IPDE patients, and Stop %.\n\n"
        "Table files:\n"
        f"{tables}\n",
        encoding="utf-8",
    )


def build_model_prior_outputs(
    output_dir: Path,
    model_label: str,
    model_description: str,
    expected_model_id: str,
    runs: tuple[PriorRun, PriorRun],
) -> list[Path]:
    """Create parallel toxicity- and efficacy-prior PDF sets for one model family."""
    truth = load_truth()
    outputs: list[Path] = []
    for run in runs:
        if not run.summary_csv.exists():
            raise FileNotFoundError(run.summary_csv)
        data = pd.read_csv(run.summary_csv)
        validate_model_prior_source(data, run, expected_model_id)
        for alpha in ALPHAS:
            for design in DESIGNS:
                outputs.append(
                    write_model_prior_both_utility_table(
                        run,
                        model_label,
                        alpha,
                        design,
                        truth,
                        extract_prior_both_utility_records(data, run, alpha, design),
                    )
                )
    write_model_prior_readme(output_dir, model_label, model_description, runs, outputs)
    return outputs


def extract_ipde_alpha_records(
    data: pd.DataFrame,
    design: dict[str, Any],
) -> dict[int, dict[float, dict[str, Any]]]:
    """Hold both model priors at Beta(0.15, 0.85) and compare IPDE alpha."""
    subset = data[
        (data["Allocation"] == design["allocation"])
        & (data["Nmax"] == 30)
        & (data["N_s1"] == design["n_s1"])
        & (data["Utility_Type"] == 2)
        & close_to(data["Lambda_T"], 0.3)
        & close_to(data["CRM_Prior_a"], 0.15)
        & close_to(data["CRM_Prior_b"], 0.85)
        & close_to(data["Efficacy_Additive_Alpha_Prior_a"], 0.15)
        & close_to(data["Efficacy_Additive_Alpha_Prior_b"], 0.85)
    ].copy()
    output: dict[int, dict[float, dict[str, Any]]] = {}
    for scenario in range(1, 38):
        output[scenario] = {}
        for alpha in ALPHAS:
            rows = subset[
                (subset["Scenario"] == scenario)
                & close_to(subset["Toxicity_IPDE_Alpha"], alpha)
                & close_to(subset["Efficacy_IPDE_Alpha"], alpha)
            ].sort_values("Dose")
            if len(rows) != 5 or rows["Dose"].astype(int).tolist() != [1, 2, 3, 4, 5]:
                raise ValueError(
                    f"Expected five dose rows for scenario {scenario}, IPDE alpha {alpha}, "
                    f"{design['allocation']}; found {len(rows)}."
                )
            output[scenario][alpha] = {
                "selected": rows["OBD_Selection_pct"].astype(float).tolist(),
                "treated": rows["Pts_Treated"].astype(float).tolist(),
                "ipde_patients": rows["IPDE_Doses"].astype(float).tolist(),
                "stop": float(rows["No_OBD_Selection_pct"].iloc[0]),
            }
    return output


def draw_ipde_alpha_table_scenario(
    pdf: canvas.Canvas,
    y: float,
    scenario: int,
    truth: dict[str, list[float]],
    records: dict[float, dict[str, Any]],
    *,
    compact: bool,
) -> float:
    scenario_height = 14 if compact else 18
    pdf.setFillColor(colors.HexColor("#EFEFEF"))
    pdf.rect(LEFT, y - scenario_height + 2, TABLE_WIDTH, scenario_height, fill=1, stroke=0)
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 10.5 if compact else 13)
    pdf.drawString(LEFT + 8, y - (9 if compact else 12), f"Scenario {scenario}")
    y -= 15 if compact else 19
    y = draw_table_row(pdf, y, "DLT rate", truth["dlt"], None, decimals=2, compact=compact)
    y = draw_table_row(pdf, y, "Efficacy rate", truth["efficacy"], None, decimals=2, compact=compact)
    y = draw_table_row(pdf, y, "Utility 2", truth["utility2"], None, compact=compact)
    y = draw_table_row(pdf, y, "Utility 3", truth["utility3"], None, compact=compact)
    y -= 1 if compact else 2
    y = draw_group_header(pdf, y, "OBD selection (%)", compact=compact)
    for alpha in ALPHAS:
        record = records[alpha]
        y = draw_table_row(
            pdf,
            y,
            f"IPDE alpha = {alpha_label(alpha)}",
            record["selected"],
            record["stop"],
            shaded=True,
            compact=compact,
        )
    y = draw_group_header(pdf, y, "Mean patients treated", compact=compact)
    for alpha in ALPHAS:
        y = draw_table_row(
            pdf,
            y,
            f"IPDE alpha = {alpha_label(alpha)}",
            records[alpha]["treated"],
            None,
            compact=compact,
        )
    y = draw_group_header(pdf, y, "Mean IPDE patients", compact=compact)
    for alpha in ALPHAS:
        y = draw_table_row(
            pdf,
            y,
            f"IPDE alpha = {alpha_label(alpha)}",
            records[alpha]["ipde_patients"],
            None,
            compact=compact,
        )
    # Preserve room for three scenarios while keeping the footer clear.
    return y - 3


def write_ipde_alpha_table(
    run: PriorRun,
    design: dict[str, Any],
    truth: dict[int, dict[str, list[float]]],
    records: dict[int, dict[float, dict[str, Any]]],
) -> Path:
    output = run.output_dir / (
        f"phase12_{run.key}_ipde_alpha_comparison_{design['allocation']}_"
        "scenarios_1_to_37_tables.pdf"
    )
    groups = table_scenario_groups()
    pdf, temporary_path = make_pdf_canvas(output)
    for page_number, scenarios in enumerate(groups, start=1):
        draw_table_header(
            pdf,
            "Phase I/II Operating Characteristics",
            f"IPDE alpha comparison - {design['label']} - N = 30",
            page_number,
            len(groups),
            (
                "Utility 2 (Lambda_T = 0.3). Fixed model priors: CRM toxicity Beta(0.15,0.85) and "
                "additive-efficacy-alpha Beta(0.15,0.85). Both IPDE alphas equal the listed row value; "
                "IPDE patients are mean recycled patients."
            ),
        )
        y = PAGE_HEIGHT - 130
        compact = len(scenarios) == 4
        for scenario in scenarios:
            y = draw_ipde_alpha_table_scenario(
                pdf, y, scenario, truth[scenario], records[scenario], compact=compact,
            )
        pdf.showPage()
    finalize_pdf(pdf, temporary_path, output)
    return output


def validate_combined_ipde_alpha_source(data: pd.DataFrame) -> None:
    """Validate the fixed-prior source used for the combined alpha comparisons."""
    required = {
        "Allocation", "Nmax", "N_s1", "Utility_Type", "Lambda_T", "Dose", "Scenario",
        "CRM_Prior_a", "CRM_Prior_b", "Efficacy_Additive_Alpha_Prior_a",
        "Efficacy_Additive_Alpha_Prior_b", "Toxicity_IPDE_Alpha", "Efficacy_IPDE_Alpha",
        "OBD_Selection_pct", "Pts_Treated", "IPDE_Doses", "No_OBD_Selection_pct",
        "n_valid", "ntrial_from_files",
    }
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(
            f"{COMBINED_IPDE_ALPHA_CSV.name} is missing required columns: {missing}"
        )
    if sorted(data["Scenario"].unique().tolist()) != list(range(1, 38)):
        raise ValueError("The combined IPDE-alpha source must contain scenarios 1 through 37.")


def extract_combined_ipde_alpha_records(
    data: pd.DataFrame,
    design: dict[str, Any],
) -> dict[int, dict[int, dict[float, dict[str, Any]]]]:
    """Compare shared IPDE alpha at fixed Beta(0.15, 0.85) priors for U2 and U3."""
    subset = data[
        (data["Allocation"] == design["allocation"])
        & (data["Nmax"] == 30)
        & (data["N_s1"] == design["n_s1"])
        & close_to(data["CRM_Prior_a"], 0.15)
        & close_to(data["CRM_Prior_b"], 0.85)
        & close_to(data["Efficacy_Additive_Alpha_Prior_a"], 0.15)
        & close_to(data["Efficacy_Additive_Alpha_Prior_b"], 0.85)
    ].copy()
    output: dict[int, dict[int, dict[float, dict[str, Any]]]] = {}
    for scenario in range(1, 38):
        output[scenario] = {}
        for utility_type, lambda_t, _ in IPDE_ALPHA_UTILITIES:
            output[scenario][utility_type] = {}
            for alpha in ALPHAS:
                rows = subset[
                    (subset["Scenario"] == scenario)
                    & (subset["Utility_Type"] == utility_type)
                    & close_to(subset["Lambda_T"], lambda_t)
                    & close_to(subset["Toxicity_IPDE_Alpha"], alpha)
                    & close_to(subset["Efficacy_IPDE_Alpha"], alpha)
                ].sort_values("Dose")
                if len(rows) != 5 or rows["Dose"].astype(int).tolist() != [1, 2, 3, 4, 5]:
                    raise ValueError(
                        f"Expected five dose rows for scenario {scenario}, utility {utility_type}, "
                        f"IPDE alpha {alpha}, {design['allocation']}; found {len(rows)}."
                    )
                if (
                    (rows["n_valid"].astype(float) <= 0).any()
                    or not (rows["n_valid"].astype(float) == rows["ntrial_from_files"].astype(float)).all()
                ):
                    raise ValueError(
                        f"Invalid replicate count for scenario {scenario}, utility {utility_type}, "
                        f"IPDE alpha {alpha}, {design['allocation']}."
                    )
                output[scenario][utility_type][alpha] = {
                    "selected": rows["OBD_Selection_pct"].astype(float).tolist(),
                    "treated": rows["Pts_Treated"].astype(float).tolist(),
                    "ipde_patients": rows["IPDE_Doses"].astype(float).tolist(),
                    "stop": float(rows["No_OBD_Selection_pct"].iloc[0]),
                }
    return output


def draw_combined_ipde_alpha_table_scenario(
    pdf: canvas.Canvas,
    y: float,
    scenario: int,
    truth: dict[str, list[float]],
    records: dict[int, dict[float, dict[str, Any]]],
) -> float:
    """Draw one scenario with both utility definitions, using a readable dense layout."""
    row_height = 10.8
    font_size = 8.2
    baseline = 8.0
    pdf.setFillColor(colors.HexColor("#EFEFEF"))
    pdf.rect(LEFT, y - 16, TABLE_WIDTH, 16, fill=1, stroke=0)
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 12)
    pdf.drawString(LEFT + 8, y - 11, f"Scenario {scenario}")
    y -= 18
    y = draw_table_row(
        pdf, y, "DLT rate", truth["dlt"], None, decimals=2, compact=False,
        row_height=row_height, font_size=font_size, baseline=baseline,
    )
    y = draw_table_row(
        pdf, y, "Efficacy rate", truth["efficacy"], None, decimals=2, compact=False,
        row_height=row_height, font_size=font_size, baseline=baseline,
    )
    y = draw_table_row(
        pdf, y, "Utility 2", truth["utility2"], None, compact=False,
        row_height=row_height, font_size=font_size, baseline=baseline,
    )
    y = draw_table_row(
        pdf, y, "Utility 3", truth["utility3"], None, compact=False,
        row_height=row_height, font_size=font_size, baseline=baseline,
    )
    y -= 1
    for group_label, metric in (
        ("OBD selection (%)", "selected"),
        ("Mean patients treated", "treated"),
        ("Mean IPDE patients", "ipde_patients"),
    ):
        y = draw_group_header(
            pdf, y, group_label, compact=False, row_height=row_height,
            font_size=font_size + 0.2, baseline=baseline,
        )
        for utility_type, _, utility_label in IPDE_ALPHA_UTILITIES:
            for alpha in ALPHAS:
                record = records[utility_type][alpha]
                y = draw_table_row(
                    pdf,
                    y,
                    f"{utility_label} - IPDE alpha = {alpha_label(alpha)}",
                    record[metric],
                    record["stop"] if metric == "selected" else None,
                    shaded=True,
                    compact=False,
                    row_height=row_height,
                    font_size=font_size,
                    baseline=baseline,
                )
    return y - 4


def combined_ipde_alpha_scenario_groups() -> list[list[int]]:
    """Two scenarios per page retain a legible label and row font with both utilities."""
    return [list(range(start, min(start + 2, 38))) for start in range(1, 38, 2)]


def write_combined_ipde_alpha_tables() -> list[Path]:
    """Create one-stage and two-stage IPDE-alpha tables with U2 and U3 together."""
    if not COMBINED_IPDE_ALPHA_CSV.exists():
        raise FileNotFoundError(COMBINED_IPDE_ALPHA_CSV)
    data = pd.read_csv(COMBINED_IPDE_ALPHA_CSV)
    validate_combined_ipde_alpha_source(data)
    truth = load_truth()
    groups = combined_ipde_alpha_scenario_groups()
    outputs: list[Path] = []
    for design in DESIGNS:
        records = extract_combined_ipde_alpha_records(data, design)
        output = N30_COMPARISON_DIR / (
            f"phase12_ipde_alpha_comparison_{design['allocation']}_"
            "scenarios_1_to_37_tables.pdf"
        )
        pdf, temporary_path = make_pdf_canvas(output)
        for page_number, scenarios in enumerate(groups, start=1):
            draw_table_header(
                pdf,
                "Phase I/II Operating Characteristics",
                f"IPDE alpha comparison - {design['label']} - both utilities - N = 30",
                page_number,
                len(groups),
                (
                    "Fixed priors: CRM toxicity Beta(0.15, 0.85) and additive efficacy Beta(0.15, 0.85). "
                    "U2: Lambda_T = 0.3; U3: Lambda_T = 1.0. Both IPDE alphas equal the listed row value; "
                    "IPDE patients are mean recycled patients."
                ),
            )
            y = PAGE_HEIGHT - 130
            for scenario in scenarios:
                y = draw_combined_ipde_alpha_table_scenario(
                    pdf, y, scenario, truth[scenario], records[scenario],
                )
            pdf.showPage()
        finalize_pdf(pdf, temporary_path, output)
        outputs.append(output)
    return outputs


PLOT_COLORS = (colors.HexColor("#0072B2"), colors.HexColor("#D55E00"), colors.HexColor("#009E73"), colors.HexColor("#CC79A7"))


def nice_axis_max(values: list[float], minimum: float) -> float:
    maximum = max([minimum, *[float(value) for value in values if pd.notna(value)]])
    if maximum <= 5:
        return 5
    if maximum <= 10:
        return 10
    if maximum <= 20:
        return 20
    return float(int((maximum + 9.999) // 10) * 10)


def draw_plot_panel(
    pdf: canvas.Canvas,
    x: float,
    y_top: float,
    width: float,
    height: float,
    title: str,
    y_label: str,
    values_by_prior: dict[tuple[float, float], list[float]],
    y_max: float,
    compact: bool,
) -> None:
    """Draw one five-dose line panel; y_top is the panel's upper coordinate."""
    title_size = 8.0 if compact else 9.5
    tick_size = 6.2 if compact else 7.5
    margin_left = 39 if compact else 45
    margin_right = 14
    margin_top = 15
    margin_bottom = 23
    bottom = y_top - height + margin_bottom
    top = y_top - margin_top
    left = x + margin_left
    right = x + width - margin_right
    pdf.setFillColor(colors.HexColor("#F4F4F4"))
    pdf.rect(x, y_top - height, width, height, fill=1, stroke=0)
    pdf.setFillColor(colors.HexColor("#202020"))
    pdf.setFont("Helvetica", title_size)
    pdf.drawCentredString((left + right) / 2, y_top - 9, title)

    for tick in range(6):
        value = y_max * tick / 5
        yy = bottom + (top - bottom) * tick / 5
        pdf.setStrokeColor(colors.white)
        pdf.setLineWidth(0.8)
        pdf.line(left, yy, right, yy)
        pdf.setFillColor(colors.HexColor("#4D4D4D"))
        pdf.setFont("Helvetica", tick_size)
        pdf.drawRightString(left - 5, yy - 2.2, display(value, 0 if y_max >= 20 else 1))

    x_values = [left + (right - left) * index / 4 for index in range(5)]
    for dose, xx in enumerate(x_values, start=1):
        pdf.setStrokeColor(colors.white)
        pdf.setLineWidth(0.55)
        pdf.line(xx, bottom, xx, top)
        pdf.setFillColor(colors.HexColor("#4D4D4D"))
        pdf.setFont("Helvetica", tick_size)
        pdf.drawCentredString(xx, bottom - 12, str(dose))

    pdf.saveState()
    pdf.setFont("Helvetica", tick_size)
    pdf.setFillColor(colors.HexColor("#333333"))
    pdf.translate(x + 11, (bottom + top) / 2)
    pdf.rotate(90)
    pdf.drawCentredString(0, 0, y_label)
    pdf.restoreState()

    for color, prior in zip(PLOT_COLORS, PRIORS, strict=True):
        values = values_by_prior[prior]
        points = [
            (xx, bottom + (max(0, value) / y_max) * (top - bottom))
            for xx, value in zip(x_values, values, strict=True)
        ]
        pdf.setStrokeColor(color)
        pdf.setLineWidth(1.25 if compact else 1.6)
        path = pdf.beginPath()
        path.moveTo(*points[0])
        for point in points[1:]:
            path.lineTo(*point)
        pdf.drawPath(path, stroke=1, fill=0)
        pdf.setFillColor(color)
        for point in points:
            pdf.circle(*point, 1.8 if compact else 2.3, fill=1, stroke=0)


def draw_plot_legend(pdf: canvas.Canvas, run: PriorRun) -> None:
    start_x = LEFT + 48
    y = PAGE_HEIGHT - 105
    pdf.setFont("Helvetica", 8.5)
    for color, prior in zip(PLOT_COLORS, PRIORS, strict=True):
        pdf.setFillColor(color)
        pdf.rect(start_x, y - 7, 13, 7, fill=1, stroke=0)
        pdf.setFillColor(colors.HexColor("#202020"))
        pdf.drawString(start_x + 18, y - 7, f"Beta({display(prior[0], 2)}, {display(prior[1], 2)})")
        start_x += 196
    pdf.setFont("Helvetica", 8)
    pdf.drawRightString(PAGE_WIDTH - LEFT, y - 7, run.display_label)


def write_prior_plots(
    run: PriorRun,
    alpha: float,
    records_by_design: dict[str, dict[int, dict[tuple[float, float], dict[str, Any]]]],
) -> Path:
    """Create one multi-page, two-design plot set for a requested IPDE alpha."""
    output = run.output_dir / (
        f"{run.plot_prefix}_alpha{alpha_label(alpha).replace('.', 'p')}_scenarios_1_to_37_plots.pdf"
    )
    groups = scenario_groups(2)
    pdf, temporary_path = make_pdf_canvas(output)
    for page_number, scenarios in enumerate(groups, start=1):
        pdf.setFillColor(colors.black)
        pdf.setFont("Helvetica-Bold", 25)
        pdf.drawCentredString(PAGE_WIDTH / 2, PAGE_HEIGHT - 54, "Phase I/II Prior-Sensitivity Plots")
        pdf.setFont("Helvetica", 13)
        pdf.drawCentredString(
            PAGE_WIDTH / 2,
            PAGE_HEIGHT - 77,
            f"{run.display_label} - IPDE alpha = {alpha_label(alpha)} - N = 30",
        )
        pdf.setStrokeColor(colors.HexColor("#999999"))
        pdf.setLineWidth(0.8)
        pdf.line(LEFT, PAGE_HEIGHT - 91, PAGE_WIDTH - LEFT, PAGE_HEIGHT - 91)
        draw_plot_legend(pdf, run)
        pdf.setFont("Helvetica", 8)
        pdf.setFillColor(colors.HexColor("#202020"))
        pdf.drawString(
            LEFT,
            34,
            "Each panel compares the four Beta priors. Left: OBD selection. Right: mean patients treated. "
            "Companion tables report Stop %.",
        )
        pdf.drawRightString(PAGE_WIDTH - LEFT, 34, f"Page {page_number} of {len(groups)}")

        usable_height = PAGE_HEIGHT - 137 - 54
        block_height = usable_height / len(scenarios)
        for index, scenario in enumerate(scenarios):
            block_top = PAGE_HEIGHT - 130 - index * block_height
            compact = len(scenarios) > 1
            pdf.setFillColor(colors.HexColor("#E7E7E7"))
            pdf.rect(LEFT, block_top - 17, PAGE_WIDTH - 2 * LEFT, 17, fill=1, stroke=0)
            pdf.setFillColor(colors.black)
            pdf.setFont("Helvetica-Bold", 10.5 if compact else 13)
            pdf.drawString(LEFT + 7, block_top - 12, f"Scenario {scenario}")

            block_bottom = block_top - block_height + 12
            panel_top = block_top - 25
            row_gap = 14 if compact else 20
            panel_height = (panel_top - block_bottom - row_gap) / 2
            panel_gap = 38
            panel_width = (PAGE_WIDTH - 2 * LEFT - panel_gap) / 2
            treated_values = [
                value
                for design in DESIGNS
                for prior in PRIORS
                for value in records_by_design[design["allocation"]][scenario][prior]["treated"]
            ]
            treated_max = nice_axis_max(treated_values, 5)

            for row, design in enumerate(DESIGNS):
                row_top = panel_top - row * (panel_height + row_gap)
                records = records_by_design[design["allocation"]][scenario]
                draw_plot_panel(
                    pdf,
                    LEFT,
                    row_top,
                    panel_width,
                    panel_height,
                    f"{design['label']} - OBD selection",
                    "Percent",
                    {prior: records[prior]["selected"] for prior in PRIORS},
                    100,
                    compact,
                )
                draw_plot_panel(
                    pdf,
                    LEFT + panel_width + panel_gap,
                    row_top,
                    panel_width,
                    panel_height,
                    f"{design['label']} - mean patients treated",
                    "Patients",
                    {prior: records[prior]["treated"] for prior in PRIORS},
                    treated_max,
                    compact,
                )
        pdf.showPage()
    finalize_pdf(pdf, temporary_path, output)
    return output


def write_prior_readme(
    run: PriorRun,
    table_outputs: list[Path],
    alpha_comparison_outputs: list[Path],
) -> None:
    tables = "\n".join(f"  - {path.name}" for path in table_outputs)
    alpha_comparison_tables = "\n".join(f"  - {path.name}" for path in alpha_comparison_outputs)
    run.output_dir.mkdir(parents=True, exist_ok=True)
    (run.output_dir / "README.txt").write_text(
        "Built from the August Phase I/II N=30 summary listed below.\n\n"
        f"Source: {run.summary_csv.name}\n\n"
        "Tables: eight 13-page landscape PDFs covering scenarios 1-37 at IPDE alpha 0, 0.3, 0.6, and 0.9. "
        "For every alpha, one table is supplied per allocation design (allocation = one_stage and allocation = two_stage). "
        "Each table uses Utility Type 2 and Lambda_T = 0.3, and shows OBD selection, mean patients treated, "
        "mean IPDE patients, and Stop %.\n\n"
        "IPDE-alpha comparison: two additional 12-page landscape PDFs compare IPDE alpha 0, 0.3, 0.6, and 0.9 "
        "while holding both the CRM toxicity prior and additive-efficacy-alpha prior at Beta(0.15,0.85). "
        "One table is supplied for allocation = one_stage and one for allocation = two_stage.\n\n"
        f"Prior varied: {run.README_prior_description}. The four settings are Beta(0.15,0.85), "
        "Beta(0.30,0.70), Beta(0.50,0.50), and Beta(1,1).\n"
        "Toxicity and efficacy IPDE alpha are both set to the listed alpha in each output.\n\n"
        "Table files:\n"
        f"{tables}\n\n"
        "IPDE-alpha comparison files:\n"
        f"{alpha_comparison_tables}\n",
        encoding="utf-8",
    )


def parse_operating_characteristics_pdf(path: Path, scenario_offset: int = 0) -> dict[int, dict[str, Any]]:
    """Parse BOIN12/U-BOIN raw result PDFs into the values used by the July table."""
    expected_labels = ("No. Pts treated", "# Pts treated")
    records: dict[int, dict[str, Any]] = {}
    with pdfplumber.open(path) as document:
        text = "\n".join(page.extract_text() or "" for page in document.pages)
    blocks = list(re.finditer(r"^Scenario\s+(\d+)\s*$", text, flags=re.MULTILINE))
    for index, match in enumerate(blocks):
        block_end = blocks[index + 1].start() if index + 1 < len(blocks) else len(text)
        block = text[match.end():block_end]
        scenario = int(match.group(1)) + scenario_offset

        def values_for(labels: tuple[str, ...], count: int) -> list[float]:
            for label in labels:
                found = re.search(re.escape(label) + r"\s+([^\n]+)", block)
                if found:
                    values = [float(value) for value in re.findall(r"-?(?:\d+\.?\d*|\.\d+)", found.group(1))]
                    if len(values) >= count:
                        return values[-count:]
            raise ValueError(f"Could not find {labels} values for scenario {scenario} in {path.name}.")

        records[scenario] = {
            "treated": values_for(expected_labels, 5),
            "selected": values_for(("Select %",), 6)[:5],
            "stop": values_for(("Select %",), 6)[5],
        }
    if not records:
        raise ValueError(f"No operating-characteristic scenarios found in {path.name}.")
    return records


def parse_efftox_html(path: Path) -> dict[int, dict[str, Any]]:
    content = path.read_text(encoding="utf-8")
    starts = list(re.finditer(r'<tr><td\s+colspan\s*=\s*"8"><b>\s*(\d+)\s*</b></tr>', content, re.IGNORECASE))
    records: dict[int, dict[str, Any]] = {}
    for index, start in enumerate(starts):
        block_end = starts[index + 1].start() if index + 1 < len(starts) else len(content)
        block = content[start.start():block_end]

        def values_after(label: str, terminator: str) -> list[float]:
            """Extract cells after an EffTox label without crossing into the next metric row.

            The raw EffTox HTML omits the closing ``</tr>`` for its ``% selected``
            rows. A conventional row-end search therefore includes the following
            ``# Patients Treated`` row. The caller supplies the first marker that
            begins the next row instead.
            """
            found = re.search(
                label + r"\s*</b>\s*</td>([\s\S]*?)" + terminator,
                block,
                re.IGNORECASE,
            )
            if not found:
                return []
            raw_cells = re.findall(r"<td[^>]*>([\s\S]*?)</td>", found.group(1), re.IGNORECASE)
            values: list[float] = []
            for cell in raw_cells:
                value = re.sub(r"\s+", " ", re.sub(r"<[^>]+>|&nbsp;", " ", html.unescape(cell))).strip()
                values.append(float("nan") if value in {"", "-", "--"} else float(value))
            return values

        selected_values = values_after(r"%\s*selected", r"<tr\b")
        treated_values = values_after(r"#\s*Patients\s*Treated", r"</tr>")
        if len(selected_values) != 6 or len(treated_values) != 6:
            raise ValueError(f"Could not parse five doses plus Stop % for EffTox scenario {start.group(1)}.")
        records[int(start.group(1))] = {
            "treated": treated_values[:5],
            "selected": selected_values[:5],
            "stop": selected_values[5],
        }
    if set(records) != set(range(1, 38)):
        raise ValueError(f"EffTox N=30 source is incomplete: found {sorted(records)}")
    return records


def validate_cycle1_all_methods_source(data: pd.DataFrame) -> None:
    required = {
        "Scenario", "Dose", "Allocation", "Nmax", "N_s1", "Model_ID",
        "CRM_r_Model", "Efficacy_Model", "Cycle_Max", "Utility_Type",
        "Lambda_T", "OBD_Selection_pct", "Pts_Treated", "No_OBD_Selection_pct",
        "n_valid", "ntrial_from_files",
    }
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"{ALL_METHODS_CYCLE1_CSV.name} is missing required columns: {missing}")
    if set(data["Model_ID"]) != {"random_carryoever", "previous_dose_additive"}:
        raise ValueError("Cycle-1 source must contain random_carryoever and previous_dose_additive models.")
    if set(data["Cycle_Max"].astype(int)) != {1}:
        raise ValueError("The all-method source must contain Cycle_Max = 1 only.")
    expected_models = {
        "random_carryoever": ("random", "dose_specific_carryover"),
        "previous_dose_additive": ("previous_dose", "previous_dose_additive"),
    }
    for model_id, (crm_model, efficacy_model) in expected_models.items():
        rows = data[data["Model_ID"] == model_id]
        observed = set(zip(rows["CRM_r_Model"], rows["Efficacy_Model"], strict=True))
        if observed != {(crm_model, efficacy_model)}:
            raise ValueError(
                f"Unexpected model configuration for {model_id}: {sorted(observed)}."
            )


def extract_cycle1_aide_records(
    data: pd.DataFrame,
    model_id: str,
    allocation: str,
    utility_type: int,
    lambda_t: float,
) -> dict[int, dict[str, Any]]:
    n_s1 = 30 if allocation == "one_stage" else 6
    subset = data[
        (data["Model_ID"] == model_id)
        & (data["Cycle_Max"] == 1)
        & (data["Allocation"] == allocation)
        & (data["Nmax"] == 30)
        & (data["N_s1"] == n_s1)
        & (data["Utility_Type"] == utility_type)
        & close_to(data["Lambda_T"], lambda_t)
    ].copy()
    output: dict[int, dict[str, Any]] = {}
    for scenario in range(1, 38):
        rows = subset[subset["Scenario"] == scenario].sort_values("Dose")
        if len(rows) != 5 or rows["Dose"].astype(int).tolist() != [1, 2, 3, 4, 5]:
            raise ValueError(
                f"Expected five {model_id} AIDE rows for scenario {scenario}, allocation={allocation}, "
                f"utility={utility_type}; found {len(rows)}."
            )
        if (
            (rows["n_valid"].astype(float) <= 0).any()
            or not (rows["n_valid"].astype(float) == rows["ntrial_from_files"].astype(float)).all()
        ):
            raise ValueError(
                f"Invalid replicate count for {model_id}, scenario {scenario}, allocation={allocation}, "
                f"utility={utility_type}."
            )
        output[scenario] = {
            "treated": rows["Pts_Treated"].astype(float).tolist(),
            "selected": rows["OBD_Selection_pct"].astype(float).tolist(),
            "stop": float(rows["No_OBD_Selection_pct"].iloc[0]),
        }
    return output


def draw_method_table_scenario(
    pdf: canvas.Canvas,
    y: float,
    scenario: int,
    truth: dict[str, list[float]],
    methods: list[tuple[str, dict[str, Any]]],
    *,
    compact: bool,
) -> float:
    scenario_height = 14 if compact else 18
    pdf.setFillColor(colors.HexColor("#EFEFEF"))
    pdf.rect(LEFT, y - scenario_height + 2, TABLE_WIDTH, scenario_height, fill=1, stroke=0)
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 10.5 if compact else 13)
    pdf.drawString(LEFT + 8, y - (9 if compact else 12), f"Scenario {scenario}")
    y -= 15 if compact else 19
    y = draw_table_row(pdf, y, "DLT rate", truth["dlt"], None, decimals=2, compact=compact)
    y = draw_table_row(pdf, y, "Efficacy rate", truth["efficacy"], None, decimals=2, compact=compact)
    y = draw_table_row(pdf, y, "Utility 2", truth["utility2"], None, compact=compact)
    y = draw_table_row(pdf, y, "Utility 3", truth["utility3"], None, compact=compact)
    y -= 1 if compact else 2

    y = draw_group_header(pdf, y, "OBD selection (%)", compact=compact)
    for method_label, record in methods:
        y = draw_table_row(
            pdf, y, method_label, record["selected"], record["stop"], compact=compact,
        )

    y = draw_group_header(pdf, y, "Mean patients treated", compact=compact)
    for method_label, record in methods:
        y = draw_table_row(pdf, y, method_label, record["treated"], None, compact=compact)
    return y - (3 if compact else 9)


def write_refreshed_n30_table() -> Path:
    """Build the N=30, Cycle_Max=1 all-method table using the random-carryover AIDE results."""
    truth = load_truth(JULY_RAW_DIR / "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv")
    if not ALL_METHODS_CYCLE1_CSV.exists():
        raise FileNotFoundError(ALL_METHODS_CYCLE1_CSV)
    aide_data = pd.read_csv(ALL_METHODS_CYCLE1_CSV)
    validate_cycle1_all_methods_source(aide_data)
    boin = parse_operating_characteristics_pdf(JULY_RAW_DIR / "BOIN12 Result N = 30.pdf")
    uboin: dict[int, dict[str, Any]] = {}
    for label, count, offset in (("Sce1-10", 10, 0), ("Sce11-20", 10, 10), ("Sce21-30", 10, 20), ("Sce31-37", 7, 30)):
        source = JULY_RAW_DIR / f"U-BOIN {label} N = 30 S1 = 9.pdf"
        parsed = parse_operating_characteristics_pdf(source, scenario_offset=offset)
        if len(parsed) != count:
            raise ValueError(f"Expected {count} scenarios in {source.name}; found {len(parsed)}.")
        uboin.update(parsed)
    efftox = parse_efftox_html(JULY_RAW_DIR / "EffTox N = 30.html")
    expected_scenarios = set(range(1, 38))
    for label, records in (("BOIN12", boin), ("U-BOIN", uboin), ("EffTox", efftox)):
        if set(records) != expected_scenarios:
            raise ValueError(f"{label} N=30 raw source is incomplete: found {sorted(records)}")

    sources = {
        "BOIN12": boin,
        "U-BOIN": uboin,
        "EffTox": efftox,
        "AIDE one_stage (U2)": extract_cycle1_aide_records(
            aide_data, "random_carryoever", "one_stage", 2, 0.3
        ),
        "AIDE one_stage (U3)": extract_cycle1_aide_records(
            aide_data, "random_carryoever", "one_stage", 3, 1.0
        ),
        "AIDE two_stage (U2)": extract_cycle1_aide_records(
            aide_data, "random_carryoever", "two_stage", 2, 0.3
        ),
        "AIDE two_stage (U3)": extract_cycle1_aide_records(
            aide_data, "random_carryoever", "two_stage", 3, 1.0
        ),
    }
    output_dir = N30_COMPARISON_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / "phase12_all_methods_N30_scenarios_1_to_37_tables.pdf"
    # Seven methods are displayed in each operating-characteristic block;
    # two scenarios per page preserves the source table's readable font size.
    groups = scenario_groups(2)
    pdf, temporary_path = make_pdf_canvas(output)
    method_order = tuple(sources)
    for page_number, scenarios in enumerate(groups, start=1):
        draw_table_header(
            pdf,
            "Phase I/II Operating Characteristics",
            "Method comparison - N = 30, Cycle_Max = 1",
            page_number,
            len(groups),
            "AIDE rows use random CRM toxicity carryover plus dose-specific efficacy carryover. "
            "U2: Type = 2, Lambda_T = 0.3; U3: Type = 3, Lambda_T = 1.",
        )
        y = PAGE_HEIGHT - 130
        compact = False
        for scenario in scenarios:
            y = draw_method_table_scenario(
                pdf,
                y,
                scenario,
                truth[scenario],
                [(method, sources[method][scenario]) for method in method_order],
                compact=compact,
            )
        pdf.showPage()
    finalize_pdf(pdf, temporary_path, output)

    (output_dir / "README.txt").write_text(
        "N=30, Cycle_Max=1 table refreshed from the newest August random-carryover AIDE summary and the July raw comparison results.\n\n"
        "AIDE rows shown: one_stage and two_stage for Utilities 2 and 3.\n\n"
        "AIDE source:\n"
        f"- Presentation 8-03-2026/{ALL_METHODS_CYCLE1_CSV.name}\n"
        "- AIDE values use the random-carryover model: random CRM toxicity carryover plus dose-specific efficacy carryover.\n"
        "- Every AIDE setting has Cycle_Max=1, 1,000 trials, and continuous enrollment.\n"
        "- One-stage and two-stage designs are shown for Utilities 2 and 3.\n\n"
        "Raw July sources retained for the comparator rows:\n"
        "- BOIN12 Result N = 30.pdf\n"
        "- U-BOIN Sce1-10 / Sce11-20 / Sce21-30 / Sce31-37 N = 30 S1 = 9.pdf\n"
        "- EffTox N = 30.html\n\n"
        "The table covers scenarios 1-37, two per page. Each scenario groups OBD selection rows together, followed by mean-patient-allocation rows. Stop % is no-dose-selection probability. "
        "U2 means Utility_Type=2 and Lambda_T=0.3; U3 means Utility_Type=3 and Lambda_T=1.\n\n"
        "IPDE-alpha comparison tables:\n"
        "- phase12_ipde_alpha_comparison_one_stage_scenarios_1_to_37_tables.pdf\n"
        "- phase12_ipde_alpha_comparison_two_stage_scenarios_1_to_37_tables.pdf\n\n"
        "These two tables use the additive previous-dose source below, fix both the CRM toxicity and additive-efficacy "
        "priors at Beta(0.15, 0.85), and show U2 and U3 together for IPDE alpha values 0, 0.3, 0.6, and 0.9:\n"
        f"- Presentation 8-03-2026/{COMBINED_IPDE_ALPHA_CSV.name}\n\n"
        "Each comparison table covers scenarios 1-37, two per page. U2 and U3 rows appear within each OBD-selection, "
        "mean-treated, and mean-IPDE-patient block.\n\n"
        "The prior July N=60 table is retained unchanged.\n",
        encoding="utf-8",
    )
    return output


def build_prior_outputs(run: PriorRun) -> list[Path]:
    if not run.summary_csv.exists():
        raise FileNotFoundError(run.summary_csv)
    run.output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(run.summary_csv)
    validate_prior_source(data, run)
    truth = load_truth()
    table_outputs: list[Path] = []
    for alpha in ALPHAS:
        for design in DESIGNS:
            table_outputs.append(
                write_prior_table(run, alpha, design, truth, extract_prior_records(data, run, alpha, design))
            )
    alpha_comparison_outputs = [
        write_ipde_alpha_table(run, design, truth, extract_ipde_alpha_records(data, design))
        for design in DESIGNS
    ]
    write_prior_readme(run, table_outputs, alpha_comparison_outputs)
    return [*table_outputs, *alpha_comparison_outputs]


def main() -> None:
    outputs = [*build_prior_outputs(TOXICITY_RUN), *build_prior_outputs(EFFICACY_RUN)]
    outputs.append(write_refreshed_n30_table())
    print("\n".join(str(path) for path in outputs))


if __name__ == "__main__":
    main()
