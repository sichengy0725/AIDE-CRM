"""Build July-27-style Phase I/II tables for the August toxicity-prior run."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
RESULTS_CSV = ROOT / "Presentation 8-03-2026" / "AIDE_phase_I_II_IDX_0001_to_1000_dose_summary.csv"
TRUTH_CSV = ROOT / "Set_5dose_adaptive_r_37_true_MTD_OBD_summary_lambda1.csv"
OUTPUT_DIR = ROOT / "Presentation 8-03-2026" / "Table and Plots" / "Toxicity Prior"

PAGE_SIZE = (18 * inch, 14 * inch)
PAGE_WIDTH, PAGE_HEIGHT = PAGE_SIZE
LEFT = 80
LABEL_WIDTH = 330
DOSE_WIDTH = 130
STOP_WIDTH = 156
TABLE_WIDTH = LABEL_WIDTH + 5 * DOSE_WIDTH + STOP_WIDTH
PRIORS = [(0.15, 0.85), (0.30, 0.70), (0.50, 0.50), (1.00, 1.00)]
ALPHAS = (0.0, 0.3, 0.6, 0.9)
DESIGNS = (
    {
        "allocation": "two_stage",
        "n_s1": 6,
        "label": "two-stage",
        "filename_prefix": "phase12_toxicity_prior",
    },
    {
        "allocation": "one_stage",
        "n_s1": 30,
        "label": "one-stage",
        "filename_prefix": "phase12_one_stage_toxicity_prior",
    },
)


def display(value: float, decimals: int = 1) -> str:
    if pd.isna(value):
        return "-"
    rounded = round(float(value), decimals)
    if rounded == 0:
        rounded = 0.0
    text = f"{rounded:.{decimals}f}"
    return text.rstrip("0").rstrip(".")


def alpha_label(alpha: float) -> str:
    return display(alpha, decimals=1)


def prior_label(a: float, b: float) -> str:
    return f"Toxicity prior Beta({display(a, 2)}, {display(b, 2)})"


def close_to(series: pd.Series, value: float) -> pd.Series:
    return (series.astype(float) - value).abs() < 1e-10


def load_truth() -> dict[int, dict[str, list[float]]]:
    truth_df = pd.read_csv(TRUTH_CSV)
    if sorted(truth_df["Scenario"].unique()) != list(range(1, 38)):
        raise ValueError("Truth file must contain scenarios 1 through 37 exactly once.")

    truth: dict[int, dict[str, list[float]]] = {}
    for _, row in truth_df.iterrows():
        scenario = int(row["Scenario"])
        truth[scenario] = {
            "dlt": [float(row[f"Tox_Dose{dose}"]) for dose in range(1, 6)],
            "efficacy": [float(row[f"Eff_Dose{dose}"]) for dose in range(1, 6)],
            "utility2": [float(row[f"Utility2_Dose{dose}"]) for dose in range(1, 6)],
            "utility3": [float(row[f"Utility3_Dose{dose}"]) for dose in range(1, 6)],
        }
    return truth


def extract_records(
    data: pd.DataFrame, alpha: float, design: dict[str, object]
) -> dict[int, dict[tuple[float, float], dict[str, object]]]:
    subset = data[
        close_to(data["Toxicity_IPDE_Alpha"], alpha)
        & close_to(data["Efficacy_IPDE_Alpha"], alpha)
        & (data["Allocation"] == design["allocation"])
        & (data["Nmax"] == 30)
        & (data["N_s1"] == design["n_s1"])
        & (data["Utility_Type"] == 2)
        & close_to(data["Lambda_T"], 0.3)
    ].copy()

    output: dict[int, dict[tuple[float, float], dict[str, object]]] = {}
    for scenario in range(1, 38):
        output[scenario] = {}
        for prior_a, prior_b in PRIORS:
            rows = subset[
                (subset["Scenario"] == scenario)
                & close_to(subset["CRM_Prior_a"], prior_a)
                & close_to(subset["CRM_Prior_b"], prior_b)
            ].sort_values("Dose")
            if len(rows) != 5 or rows["Dose"].tolist() != [1, 2, 3, 4, 5]:
                raise ValueError(
                    f"Expected five ordered rows for scenario {scenario}, alpha {alpha}, "
                    f"and prior Beta({prior_a}, {prior_b}); found {len(rows)}."
                )
            if not (rows["n_valid"] == rows["ntrial_from_files"]).all() or (rows["n_valid"] <= 0).any():
                raise ValueError(
                    f"Scenario {scenario}, alpha {alpha}, prior Beta({prior_a}, {prior_b}) "
                    "has an invalid successful-replicate count."
                )
            output[scenario][(prior_a, prior_b)] = {
                "treated": rows["Pts_Treated"].astype(float).tolist(),
                "selected": rows["OBD_Selection_pct"].astype(float).tolist(),
                "stop": float(rows["No_OBD_Selection_pct"].iloc[0]),
            }
    return output


def draw_header(
    pdf: canvas.Canvas,
    alpha: float,
    design: dict[str, object],
    page_number: int,
    page_count: int,
) -> None:
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 28)
    pdf.drawCentredString(PAGE_WIDTH / 2, PAGE_HEIGHT - 56, "Phase I/II Operating Characteristics")
    pdf.setFont("Helvetica", 14)
    pdf.drawCentredString(
        PAGE_WIDTH / 2,
        PAGE_HEIGHT - 78,
        f"Toxicity-prior comparison - {design['label']} - IPDE alpha = {alpha_label(alpha)} - N = 30",
    )
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

    pdf.setFont("Helvetica", 8.5)
    footer = (
        f"AIDE Phase I/II: {design['label']}, Utility 2 (Lambda_T = 0.3). "
        "Selection % is OBD selection; Stop % is no OBD selection. "
        "Toxicity and efficacy IPDE alphas are both set to the listed alpha; available replicates vary by setting."
    )
    pdf.drawString(LEFT, 38, footer)
    pdf.drawRightString(PAGE_WIDTH - LEFT, 38, f"Page {page_number} of {page_count}")


def draw_row(
    pdf: canvas.Canvas,
    y: float,
    label: str,
    values: list[float],
    stop: float | None,
    *,
    shaded: bool = False,
    decimals: int = 1,
    compact: bool = False,
) -> float:
    row_height = 9.8 if compact else 13.2
    font_size = 7.4 if compact else 9.5
    baseline = 7.2 if compact else 9.6
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


def draw_group_header(pdf: canvas.Canvas, y: float, label: str, *, compact: bool = False) -> float:
    row_height = 9.8 if compact else 13.2
    font_size = 7.6 if compact else 9.6
    baseline = 7.2 if compact else 9.6
    pdf.setFillColor(colors.HexColor("#E7E7E7"))
    pdf.rect(LEFT, y - row_height + 2, TABLE_WIDTH, row_height, fill=1, stroke=0)
    pdf.setFillColor(colors.HexColor("#202020"))
    pdf.setFont("Helvetica-Bold", font_size)
    pdf.drawString(LEFT + 8, y - baseline, label)
    pdf.setStrokeColor(colors.HexColor("#D0D0D0"))
    pdf.setLineWidth(0.45)
    pdf.line(LEFT, y - row_height, LEFT + TABLE_WIDTH, y - row_height)
    return y - row_height


def draw_scenario(
    pdf: canvas.Canvas,
    y: float,
    scenario: int,
    truth: dict[str, list[float]],
    records: dict[tuple[float, float], dict[str, object]],
    *,
    compact: bool = False,
) -> float:
    scenario_height = 14 if compact else 18
    pdf.setFillColor(colors.HexColor("#EFEFEF"))
    pdf.rect(LEFT, y - scenario_height + 2, TABLE_WIDTH, scenario_height, fill=1, stroke=0)
    pdf.setFillColor(colors.black)
    pdf.setFont("Helvetica-Bold", 10.5 if compact else 13)
    pdf.drawString(LEFT + 8, y - (9 if compact else 12), f"Scenario {scenario}")
    y -= 15 if compact else 19

    y = draw_row(pdf, y, "DLT rate", truth["dlt"], None, decimals=2, compact=compact)
    y = draw_row(pdf, y, "Efficacy rate", truth["efficacy"], None, decimals=2, compact=compact)
    y = draw_row(pdf, y, "Utility 2", truth["utility2"], None, compact=compact)
    y = draw_row(pdf, y, "Utility 3", truth["utility3"], None, compact=compact)
    y -= 1 if compact else 2

    y = draw_group_header(pdf, y, "OBD selection (%)", compact=compact)
    for prior_a, prior_b in PRIORS:
        record = records[(prior_a, prior_b)]
        label = prior_label(prior_a, prior_b)
        y = draw_row(
            pdf,
            y,
            label,
            record["selected"],
            float(record["stop"]),
            shaded=True,
            compact=compact,
        )

    y = draw_group_header(pdf, y, "Mean patients treated", compact=compact)
    for prior_a, prior_b in PRIORS:
        record = records[(prior_a, prior_b)]
        label = prior_label(prior_a, prior_b)
        y = draw_row(
            pdf,
            y,
            label,
            record["treated"],
            None,
            compact=compact,
        )
    return y - (3 if compact else 9)


def scenario_groups() -> list[list[int]]:
    groups = [list(range(start, start + 3)) for start in range(1, 34, 3)]
    groups.append([34, 35, 36, 37])
    return groups


def write_pdf(
    alpha: float,
    design: dict[str, object],
    truth: dict[int, dict[str, list[float]]],
    records: dict[int, dict[tuple[float, float], dict[str, object]]],
) -> Path:
    groups = scenario_groups()
    filename = (
        f"{design['filename_prefix']}_alpha{alpha_label(alpha).replace('.', 'p')}_"
        "scenarios_1_to_37_tables.pdf"
    )
    output = OUTPUT_DIR / filename
    pdf = canvas.Canvas(str(output), pagesize=PAGE_SIZE)
    for page_number, scenarios in enumerate(groups, start=1):
        draw_header(pdf, alpha, design, page_number, len(groups))
        y = PAGE_HEIGHT - 130
        compact = len(scenarios) == 4
        for scenario in scenarios:
            y = draw_scenario(pdf, y, scenario, truth[scenario], records[scenario], compact=compact)
        pdf.showPage()
    pdf.save()
    return output


def write_readme(outputs: dict[str, list[Path]]) -> None:
    lines = [
        "Generated 2026-07-30 from Presentation 8-03-2026/AIDE_phase_I_II_IDX_0001_to_1000_dose_summary.csv.",
        "",
        "Eight 12-page landscape PDFs cover scenarios 1-37 at IPDE alpha 0, 0.3, 0.6, and 0.9.",
        "Each table compares toxicity priors Beta(0.15,0.85), Beta(0.30,0.70), Beta(0.50,0.50), and Beta(1,1).",
        "Separate files are provided for AIDE Phase I/II two-stage and one-stage allocation, N=30, Utility Type 2, and Lambda_T=0.3.",
        "Within every scenario, all OBD-selection rows are grouped together, followed by all mean patient-allocation rows.",
        "Selection percentages are OBD-selection probabilities; Stop % is No_OBD_Selection_pct.",
        "The source data pairs toxicity and efficacy IPDE alpha at the listed alpha value.",
        "Successful-replicate counts vary by scenario and setting (811 to 1,000); values use the source denominators.",
        "",
        "Files:",
    ]
    for design_label, paths in outputs.items():
        lines.append(f"- {design_label}:")
        lines.extend(f"  - {path.name}" for path in paths)
    (OUTPUT_DIR / "README.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(RESULTS_CSV)
    truth = load_truth()
    outputs: dict[str, list[Path]] = {}
    for design in DESIGNS:
        label = str(design["label"])
        outputs[label] = [
            write_pdf(alpha, design, truth, extract_records(data, alpha, design))
            for alpha in ALPHAS
        ]
    write_readme(outputs)
    print("\n".join(str(path) for paths in outputs.values() for path in paths))


if __name__ == "__main__":
    main()
