"""Copy the July IPDE design PDF and append the August additive-model section."""

from __future__ import annotations

import io
import shutil
from pathlib import Path
from textwrap import wrap

from pypdf import PdfReader, PdfWriter
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
JULY_PDF = ROOT / "Presentation 7-27-2026" / "Current_IPDE_Design_with_Efficacy_Revised.pdf"
OUTPUT_PDF = ROOT / "Presentation 8-03-2026" / "Current_IPDE_Design_with_Efficacy_Revised.pdf"
TMP_DIR = ROOT / "tmp" / "pdfs" / "current_ipde_design_additive"
PAGE_WIDTH, PAGE_HEIGHT = letter
LEFT = 54
RIGHT = PAGE_WIDTH - LEFT
TOTAL_PAGES = 11


def draw_footer(pdf: canvas.Canvas, page_number: int) -> None:
    pdf.setFillColor(colors.HexColor("#333333"))
    pdf.setFont("Times-Roman", 8)
    pdf.drawString(LEFT, 31, "Technical design specification | August 2026")
    pdf.drawRightString(RIGHT, 31, f"Page {page_number} of {TOTAL_PAGES}")


def footer_overlay(page_number: int) -> PdfReader:
    """Cover the inherited footer so the copied document has consistent pagination."""
    stream = io.BytesIO()
    pdf = canvas.Canvas(stream, pagesize=letter)
    pdf.setFillColor(colors.white)
    pdf.rect(0, 0, PAGE_WIDTH, 52, fill=1, stroke=0)
    draw_footer(pdf, page_number)
    pdf.save()
    stream.seek(0)
    return PdfReader(stream)


def draw_chrome(pdf: canvas.Canvas, page_number: int) -> None:
    pdf.setFillColor(colors.black)
    pdf.setFont("Times-Bold", 14)
    pdf.drawString(LEFT, 752, "CURRENT IPDE DESIGN WITH EFFICACY")
    pdf.setStrokeColor(colors.HexColor("#888888"))
    pdf.setLineWidth(0.45)
    pdf.line(LEFT, 744, RIGHT, 744)
    draw_footer(pdf, page_number)


def draw_wrapped(
    pdf: canvas.Canvas,
    text: str,
    x: float,
    y: float,
    *,
    width_chars: int = 96,
    font: str = "Times-Roman",
    size: float = 10,
    leading: float = 13,
) -> float:
    pdf.setFillColor(colors.black)
    pdf.setFont(font, size)
    for line in wrap(text, width=width_chars, break_long_words=False, break_on_hyphens=False):
        pdf.drawString(x, y, line)
        y -= leading
    return y


def draw_heading(pdf: canvas.Canvas, text: str, y: float, *, size: float = 13) -> float:
    pdf.setFillColor(colors.black)
    pdf.setFont("Times-Bold", size)
    pdf.drawString(LEFT, y, text)
    return y - (size + 7)


def draw_formula(pdf: canvas.Canvas, text: str, y: float) -> float:
    pdf.setFillColor(colors.black)
    pdf.setFont("Courier", 9.8)
    pdf.drawCentredString(PAGE_WIDTH / 2, y, text)
    return y - 18


def draw_bullet(pdf: canvas.Canvas, text: str, y: float) -> float:
    pdf.setFillColor(colors.black)
    pdf.setFont("Times-Roman", 10)
    lines = wrap(text, width=88, break_long_words=False, break_on_hyphens=False)
    pdf.drawString(LEFT + 10, y, "- " + lines[0])
    for line in lines[1:]:
        y -= 13
        pdf.drawString(LEFT + 22, y, line)
    return y - 17


def draw_additive_pages(path: Path) -> None:
    pdf = canvas.Canvas(str(path), pagesize=letter)

    draw_chrome(pdf, 10)
    y = 711
    y = draw_heading(pdf, "6. Additive Previous-Dose Outcome Models", y, size=18)
    y = draw_wrapped(
        pdf,
        "The current implementation adds the immediately preceding dose contribution directly to the "
        "regular-patient probability for an IPDE administration. This section applies when the "
        "previous_dose_additive toxicity and efficacy models are selected; it replaces the random "
        "carryover specifications in Sections 1.1-1.3 for that analysis.",
        LEFT,
        y,
        width_chars=95,
    ) - 9

    y = draw_heading(pdf, "6.1 Toxicity: additive previous-dose CRM", y)
    y = draw_wrapped(
        pdf,
        "Let j be the current dose and l be the same patient's immediately preceding dose. The regular "
        "toxicity model remains the power CRM p_T,j = p_0,j ^ exp(beta), with beta given its CRM prior. "
        "For a regular administration, Y_T,i follows Bernoulli(p_T,j).",
        LEFT,
        y,
        width_chars=95,
    ) - 3
    y = draw_formula(pdf, "For an IPDE administration:  q_T,i = min{1, p_T,j + alpha_T * p_T,l}", y)
    y = draw_wrapped(
        pdf,
        "The shared additive coefficient alpha_T has a Beta(a_T, b_T) prior, implemented as a ratio "
        "of independent Gamma variables. Thus, alpha_T = 0 makes the IPDE toxicity probability equal "
        "to the regular toxicity probability at the current dose. Larger values add the regular "
        "toxicity risk at the patient's prior dose, with the probability capped at one.",
        LEFT,
        y,
        width_chars=95,
    ) - 9

    y = draw_heading(pdf, "6.2 Efficacy: additive previous-dose beta-binomial model", y)
    y = draw_wrapped(
        pdf,
        "Each dose has a regular-patient efficacy probability p_E,j with its dose-specific Beta(a_E,j, "
        "b_E,j) prior. For a regular administration, Y_E,i follows Bernoulli(p_E,j). The additive "
        "model assigns a separate shared efficacy coefficient alpha_E to IPDE administrations.",
        LEFT,
        y,
        width_chars=95,
    ) - 3
    y = draw_formula(pdf, "For an IPDE administration:  q_E,i = min{1, p_E,j + alpha_E * p_E,l}", y)
    draw_wrapped(
        pdf,
        "Here alpha_E has its own Beta(a_E, b_E) prior and is distinct from alpha_T. The efficacy "
        "model uses the same current/previous dose pair as the toxicity model, rather than assigning a "
        "separate carryover parameter to every destination dose.",
        LEFT,
        y,
        width_chars=95,
    )
    pdf.showPage()

    draw_chrome(pdf, 11)
    y = 711
    y = draw_heading(pdf, "6.3 Likelihood, posterior use, and operating rules", y, size=18)
    y = draw_wrapped(
        pdf,
        "The likelihood is evaluated at the individual-administration level. Each IPDE row retains its "
        "current dose j, its preceding dose l, and its IPDE indicator. This is essential because two "
        "IPDE patients treated at the same destination dose can have different preceding doses and, "
        "therefore, different additive probabilities.",
        LEFT,
        y,
        width_chars=95,
    ) - 8

    y = draw_heading(pdf, "Posterior quantities", y)
    y = draw_bullet(
        pdf,
        "The toxicity fit estimates beta, alpha_T, and the regular-dose probabilities p_T,j from all regular and IPDE toxicity outcomes.",
        y,
    )
    y = draw_bullet(
        pdf,
        "The efficacy fit estimates alpha_E and every regular-dose probability p_E,j from all regular and IPDE efficacy outcomes.",
        y,
    )
    y = draw_bullet(
        pdf,
        "For a candidate IPDE treatment, q_T,i and q_E,i are recomputed from posterior draws using that patient's actual pair (l, j).",
        y,
    ) - 4

    y = draw_heading(pdf, "Decision interpretation", y)
    y = draw_bullet(
        pdf,
        "Use p_T,j for dose-level safety, MTD selection, and utility calculations; use q_T,i for the individual IPDE safety criterion.",
        y,
    )
    y = draw_bullet(
        pdf,
        "Use p_E,j for efficacy ranking, futility monitoring, and utility calculations; use q_E,i for the expected efficacy improvement required for recycling.",
        y,
    )
    y = draw_bullet(
        pdf,
        "The one-stage and two-stage allocation rules, enrollment schemes, overtoxicity rule, and futility rule are otherwise unchanged.",
        y,
    ) - 4

    y = draw_heading(pdf, "Implementation note", y)
    draw_wrapped(
        pdf,
        "The additive coefficient can be fixed or estimated with a prespecified Beta prior. Toxicity and "
        "efficacy coefficients may be assigned separately. The min{1, ...} cap is applied to every IPDE "
        "probability so the additive contribution remains a valid probability.",
        LEFT,
        y,
        width_chars=95,
    )
    pdf.save()


def main() -> None:
    if not JULY_PDF.exists():
        raise FileNotFoundError(JULY_PDF)
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    addendum_path = TMP_DIR / "additive_previous_dose_models.pdf"
    temporary_output = TMP_DIR / OUTPUT_PDF.name
    draw_additive_pages(addendum_path)

    source = PdfReader(JULY_PDF)
    addendum = PdfReader(addendum_path)
    if len(source.pages) != 9 or len(addendum.pages) != 2:
        raise ValueError("Expected a nine-page source PDF and a two-page additive-model section.")

    writer = PdfWriter()
    for page_number, page in enumerate(source.pages, start=1):
        if page_number > 1:
            page.merge_page(footer_overlay(page_number).pages[0])
        writer.add_page(page)
    for page in addendum.pages:
        writer.add_page(page)
    writer.add_metadata(
        {
            "/Title": "Current IPDE Design with Efficacy - Additive Previous-Dose Models",
            "/Subject": "Technical design specification updated with additive toxicity and efficacy models",
        }
    )
    with temporary_output.open("wb") as stream:
        writer.write(stream)

    OUTPUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(temporary_output, OUTPUT_PDF)
    print(OUTPUT_PDF)


if __name__ == "__main__":
    main()
