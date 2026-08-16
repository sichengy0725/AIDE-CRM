from pathlib import Path
import re

from PIL import Image


root = Path("tmp/pdfs/render_tite_56")
out = root / "contacts"
out.mkdir(exist_ok=True)


def page_number(path: Path) -> int:
    return int(re.search(r"(\d+)$", path.stem).group(1))


for folder in sorted(path for path in root.iterdir() if path.is_dir() and path.name != "contacts"):
    pages = sorted(folder.glob("page-*.png"), key=page_number)
    for start in range(0, len(pages), 20):
        sheet = Image.new("RGB", (1296, 1260), "white")
        for index, page in enumerate(pages[start:start + 20]):
            image = Image.open(page).convert("RGB").resize((324, 252))
            sheet.paste(image, ((index % 4) * 324, (index // 4) * 252))
        sheet.save(out / f"{folder.name}-{start // 20 + 1}.png")

print(f"Created {len(list(out.glob('*.png')))} contact sheets.")
