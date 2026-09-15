"""Extract the paper evidence figures used by the standalone Level 5 report."""
from pathlib import Path

import pypdfium2 as pdfium


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / ".build" / "paper"
OUT.mkdir(parents=True, exist_ok=True)
doc = pdfium.PdfDocument(ROOT / "2403.12021v4.pdf")
# Rectangles refer to the 1.7x rendered pages. Preserve original scientific art.
for page, name, rect in (
    (6, "fig5.png", (534, 84, 983, 391)),
    (7, "fig6abc.png", (159, 85, 657, 422)),
    (19, "ext10e.png", (146, 552, 675, 880)),
):
    image = doc[page - 1].render(scale=1.7).to_pil()
    image.crop(rect).save(OUT / name)
