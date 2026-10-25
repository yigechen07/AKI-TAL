#!/usr/bin/env python3
"""Assemble journal-facing main and supplementary figure files.

The source figures remain unchanged. Figure 3A-D and Figure 3E are combined
into one submission asset. The exploratory Figure 4F bridge is copied to the
supplement as Figure S2; Figure S1 is the submission sensitivity audit.
"""

from pathlib import Path
from shutil import copy2

from PIL import Image, ImageOps


PROJECT = Path(__file__).resolve().parents[2]
EXT = PROJECT / "figures" / "extension"
MAIN = PROJECT / "figures" / "submission"
SUPP = PROJECT / "figures" / "supplement"
MAIN.mkdir(parents=True, exist_ok=True)
SUPP.mkdir(parents=True, exist_ok=True)


def rgb(path: Path) -> Image.Image:
    return Image.open(path).convert("RGB")


top = rgb(EXT / "figure3_tal_pseudobulk.png")
bottom = rgb(EXT / "figure3e_consensus_external_validation.png")
if bottom.width != top.width:
    new_h = round(bottom.height * top.width / bottom.width)
    bottom = bottom.resize((top.width, new_h), Image.Resampling.LANCZOS)

gap = 36
canvas = Image.new("RGB", (top.width, top.height + gap + bottom.height), "white")
canvas.paste(top, (0, 0))
canvas.paste(bottom, (0, top.height + gap))
canvas = ImageOps.expand(canvas, border=(0, 0, 0, 0), fill="white")
canvas.save(MAIN / "figure3_complete.png", dpi=(300, 300), optimize=True)
canvas.save(MAIN / "figure3_complete.pdf", resolution=300.0)

for number, stem in (
    (1, "figure1_atlas"),
    (2, "figure2_frozen"),
    (4, "figure4_tal_states"),
    (5, "figure5_nephron_context"),
    (6, "figure6_paired_iri"),
):
    for suffix in (".png", ".pdf"):
        copy2(EXT / f"{stem}{suffix}", MAIN / f"figure{number}{suffix}")

for suffix in (".png", ".pdf"):
    copy2(EXT / f"figure4f_state_programme_bridge{suffix}",
          SUPP / f"figureS2_state_programme_bridge{suffix}")

print("Submission figure assembly complete")
