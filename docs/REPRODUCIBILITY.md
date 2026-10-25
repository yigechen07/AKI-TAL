# Figures and poster

The poster uses the approved BioRender schematics in Figure 1a (study design)
and Figure 2a (TAL programme). Numerical plots use frozen results. Biological
replicates are donors, samples or paired kidney allografts.

## Files

- `poster/AKI_TAL_poster.pptx`: updated original landscape poster.
- `poster/Poster.png`: README preview of the final poster.
- `figures/`: updated Figures 1 and 2 in PDF, SVG and 400 dpi PNG.
- `assets/biorender/`: the two approved schematic images at source resolution.
- `data/README.md`: required input-file inventory. Plotting tables and atlas coordinates remain local and are not distributed in this repository.
- `scripts/`: portable figure renderer and poster-image replacement script.
- `analysis/code/`: original scientific analysis scripts, retained for reference.
- `SHA256SUMS.txt`: checksums of the shared files.

## Replot Figure 1 and Figure 2

From the repository directory, using Python 3.12:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/replot.py --project /path/to/local/AKI_TAL_severity_analysis
```

Outputs go to `rebuilt/figures`. Supply the local original analysis directory
containing the required saved result tables, listed in `data/README.md`.
The command does not refit models. Arial is preferred, with DejaVu Sans as a fallback.
The BioRender artwork is embedded as an image. Other figure labels and plots
are vector objects in PDF/SVG; dense UMAP points are rasterised at 400 dpi.

## Rebuild the poster

Use Node.js 20.11 or later:

```bash
npm install
npm run poster
```

This inserts the supplied figures into a copy of the original landscape
poster, producing `rebuilt/AKI_TAL_poster.pptx`. To insert freshly replotted
figures instead:

```bash
node scripts/update_poster.mjs --figures rebuilt/figures --out rebuilt/AKI_TAL_replotted.pptx
```

Only the two embedded full-figure image parts change. The script verifies
aspect ratios, preserved package parts and source integrity. Existing output
files are not overwritten.

## BioRender sources

- [Study design and evidence framework](https://app.biorender.com/illustrations/698e3f6b10213ffb2a99a3bb?slideId=160a3a84-3999-46b0-c0fa-7be5a367d5ff)
- [TAL programme schematic](https://app.biorender.com/illustrations/aaedd32e0bdab433d58b9a20?slideId=98d4bb6f-516c-6bbe-3997-4151b2d4a853)

Editing access to those source illustrations is controlled by BioRender.
Third-party illustration assets retain their respective terms.

## Scientific scope

The frozen six genes are SLC12A1, UMOD, CLDN16, KCNJ10, PTH1R and CLDN10.
The evidence supports lower TAL programme detectability across the studied
contexts. It does not establish TAL selectivity, PT-injury independence,
biological abundance change, cell-intrinsic regulation or causality.

The current cell-resolved primary comparisons use 10 AKI / 17 reference donors
in GSE183276 and 8 AKI / 4 non-AKI donors in GSE210622. Atlas coverage and
primary-analysis eligibility are different. In GSE183276 the model P value
(5.027 × 10⁻⁵) differs from the permutation P value (0.0001).

The poster's author/affiliation placeholders are retained from the source.
This repository presents a research output, without claiming publication or
conference acceptance.
