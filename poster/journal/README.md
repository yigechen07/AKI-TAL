# Journal-style poster (A0 landscape)

- `AKI_TAL_poster_journal.pdf` — print file, 1189 × 841 mm, fonts embedded.
- `AKI_TAL_poster_journal.html` — source; uses the macOS fonts Avenir Next and Charter.
- `AKI_TAL_poster_journal_preview.png` — preview image.
- `figures/` — poster Figures 1–4 (files `fig2`–`fig5`, named after the manuscript figures they replace), drawn by `scripts/poster_redesign/make_poster_figures.py`.

The figure script and `scripts/poster_redesign/check_poster_numbers.py` read the
result tables of the local analysis directory, which are not in this repository
(see `data/README.md`). Run both from that directory:

```bash
python3 -I scripts/poster_redesign/make_poster_figures.py /path/to/AKI_TAL_severity_analysis
python3 -I scripts/poster_redesign/check_poster_numbers.py /path/to/AKI_TAL_severity_analysis
```

To export the PDF from the HTML, print it with Chrome (no headers/footers, page
size taken from the file).
