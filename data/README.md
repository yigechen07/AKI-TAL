# Required local plotting inputs

The result tables and atlas coordinates are not distributed in this repository.
Pass the original local analysis directory to `scripts/replot.py --project`.
It must contain:

```text
results/formal_evidence.csv
results/formal_effects.csv
results/formal_genes.csv
results/bulk_effects.csv
results/bulk_genes.csv
results/extension/atlas_GSE183276_umap.csv
results/extension/atlas_GSE183276_markers.csv
results/extension/atlas_GSE183276_counts.csv
results/extension/atlas_GSE210622_umap.csv
results/extension/atlas_GSE210622_markers.csv
results/extension/atlas_GSE210622_counts.csv
```

`formal_*` and `bulk_*` are the saved numerical inputs for Figure 2.
`atlas_*` contains saved coordinates, annotation-marker summaries and
recovered-cell counts for Figure 1.

Public cell-resolved cohorts: GSE183276 and GSE210622. Cross-sectional tissue
support: GSE139061 and GSE30718. The poster additionally includes the existing
GSE126805 paired-allograft results. Accession IDs identify the public source
studies; the required local tables are processed research outputs.

The atlas coordinates are not recomputed, and cells are not independent
patient replicates. Recovered-cell counts describe sampling coverage and
must not be interpreted as biological TAL abundance.
