# Original analysis programs

`code/` preserves the original R and Python analysis scripts and their
extension directory. These are the scientific analysis provenance. The
portable, tested figure-replot entry point is `../scripts/replot.py`.

The full original pipeline requires the original project directory layout,
large processed expression inputs, saved pseudobulk matrices and the sibling
`AKI_DCT_analysis` resources. Those large inputs and historical R environments
are not included here. Do not run this archive as though it were a fully
provisioned independent pipeline.

| Entry | Responsibility |
|---|---|
| `01_bulk_validation.R` | Cross-sectional tissue support |
| `02_primary_cell_resolved.R` | Donor-level TAL programme analyses |
| `03_longitudinal_validation.R` | Paired transplant reperfusion |
| `04_tubular_coordination.R` | Cross-tubular identity context |
| `05_mitochondrial_audit.R` | Documented negative mitochondrial audit |
| `06_evidence_freeze_check.R` | Frozen numerical anchors |
| `extension/` | Atlas, pseudobulk, state, external-response and sensitivity programs |

Packages referenced by these R scripts include Matrix, ggplot2, tidyverse,
edgeR, limma, digest, data.table, Seurat, UCell, jsonlite, patchwork, scales
and gridExtra. This is an inventory, not a newly reconstructed lockfile.
Full model reruns were not performed for the poster update.
