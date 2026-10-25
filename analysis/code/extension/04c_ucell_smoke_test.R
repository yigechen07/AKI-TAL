#!/usr/bin/env Rscript
# =============================================================================
# 04c_ucell_smoke_test.R — diagnosis-blinded UCell smoke test (Figure 4)
# =============================================================================
# Verifies that UCell (project-local, metadata/extension/r-lib) loads and can
# score the frozen state programs on a small strict author-annotated TAL
# subset, WITHOUT using diagnosis in any way (no diagnosis comparison, no
# outcome-informed decision).
#
# On success it writes the PASS marker:
#   metadata/extension/r-lib/UCell_SMOKE_PASS.txt
# which figure4_tal_states_scaffold.R requires before it will report READY.
#
# Usage:
#   cd code/extension && Rscript 04c_ucell_smoke_test.R
#   Optional first CLI arg: project root (default: parent of code/extension)
#
# Requirements: Seurat, Matrix, data.table, UCell (project-local r-lib).
# Memory: loads the 2.87 GB GSE183276 object; needs ~8 GB free RAM.
# =============================================================================

suppressMessages({
  library(Seurat); library(data.table); library(Matrix)
})

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1 && nzchar(cli_args[1])) {
  proj_root <- normalizePath(cli_args[1])
} else {
  proj_root <- normalizePath(file.path(getwd(), "..", ".."))
}

rl <- file.path(proj_root, "metadata", "extension", "r-lib")
.libPaths(c(rl, .libPaths()))

if (!requireNamespace("UCell", quietly = TRUE)) {
  stop("UCell not found in project-local library: ", rl,
       ". Install via BiocManager::install(\"UCell\", lib=...) first.")
}
suppressMessages(library(UCell))
cat("UCell version:", as.character(packageVersion("UCell")), "\n")

obj_path <- file.path(proj_root, "results", "extension", "atlas_GSE183276_seurat.rds")
ref_path <- file.path(proj_root, "metadata", "extension", "tal_state_reference.csv")
stopifnot(file.exists(obj_path), file.exists(ref_path))

obj <- readRDS(obj_path)
md <- as.data.table(obj@meta.data)
# IMPORTANT: rownames(md) are integer positions in this object; map via colnames()
strict_pos <- which(md$author_label == "TAL")
strict_names <- colnames(obj)[strict_pos]
cat("strict TAL cells:", length(strict_names), "\n")

set.seed(20260806)  # deterministic subset
sub_names <- sample(strict_names, min(500, length(strict_names)))

counts <- GetAssayData(obj, layer = "counts")
sub_counts <- counts[, sub_names, drop = FALSE]

ref <- fread(ref_path, na.strings = c("NA", ""))
pair <- ref[, .(has_none = any(flag == "none"),
                any_non_none = any(flag != "none"),
                any_frozen = any(exclusion_reason == "frozen_primary_endpoint_circularity")),
            by = .(state, gene)]
primary <- pair[has_none == TRUE & any_non_none == FALSE & any_frozen == FALSE, .(state, gene)]
sensitivity <- unique(ref[flag == "shared_adaptive_epithelial", .(state, gene)])

keep_genes <- intersect(unique(ref$gene), rownames(sub_counts))
sub_counts <- sub_counts[keep_genes, , drop = FALSE]
cat("subset matrix dims:", paste(dim(sub_counts), collapse = " x "), "\n")
sub_obj <- CreateSeuratObject(counts = sub_counts, min.cells = 0, min.features = 0)

sets <- list(
  reference_TAL = intersect(primary[state == "reference_TAL", gene], rownames(sub_obj)),
  adaptive_TAL = intersect(primary[state == "adaptive_TAL", gene], rownames(sub_obj)),
  degenerative_TAL = intersect(primary[state == "degenerative_TAL", gene], rownames(sub_obj)),
  adaptive_sensitivity = intersect(sensitivity[state == "adaptive_TAL", gene], rownames(sub_obj))
)
cat("genes scored:", paste(sapply(sets, length), collapse = " / "), "\n")

# Use serial (single-process) execution: the sandbox blocks BiocParallel
# socket ports (SnowParam), so force SerialParam for the smoke test.
bp <- BiocParallel::SerialParam()
scores <- AddModuleScore_UCell(sub_obj, features = sets, BPPARAM = bp)
nm <- grep("_UCell$", colnames(scores[[]]), value = TRUE)
cat("UCell columns:", paste(nm, collapse = ", "), "\n")
print(apply(scores[[nm]], 2, summary))

stopifnot(length(nm) == 4)
stopifnot(!anyNA(scores[[nm]]))

marker <- file.path(rl, "UCell_SMOKE_PASS.txt")
writeLines(c(
  "UCell smoke test PASS (diagnosis-blinded, strict-TAL subset)",
  sprintf("date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")),
  sprintf("ucell_version: %s", as.character(packageVersion("UCell"))),
  sprintf("n_cells: %d", length(sub_names)),
  sprintf("genes_scored: %s", paste(sapply(sets, length), collapse = "/"))
), marker)
cat("SMOKE TEST PASS — marker written:", marker, "\n")
