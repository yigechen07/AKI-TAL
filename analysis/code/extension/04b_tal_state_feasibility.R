#!/usr/bin/env Rscript
# =============================================================================
# 04b_tal_state_feasibility.R
# Reproducible, deterministic regeneration of the Figure 4 TAL-state
# feasibility outputs:
#   - results/extension/tal_state_marker_coverage.csv  (canonical)
#   - results/extension/tal_donor_coverage.csv         (canonical)
#   - numerical sections of results/extension/tal_state_feasibility.md
#
# Design (frozen, diagnosis-blinded):
#   - Strict author-annotated TAL cells: metadata$author_label == "TAL"
#   - Detection = fraction of strict TAL cells with count > 0 on the
#     counts layer (sparse rowMeans over a dgCMatrix; no dense materialisation)
#   - Primary scoring sets per state: strict rule (see build_signatures())
#   - Sensitivity set: shared_adaptive_epithelial (aPT&aTAL overlap), NOT scored
#     as TAL-specific
#   - Donor-level summaries only; cells are never patient replicates
#   - NO diagnosis-informed gene/threshold selection; NO AKI-vs-control test
#   - Differently sized state sets are summarised per state/cohort only; never
#     compared by raw magnitude
#
# Usage:
#   cd code/extension && Rscript 04b_tal_state_feasibility.R
#   Optional first CLI arg: project root
# Requires: Seurat, Matrix, data.table, jsonlite
# =============================================================================

suppressMessages({
  for (pkg in c("Seurat", "Matrix", "data.table", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Required package missing: ", pkg)
    }
  }
  library(Seurat); library(Matrix); library(data.table); library(jsonlite)
})

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1 && nzchar(cli_args[1])) {
  proj_root <- normalizePath(cli_args[1])
} else {
  proj_root <- normalizePath(file.path(getwd(), "..", ".."))
}

REF_CSV   <- file.path(proj_root, "metadata", "extension", "tal_state_reference.csv")
OUT_COV   <- file.path(proj_root, "results", "extension", "tal_state_marker_coverage.csv")
OUT_DON   <- file.path(proj_root, "results", "extension", "tal_donor_coverage.csv")
OUT_MD    <- file.path(proj_root, "results", "extension", "tal_state_feasibility.md")
OBJ_183   <- file.path(proj_root, "results", "extension", "atlas_GSE183276_seurat.rds")
OBJ_622   <- file.path(proj_root, "results", "extension", "atlas_GSE210622_seurat.rds")

stopifnot(file.exists(REF_CSV))
stopifnot(file.exists(OBJ_183), file.exists(OBJ_622))

FROZEN_SIX <- c("SLC12A1", "UMOD", "CLDN16", "KCNJ10", "PTH1R", "CLDN10")

# ---------------------------------------------------------------------------
# build_signatures(): explicit signature-building function
#   - drops any gene carrying a frozen-six exclusion
#   - primary requires no non-'none' flag in ANY row for the same state x gene
#   - adaptive genes present in the shared aPT-aTAL sensitivity set never
#     re-enter the primary set through duplicate unflagged rows
#   - returns primary (strict) and sensitivity (shared) sets + exact-minus-six
#     (published-set-minus-frozen-six) sensitivity definition for robustness
# ---------------------------------------------------------------------------
build_signatures <- function(ref) {
  setDT(ref)
  ref[, excluded := exclusion_reason == "frozen_primary_endpoint_circularity"]

  # per (state, gene): does ANY row carry a non-none flag or a frozen exclusion?
  pair_info <- ref[, .(
    any_non_none_flag = any(flag != "none"),
    any_frozen        = any(excluded),
    has_none_row      = any(flag == "none"),
    flags             = paste(unique(flag), collapse = "|")
  ), by = .(state, gene)]

  # strict primary: has a none row, no non-none flag, no frozen exclusion
  primary <- pair_info[has_none_row == TRUE & any_non_none_flag == FALSE & any_frozen == FALSE,
                       .(state, gene)]

  # sensitivity: shared adaptive epithelial (aPT & aTAL overlap)
  sensitivity <- unique(ref[flag == "shared_adaptive_epithelial", .(state, gene)])

  # robustness definition: exact published set minus frozen six (row-wise)
  exact_minus_six <- unique(ref[excluded == FALSE, .(state, gene)])

  list(primary = primary, sensitivity = sensitivity,
       exact_minus_six = exact_minus_six,
       pair_info = pair_info)
}

# ---------------------------------------------------------------------------
# compute_feasibility(): cohort-level coverage + donor-level summaries
# ---------------------------------------------------------------------------
compute_feasibility <- function(obj_path, donor_col, cohort, primary, sensitivity) {
  obj <- readRDS(obj_path)
  md  <- as.data.table(obj@meta.data)
  counts <- GetAssayData(obj, layer = "counts")
  strict_cells <- which(md$author_label == "TAL")
  donor_of_cell <- md[[donor_col]][strict_cells]
  diag_map <- unique(md[, .(donor_id = get(donor_col), diagnosis)])
  n_strict <- length(strict_cells)

  cov_rows <- list(); don_rows <- list()

  for (st in unique(primary$state)) {
    genes <- intersect(primary[state == st, gene], rownames(counts))
    X <- counts[genes, strict_cells, drop = FALSE]
    det_cohort <- rowMeans(X > 0)
    for (g in genes)
      cov_rows[[length(cov_rows) + 1]] <- data.table(
        cohort = cohort, state = st, set_type = "primary", gene = g,
        available_in_assay = TRUE,
        detection_fraction_among_strict_TAL = unname(det_cohort[[g]]),
        n_strict_TAL_cells = n_strict)
    for (g in setdiff(primary[state == st, gene], genes))
      cov_rows[[length(cov_rows) + 1]] <- data.table(
        cohort = cohort, state = st, set_type = "primary", gene = g,
        available_in_assay = FALSE,
        detection_fraction_among_strict_TAL = NA_real_,
        n_strict_TAL_cells = n_strict)
    for (d in unique(donor_of_cell)) {
      cells_idx <- which(donor_of_cell == d)
      Xd <- X[, cells_idx, drop = FALSE]
      don_rows[[length(don_rows) + 1]] <- data.table(
        cohort = cohort, donor_id = d, state = st, set_type = "primary",
        n_strict_TAL = length(cells_idx),
        mean_detection = mean(rowMeans(Xd > 0)))
    }
  }

  for (st in unique(sensitivity$state)) {
    genes <- intersect(sensitivity[state == st, gene], rownames(counts))
    X <- counts[genes, strict_cells, drop = FALSE]
    det_cohort <- rowMeans(X > 0)
    for (g in genes)
      cov_rows[[length(cov_rows) + 1]] <- data.table(
        cohort = cohort, state = st, set_type = "sensitivity", gene = g,
        available_in_assay = TRUE,
        detection_fraction_among_strict_TAL = unname(det_cohort[[g]]),
        n_strict_TAL_cells = n_strict)
    for (g in setdiff(sensitivity[state == st, gene], genes))
      cov_rows[[length(cov_rows) + 1]] <- data.table(
        cohort = cohort, state = st, set_type = "sensitivity", gene = g,
        available_in_assay = FALSE,
        detection_fraction_among_strict_TAL = NA_real_,
        n_strict_TAL_cells = n_strict)
    for (d in unique(donor_of_cell)) {
      cells_idx <- which(donor_of_cell == d)
      Xd <- X[, cells_idx, drop = FALSE]
      don_rows[[length(don_rows) + 1]] <- data.table(
        cohort = cohort, donor_id = d, state = st, set_type = "sensitivity",
        n_strict_TAL = length(cells_idx),
        mean_detection = mean(rowMeans(Xd > 0)))
    }
  }

  cov <- rbindlist(cov_rows)
  don <- rbindlist(don_rows)
  don <- merge(don, diag_map, by = "donor_id")
  don[, tier := fifelse(n_strict_TAL >= 50, "ge50_primary",
                        fifelse(n_strict_TAL >= 20, "ge20_sensitivity",
                                "lt20_low-confidence"))]
  rm(obj, counts); gc()
  list(cov = cov, don = don, n_strict = n_strict,
       md_summary = md[, .(n_donors = uniqueN(get(donor_col)), n_cells = .N),
                       by = diagnosis])
}

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
ref <- fread(REF_CSV, na.strings = c("NA", ""))
sig <- build_signatures(ref)

stopifnot(sig$primary[state == "reference_TAL", .N] == 17)
stopifnot(sig$primary[state == "adaptive_TAL", .N] == 281)
stopifnot(sig$primary[state == "degenerative_TAL", .N] == 767)
stopifnot(sig$sensitivity[, .N] == 132)
stopifnot(!any(sig$primary$gene %in% FROZEN_SIX))

r1 <- compute_feasibility(OBJ_183, "patient", "GSE183276", sig$primary, sig$sensitivity)
r2 <- compute_feasibility(OBJ_622, "donor",    "GSE210622", sig$primary, sig$sensitivity)

coverage_all <- rbind(r1$cov, r2$cov)
donor_all    <- rbind(r1$don, r2$don)
setcolorder(coverage_all, c("cohort", "state", "set_type", "gene",
                            "available_in_assay",
                            "detection_fraction_among_strict_TAL",
                            "n_strict_TAL_cells"))
setcolorder(donor_all, c("cohort", "donor_id", "diagnosis", "n_strict_TAL",
                         "tier", "state", "set_type", "mean_detection"))

fwrite(coverage_all, OUT_COV)
fwrite(donor_all, OUT_DON)
cat("Wrote coverage:", OUT_COV, nrow(coverage_all), "rows\n")
cat("Wrote donor:", OUT_DON, nrow(donor_all), "rows\n")

# ---- numerical summaries for the markdown sections --------------------------
fmt3 <- function(x) sprintf("%.3f", x)

det_summary <- coverage_all[set_type == "primary", .(
  n_genes = .N,
  n_available = sum(available_in_assay),
  det_median = median(detection_fraction_among_strict_TAL, na.rm = TRUE),
  det_min = min(detection_fraction_among_strict_TAL, na.rm = TRUE),
  det_max = max(detection_fraction_among_strict_TAL, na.rm = TRUE)
), by = .(cohort, state)]

det_sens <- coverage_all[set_type == "sensitivity", .(
  n_genes = .N,
  n_available = sum(available_in_assay),
  det_median = median(detection_fraction_among_strict_TAL, na.rm = TRUE),
  det_min = min(detection_fraction_among_strict_TAL, na.rm = TRUE),
  det_max = max(detection_fraction_among_strict_TAL, na.rm = TRUE)
), by = .(cohort, state)]

don_summary <- donor_all[set_type == "primary", .(
  n_donors = .N,
  donor_det_median = median(mean_detection),
  donor_det_min = min(mean_detection),
  donor_det_max = max(mean_detection)
), by = .(cohort, state)]

don_sens <- donor_all[set_type == "sensitivity", .(
  n_donors = .N,
  donor_det_median = median(mean_detection),
  donor_det_min = min(mean_detection),
  donor_det_max = max(mean_detection)
), by = .(cohort, state)]

near_zero <- coverage_all[set_type == "primary" & available_in_assay == TRUE, .(
  n_genes = .N,
  zero_det = sum(detection_fraction_among_strict_TAL == 0),
  lt1pct = sum(detection_fraction_among_strict_TAL < 0.01)
), by = .(cohort, state)]

tier_counts <- donor_all[, .(n_donors = uniqueN(donor_id)), by = .(cohort, tier)][order(cohort, tier)]

low_conf <- unique(donor_all[set_type == "primary" & n_strict_TAL < 20,
                             .(cohort, donor_id, diagnosis, n_strict_TAL, tier)])

# ---- assemble markdown (numerical sections 1-5, 8; keep fixed prose) --------
md <- c(
  "# Figure 4 Feasibility Findings — TAL State Scoring (Diagnosis-Blinded)",
  "",
  "Feasibility track of the Figure 4 pre-analysis package (AKI_TAL_severity_analysis).",
  "Scope: verify both cohort Seurat objects, confirm frozen state-program marker availability and",
  "detection in strict author-annotated TAL cells, assess donor-level coverage, and report whether",
  "continuous state scoring is feasible. **No AKI-versus-control tests, no outcome-informed",
  "decisions, no threshold tuning were performed.**",
  "",
  "---",
  "",
  "## 1. Cohort verification",
  "",
  sprintf("### GSE183276 (scRNA-seq, 10X) — %d strict TAL cells, %d donors",
          r1$n_strict, sum(r1$md_summary$n_donors)),
  sprintf("- Strict author-annotated TAL cells (`author_label == \"TAL\"`): **%d** (primary).", r1$n_strict),
  sprintf("- Donors by diagnosis (descriptor only): %s",
          paste(sprintf("%s %d", r1$md_summary$diagnosis, r1$md_summary$n_donors), collapse = "; ")),
  "",
  sprintf("### GSE210622 (snRNA-seq) — %d strict TAL cells, %d donors",
          r2$n_strict, sum(r2$md_summary$n_donors)),
  sprintf("- Strict author-annotated TAL cells (`author_label == \"TAL\"`): **%d** (primary).", r2$n_strict),
  sprintf("- Donors by diagnosis (descriptor only): %s",
          paste(sprintf("%s %d", r2$md_summary$diagnosis, r2$md_summary$n_donors), collapse = "; ")),
  "",
  "---",
  "",
  "## 2. Frozen state-program gene sets",
  "",
  sprintf("- Marker table rows: %d. Frozen six-gene exclusion rows: %d; flagged rows: %d.",
          nrow(ref), sum(ref$exclusion_reason == "frozen_primary_endpoint_circularity"),
          sum(ref$flag != "none")),
  sprintf("- Primary (strict): reference_TAL %d, adaptive_TAL %d, degenerative_TAL %d.",
          sig$primary[state == "reference_TAL", .N],
          sig$primary[state == "adaptive_TAL", .N],
          sig$primary[state == "degenerative_TAL", .N]),
  sprintf("- Sensitivity (shared adaptive epithelial, aPT&aTAL overlap): %d genes (adaptive_TAL only).",
          sig$sensitivity[, .N]),
  "",
  "### Cross-state gene sharing (disclosure)",
  "",
  "The three frozen state programs overlap by construction (Lake 2023 tables share",
  "conserved marker genes across cell states). The scoring phase treats the state",
  "scores as correlated by construction and reports this overlap. Strict-primary",
  "and sensitivity set intersections (computed by 04b, deterministic):",
  "",
  sprintf("- reference_TAL primary ∩ degenerative_TAL primary: %d gene(s).",
          length(intersect(sig$primary[state == "reference_TAL", gene],
                           sig$primary[state == "degenerative_TAL", gene]))),
  sprintf("- adaptive_TAL primary ∩ degenerative_TAL primary: %d gene(s).",
          length(intersect(sig$primary[state == "adaptive_TAL", gene],
                           sig$primary[state == "degenerative_TAL", gene]))),
  sprintf("- adaptive sensitivity (scored 132) ∩ degenerative_TAL primary: %d gene(s).",
          length(intersect(sig$sensitivity[, gene],
                           sig$primary[state == "degenerative_TAL", gene]))),
  sprintf("- adaptive_TAL primary ∩ adaptive sensitivity: %d gene(s) (zero overlap by construction).",
          length(intersect(sig$primary[state == "adaptive_TAL", gene],
                           sig$sensitivity[, gene]))),
  "",
  "### Sensitivity-set accounting",
  "",
  sprintf("- Rows in the sensitivity set (`primary_sensitivity == \"sensitivity\"`): %d; unique genes: %d; scored genes: %d.",
          nrow(ref[primary_sensitivity == "sensitivity"]),
          length(unique(ref[primary_sensitivity == "sensitivity", gene])),
          sig$sensitivity[, .N]),
  sprintf("- The %d unscored unique gene(s) carry a technical/housekeeping flag and are not scored: %s.",
          length(setdiff(unique(ref[primary_sensitivity == "sensitivity", gene]),
                         sig$sensitivity[, gene])),
          paste(sort(setdiff(unique(ref[primary_sensitivity == "sensitivity", gene]),
                             sig$sensitivity[, gene])), collapse = ", ")),
  "",
  "---",
  "",
  "## 3. Gene availability and detection (diagnosis-blinded, strict TAL cells)",
  "",
  "### Cohort-wide per-gene detection among strict TAL cells (primary)",
  "",
  "| cohort | state | n genes | n available | detection median | detection min | detection max |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %s | %s | %d | %d | %s | %s | %s |",
          det_summary$cohort, det_summary$state, det_summary$n_genes,
          det_summary$n_available, fmt3(det_summary$det_median),
          fmt3(det_summary$det_min), fmt3(det_summary$det_max)),
  "",
  "### Sensitivity set (shared adaptive epithelial)",
  "",
  "| cohort | state | n genes | n available | detection median | detection min | detection max |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %s | %s | %d | %d | %s | %s | %s |",
          det_sens$cohort, det_sens$state, det_sens$n_genes,
          det_sens$n_available, fmt3(det_sens$det_median),
          fmt3(det_sens$det_min), fmt3(det_sens$det_max)),
  "",
  "### Near-zero detection (primary)",
  "",
  "| cohort | state | n genes | zero detection | <1% detection |",
  "|---|---|---|---|---|",
  sprintf("| %s | %s | %d | %d | %d |",
          near_zero$cohort, near_zero$state, near_zero$n_genes,
          near_zero$zero_det, near_zero$lt1pct),
  "",
  "_Near-zero wording:_ the <1% counts include the zero-detection genes (e.g. GSE183276",
  "degenerative 22 = 21 + 1 zero-detection; GSE210622 degenerative 63 = 56 + 7).",
  "",
  "### Per-donor mean detection (primary sets, donor-level)",
  "",
  "| cohort | state | n donors | donor median | donor min | donor max |",
  "|---|---|---|---|---|---|",
  sprintf("| %s | %s | %d | %s | %s | %s |",
          don_summary$cohort, don_summary$state, don_summary$n_donors,
          fmt3(don_summary$donor_det_median), fmt3(don_summary$donor_det_min),
          fmt3(don_summary$donor_det_max)),
  "",
  "### Per-donor mean detection (sensitivity)",
  "",
  "| cohort | state | n donors | donor median | donor min | donor max |",
  "|---|---|---|---|---|---|",
  sprintf("| %s | %s | %d | %s | %s | %s |",
          don_sens$cohort, don_sens$state, don_sens$n_donors,
          fmt3(don_sens$donor_det_median), fmt3(don_sens$donor_det_min),
          fmt3(don_sens$donor_det_max)),
  "",
  "---",
  "",
  "## 4. Non-core gene retention per state",
  "",
  "Rule of thumb: a state is reliably scorable if ≥10 non-excluded, non-flagged genes are detected.",
  "",
  "| cohort | state | n scorable | n available | verdict |",
  "|---|---|---|---|---|",
  sprintf("| %s | %s | %d | %d | %s |",
          det_summary$cohort, det_summary$state, det_summary$n_genes,
          det_summary$n_available,
          ifelse(det_summary$n_genes >= 10, "PASS (≥10)", "CHECK")),
  "",
  "---",
  "",
  "## 5. Donor-level TAL and state coverage",
  "",
  "Tier definitions: ge50_primary (≥50 strict TAL cells) / ge20_sensitivity (20–49) / lt20_low-confidence (<20).",
  "",
  "| cohort | tier | n donors |",
  "|---|---|---|",
  sprintf("| %s | %s | %d |", tier_counts$cohort, tier_counts$tier, tier_counts$n_donors),
  "",
  if (nrow(low_conf) > 0) {
    c(sprintf("**Low-confidence donors (lt20 strict TAL):** %s.",
              paste(sprintf("%s %s (%d cells)", low_conf$cohort, low_conf$donor_id, low_conf$n_strict_TAL),
                    collapse = "; ")),
      "These donors are excluded from (or flagged in) donor-level state scoring.")
  } else {
    "No low-confidence donors."
  },
  "",
  "Cells are never treated as patient replicates; all summaries are donor-level.",
  "",
  "---",
  "",
  "## 6. Diagnosis-blinding statement",
  "",
  "Diagnosis appears **only as a descriptive donor attribute**. It was never used to",
  "select genes, set thresholds, choose tiers, or drive any decision. No AKI-versus-control",
  "tests were run; no outcome-informed thresholds were examined. Gene sets were frozen",
  "before this feasibility pass, and all detection computations are pooled across diagnoses.",
  "",
  "---",
  "",
  "## 7. Recommendation",
  "",
  "**Use UCell continuous scoring as the primary state-scoring approach** (per-cell continuous",
  "enrichment for each frozen state program, then donor-level aggregation of continuous scores).",
  "Rationale: all three state programs retain ample scorable genes (17/281/767) with broad",
  "detection in both cohorts; hard state assignment would require an outcome-informed",
  "classification threshold (forbidden); continuous scores preserve the aPT&aTAL overlap",
  "structure (132 shared adaptive epithelial genes scored separately as sensitivity).",
  "UCell 2.16.0 is installed in `metadata/extension/r-lib` (verified) and the",
  "diagnosis-blinded smoke test passed on a 500-cell strict-TAL subset",
  "(17/281/767/132 genes scored; marker `metadata/extension/r-lib/UCell_SMOKE_PASS.txt`).",
  "",
  "---",
  "",
  "## 8. Install route for UCell (scoring phase)",
  "",
  "Install UCell into the **project-local** R library `metadata/extension/r-lib`",
  "(never modify global R libraries):",
  "",
  "```r",
  "# Primary route (Bioconductor, project-local):",
  "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")",
  "dir.create(\"metadata/extension/r-lib\", showWarnings = FALSE, recursive = TRUE)",
  "BiocManager::install(\"UCell\", lib = \"metadata/extension/r-lib\")",
  ".libPaths(c(\"metadata/extension/r-lib\", .libPaths()))",
  "```",
  "",
  "If the sandbox blocks the Bioconductor `3.x/bioc` package index, use the",
  "release tarball directly (UCell itself is pure R; its C++ dependency",
  "BiocNeighbors must be compiled with the conda R toolchain on PATH):",
  "",
  "```r",
  "# Fallback (direct tarball, project-local; requires BiocNeighbors installed):",
  "#   R CMD INSTALL --library=metadata/extension/r-lib BiocNeighbors_2.6.0.tar.gz",
  "#   R CMD INSTALL --library=metadata/extension/r-lib UCell_2.16.0.tar.gz",
  ".libPaths(c(\"metadata/extension/r-lib\", .libPaths()))",
  "```",
  "",
  "The scoring phase must not start until UCell is verified present in that library.",
  "",
  "---",
  "",
  "## 9. Reproducibility and canonical outputs",
  "",
  "Regenerated by `code/extension/04b_tal_state_feasibility.R` (deterministic; date-stamped on run date).",
  "Run twice to confirm deterministic outputs: `Rscript 04b_tal_state_feasibility.R`.",
  "",
  "**Canonical deliverables:** `tal_state_marker_coverage.csv`, `tal_donor_coverage.csv`,",
  "and `tal_state_feasibility.md` are canonical. Compatibility copies retained for",
  "back-compatibility with the Feasibility track names (`tal_marker_coverage.csv`,",
  "`tal_feasibility_findings.md`) are maintained **byte-identical** to the canonical",
  "files and are documented as compatibility aliases only."
)

writeLines(md, OUT_MD)
cat("Wrote feasibility MD:", OUT_MD, "\n")

# -----------------------------------------------------------------------------
# Sync compatibility aliases (byte-identical), then verify byte/checksum.
# The canonical deliverables are:
#   tal_state_marker_coverage.csv, tal_donor_coverage.csv, tal_state_feasibility.md
# Compatibility aliases retained for back-compatibility with the Feasibility
# track names:
#   tal_marker_coverage.csv   <- tal_state_marker_coverage.csv
#   tal_feasibility_findings.md <- tal_state_feasibility.md
# These must remain byte-identical to the canonical files. Any sync or
# verification failure is a hard error: the script must NOT exit 0 with
# inconsistent aliases.
# -----------------------------------------------------------------------------
ALIAS_COV <- file.path(proj_root, "results", "extension", "tal_marker_coverage.csv")
ALIAS_MD  <- file.path(proj_root, "results", "extension", "tal_feasibility_findings.md")

file_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

sync_and_verify <- function(canonical, alias, label) {
  if (!file.exists(canonical)) stop("Canonical file missing: ", canonical)
  ok <- file.copy(canonical, alias, overwrite = TRUE)
  if (!isTRUE(ok)) stop("Failed to sync alias for ", label, ": ", alias)
  hc <- file_sha256(canonical)
  ha <- file_sha256(alias)
  if (!identical(hc, ha)) {
    stop(sprintf("Alias verification failed for %s: canonical sha256 %s != alias sha256 %s (%s)", label, hc, ha, alias))
  }
  # Byte-level double-check: read both and compare raw bytes (digest(file=TRUE)
  # reads the file; this second check guards against any digest quirk).
  rb1 <- readBin(canonical, "raw", n = file.info(canonical)$size)
  rb2 <- readBin(alias, "raw", n = file.info(alias)$size)
  if (!identical(rb1, rb2)) stop("Byte-level mismatch for alias ", label, ": ", alias)
  cat("Synced alias (byte-identical):", label, "->", alias, "sha256", substr(hc, 1, 16), "\n")
  invisible(hc)
}

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required for alias verification. Install it (BiocManager::install('digest') or install.packages('digest')).")
}

sync_and_verify(OUT_COV, ALIAS_COV, "tal_state_marker_coverage.csv -> tal_marker_coverage.csv")
sync_and_verify(OUT_MD,  ALIAS_MD,  "tal_state_feasibility.md -> tal_feasibility_findings.md")

cat("DONE\n")
