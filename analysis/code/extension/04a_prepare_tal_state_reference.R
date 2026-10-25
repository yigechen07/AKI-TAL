#!/usr/bin/env Rscript
# =============================================================================
# 04a_prepare_tal_state_reference.R
# Reproducible rebuild of the frozen TAL state reference for Figure 4.
#
# Input  : metadata/extension/reference_sources/raw_extracted_marker_sets.json
#          (verbatim published gene sets extracted from Lake et al. 2023,
#          DOI 10.1038/s41586-023-05769-3, supplementary Tables 5, 10, 12, 13;
#          produced once by the Markers track via python openpyxl and verified
#          by the Audit track — see reference_sources/provenance_audit_report.md)
# Output : metadata/extension/tal_state_reference.csv (long format, frozen)
#
# This script is the deterministic RULE ENGINE: it applies the six-gene
# anti-circularity exclusion and the flag vocabulary to the verbatim
# published sets, and assembles the long-format reference. It does not
# hand-pick or trim genes.
#
# Usage  : Rscript 04a_prepare_tal_state_reference.R
#          (base R + jsonlite; runs on any R >= 4.0)
# =============================================================================

suppressMessages({
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required: install.packages('jsonlite')")
  }
})

# --- resolve project root ---------------------------------------------------
# Robust to spaces in the repo path: resolve from the working directory.
# Documented usage: run from code/extension/ (cd code/extension && Rscript 04a_...).
# Optional first CLI arg overrides the project root.
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1 && nzchar(cli_args[1])) {
  proj_root <- normalizePath(cli_args[1])
} else {
  proj_root <- normalizePath(file.path(getwd(), "..", ".."))
}

json_path <- file.path(proj_root, "metadata", "extension", "reference_sources",
                       "raw_extracted_marker_sets.json")
out_path  <- file.path(proj_root, "metadata", "extension",
                       "tal_state_reference.csv")

stopifnot(file.exists(json_path))
raw <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)

# --- frozen six-gene primary-endpoint exclusion -----------------------------
FROZEN_SIX <- c("SLC12A1", "UMOD", "CLDN16", "KCNJ10", "PTH1R", "CLDN10")

# --- flag vocabularies (consistent, deterministic) --------------------------
HOUSEKEEPING <- c("ACTB","ACTG1","B2M","EEF1A1","FTH1","FTL","GAPDH","PTMA",
                  "TMSB10","TMSB4X","TPT1","UBA52","UBB","YWHAZ")
LNCRNA <- c("AC005699.1","AC010894.3","AC012593.1","AC016831.7","AL445472.1",
            "CASC15","FP236383.1","FP671120.1","LINC01320","MALAT1",
            "MTRNR2L8","NEAT1","SPINT1-AS1")
RIBOSOMAL_RE <- "^RPL|^RPS|^MRPL|^MRPS|^RPLP[0-9]|^RPSA"

flag_gene <- function(g) {
  if (grepl("^MT-", g))                          return("mitochondrial")
  if (g %in% HOUSEKEEPING)                       return("generic_housekeeping")
  if (g %in% LNCRNA)                             return("technical_lncRNA")
  if (grepl(RIBOSOMAL_RE, g))                    return("ribosomal")
  return("none")
}

# --- selection-rule strings (verbatim from the extraction log) --------------
RULE_T5 <- list(
  "Positive Markers (M-TAL)"     = "Table 5 row Subclass Level 3 = 'M-TAL' (Subclass Level 1 = 'TAL'); Positive Markers column",
  "Positive Markers (C-TAL)"     = "Table 5 row Subclass Level 3 = 'C-TAL' (Subclass Level 1 = 'TAL'); Positive Markers column",
  "Positive Markers (generic TAL)" = "Table 5 generic row Subclass Level 1 = 'TAL'; Positive Markers column"
)
RULE_T13 <- list(
  "Expanded Adaptive State Gene Set" = "Table 13 Expanded Adaptive State Gene Set (log2fc > 0.6 & pval_adj < 0.05 in Ref or AKI or CKD); State/Cluster = aTAL",
  "Full Conserved Adaptive State Gene Set" = "Table 13 Full Conserved Adaptive State Gene Set (p val < 0.05 for Ref, AKI and CKD); State/Cluster = aTAL",
  "Minimal Overlapping (aPT & aTAL) Set" = "Table 13 Minimal Overlapping (aPT & aTAL) Set (aTUB enriched pval < 0.05); shared adaptive epithelial, NOT TAL-specific"
)
RULE_T12 <- list(
  "Minimal Conserved Gene Set" = "Table 12 Minimal Conserved Gene Set (log2fc > 0.25, pval < 0.05 for Ref & AKI & CKD; >= 4 cell types)",
  "Expanded Gene Set"          = "Table 12 Expanded Gene Set (log2fc > 0.6 & pval_adj < 0.05 in Ref or AKI or CKD)"
)

# --- set-name normalisation (modality suffix -> source_set + modality) -------
normalise_set <- function(set_name, table) {
  if (table == "Table 5") {
    return(list(source_set = set_name, modality = "snCv3"))
  }
  mod <- ifelse(grepl("scCv3", set_name), "scCv3", "snCv3")
  ss  <- sub(" \\((snCv3|scCv3)\\)$", "", set_name)
  if (table == "Table 13") {
    ss <- sub(" aTAL$", "", ss)                  # drop state token
    if (grepl("Minimal Overlapping", ss)) {
      ss <- "Minimal Overlapping (aPT & aTAL) Set"
    }
  }
  list(source_set = ss, modality = mod)
}

rows <- list()

add_rows <- function(state, program, table) {
  sets <- program$sets
  dirn <- program$direction
  for (set_name in names(sets)) {
    # Negative Markers and Macula Densa positive markers are recorded for
    # reference/exclusion only and are NOT scored (documented in the extraction log)
    if (grepl("^Negative Markers", set_name)) next
    if (grepl("reference only|not scored", set_name)) next
    norm <- normalise_set(set_name, table)
    ss   <- norm$source_set
    mod  <- norm$modality
    rule <- switch(table,
                   "Table 5"  = RULE_T5[[ss]],
                   "Table 13" = RULE_T13[[ss]],
                   "Table 12" = RULE_T12[[ss]])
    if (is.null(rule)) stop("No selection rule for ", table, " / ", ss)
    sensitivity <- (ss == "Minimal Overlapping (aPT & aTAL) Set")
    for (g in sets[[set_name]]) {
      fl <- flag_gene(g)
      if (sensitivity && fl == "none") fl <- "shared_adaptive_epithelial"
      excl <- if (g %in% FROZEN_SIX) "frozen_primary_endpoint_circularity" else ""
      rows[[length(rows) + 1]] <<- data.frame(
        state              = state,
        gene               = g,
        direction          = dirn,
        source_table       = table,
        source_set         = ss,
        modality           = mod,
        selection_rule     = rule,
        flag               = fl,
        primary_sensitivity = ifelse(sensitivity, "sensitivity", "primary"),
        exclusion_reason   = excl,
        stringsAsFactors   = FALSE
      )
    }
  }
}

add_rows("reference_TAL",    raw$program1_reference_TAL,    "Table 5")
add_rows("adaptive_TAL",     raw$program2_adaptive_TAL,     "Table 13")
add_rows("degenerative_TAL", raw$program3_degenerative_TAL, "Table 12")

ref <- do.call(rbind, rows)
rownames(ref) <- NULL

write.csv(ref, out_path, row.names = FALSE)

# --- verification summary ----------------------------------------------------
cat("Wrote:", out_path, "\n")
cat("Rows:", nrow(ref), "\n")
cat("Rows by state:\n")
print(table(ref$state))
cat("Excluded rows (frozen six):", sum(ref$exclusion_reason != ""), "\n")
cat("Flagged rows (flag != none):", sum(ref$flag != "none"), "\n")
cat("Flag counts:\n")
print(table(ref$flag))
cat("Six-gene check:\n")
six <- ref[ref$gene %in% FROZEN_SIX, c("state", "gene", "exclusion_reason")]
print(six)
cat("Sensitivity rows:", sum(ref$primary_sensitivity == "sensitivity"), "\n")
