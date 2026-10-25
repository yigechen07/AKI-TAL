#!/usr/bin/env Rscript
# 06_evidence_freeze_check.R — deterministic cross-file numerical freeze check
# Statistical units and claim limits are documented in manuscript/evidence_freeze.md.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_arg) == 1L) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
} else {
  normalizePath("projects/AKI_TAL_severity_analysis")
}
result_dir <- file.path(project_dir, "results")
figure_dir <- file.path(project_dir, "figures")

formal <- read.csv(file.path(result_dir, "formal_effects.csv"), stringsAsFactors = FALSE)
bulk <- read.csv(file.path(result_dir, "bulk_effects.csv"), stringsAsFactors = FALSE)
longitudinal <- read.csv(file.path(result_dir, "gse126805_results.csv"), stringsAsFactors = FALSE)
coordination <- read.csv(file.path(result_dir, "tal_pt_coupling_results.csv"), stringsAsFactors = FALSE)
mito <- read.csv(file.path(result_dir, "tal_mito_results.csv"), stringsAsFactors = FALSE)

get_formal <- function(dataset, comparison, model) {
  x <- formal[formal$dataset == dataset & formal$comparison == comparison & formal$model == model, ]
  stopifnot(nrow(x) == 1L)
  x
}
get_bulk <- function(dataset, model) {
  x <- bulk[bulk$dataset == dataset & bulk$model == model, ]
  stopifnot(nrow(x) == 1L)
  x
}
get_longitudinal <- function(label) {
  x <- longitudinal[longitudinal$label == label, ]
  stopifnot(nrow(x) == 1L)
  x
}
get_coordination <- function(dataset, analysis) {
  x <- coordination[coordination$dataset == dataset & coordination$analysis == analysis, ]
  stopifnot(nrow(x) == 1L)
  x
}
get_mito <- function(dataset, analysis) {
  x <- mito[mito$dataset == dataset & mito$tier == "primary_ge50" & mito$analysis == analysis, ]
  stopifnot(nrow(x) == 1L)
  x
}

checks <- list(
  c("GSE183276 primary donor-pseudobulk", -0.75053196158748,
    get_formal("GSE183276", "AKI_vs_Reference", "unadjusted")$estimate,
    "donor"),
  c("GSE210622 primary donor-pseudobulk", -0.983086775217817,
    get_formal("GSE210622", "AKI_vs_nonAKI", "unadjusted")$estimate,
    "donor"),
  c("GSE139061 tissue-level bulk", -0.858140165162827,
    get_bulk("GSE139061", "unadjusted")$estimate,
    "sample/donor"),
  c("GSE30718 tissue-level microarray", -0.649053618699422,
    get_bulk("GSE30718", "unadjusted")$estimate,
    "unique donor (primary biopsy)"),
  c("GSE126805 paired post-minus-pre", -0.222835359780525,
    get_longitudinal("baseline_post_vs_pre_delta")$estimate,
    "paired kidney allograft"),
  c("Direction B pooled TAL-PT identity", 0.39870162270799,
    get_coordination("pooled", "TAL_identity_vs_PT_identity")$partial_r,
    "effect-level meta-analysis; donors remain cohort-specific"),
  c("Direction C GSE183276 TAL-OXPHOS partial r", -0.279155902762716,
    get_mito("GSE183276", "TAL_identity_vs_OXPHOS")$partial_r,
    "donor"),
  c("Direction C GSE210622 TAL-OXPHOS partial r", 0.584842042818858,
    get_mito("GSE210622", "TAL_identity_vs_OXPHOS")$partial_r,
    "donor")
)

out <- data.frame(
  item = vapply(checks, `[[`, character(1), 1),
  expected = as.numeric(vapply(checks, `[[`, character(1), 2)),
  observed = as.numeric(vapply(checks, `[[`, character(1), 3)),
  biological_unit = vapply(checks, `[[`, character(1), 4),
  stringsAsFactors = FALSE
)
out$absolute_difference <- abs(out$expected - out$observed)
out$status <- ifelse(out$absolute_difference < 1e-12, "PASS", "FAIL")

required_figures <- c("formal_1_cohort.pdf", "formal_2_scores.pdf", "formal_3_specificity.pdf",
                      "bulk_scores.png", "bulk_forest.png", "bulk_genes.png",
                      "longitudinal_gse126805.pdf", "longitudinal_gse126805.png")
figure_check <- data.frame(
  item = paste0("Required figure: ", required_figures), expected = NA_real_, observed = NA_real_,
  biological_unit = "not applicable", absolute_difference = NA_real_,
  status = ifelse(file.exists(file.path(figure_dir, required_figures)), "PASS", "FAIL"),
  stringsAsFactors = FALSE
)
out <- rbind(out, figure_check)

write.csv(out, file.path(result_dir, "evidence_freeze_check.csv"), row.names = FALSE)
if (any(out$status != "PASS")) stop("Evidence freeze check failed; see results/evidence_freeze_check.csv")
message("Evidence freeze check PASS: all numerical anchors and required figures are present.")
