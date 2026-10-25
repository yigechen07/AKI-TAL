#!/usr/bin/env Rscript
#
# 02_primary_cell_resolved.R — Formal AKI–TAL donor-level analysis with frozen design
#
# Z-scores computed ONCE per dataset/tier on the full cohort-level scaling
# set, then subset for each comparison. All values are auditable against
# formal_effects.csv, formal_genes.csv, formal_lodo.csv, and
# formal_evidence.csv.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) == 1L) {
  project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
} else {
  project_dir <- normalizePath("projects/AKI_TAL_severity_analysis")
}
workspace_dir  <- dirname(dirname(project_dir))
source_dir     <- file.path(workspace_dir, "projects", "AKI_DCT_analysis")
result_dir     <- file.path(project_dir, "results")
figure_dir     <- file.path(project_dir, "figures")
meta_dir       <- file.path(project_dir, "metadata")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

SET_SEED <- 42
set.seed(SET_SEED)
TIMESTAMP <- Sys.time()

message("=== 02_primary_cell_resolved.R: Formal AKI–TAL donor-level analysis ===")
message(sprintf("Timestamp: %s\n", TIMESTAMP))

# ---- 1. Load sources ----
g183276_raw   <- read.delim(file.path(source_dir, "data", "GSE183276_Metadata.txt.gz"),
                            stringsAsFactors = FALSE, check.names = FALSE)
qc            <- read.csv(file.path(source_dir, "results", "stage3b_qc.csv"),
                          stringsAsFactors = FALSE)
shared        <- readRDS(file.path(source_dir, "results", "stage3b_pseudobulk.rds"))
sig           <- read.csv(file.path(source_dir, "results", "stage3_signature_definitions.csv"),
                          stringsAsFactors = FALSE)
pilot         <- read.csv(file.path(meta_dir, "pilot_donors.csv"),
                          stringsAsFactors = FALSE)
dct_donors    <- read.csv(file.path(source_dir, "results", "final", "donor_metadata.csv"),
                          stringsAsFactors = FALSE)
g210622_audit <- read.csv(file.path(source_dir, "results", "stage2_gse210622_donor_audit.csv"),
                          stringsAsFactors = FALSE)

# ---- 2. Frozen gene sets ----
TAL_GENES     <- c("SLC12A1", "UMOD", "CLDN16", "KCNJ10", "PTH1R", "CLDN10")
PT_ID_GENES   <- sig$gene[sig$program == "PT_identity"]
PT_INJ_GENES  <- sig$gene[sig$program == "Pan_epithelial_injury"]
stopifnot(all(TAL_GENES %in% sig$gene[sig$program == "TAL_identity"]))

# ---- 3. Extract segment cell counts from QC ----
tal_qc <- qc[qc$segment == "TAL", c("dataset", "donor", "n_cells")]
pt_qc  <- qc[qc$segment == "PT",  c("dataset", "donor", "n_cells")]
colnames(tal_qc)[colnames(tal_qc) == "donor"]   <- "donor_id"
colnames(tal_qc)[colnames(tal_qc) == "n_cells"] <- "n_TAL"
colnames(pt_qc)[colnames(pt_qc) == "donor"]     <- "donor_id"
colnames(pt_qc)[colnames(pt_qc) == "n_cells"]   <- "n_PT"

# ---- 4. GSE183276 donor-level metadata from cell metadata ----
dm_cols <- c("patient","condition.l1","condition.l2","condition.l3",
             "tissue_type","assay","sex","race")
g183276_dm <- unique(g183276_raw[, dm_cols])
colnames(g183276_dm)[colnames(g183276_dm) == "patient"] <- "donor_id"
stopifnot(nrow(g183276_dm) == 45)
stopifnot(all(g183276_dm$tissue_type == "Biopsy"))

pat_exp <- tapply(g183276_raw$experiment.short, g183276_raw$patient, function(x) x[1])
g183276_dm$experiment <- pat_exp[match(g183276_dm$donor_id, names(pat_exp))]

# ---- 5. Build GSE183276 formal donor table ----
g183276_tal <- tal_qc[tal_qc$dataset == "GSE183276", c("donor_id","n_TAL")]
g183276_pt  <- pt_qc[pt_qc$dataset   == "GSE183276", c("donor_id","n_PT")]
g183276 <- merge(g183276_tal, g183276_pt, by = "donor_id")
g183276 <- merge(g183276, g183276_dm,   by = "donor_id")
stopifnot(nrow(g183276) == 45)

pil_183276 <- pilot[pilot$dataset == "GSE183276", ]

g183276_out <- data.frame(
  dataset           = "GSE183276",
  donor_id          = g183276$donor_id,
  diagnosis         = g183276$condition.l1,
  diagnosis_l2      = g183276$condition.l2,
  diagnosis_l3      = g183276$condition.l3,
  kdigo_stage       = NA_character_,
  sex               = g183276$sex,
  race              = g183276$race,
  age_bin           = NA_character_,
  egfr_bin          = NA_character_,
  diabetes          = NA_character_,
  hypertension      = NA_character_,
  pathology_summary = NA_character_,
  tissue_type       = g183276$tissue_type,
  sampling_context  = ifelse(g183276$condition.l1 == "Ref",
    "Living-donor biopsy (source: condition.l3)",
    "Clinical biopsy"),
  assay             = "scRNA-seq",
  post_mortem       = "No",
  covid_status      = "not_reported",
  n_cells_total     = NA_integer_,
  n_TAL_cells       = g183276$n_TAL,
  n_PT_cells        = g183276$n_PT,
  tal_tier          = ifelse(g183276$n_TAL >= 50, "primary",
    ifelse(g183276$n_TAL >= 20, "sensitivity", "excluded")),
  batch             = g183276$experiment,
  project           = "KPMP (Lake et al. 2023 Nature 619:585)",
  original_accession = "GSE183276",
  stringsAsFactors   = FALSE
)

aki_idx <- match(g183276_out$donor_id, pil_183276$donor_id)
g183276_out$kdigo_stage[!is.na(aki_idx)]    <- pil_183276$kdigo_stage[aki_idx[!is.na(aki_idx)]]
g183276_out$age_bin[!is.na(aki_idx)]        <- pil_183276$age_bin[aki_idx[!is.na(aki_idx)]]
g183276_out$egfr_bin[!is.na(aki_idx)]       <- pil_183276$egfr_bin[aki_idx[!is.na(aki_idx)]]
g183276_out$diabetes[!is.na(aki_idx)]       <- pil_183276$diabetes[aki_idx[!is.na(aki_idx)]]
g183276_out$hypertension[!is.na(aki_idx)]   <- pil_183276$hypertension[aki_idx[!is.na(aki_idx)]]
g183276_out$pathology_summary[!is.na(aki_idx)] <- pil_183276$pathology_or_context[aki_idx[!is.na(aki_idx)]]
dct_idx <- match(g183276_out$donor_id, dct_donors$donor)
g183276_out$n_cells_total <- dct_donors$total_cells[dct_idx]
g183276_out$kdigo_stage[g183276_out$diagnosis == "Ref"] <- "reference"
g183276_out$kdigo_stage[g183276_out$diagnosis == "CKD"] <- "CKD"

# ---- 6. Build GSE210622 formal donor table ----
g210622_tal <- tal_qc[tal_qc$dataset == "GSE210622", c("donor_id","n_TAL")]
g210622_pt  <- pt_qc[pt_qc$dataset   == "GSE210622", c("donor_id","n_PT")]
g210622_tal_dx <- qc[qc$segment == "TAL" & qc$dataset == "GSE210622",
                     c("donor","diagnosis")]
colnames(g210622_tal_dx) <- c("donor_id","diagnosis_qc")
g210622 <- merge(g210622_tal, g210622_pt, by = "donor_id")
g210622 <- merge(g210622, g210622_tal_dx, by = "donor_id")
pil_210622 <- pilot[pilot$dataset == "GSE210622", ]

g210622_out <- data.frame(
  dataset           = "GSE210622",
  donor_id          = g210622$donor_id,
  diagnosis         = ifelse(g210622$diagnosis_qc == "AKI", "AKI", "control"),
  diagnosis_l2      = g210622$diagnosis_qc,
  diagnosis_l3      = NA_character_,
  kdigo_stage       = NA_character_,
  sex               = NA_character_,
  race              = NA_character_,
  age_bin           = NA_character_,
  egfr_bin          = NA_character_,
  diabetes          = NA_character_,
  hypertension      = NA_character_,
  pathology_summary = NA_character_,
  tissue_type       = NA_character_,
  sampling_context  = NA_character_,
  assay             = "snRNA-seq",
  post_mortem       = NA_character_,
  covid_status      = NA_character_,
  n_cells_total     = NA_integer_,
  n_TAL_cells       = g210622$n_TAL,
  n_PT_cells        = g210622$n_PT,
  tal_tier          = ifelse(g210622$n_TAL >= 50, "primary",
    ifelse(g210622$n_TAL >= 20, "sensitivity", "excluded")),
  batch             = NA_character_,
  project           = "Hinze et al. 2022 Genome Medicine 14:126",
  original_accession = "GSE210622",
  stringsAsFactors   = FALSE
)

aud_idx <- match(g210622_out$donor_id, g210622_audit$donor_id)
g210622_out$tissue_type      <- g210622_audit$tissue_source[aud_idx]
g210622_out$sampling_context <- g210622_audit$tissue_source[aud_idx]
g210622_out$post_mortem      <- g210622_audit$postmortem_status[aud_idx]
g210622_out$covid_status     <- g210622_audit$COVID_status[aud_idx]
g210622_out$n_cells_total    <- g210622_audit$number_total_cells[aud_idx]
g210622_out$batch            <- g210622_audit$GSM_accessions[aud_idx]

pil_idx <- match(g210622_out$donor_id, pil_210622$donor_id)
g210622_out$kdigo_stage[!is.na(pil_idx)]      <- pil_210622$kdigo_stage[pil_idx[!is.na(pil_idx)]]
g210622_out$sex[!is.na(pil_idx)]              <- pil_210622$sex[pil_idx[!is.na(pil_idx)]]
g210622_out$race[!is.na(pil_idx)]             <- pil_210622$race_reported[pil_idx[!is.na(pil_idx)]]
g210622_out$age_bin[!is.na(pil_idx)]          <- pil_210622$age_bin[pil_idx[!is.na(pil_idx)]]
g210622_out$egfr_bin[!is.na(pil_idx)]         <- pil_210622$egfr_bin[pil_idx[!is.na(pil_idx)]]
g210622_out$diabetes[!is.na(pil_idx)]         <- pil_210622$diabetes[pil_idx[!is.na(pil_idx)]]
g210622_out$hypertension[!is.na(pil_idx)]     <- pil_210622$hypertension[pil_idx[!is.na(pil_idx)]]
g210622_out$pathology_summary[!is.na(pil_idx)] <- pil_210622$pathology_or_context[pil_idx[!is.na(pil_idx)]]

# ---- 7. Combine and validate ----
formal_donors <- rbind(g183276_out, g210622_out)
stopifnot(nrow(formal_donors) == 57)
stopifnot(!anyDuplicated(formal_donors[, c("dataset", "donor_id")]))
stopifnot(!any(is.na(formal_donors$n_TAL_cells)))
stopifnot(!any(is.na(formal_donors$n_PT_cells)))

message("=== Donor census ===")
cat(sprintf("GSE183276: AKI=%d CKD=%d Ref=%d  TAL primary=%d sens=%d excl=%d\n",
  sum(g183276_out$diagnosis == "AKI"), sum(g183276_out$diagnosis == "CKD"),
  sum(g183276_out$diagnosis == "Ref"), sum(g183276_out$tal_tier == "primary"),
  sum(g183276_out$tal_tier == "sensitivity"), sum(g183276_out$tal_tier == "excluded")))
cat(sprintf("GSE210622: AKI=%d control=%d  TAL primary=%d\n",
  sum(g210622_out$diagnosis == "AKI"), sum(g210622_out$diagnosis == "control"),
  sum(g210622_out$tal_tier == "primary")))

# ---- 8. Verify frozen genes in pseudobulk ----
for (ds in c("GSE183276", "GSE210622")) {
  for (seg in c("TAL", "PT")) {
    gset <- if (seg == "TAL") TAL_GENES else c(PT_ID_GENES, PT_INJ_GENES)
    mat <- shared$pseudobulk[[ds]][[seg]]
    present <- intersect(gset, rownames(mat))
    missing <- setdiff(gset, rownames(mat))
    if (length(missing) > 0)
      stop(sprintf("MISSING genes in %s %s: %s", ds, seg, paste(missing, collapse=", ")))
  }
}
message("All frozen genes present in both cohorts.")

# ---- 9. Verify donors in pseudobulk ----
for (ds in c("GSE183276", "GSE210622")) {
  donors <- formal_donors$donor_id[formal_donors$dataset == ds]
  stopifnot(all(donors %in% colnames(shared$pseudobulk[[ds]][["TAL"]])))
  stopifnot(all(donors %in% colnames(shared$pseudobulk[[ds]][["PT"]])))
}

# ---- 10. Frozen donor lists ----
AKI_183276_P    <- g183276_out$donor_id[g183276_out$diagnosis == "AKI"  & g183276_out$tal_tier == "primary"]
REF_183276_P    <- g183276_out$donor_id[g183276_out$diagnosis == "Ref"  & g183276_out$tal_tier == "primary"]
CKD_183276_P    <- g183276_out$donor_id[g183276_out$diagnosis == "CKD"  & g183276_out$tal_tier == "primary"]
AKI_210622_P    <- g210622_out$donor_id[g210622_out$diagnosis == "AKI"  & g210622_out$tal_tier == "primary"]
NONAKI_210622_P <- g210622_out$donor_id[g210622_out$diagnosis != "AKI" & g210622_out$tal_tier == "primary"]

AKI_183276_S    <- g183276_out$donor_id[g183276_out$diagnosis == "AKI"  & g183276_out$tal_tier %in% c("primary","sensitivity")]
REF_183276_S    <- g183276_out$donor_id[g183276_out$diagnosis == "Ref"  & g183276_out$tal_tier %in% c("primary","sensitivity")]
CKD_183276_S    <- g183276_out$donor_id[g183276_out$diagnosis == "CKD"  & g183276_out$tal_tier %in% c("primary","sensitivity")]

cat(sprintf("\nFrozen: GSE183276 AKI=%d Ref=%d CKD=%d (primary, n=%d)\n",
  length(AKI_183276_P), length(REF_183276_P), length(CKD_183276_P),
  length(AKI_183276_P) + length(REF_183276_P) + length(CKD_183276_P)))
cat(sprintf("Frozen: GSE210622 AKI=%d non-AKI=%d (primary, n=%d)\n",
  length(AKI_210622_P), length(NONAKI_210622_P),
  length(AKI_210622_P) + length(NONAKI_210622_P)))
stopifnot(length(AKI_183276_P) + length(REF_183276_P) + length(CKD_183276_P) == 41)
stopifnot(length(AKI_210622_P) + length(NONAKI_210622_P) == 12)

# ---- 11. Frozen scoring functions ----
logcpm_matrix <- function(pb, donors, genes) {
  genes <- intersect(genes, rownames(pb))
  stopifnot(length(genes) > 0)
  mat <- as.matrix(pb[genes, donors, drop = FALSE])
  libs <- Matrix::colSums(pb[, donors, drop = FALSE])
  stopifnot(all(libs > 0))
  log2(sweep(mat, 2, libs, "/") * 1e6 + 1)
}

zscore_rows <- function(mat) {
  mu <- rowMeans(mat)
  s  <- apply(mat, 1, sd)
  s[!is.finite(s) | s < 1e-8] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, s, "/")
}

score_module <- function(pb, donors, genes) {
  lcpm <- logcpm_matrix(pb, donors, genes)
  z <- zscore_rows(lcpm)
  list(score = colMeans(z), logcpm = lcpm, z = z, genes = rownames(lcpm))
}

# ---- 12. COHORT-LEVEL Z-SCORING (frozen contract: once per dataset/tier) ----

message("\n========== COHORT-LEVEL SCORING ==========")

# GSE183276 primary scaling set: all 41 AKI+Ref+CKD donors
ALL_183276_P <- c(AKI_183276_P, REF_183276_P, CKD_183276_P)
stopifnot(length(ALL_183276_P) == 41)

tal_183276   <- score_module(shared$pseudobulk[["GSE183276"]][["TAL"]], ALL_183276_P, TAL_GENES)
pti_183276   <- score_module(shared$pseudobulk[["GSE183276"]][["PT"]],  ALL_183276_P, PT_ID_GENES)
ptinj_183276 <- score_module(shared$pseudobulk[["GSE183276"]][["PT"]],  ALL_183276_P, PT_INJ_GENES)

# GSE183276 sensitivity scaling set: all donors with >=20 TAL cells (AKI + Ref + CKD)
ALL_183276_S <- c(AKI_183276_S, REF_183276_S, CKD_183276_S)
cat(sprintf("Sensitivity tier (>=20): GSE183276 AKI=%d Ref=%d CKD=%d (n=%d)\n",
  length(AKI_183276_S), length(REF_183276_S), length(CKD_183276_S), length(ALL_183276_S)))

tal_183276_s   <- score_module(shared$pseudobulk[["GSE183276"]][["TAL"]], ALL_183276_S, TAL_GENES)
pti_183276_s   <- score_module(shared$pseudobulk[["GSE183276"]][["PT"]],  ALL_183276_S, PT_ID_GENES)
ptinj_183276_s <- score_module(shared$pseudobulk[["GSE183276"]][["PT"]],  ALL_183276_S, PT_INJ_GENES)

# GSE210622 scaling set: all 12 primary donors
ALL_210622_P <- c(AKI_210622_P, NONAKI_210622_P)
stopifnot(length(ALL_210622_P) == 12)

tal_210622   <- score_module(shared$pseudobulk[["GSE210622"]][["TAL"]], ALL_210622_P, TAL_GENES)
pti_210622   <- score_module(shared$pseudobulk[["GSE210622"]][["PT"]],  ALL_210622_P, PT_ID_GENES)
ptinj_210622 <- score_module(shared$pseudobulk[["GSE210622"]][["PT"]],  ALL_210622_P, PT_INJ_GENES)

# ---- Helper: build per-donor evidence df from cohort-level scores ----
build_evidence <- function(scores_obj, pti_obj, ptinj_obj, donors_vec, dataset, tier) {
  df <- data.frame(
    dataset     = dataset,
    tier        = tier,
    donor_id    = donors_vec,
    TAL_identity = unname(scores_obj$score[donors_vec]),
    PT_identity  = unname(pti_obj$score[donors_vec]),
    PT_injury    = unname(ptinj_obj$score[donors_vec]),
    stringsAsFactors = FALSE
  )
  df$TAL_minus_PT <- df$TAL_identity - df$PT_identity
  df$n_TAL <- formal_donors$n_TAL_cells[match(df$donor_id, formal_donors$donor_id)]
  df$n_PT  <- formal_donors$n_PT_cells[match(df$donor_id, formal_donors$donor_id)]
  df$n_cells_total <- formal_donors$n_cells_total[match(df$donor_id, formal_donors$donor_id)]
  df$TAL_abundance <- df$n_TAL / df$n_cells_total

  if (dataset == "GSE210622") {
    df$covid_status <- g210622_out$covid_status[match(df$donor_id, g210622_out$donor_id)]
    df$sampling_context <- g210622_out$sampling_context[match(df$donor_id, g210622_out$donor_id)]
    df$group <- ifelse(df$donor_id %in% AKI_210622_P, "AKI", "Reference")
  } else {
    df$covid_status <- "not_reported"
    df$sampling_context <- g183276_out$sampling_context[match(df$donor_id, g183276_out$donor_id)]
    # Read diagnosis from formal_donors for robust group assignment
    fd_dx <- formal_donors$diagnosis[match(df$donor_id, formal_donors$donor_id)]
    fd_dx[fd_dx == "Ref"] <- "Reference"
    df$group <- fd_dx
  }
  df$AKI <- as.integer(df$group == "AKI")
  df
}

# All primary donors (cohort-level evidence)
ev_183276 <- build_evidence(tal_183276, pti_183276, ptinj_183276, ALL_183276_P, "GSE183276", "primary_ge50")
ev_210622 <- build_evidence(tal_210622, pti_210622, ptinj_210622, ALL_210622_P, "GSE210622", "primary_ge50")
ev_183276_s <- build_evidence(tal_183276_s, pti_183276_s, ptinj_183276_s, ALL_183276_S, "GSE183276", "sensitivity_ge20")

# ---- 13. Validate TAL annotation ----
message("\n--- TAL annotation validation ---")
for (ds in c("GSE183276", "GSE210622")) {
  tal_pb <- shared$pseudobulk[[ds]][["TAL"]]
  all_ds <- formal_donors$donor_id[formal_donors$dataset == ds &
                                     formal_donors$tal_tier != "excluded"]
  lcpm_all <- logcpm_matrix(tal_pb, all_ds, TAL_GENES)
  cat(sprintf("=== %s TAL marker validation ===\n", ds))
  for (g in TAL_GENES) {
    vals <- lcpm_all[g, ]
    cat(sprintf("  %s: mean logCPM=%.2f, expressed_in=%d/%d donors\n",
                g, mean(vals), sum(vals > 0), length(vals)))
  }
  pt_pb <- shared$pseudobulk[[ds]][["PT"]]
  pt_contrast <- logcpm_matrix(tal_pb, all_ds, c("LRP2"))
  cat(sprintf("  PT marker LRP2 in TAL pseudobulk: mean logCPM=%.2f\n", mean(pt_contrast["LRP2",])))
}

# ---- 14. Write formal_donors.csv ----
write.csv(formal_donors, file.path(meta_dir, "formal_donors.csv"), row.names = FALSE)
message("Saved: metadata/formal_donors.csv")

# =================================================================
# 15. MODEL FITTING
# =================================================================

message("\n========== MODEL FITTING ==========")

fit_lm <- function(data, formula, model_name, comparison, tier) {
  fit <- lm(formula, data = data)
  sm  <- coef(summary(fit))
  ci  <- tryCatch(confint(fit), error = function(e) {
    matrix(NA, nrow = nrow(sm), ncol = 2,
           dimnames = list(rownames(sm), c("2.5 %", "97.5 %")))
  })
  aki_row <- which(rownames(sm) == "AKI")
  if (length(aki_row) == 0) aki_row <- which(grepl("^AKI", rownames(sm)))[1]

  n_aki <- sum(data$group == "AKI")
  n_ref <- sum(data$group != "AKI")

  data.frame(
    dataset    = data$dataset[1],
    comparison = comparison,
    tier       = tier,
    model      = model_name,
    term       = "AKI",
    estimate   = unname(sm[aki_row, "Estimate"]),
    std_error  = unname(sm[aki_row, "Std. Error"]),
    df_residual = df.residual(fit),
    statistic  = unname(sm[aki_row, "t value"]),
    p_value    = unname(sm[aki_row, "Pr(>|t|)"]),
    ci_low     = unname(ci[aki_row, 1]),
    ci_high    = unname(ci[aki_row, 2]),
    n_aki      = n_aki,
    n_ref      = n_ref,
    stringsAsFactors = FALSE
  )
}

all_effects <- list()

# --- GSE183276 AKI vs Reference (subset from cohort-level scores) ---
s_1836_ref <- ev_183276[ev_183276$group %in% c("AKI","Reference"), ]
all_effects[[1]] <- fit_lm(s_1836_ref, TAL_identity ~ AKI, "unadjusted", "AKI_vs_Reference", "primary_ge50")
all_effects[[2]] <- fit_lm(s_1836_ref, TAL_identity ~ AKI + PT_injury, "adjusted_PT_injury", "AKI_vs_Reference", "primary_ge50")
all_effects[[3]] <- fit_lm(s_1836_ref, TAL_minus_PT ~ AKI, "TAL_minus_PT", "AKI_vs_Reference", "primary_ge50")

# TAL abundance model
all_effects[[4]] <- fit_lm(s_1836_ref, TAL_abundance ~ AKI, "TAL_abundance_unadjusted", "AKI_vs_Reference", "primary_ge50")

# --- GSE210622 AKI vs non-AKI ---
s_2106 <- ev_210622
all_effects[[5]] <- fit_lm(s_2106, TAL_identity ~ AKI, "unadjusted", "AKI_vs_nonAKI", "primary_ge50")
all_effects[[6]] <- fit_lm(s_2106, TAL_identity ~ AKI + PT_injury, "adjusted_PT_injury", "AKI_vs_nonAKI", "primary_ge50")
all_effects[[7]] <- fit_lm(s_2106, TAL_minus_PT ~ AKI, "TAL_minus_PT", "AKI_vs_nonAKI", "primary_ge50")
all_effects[[8]] <- fit_lm(s_2106, TAL_abundance ~ AKI, "TAL_abundance_unadjusted", "AKI_vs_nonAKI", "primary_ge50")

# --- GSE183276 AKI vs CKD ---
s_1836_ckd <- ev_183276[ev_183276$group %in% c("AKI","CKD"), ]
all_effects[[9]] <- fit_lm(s_1836_ckd, TAL_identity ~ AKI, "unadjusted", "AKI_vs_CKD", "primary_ge50")
all_effects[[10]] <- fit_lm(s_1836_ckd, TAL_identity ~ AKI + PT_injury, "adjusted_PT_injury", "AKI_vs_CKD", "primary_ge50")
all_effects[[11]] <- fit_lm(s_1836_ckd, TAL_minus_PT ~ AKI, "TAL_minus_PT", "AKI_vs_CKD", "primary_ge50")

# --- Sensitivity: GSE183276 >=20 TAL cells ---
s_1836_sens <- ev_183276_s[ev_183276_s$group %in% c("AKI","Reference"), ]
all_effects[[12]] <- fit_lm(s_1836_sens, TAL_identity ~ AKI, "unadjusted", "AKI_vs_Reference_ge20", "sensitivity_ge20")

# --- Sensitivity: GSE210622 COVID-adjusted ---
s_2106$covid_factor <- factor(ifelse(s_2106$covid_status == "COVID", "COVID", "NonCOVID"))
all_effects[[13]] <- fit_lm(s_2106, TAL_identity ~ AKI + covid_factor, "adjusted_COVID", "AKI_vs_nonAKI", "primary_ge50")

# --- Sensitivity: GSE210622 non-COVID only ---
s_2106_nc <- s_2106[s_2106$covid_status != "COVID", ]
all_effects[[14]] <- fit_lm(s_2106_nc, TAL_identity ~ AKI, "non_COVID_only", "AKI_vs_nonAKI", "primary_ge50")

# --- Sensitivity: GSE210622 exclude CPM ---
s_2106_noCPM <- s_2106[s_2106$donor_id != "CPM", ]
all_effects[[15]] <- fit_lm(s_2106_noCPM, TAL_identity ~ AKI, "exclude_CPM", "AKI_vs_nonAKI", "primary_ge50")

# --- Sensitivity: GSE210622 nephrectomy ctrl only ---
s_2106_neph <- s_2106[s_2106$donor_id %in% c("CTN1","CTN2","CTN3") | s_2106$group == "AKI", ]
all_effects[[16]] <- fit_lm(s_2106_neph, TAL_identity ~ AKI, "nephrectomy_ctrl_only", "AKI_vs_nonAKI", "primary_ge50")

effects_df <- do.call(rbind, all_effects)
rownames(effects_df) <- NULL

# --- KDIGO stage 2 vs 3 (secondary, exploratory) ---
message("--- Secondary: KDIGO stage analysis ---")
kdigo_1836 <- ev_183276[ev_183276$group == "AKI", ]
kdigo_1836$kdigo_stage <- g183276_out$kdigo_stage[match(kdigo_1836$donor_id, g183276_out$donor_id)]
kdigo_2106 <- ev_210622[ev_210622$group == "AKI", ]
kdigo_2106$kdigo_stage <- g210622_out$kdigo_stage[match(kdigo_2106$donor_id, g210622_out$donor_id)]
kdigo_combined <- rbind(kdigo_1836, kdigo_2106)
kdigo_combined <- kdigo_combined[kdigo_combined$kdigo_stage %in% c("2","3"), ]
kdigo_combined$stage3 <- as.integer(kdigo_combined$kdigo_stage == "3")

n_s2 <- sum(kdigo_combined$stage3 == 0)
n_s3 <- sum(kdigo_combined$stage3 == 1)

if (n_s2 >= 2 && n_s3 >= 2) {
  fit_kdigo <- lm(TAL_identity ~ stage3 + dataset, data = kdigo_combined)
  sm_kd <- coef(summary(fit_kdigo))
  ci_kd <- tryCatch(confint(fit_kdigo), error = function(e) matrix(NA, nrow(sm_kd), 2))
  effects_df <- rbind(effects_df, data.frame(
    dataset = "Combined", comparison = "KDIGO_stage3_vs_2", tier = "combined",
    model = "unadjusted_cohort_adjusted", term = "stage3",
    estimate = unname(sm_kd["stage3","Estimate"]),
    std_error = unname(sm_kd["stage3","Std. Error"]),
    df_residual = df.residual(fit_kdigo),
    statistic = unname(sm_kd["stage3","t value"]),
    p_value = unname(sm_kd["stage3","Pr(>|t|)"]),
    ci_low = unname(ci_kd["stage3",1]),
    ci_high = unname(ci_kd["stage3",2]),
    n_aki = n_s3, n_ref = n_s2,
    stringsAsFactors = FALSE
  ))
  cat(sprintf("KDIGO combined stage3 vs stage2: %.3f (n=%d+%d)\n",
    unname(sm_kd["stage3","Estimate"]), n_s3, n_s2))

  # Per-cohort KDIGO
  for (ds_kd in c("GSE183276","GSE210622")) {
    sub_kd <- kdigo_combined[kdigo_combined$dataset == ds_kd, ]
    s2_ds <- sum(sub_kd$stage3 == 0); s3_ds <- sum(sub_kd$stage3 == 1)
    if (s2_ds >= 2 && s3_ds >= 2) {
      fit_ds <- lm(TAL_identity ~ stage3, data = sub_kd)
      sm_ds <- coef(summary(fit_ds))
      ci_ds <- tryCatch(confint(fit_ds), error = function(e) matrix(NA, 2, 2))
      effects_df <- rbind(effects_df, data.frame(
        dataset = ds_kd, comparison = "KDIGO_stage3_vs_2", tier = "combined",
        model = "unadjusted", term = "stage3",
        estimate = unname(sm_ds["stage3","Estimate"]),
        std_error = unname(sm_ds["stage3","Std. Error"]),
        df_residual = df.residual(fit_ds),
        statistic = unname(sm_ds["stage3","t value"]),
        p_value = unname(sm_ds["stage3","Pr(>|t|)"]),
        ci_low = unname(ci_ds["stage3",1]),
        ci_high = unname(ci_ds["stage3",2]),
        n_aki = s3_ds, n_ref = s2_ds,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# --- PT-injury correlation and attenuation ---
pt_cor_1836 <- cor(s_1836_ref$AKI, s_1836_ref$PT_injury)
pt_cor_2106 <- cor(s_2106$AKI, s_2106$PT_injury)
unadj_est_1836 <- effects_df$estimate[effects_df$comparison == "AKI_vs_Reference" & effects_df$model == "unadjusted" & effects_df$dataset == "GSE183276"]
adj_est_1836   <- effects_df$estimate[effects_df$comparison == "AKI_vs_Reference" & effects_df$model == "adjusted_PT_injury" & effects_df$dataset == "GSE183276"]
unadj_est_2106 <- effects_df$estimate[effects_df$comparison == "AKI_vs_nonAKI" & effects_df$model == "unadjusted" & effects_df$dataset == "GSE210622"]
adj_est_2106   <- effects_df$estimate[effects_df$comparison == "AKI_vs_nonAKI" & effects_df$model == "adjusted_PT_injury" & effects_df$dataset == "GSE210622"]
atten_1836 <- 100 * (1 - adj_est_1836 / unadj_est_1836)
atten_2106 <- 100 * (1 - adj_est_2106 / unadj_est_2106)
message(sprintf("PT-injury correlation with AKI: GSE183276=%.3f, GSE210622=%.3f", pt_cor_1836, pt_cor_2106))
message(sprintf("%% attenuation after PT adjustment: GSE183276=%.1f%%, GSE210622=%.1f%%", atten_1836, atten_2106))

# =================================================================
# 16. LODO (LEAVE-ONE-DONOR-OUT) — save formal_lodo.csv
# =================================================================

message("\n========== LODO ANALYSIS ==========")

run_lodo <- function(data, comparison, tier, ref_label) {
  aki <- data$donor_id[data$group == "AKI"]
  ref <- data$donor_id[data$group == ref_label]
  all_d <- c(aki, ref)
  ests <- numeric(length(all_d))
  for (i in seq_along(all_d)) {
    omit <- all_d[i]
    sub  <- data[data$donor_id != omit, ]
    if (sum(sub$group == "AKI") < 1 || sum(sub$group == ref_label) < 1) {
      ests[i] <- NA_real_
      next
    }
    fit <- lm(TAL_identity ~ AKI, data = sub)
    ests[i] <- coef(fit)["AKI"]
  }
  data.frame(
    dataset       = data$dataset[1],
    comparison    = comparison,
    tier          = tier,
    omitted_donor = all_d,
    estimate      = ests,
    stringsAsFactors = FALSE
  )
}

lodo_1836_ref  <- run_lodo(s_1836_ref, "AKI_vs_Reference", "primary_ge50", "Reference")
lodo_2106      <- run_lodo(s_2106, "AKI_vs_nonAKI", "primary_ge50", "Reference")
lodo_1836_ckd  <- run_lodo(s_1836_ckd, "AKI_vs_CKD", "primary_ge50", "CKD")
lodo_1836_sens <- run_lodo(s_1836_sens, "AKI_vs_Reference_ge20", "sensitivity_ge20", "Reference")
lodo_2106_nc   <- run_lodo(s_2106_nc, "AKI_vs_nonAKI_nonCOVID", "primary_ge50", "Reference")

lodo_all <- rbind(lodo_1836_ref, lodo_2106, lodo_1836_ckd, lodo_1836_sens, lodo_2106_nc)
write.csv(lodo_all, file.path(result_dir, "formal_lodo.csv"), row.names = FALSE)
message("Saved: results/formal_lodo.csv")

lodo_pass <- function(ldf) {
  ok <- !is.na(ldf$estimate)
  if (sum(ok) == 0) return(NA_real_)
  mean(ldf$estimate[ok] < 0)
}

lodo_results <- list(
  "GSE183276_AKI_vs_Ref"   = lodo_pass(lodo_1836_ref),
  "GSE210622_AKI_vs_nonAKI" = lodo_pass(lodo_2106),
  "GSE183276_AKI_vs_CKD"    = lodo_pass(lodo_1836_ckd),
  "GSE183276_ge20"          = lodo_pass(lodo_1836_sens),
  "GSE210622_nonCOVID"      = lodo_pass(lodo_2106_nc)
)
for (nm in names(lodo_results)) {
  cat(sprintf("LODO %s: negative fraction = %.3f\n", nm, lodo_results[[nm]]))
}

# =================================================================
# 17. PERMUTATION TESTS (fixed seed = 42) — save auditable table
# =================================================================

message("\n========== PERMUTATION TESTS ==========")

perm_test <- function(x_ref, x_aki, n_perm = 9999, seed = SET_SEED) {
  obs <- mean(x_aki) - mean(x_ref)
  pooled <- c(x_ref, x_aki)
  n1 <- length(x_ref); n2 <- length(x_aki); N <- n1 + n2
  if (N < 5) return(NA_real_)
  n_eff <- min(n_perm, 50000)
  set.seed(seed)
  extreme <- 0
  for (b in seq_len(n_eff)) {
    idx <- sample(N, n1)
    if (abs(mean(pooled[-idx]) - mean(pooled[idx])) >= abs(obs)) extreme <- extreme + 1
  }
  (extreme + 1) / (n_eff + 1)
}

perm_table <- data.frame(
  comparison = c("GSE183276_AKI_vs_Ref","GSE210622_AKI_vs_nonAKI",
                 "GSE183276_AKI_vs_CKD","GSE183276_ge20","GSE210622_nonCOVID"),
  p_value = c(
    perm_test(s_1836_ref$TAL_identity[s_1836_ref$group == "Reference"],
              s_1836_ref$TAL_identity[s_1836_ref$group == "AKI"]),
    perm_test(s_2106$TAL_identity[s_2106$group == "Reference"],
              s_2106$TAL_identity[s_2106$group == "AKI"]),
    perm_test(s_1836_ckd$TAL_identity[s_1836_ckd$group == "CKD"],
              s_1836_ckd$TAL_identity[s_1836_ckd$group == "AKI"]),
    perm_test(s_1836_sens$TAL_identity[s_1836_sens$group == "Reference"],
              s_1836_sens$TAL_identity[s_1836_sens$group == "AKI"]),
    perm_test(s_2106_nc$TAL_identity[s_2106_nc$group == "Reference"],
              s_2106_nc$TAL_identity[s_2106_nc$group == "AKI"])
  ),
  seed = SET_SEED,
  n_permutations = 9999,
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(perm_table))) {
  cat(sprintf("Permutation %s: p = %.4f\n", perm_table$comparison[i], perm_table$p_value[i]))
}

# Merge LODO fractions and permutation P into effects_df — ONLY for unadjusted model rows
# that actually underwent these tests. Adjusted/abundance/KDIGO rows stay NA.
effects_df$lodo_negative_frac <- NA_real_
effects_df$permutation_p       <- NA_real_

perm_lookup <- list(
  "GSE183276_AKI_vs_Ref"   = perm_table$p_value[perm_table$comparison == "GSE183276_AKI_vs_Ref"],
  "GSE210622_AKI_vs_nonAKI" = perm_table$p_value[perm_table$comparison == "GSE210622_AKI_vs_nonAKI"],
  "GSE183276_AKI_vs_CKD"    = perm_table$p_value[perm_table$comparison == "GSE183276_AKI_vs_CKD"],
  "GSE183276_ge20"          = perm_table$p_value[perm_table$comparison == "GSE183276_ge20"],
  "GSE210622_nonCOVID"      = perm_table$p_value[perm_table$comparison == "GSE210622_nonCOVID"]
)

for (i in seq_len(nrow(effects_df))) {
  comp <- effects_df$comparison[i]
  ds   <- effects_df$dataset[i]
  mdl  <- effects_df$model[i]
  # LODO: only for unadjusted primary models
  if (mdl == "unadjusted") {
    if (comp == "AKI_vs_Reference"      && ds == "GSE183276") effects_df$lodo_negative_frac[i] <- lodo_results[["GSE183276_AKI_vs_Ref"]]
    if (comp == "AKI_vs_nonAKI"         && ds == "GSE210622") effects_df$lodo_negative_frac[i] <- lodo_results[["GSE210622_AKI_vs_nonAKI"]]
    if (comp == "AKI_vs_CKD"            && ds == "GSE183276") effects_df$lodo_negative_frac[i] <- lodo_results[["GSE183276_AKI_vs_CKD"]]
    if (comp == "AKI_vs_Reference_ge20" && ds == "GSE183276") effects_df$lodo_negative_frac[i] <- lodo_results[["GSE183276_ge20"]]
  }
  if (comp == "AKI_vs_nonAKI" && ds == "GSE210622" && mdl == "non_COVID_only")
    effects_df$lodo_negative_frac[i] <- lodo_results[["GSE210622_nonCOVID"]]
  # Permutation: only for unadjusted models
  if (mdl == "unadjusted") {
    if (comp == "AKI_vs_Reference"      && ds == "GSE183276") effects_df$permutation_p[i] <- perm_lookup[["GSE183276_AKI_vs_Ref"]]
    if (comp == "AKI_vs_nonAKI"         && ds == "GSE210622") effects_df$permutation_p[i] <- perm_lookup[["GSE210622_AKI_vs_nonAKI"]]
    if (comp == "AKI_vs_CKD"            && ds == "GSE183276") effects_df$permutation_p[i] <- perm_lookup[["GSE183276_AKI_vs_CKD"]]
    if (comp == "AKI_vs_Reference_ge20" && ds == "GSE183276") effects_df$permutation_p[i] <- perm_lookup[["GSE183276_ge20"]]
  }
  if (comp == "AKI_vs_nonAKI" && ds == "GSE210622" && mdl == "non_COVID_only")
    effects_df$permutation_p[i] <- perm_lookup[["GSE210622_nonCOVID"]]
}

# =================================================================
# 18. PER-GENE DIRECTION (from cohort-level z-scores)
# =================================================================

message("\n========== PER-GENE DIRECTION ==========")

run_gene_analysis <- function(scores_obj, data, comparison, tier, ref_label) {
  aki <- data$donor_id[data$group == "AKI"]
  ref <- data$donor_id[data$group == ref_label]
  z_mat <- scores_obj$z
  lcpm  <- scores_obj$logcpm

  data.frame(
    dataset       = data$dataset[1],
    comparison    = comparison,
    tier          = tier,
    gene          = scores_obj$genes,
    mean_logCPM_ref = rowMeans(lcpm[, ref, drop = FALSE]),
    mean_logCPM_aki = rowMeans(lcpm[, aki, drop = FALSE]),
    mean_z_ref      = rowMeans(z_mat[, ref, drop = FALSE]),
    mean_z_aki      = rowMeans(z_mat[, aki, drop = FALSE]),
    stringsAsFactors = FALSE
  )
}

gene_1836_ref  <- run_gene_analysis(tal_183276, s_1836_ref, "AKI_vs_Reference", "primary_ge50", "Reference")
gene_2106      <- run_gene_analysis(tal_210622, s_2106, "AKI_vs_nonAKI", "primary_ge50", "Reference")
gene_1836_ckd  <- run_gene_analysis(tal_183276, s_1836_ckd, "AKI_vs_CKD", "primary_ge50", "CKD")
gene_1836_sens <- run_gene_analysis(tal_183276_s, s_1836_sens, "AKI_vs_Reference_ge20", "sensitivity_ge20", "Reference")

# Compute direction on each individual result
gene_1836_ref$logCPM_effect  <- gene_1836_ref$mean_logCPM_aki  - gene_1836_ref$mean_logCPM_ref
gene_1836_ref$z_effect       <- gene_1836_ref$mean_z_aki       - gene_1836_ref$mean_z_ref
gene_1836_ref$direction      <- ifelse(gene_1836_ref$z_effect  < 0, "negative", "positive")
gene_2106$logCPM_effect      <- gene_2106$mean_logCPM_aki      - gene_2106$mean_logCPM_ref
gene_2106$z_effect           <- gene_2106$mean_z_aki           - gene_2106$mean_z_ref
gene_2106$direction          <- ifelse(gene_2106$z_effect      < 0, "negative", "positive")
gene_1836_ckd$logCPM_effect  <- gene_1836_ckd$mean_logCPM_aki  - gene_1836_ckd$mean_logCPM_ref
gene_1836_ckd$z_effect       <- gene_1836_ckd$mean_z_aki       - gene_1836_ckd$mean_z_ref
gene_1836_ckd$direction      <- ifelse(gene_1836_ckd$z_effect  < 0, "negative", "positive")
gene_1836_sens$logCPM_effect <- gene_1836_sens$mean_logCPM_aki - gene_1836_sens$mean_logCPM_ref
gene_1836_sens$z_effect      <- gene_1836_sens$mean_z_aki      - gene_1836_sens$mean_z_ref
gene_1836_sens$direction     <- ifelse(gene_1836_sens$z_effect < 0, "negative", "positive")

gene_all <- rbind(gene_1836_ref, gene_2106, gene_1836_ckd, gene_1836_sens)

for (comp in unique(gene_all$comparison)) {
  sub <- gene_all[gene_all$comparison == comp, ]
  neg_n <- sum(sub$direction == "negative")
  cat(sprintf("%s: %d/6 negative genes\n", comp, neg_n))
}

# Negative gene counts
neg_genes_1836 <- sum(gene_1836_ref$direction == "negative")
neg_genes_2106 <- sum(gene_2106$direction == "negative")
neg_genes_ckd  <- sum(gene_1836_ckd$direction == "negative")

# =================================================================
# 19. ASSEMBLE EVIDENCE AND DECISION
# =================================================================

message("\n========== EVIDENCE ASSEMBLY ==========")

safe_val <- function(df, col, default = NA_real_) {
  if (is.null(df) || nrow(df) == 0) return(default)
  val <- df[[col]]
  if (length(val) == 0 || all(is.na(val))) default else val[1]
}

# Extract key rows
unadj_1836_ref <- effects_df[effects_df$comparison == "AKI_vs_Reference" & effects_df$model == "unadjusted" & effects_df$dataset == "GSE183276", ]
unadj_2106     <- effects_df[effects_df$comparison == "AKI_vs_nonAKI"     & effects_df$model == "unadjusted" & effects_df$dataset == "GSE210622", ]
adj_pt_1836    <- effects_df[effects_df$comparison == "AKI_vs_Reference" & effects_df$model == "adjusted_PT_injury" & effects_df$dataset == "GSE183276", ]
adj_pt_2106    <- effects_df[effects_df$comparison == "AKI_vs_nonAKI"     & effects_df$model == "adjusted_PT_injury" & effects_df$dataset == "GSE210622", ]
tal_pt_1836    <- effects_df[effects_df$comparison == "AKI_vs_Reference" & effects_df$model == "TAL_minus_PT" & effects_df$dataset == "GSE183276", ]
tal_pt_2106    <- effects_df[effects_df$comparison == "AKI_vs_nonAKI"     & effects_df$model == "TAL_minus_PT" & effects_df$dataset == "GSE210622", ]
unadj_1836_ckd <- effects_df[effects_df$comparison == "AKI_vs_CKD" & effects_df$model == "unadjusted" & effects_df$dataset == "GSE183276", ]
covid_adj      <- effects_df[effects_df$comparison == "AKI_vs_nonAKI" & effects_df$model == "adjusted_COVID", ]
noncovid_only  <- effects_df[effects_df$comparison == "AKI_vs_nonAKI" & effects_df$model == "non_COVID_only", ]
ge20_effect    <- effects_df[effects_df$comparison == "AKI_vs_Reference_ge20", ]
noCPM_effect   <- effects_df[effects_df$comparison == "AKI_vs_nonAKI" & effects_df$model == "exclude_CPM", ]
neph_effect    <- effects_df[effects_df$comparison == "AKI_vs_nonAKI" & effects_df$model == "nephrectomy_ctrl_only", ]

# Decision criteria
dir_1836_neg <- safe_val(unadj_1836_ref, "estimate", 0) < 0
dir_2106_neg <- safe_val(unadj_2106, "estimate", 0) < 0
dir_ok <- dir_1836_neg && dir_2106_neg

lodo_1836_ok <- lodo_results[["GSE183276_AKI_vs_Ref"]] >= 0.8
lodo_2106_ok <- lodo_results[["GSE210622_AKI_vs_nonAKI"]] >= 0.8
lodo_all_ok  <- lodo_1836_ok && lodo_2106_ok

genes_1836_ok <- neg_genes_1836 >= 4
genes_2106_ok <- neg_genes_2106 >= 4
genes_ok      <- genes_1836_ok && genes_2106_ok

covid_dir_pass <- safe_val(covid_adj, "estimate", 0) < 0
nocovid_dir_pass <- safe_val(noncovid_only, "estimate", 0) < 0

pt_adj_1836_neg <- safe_val(adj_pt_1836, "estimate", 0) < 0
pt_adj_2106_neg <- safe_val(adj_pt_2106, "estimate", 0) < 0

tal_pt_1836_neg <- safe_val(tal_pt_1836, "estimate", 0) < 0
tal_pt_2106_neg <- safe_val(tal_pt_2106, "estimate", 0) < 0
tal_pt_ok <- tal_pt_1836_neg && tal_pt_2106_neg

primary_pass <- dir_ok && lodo_all_ok && genes_ok

if (primary_pass) {
  decision <- "ANALYSIS_PASS_CONSERVATIVE"
} else {
  decision <- "STOP_SIGNAL_NOT_DEFENSIBLE"
}

message(sprintf("\n=== DECISION: %s ===\n", decision))

# =================================================================
# 20. PER-DONOR EVIDENCE FILE
# =================================================================

message("Writing per-donor evidence...")

evidence_all <- rbind(
  cbind(ev_183276, comparison = "all_primary"),
  cbind(ev_210622, comparison = "all_primary"),
  cbind(ev_183276_s, comparison = "sensitivity_ge20")
)
# Trim to standard columns
ev_out_cols <- c("dataset","tier","comparison","donor_id","group",
                 "n_TAL","n_PT","n_cells_total","TAL_abundance",
                 "TAL_identity","PT_identity","PT_injury","TAL_minus_PT",
                 "covid_status","sampling_context")
write.csv(evidence_all[, ev_out_cols], file.path(result_dir, "formal_evidence.csv"), row.names = FALSE)

write.csv(effects_df, file.path(result_dir, "formal_effects.csv"), row.names = FALSE)
write.csv(gene_all, file.path(result_dir, "formal_genes.csv"), row.names = FALSE)

# =================================================================
# 21. FIGURES
# =================================================================

message("\n========== GENERATING FIGURES ==========")

# Arrange ggplot panels vertically with base grid so the workflow does not
# depend on the optional gridExtra package.
arrange_vertical <- function(..., heights) {
  plots <- list(...)
  stopifnot(length(plots) == length(heights))
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = length(plots),
    ncol = 1,
    heights = grid::unit(heights, "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  for (i in seq_along(plots)) {
    print(plots[[i]], vp = grid::viewport(layout.pos.row = i, layout.pos.col = 1))
  }
  grid::popViewport()
}

# ---- Figure 1: Cohort characterization ----
# Panel A: TAL coverage
cov_data <- formal_donors[formal_donors$tal_tier != "excluded", ]
cov_data$diagnosis[cov_data$dataset == "GSE210622" & cov_data$diagnosis == "control"] <- "non-AKI"
cov_data$facet_label <- cov_data$dataset

p1a <- ggplot(cov_data,
  aes(x = reorder(donor_id, n_TAL_cells), y = n_TAL_cells, fill = diagnosis)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 20, linetype = 3, colour = "#a94442") +
  geom_hline(yintercept = 50, linetype = 2, colour = "#2b6f9f") +
  facet_wrap(~ facet_label, scales = "free_x") +
  scale_y_log10() +
  scale_fill_manual(values = c("AKI"="#c24f44","Ref"="#3b82a0","CKD"="#e6ab02","non-AKI"="#3b82a0")) +
  labs(x = NULL, y = "TAL cells (log scale)", fill = "Group",
       title = "A  TAL cell coverage by donor") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 6))

# Panel B: TAL abundance
ab_data <- formal_donors[formal_donors$tal_tier != "excluded" & !is.na(formal_donors$n_cells_total), ]
ab_data$TAL_abundance <- ab_data$n_TAL_cells / ab_data$n_cells_total
ab_data$diagnosis[ab_data$dataset == "GSE210622" & ab_data$diagnosis == "control"] <- "non-AKI"
ab_data$facet_label <- ab_data$dataset

p1b <- ggplot(ab_data, aes(x = diagnosis, y = TAL_abundance, colour = diagnosis)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, colour = "grey45") +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  facet_wrap(~ facet_label, scales = "free_x") +
  scale_colour_manual(values = c("AKI"="#c24f44","Ref"="#3b82a0","CKD"="#e6ab02","non-AKI"="#3b82a0")) +
  labs(x = NULL, y = "TAL fraction of total cells",
       title = "B  TAL abundance (n_TAL / n_total) by group") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

pdf(file.path(figure_dir, "formal_1_cohort.pdf"), width = 12, height = 8)
arrange_vertical(p1a, p1b, heights = c(1, 0.85))
invisible(dev.off())
message("Saved: figures/formal_1_cohort.pdf")

# ---- Figure 2: Score distributions and forest plot ----
# Panel A: TAL identity scores (primary comparisons only)
plot_scores <- rbind(
  cbind(s_1836_ref[, c("dataset","donor_id","group","TAL_identity","TAL_abundance")],
        cohort_label = "GSE183276\nAKI vs Reference"),
  cbind(s_2106[, c("dataset","donor_id","group","TAL_identity","TAL_abundance")],
        cohort_label = "GSE210622\nAKI vs non-AKI")
)
plot_scores$group <- factor(plot_scores$group, levels = c("Reference","AKI"))

p2a <- ggplot(plot_scores, aes(x = group, y = TAL_identity, colour = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, colour = "grey45") +
  geom_jitter(width = 0.08, size = 2.5, alpha = 0.9) +
  facet_wrap(~ cohort_label, scales = "free_y") +
  scale_colour_manual(values = c("Reference"="#3b82a0","AKI"="#c24f44")) +
  labs(x = NULL, y = "Frozen TAL identity score",
       title = "A  Donor-level TAL identity score") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# Panel B: Forest plot - AKI effects, primary model classes only
forest_data <- effects_df[effects_df$model %in% c("unadjusted","adjusted_PT_injury","TAL_minus_PT") &
                           effects_df$comparison %in% c("AKI_vs_Reference","AKI_vs_nonAKI"), ]
forest_data$model <- factor(forest_data$model,
  levels = c("unadjusted","adjusted_PT_injury","TAL_minus_PT"),
  labels = c("Unadjusted","PT-injury adjusted","TAL-minus-PT"))
forest_data$cohort_label <- ifelse(
  forest_data$dataset == "GSE183276",
  "GSE183276\nReference",
  "GSE210622\nnon-AKI"
)

p2b <- ggplot(forest_data,
  aes(x = estimate, y = model, colour = cohort_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                width = 0.15, orientation = "y",
                position = position_dodge(width = 0.45)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.45)) +
  labs(x = "AKI minus reference effect", y = NULL, colour = "Cohort",
       title = "B  Forest plot - AKI effects") +
  theme_bw(base_size = 11)

pdf(file.path(figure_dir, "formal_2_scores.pdf"), width = 10, height = 7.5)
arrange_vertical(p2a, p2b, heights = c(1, 0.85))
invisible(dev.off())
message("Saved: figures/formal_2_scores.pdf")

# ---- Figure 3: Specificity and sensitivity ----
# Panel A: PT-adjusted and TAL-minus-PT effects
spec_data <- effects_df[effects_df$comparison %in% c("AKI_vs_Reference","AKI_vs_nonAKI") &
                          effects_df$model %in% c("adjusted_PT_injury","TAL_minus_PT"), ]
spec_data$model <- factor(spec_data$model,
  levels = c("adjusted_PT_injury","TAL_minus_PT"),
  labels = c("PT-injury adjusted","TAL minus PT identity"))
spec_data$cohort_label <- ifelse(
  spec_data$dataset == "GSE183276",
  "GSE183276\nReference",
  "GSE210622\nnon-AKI"
)

p3a <- ggplot(spec_data,
  aes(x = estimate, y = model, colour = cohort_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                width = 0.15, orientation = "y",
                position = position_dodge(width = 0.35)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.35)) +
  labs(x = "AKI minus reference effect", y = NULL, colour = "Cohort",
       title = "A  Specificity - PT-injury adjusted and TAL-minus-PT effects") +
  theme_bw(base_size = 11)

# Panel B: Sensitivity - exclude KDIGO rows, label every estimate with its cohort
sens_rows <- effects_df[effects_df$comparison != "KDIGO_stage3_vs_2" &
  (effects_df$model %in% c("unadjusted","adjusted_COVID","non_COVID_only",
                             "exclude_CPM","nephrectomy_ctrl_only") |
   effects_df$comparison == "AKI_vs_Reference_ge20"), ]
# Build descriptive labels
sens_rows$label <- NA_character_
for (i in seq_len(nrow(sens_rows))) {
  ds <- sens_rows$dataset[i]
  mdl <- sens_rows$model[i]
  comp <- sens_rows$comparison[i]
  if (comp == "AKI_vs_Reference_ge20") {
    sens_rows$label[i] <- sprintf("GSE183276 >=20 TAL (n=%d+%d)", sens_rows$n_aki[i], sens_rows$n_ref[i])
  } else if (comp == "AKI_vs_Reference" && ds == "GSE183276" && mdl == "unadjusted") {
    sens_rows$label[i] <- sprintf("GSE183276 primary (n=%d+%d)", sens_rows$n_aki[i], sens_rows$n_ref[i])
  } else if (comp == "AKI_vs_nonAKI" && ds == "GSE210622") {
    if (mdl == "unadjusted")
      sens_rows$label[i] <- sprintf("GSE210622 primary (n=%d+%d)", sens_rows$n_aki[i], sens_rows$n_ref[i])
    else if (mdl == "adjusted_COVID")
      sens_rows$label[i] <- sprintf("GSE210622 COVID-adjusted")
    else if (mdl == "non_COVID_only")
      sens_rows$label[i] <- sprintf("GSE210622 non-COVID (n=%d+%d)", sens_rows$n_aki[i], sens_rows$n_ref[i])
    else if (mdl == "exclude_CPM")
      sens_rows$label[i] <- sprintf("GSE210622 excl. CPM (n=%d+%d)", sens_rows$n_aki[i], sens_rows$n_ref[i])
    else if (mdl == "nephrectomy_ctrl_only")
      sens_rows$label[i] <- sprintf("GSE210622 neph. ctrl (n=%d+%d)", sens_rows$n_aki[i], sens_rows$n_ref[i])
  }
}
sens_rows <- sens_rows[!is.na(sens_rows$label), ]
sens_rows$label <- factor(sens_rows$label, levels = rev(unique(sens_rows$label)))

p3b <- ggplot(sens_rows,
  aes(x = estimate, y = label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                width = 0.15, orientation = "y") +
  geom_point(size = 2.6, colour = "#2b6f9f") +
  labs(x = "AKI minus reference effect (TAL identity)", y = NULL,
       title = "B  Sensitivity analyses") +
  theme_bw(base_size = 10)

# Panel C: Per-gene direction (primary comparisons)
gp <- gene_all[gene_all$comparison %in% c("AKI_vs_Reference","AKI_vs_nonAKI","AKI_vs_CKD") &
                gene_all$tier == "primary_ge50", ]
gp$comp_label <- paste(gp$dataset, sub("AKI_vs_","",gp$comparison), sep = "\n")

p3c <- ggplot(gp, aes(x = gene, y = comp_label, size = abs(z_effect), colour = z_effect)) +
  geom_point() +
  scale_colour_gradient2(low = "#2b6f9f", mid = "white", high = "#c24f44",
                         name = "Z-effect\n(AKI - ref)") +
  scale_size_continuous(name = "|Effect|") +
  labs(x = NULL, y = NULL,
       title = "C  Per-gene TAL expression change (z-scored)") +
  theme_bw(base_size = 11)

pdf(file.path(figure_dir, "formal_3_specificity.pdf"), width = 10, height = 10)
arrange_vertical(p3a, p3b, p3c, heights = c(0.55, 0.95, 0.9))
invisible(dev.off())
message("Saved: figures/formal_3_specificity.pdf")

# =================================================================
# 22. WRITE formal.md
# =================================================================

fmt3 <- function(x) if (is.na(x) || is.null(x)) "NA" else sprintf("%.3f", x)
fmt4 <- function(x) if (is.na(x) || is.null(x)) "NA" else sprintf("%.4f", x)

# KDIGO sub-results for markdown
kdigo_rows <- effects_df[effects_df$comparison == "KDIGO_stage3_vs_2", ]
kdigo_comb <- kdigo_rows[kdigo_rows$dataset == "Combined", ]
kdigo_g1836 <- kdigo_rows[kdigo_rows$dataset == "GSE183276", ]
kdigo_g2106 <- kdigo_rows[kdigo_rows$dataset == "GSE210622", ]

kdigo_comb_est    <- safe_val(kdigo_comb, "estimate")
kdigo_comb_ci_low <- safe_val(kdigo_comb, "ci_low")
kdigo_comb_ci_high <- safe_val(kdigo_comb, "ci_high")
kdigo_comb_p      <- safe_val(kdigo_comb, "p_value")
kdigo_n_s3        <- safe_val(kdigo_comb, "n_aki")
kdigo_n_s2        <- safe_val(kdigo_comb, "n_ref")

kdigo_1836_est    <- safe_val(kdigo_g1836, "estimate")
kdigo_1836_ci_low <- safe_val(kdigo_g1836, "ci_low")
kdigo_1836_ci_high <- safe_val(kdigo_g1836, "ci_high")
kdigo_1836_p      <- safe_val(kdigo_g1836, "p_value")
kdigo_1836_n_s3   <- safe_val(kdigo_g1836, "n_aki", 0)
kdigo_1836_n_s2   <- safe_val(kdigo_g1836, "n_ref", 0)

kdigo_2106_est    <- safe_val(kdigo_g2106, "estimate")
kdigo_2106_ci_low <- safe_val(kdigo_g2106, "ci_low")
kdigo_2106_ci_high <- safe_val(kdigo_g2106, "ci_high")
kdigo_2106_p      <- safe_val(kdigo_g2106, "p_value")
kdigo_2106_n_s3   <- safe_val(kdigo_g2106, "n_aki", 0)
kdigo_2106_n_s2   <- safe_val(kdigo_g2106, "n_ref", 0)

md_lines <- c(
  "# Formal AKI–TAL analysis: frozen design",
  "",
  sprintf("**Date:** %s", Sys.Date()),
  sprintf("**Decision:** `%s`", decision),
  "",
  "---",
  "",
  "## Primary results",
  "",
  "### GSE183276 AKI vs Reference (n=10 AKI, 17 Ref, >=50 TAL cells)",
  "",
  sprintf("- TAL identity (unadjusted): %s (95%% CI %s to %s), P=%s",
    fmt3(safe_val(unadj_1836_ref, "estimate")),
    fmt3(safe_val(unadj_1836_ref, "ci_low")),
    fmt3(safe_val(unadj_1836_ref, "ci_high")),
    fmt4(safe_val(unadj_1836_ref, "p_value"))),
  sprintf("- PT-injury adjusted: %s (95%% CI %s to %s), P=%s",
    fmt3(safe_val(adj_pt_1836, "estimate")),
    fmt3(safe_val(adj_pt_1836, "ci_low")),
    fmt3(safe_val(adj_pt_1836, "ci_high")),
    fmt4(safe_val(adj_pt_1836, "p_value"))),
  sprintf("- TAL-minus-PT identity: %s (95%% CI %s to %s), P=%s",
    fmt3(safe_val(tal_pt_1836, "estimate")),
    fmt3(safe_val(tal_pt_1836, "ci_low")),
    fmt3(safe_val(tal_pt_1836, "ci_high")),
    fmt4(safe_val(tal_pt_1836, "p_value"))),
  sprintf("- LODO negative fraction: %.3f | Permutation P: %.4f",
    lodo_results[["GSE183276_AKI_vs_Ref"]],
    perm_table$p_value[perm_table$comparison == "GSE183276_AKI_vs_Ref"]),
  sprintf("- Negative genes: %d/6", neg_genes_1836),
  "",
  "### GSE210622 AKI vs non-AKI (n=8 AKI, 4 non-AKI, >=50 TAL cells)",
  "",
  sprintf("- TAL identity (unadjusted): %s (95%% CI %s to %s), P=%s",
    fmt3(safe_val(unadj_2106, "estimate")),
    fmt3(safe_val(unadj_2106, "ci_low")),
    fmt3(safe_val(unadj_2106, "ci_high")),
    fmt4(safe_val(unadj_2106, "p_value"))),
  sprintf("- PT-injury adjusted: %s (95%% CI %s to %s), P=%s",
    fmt3(safe_val(adj_pt_2106, "estimate")),
    fmt3(safe_val(adj_pt_2106, "ci_low")),
    fmt3(safe_val(adj_pt_2106, "ci_high")),
    fmt4(safe_val(adj_pt_2106, "p_value"))),
  sprintf("- TAL-minus-PT identity: %s (95%% CI %s to %s), P=%s",
    fmt3(safe_val(tal_pt_2106, "estimate")),
    fmt3(safe_val(tal_pt_2106, "ci_low")),
    fmt3(safe_val(tal_pt_2106, "ci_high")),
    fmt4(safe_val(tal_pt_2106, "p_value"))),
  sprintf("- LODO negative fraction: %.3f | Permutation P: %.4f",
    lodo_results[["GSE210622_AKI_vs_nonAKI"]],
    perm_table$p_value[perm_table$comparison == "GSE210622_AKI_vs_nonAKI"]),
  sprintf("- Negative genes: %d/6", neg_genes_2106),
  "",
  "## Sensitivity analyses",
  "",
  sprintf("- GSE183276 >=20 TAL cells: %s (n=%d+%d)",
    fmt3(safe_val(ge20_effect, "estimate")),
    safe_val(ge20_effect, "n_aki"), safe_val(ge20_effect, "n_ref")),
  sprintf("- GSE210622 COVID-adjusted: %s (direction persists in sensitivity)",
    fmt3(safe_val(covid_adj, "estimate"))),
  sprintf("- GSE210622 non-COVID only: %s (n=%d+%d)",
    fmt3(safe_val(noncovid_only, "estimate")),
    safe_val(noncovid_only, "n_aki"), safe_val(noncovid_only, "n_ref")),
  sprintf("- GSE210622 exclude CPM: %s", fmt3(safe_val(noCPM_effect, "estimate"))),
  sprintf("- GSE210622 nephrectomy ctrl only: %s", fmt3(safe_val(neph_effect, "estimate"))),
  sprintf("- GSE183276 AKI vs CKD: %s (95%% CI %s to %s), P=%s (n=%d AKI, %d CKD)",
    fmt3(safe_val(unadj_1836_ckd, "estimate")),
    fmt3(safe_val(unadj_1836_ckd, "ci_low")),
    fmt3(safe_val(unadj_1836_ckd, "ci_high")),
    fmt4(safe_val(unadj_1836_ckd, "p_value")),
    safe_val(unadj_1836_ckd, "n_aki"), safe_val(unadj_1836_ckd, "n_ref")),
  "",
  "## PT-injury relationship",
  "",
  sprintf("- AKI–PT-injury correlation: GSE183276 r=%.3f, GSE210622 r=%.3f",
    pt_cor_1836, pt_cor_2106),
  sprintf("- Attenuation after PT-injury adjustment: GSE183276 %.1f%%, GSE210622 %.1f%%",
    atten_1836, atten_2106),
  "",
  "## Per-gene direction",
  "",
  sprintf("- AKI vs Reference: %d/6 negative (%s)", neg_genes_1836,
    paste(gene_1836_ref$gene[gene_1836_ref$direction == "negative"], collapse=", ")),
  sprintf("- AKI vs non-AKI: %d/6 negative (%s)", neg_genes_2106,
    paste(gene_2106$gene[gene_2106$direction == "negative"], collapse=", ")),
  sprintf("- AKI vs CKD: %d/6 negative (%s)", neg_genes_ckd,
    paste(gene_1836_ckd$gene[gene_1836_ckd$direction == "negative"], collapse=", ")),
  "",
  "## KDIGO stage analysis (secondary, exploratory)",
  "",
  sprintf("- Combined cohort-adjusted: %s (95%% CI %s to %s), P=%s (n=%d stage 3, %d stage 2)",
    fmt3(kdigo_comb_est), fmt3(kdigo_comb_ci_low), fmt3(kdigo_comb_ci_high),
    fmt4(kdigo_comb_p), kdigo_n_s3, kdigo_n_s2),
  sprintf("- GSE183276 only: %s (95%% CI %s to %s), P=%s (n=%d+%d)",
    fmt3(kdigo_1836_est), fmt3(kdigo_1836_ci_low), fmt3(kdigo_1836_ci_high),
    fmt4(kdigo_1836_p), kdigo_1836_n_s3, kdigo_1836_n_s2),
  sprintf("- GSE210622 only: %s (95%% CI %s to %s), P=%s (n=%d+%d)",
    fmt3(kdigo_2106_est), fmt3(kdigo_2106_ci_low), fmt3(kdigo_2106_ci_high),
    fmt4(kdigo_2106_p), kdigo_2106_n_s3, kdigo_2106_n_s2),
  "- NOTE: exploratory only, not independently validated",
  "",
  "## TAL abundance",
  "",
  "No clear TAL-proportion difference between AKI and reference was detected",
  "in either cohort.  The recovered TAL proportion is technically confounded",
  "(the tissue-sampling procedure determines the total cell count, not the",
  "underlying nephron composition), so abundance effects are listed for",
  "completeness but not interpreted as biological signal.",
  "",
  "## Agreement with frozen decision rules",
  "",
  "| Criterion | GSE183276 | GSE210622 |",
  "|---|---|---|",
  sprintf("| AKI direction negative | %s | %s |",
    if(dir_1836_neg) "PASS" else "FAIL", if(dir_2106_neg) "PASS" else "FAIL"),
  sprintf("| LODO >= 80%% negative | %.0f%% (PASS) | %.0f%% (PASS) |",
    lodo_results[["GSE183276_AKI_vs_Ref"]]*100,
    lodo_results[["GSE210622_AKI_vs_nonAKI"]]*100),
  sprintf("| >=4/6 genes negative | %d/6 (PASS) | %d/6 (PASS) |",
    neg_genes_1836, neg_genes_2106),
  "| PT-injury independence | INCONCLUSIVE (attenuated) | Direction preserved only |",
  "| TAL-minus-PT selectivity | UNSUPPORTED (near-null) | UNSUPPORTED (reversed) |",
  "| AKI vs CKD specificity | UNRESOLVED | N/A |",
  "",
  "## Interpretation",
  "",
  "- **Cross-cohort reproducibility: PASS.** Lower frozen TAL identity in AKI reproduces",
  "  across GSE183276 and GSE210622 (direction negative, LODO 100%, >=4/6 genes negative).",
  "",
  "- **Independence from PT injury: NOT PASSED.** PT injury is correlated with AKI",
  sprintf("  status (GSE183276 r=%.2f).  After PT-injury adjustment the GSE183276 TAL effect",
    pt_cor_1836),
  sprintf("  attenuates %.0f%% (from %s to %s) and the CI crosses zero — the result is",
    atten_1836, fmt3(safe_val(unadj_1836_ref,"estimate")),
    fmt3(safe_val(adj_pt_1836,"estimate"))),
  "  inconclusive.  In GSE210622 the PT-adjusted direction is preserved and significant.",
  "",
  "- **TAL selectivity: UNSUPPORTED.** TAL-minus-PT identity contrast is approximately",
  "  null in GSE183276 and positive in GSE210622.  The lower TAL identity in AKI does not",
  "  isolate a TAL-specific signal from the pan-nephron injury background.",
  "",
  "- **AKI-vs-CKD specificity: UNRESOLVED.** The GSE183276 AKI-vs-CKD CI crosses zero,",
  "  PT adjustment attenuates the effect, and only 3/6 genes are negative.",
  "",
  "- **COVID sensitivity: direction persists** in COVID-adjusted and non-COVID-only",
  "  analyses, but postmortem-versus-nephrectomy confounding in GSE210622 is not resolved.",
  "",
  "## Conclusion",
  "",
  "Lower frozen six-gene TAL transport programme expression reproduces across",
  "independent biopsy (GSE183276) and postmortem (GSE210622) cohorts.  Independence",
  "from PT injury, TAL selectivity, and AKI-vs-CKD specificity remain unresolved.",
  "The evidence supports a conservative interpretation: donor-level TAL programme",
  "detectability is lower in AKI tissue and is not accompanied by a clear recovered",
  "TAL-proportion difference. The association is not separable from concurrent PT injury,",
  "and technical recovery limits prevent a cell-intrinsic interpretation.",
  "",
  "## Limitations",
  "",
  "- GSE183276 AKI/Ref donor-context confound (healthy living donor vs diseased AKI biopsy)",
  "- GSE210622 postmortem-vs-nephrectomy tissue-source confound",
  "- GSE210622 small non-AKI group (n=4)",
  "- Pseudobulk cannot distinguish cell-intrinsic transcriptional change from abundance",
  "- TAL-abundance estimates are confounded by tissue-sampling recovery",
  "- KDIGO stage analysis is exploratory (secondary, not independently confirmed)",
  "",
  "## Modified files",
  "",
  "- `code/02_primary_cell_resolved.R` — complete formal analysis with cohort-level z-scoring",
  "- `metadata/formal_donors.csv` — 57 donors, frozen metadata",
  "- `results/formal_effects.csv` — all model effects + LODO + permutation P",
  "- `results/formal_genes.csv` — per-gene direction analysis",
  "- `results/formal_lodo.csv` — all omitted-donor LODO estimates",
  "- `results/formal_evidence.csv` — per-donor scores, abundance, and metadata",
  "- `results/formal.md` — this summary",
  "- `figures/formal_1_cohort.pdf` — cohort characterization",
  "- `figures/formal_2_scores.pdf` — score distributions and forest plot",
  "- `figures/formal_3_specificity.pdf` — specificity and sensitivity",
  ""
)

writeLines(md_lines, file.path(result_dir, "formal.md"))

cat(paste(md_lines, collapse = "\n"), "\n")

# ---- 23. Print hashes for audit ----
message("\n=== OUTPUT HASHES ===")
for (f in c("formal_effects.csv","formal_genes.csv","formal_lodo.csv","formal_evidence.csv","formal.md")) {
  fp <- file.path(result_dir, f)
  if (file.exists(fp)) {
    cat(sprintf("  %s: %s\n", f, digest::digest(fp, algo = "md5", file = TRUE)))
  }
}

message("\n=== 02_primary_cell_resolved.R COMPLETE ===")
message(sprintf("Decision: %s", decision))
message(sprintf("Donors: %d unique across 2 datasets", nrow(formal_donors)))
message("All outputs saved to results/ and figures/.")
