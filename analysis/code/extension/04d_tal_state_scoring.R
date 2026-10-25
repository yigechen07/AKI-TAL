#!/usr/bin/env Rscript
# =============================================================================
# 04d_tal_state_scoring.R — Figure 4 (extension): TAL state context and score
#                            decomposition (formal analysis)
# =============================================================================
# Implements the frozen Figure 4 analysis (analysis_contract_figure4.md §6-§7)
# with definitions frozen BEFORE any diagnosis comparison:
#
#   STATE PROGRAMS (frozen; metadata/extension/tal_state_reference.csv)
#     primary   : reference_TAL 17, adaptive_TAL 281, degenerative_TAL 767
#                 (strict build_signatures rule; six frozen genes excluded)
#     sensitivity: adaptive_TAL 132-gene shared-adaptive-epithelial set
#                 (aPT&aTAL overlap), labelled NOT TAL-specific
#
#   CELLS
#     strict author-annotated TAL cells only (author_label == "TAL");
#     GSE183276 and GSE210622 are scored separately, never merged.
#
#   DONOR AGGREGATION (frozen)
#     PRIMARY        : donor mean UCell on >=50-cell donors
#     SENSITIVITY-A  : donor median UCell on the SAME >=50-cell donors
#     SENSITIVITY-B  : donor mean UCell on the expanded >=20-cell tier
#     <20 donors excluded (reported separately)
#
#   CONTRASTS (frozen)
#     GSE183276: 10 AKI vs 17 Ref (primary); 10 AKI vs 14 CKD (secondary)
#     GSE210622:  8 AKI vs  4 non-AKI
#
#   MODELS (donor is the unit; per cohort, per state, per aggregation)
#     score ~ AKI                                  (effect + 95% CI + raw P)
#     BH-FDR across the SIX primary cohort x state tests
#     LODO: leave-one-donor-out; sign concordance with full-data effect;
#           "negative fraction" ONLY when the full effect is negative
#
#   DECOMPOSITION (frozen before outcome testing)
#     PRIMARY : aggregate six-gene detection fraction (donor-level mean over
#               the six genes of the fraction of strict TAL cells with
#               counts > 0) ~ AKI
#     SECONDARY: per-gene detection fraction and per-gene positive-expression
#               (mean data-layer value among expressing cells; NO zero
#               imputation when no positive cells exist; valid donor n;
#               BH-FDR across the six genes per cohort)
#     six-gene TAL_identity (frozen, results/formal_evidence.csv,
#               tier == "primary_ge50" & comparison == "all_primary",
#               asserted one row per cohort-donor) ~ AKI, and
#               ~ AKI + state_score (change in AKI effect reported; never
#               called mediation or mechanism)
#     GSE210622 sensitivity: six-gene ~ AKI + COVID_status;
#               six-gene ~ AKI + sampling_source (full-rank here but AKI is
#               identified among postmortem samples only - flagged limited)
#
#   BOUNDARIES: no max-of-three state assignment, no state fractions, no
#   pseudotime, no cell-level P values, no cross-cohort pooling, no raw
#   cross-program magnitude comparison. All expression analyses are described
#   as within recovered strict TAL (never "within state").
#
# Usage:
#   cd code/extension && Rscript 04d_tal_state_scoring.R
#   Optional first CLI arg: project root
# Requires: Seurat, Matrix, data.table, UCell (project-local r-lib)
# =============================================================================

cli_args <- commandArgs(trailingOnly = TRUE)
# project root: explicit CLI arg, or robust auto-detection (works when run
# from the project root OR from code/extension)
if (length(cli_args) >= 1 && nzchar(cli_args[1])) {
  proj_root <- normalizePath(cli_args[1])
} else {
  cand <- c(getwd(), file.path(getwd(), "..", ".."))
  hit <- cand[vapply(cand, function(p) file.exists(file.path(p, "results", "extension")), logical(1))]
  if (length(hit) == 0) stop("Could not locate project root (no results/extension). Pass project root as first CLI arg.")
  proj_root <- normalizePath(hit[1])
}
# project-local R library must be on the search path BEFORE the package check
rl <- file.path(proj_root, "metadata", "extension", "r-lib")
.libPaths(c(rl, .libPaths()))

suppressMessages({
  for (pkg in c("Seurat", "Matrix", "data.table", "UCell")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package missing: ", pkg)
  }
  library(Seurat); library(Matrix); library(data.table); library(UCell)
})

REF_CSV  <- file.path(proj_root, "metadata", "extension", "tal_state_reference.csv")
FE_CSV   <- file.path(proj_root, "results", "formal_evidence.csv")
OBJ_183  <- file.path(proj_root, "results", "extension", "atlas_GSE183276_seurat.rds")
OBJ_622  <- file.path(proj_root, "results", "extension", "atlas_GSE210622_seurat.rds")
OUT_DIR  <- file.path(proj_root, "results", "extension")
OUT_DON  <- file.path(OUT_DIR, "tal_state_donor_scores.csv")
OUT_EFF  <- file.path(OUT_DIR, "tal_state_effects.csv")
OUT_LODO <- file.path(OUT_DIR, "tal_state_lodo.csv")
OUT_DEC  <- file.path(OUT_DIR, "tal_state_decomposition.csv")
OUT_PLOT <- file.path(OUT_DIR, "figure4_plotting_data.rds")

stopifnot(file.exists(REF_CSV), file.exists(FE_CSV), file.exists(OBJ_183), file.exists(OBJ_622))

FROZEN_SIX <- c("SLC12A1", "UMOD", "CLDN16", "KCNJ10", "PTH1R", "CLDN10")

cat("Software versions:\n")
cat("  R:", R.version.string, "\n")
cat("  Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("  UCell:", as.character(packageVersion("UCell")), "\n")
cat("  data.table:", as.character(packageVersion("data.table")), "\n")

# ---------------------------------------------------------------------------
# build_signatures(): identical rule to 04b (frozen)
# ---------------------------------------------------------------------------
build_signatures <- function(ref) {
  setDT(ref)
  ref[, excluded := exclusion_reason == "frozen_primary_endpoint_circularity"]
  pair_info <- ref[, .(
    any_non_none_flag = any(flag != "none"),
    any_frozen        = any(excluded),
    has_none_row      = any(flag == "none")
  ), by = .(state, gene)]
  primary <- pair_info[has_none_row == TRUE & any_non_none_flag == FALSE & any_frozen == FALSE,
                       .(state, gene)]
  sensitivity <- unique(ref[flag == "shared_adaptive_epithelial", .(state, gene)])
  list(primary = primary, sensitivity = sensitivity)
}

# ---------------------------------------------------------------------------
# Per-cell UCell scoring on strict author TAL cells (per cohort)
# ---------------------------------------------------------------------------
score_cohort <- function(obj_path, donor_col, cohort, primary, sensitivity) {
  obj <- readRDS(obj_path)
  md  <- as.data.table(obj@meta.data)
  strict_pos <- which(md$author_label == "TAL")
  strict_names <- colnames(obj)[strict_pos]
  cat(sprintf("[%s] strict TAL cells: %d\n", cohort, length(strict_names)))

  counts <- GetAssayData(obj, layer = "counts")
  sub_counts <- counts[, strict_names, drop = FALSE]
  ref_genes <- unique(c(primary$gene, sensitivity$gene))
  keep <- intersect(ref_genes, rownames(sub_counts))
  sub_counts <- sub_counts[keep, , drop = FALSE]

  sub_meta <- md[strict_pos]
  rownames(sub_meta) <- strict_names
  sub_obj <- CreateSeuratObject(counts = sub_counts, meta.data = sub_meta,
                                min.cells = 0, min.features = 0)

  sets <- list(
    reference_TAL    = intersect(primary[state == "reference_TAL", gene], rownames(sub_obj)),
    adaptive_TAL     = intersect(primary[state == "adaptive_TAL", gene], rownames(sub_obj)),
    degenerative_TAL = intersect(primary[state == "degenerative_TAL", gene], rownames(sub_obj)),
    adaptive_sensitivity = intersect(sensitivity[state == "adaptive_TAL", gene], rownames(sub_obj))
  )
  cat(sprintf("[%s] genes scored: %s\n", cohort,
              paste(sapply(sets, length), collapse = " / ")))
  stopifnot(all(sapply(sets, length) > 0))

  # SerialParam: the sandbox blocks BiocParallel socket ports
  bp <- BiocParallel::SerialParam()
  scored <- AddModuleScore_UCell(sub_obj, features = sets, BPPARAM = bp)
  sc <- as.data.table(scored[[]])
  nm <- grep("_UCell$", colnames(sc), value = TRUE)
  cat(sprintf("[%s] UCell columns: %s\n", cohort, paste(nm, collapse = ", ")))
  stopifnot(length(nm) == 4, !anyNA(sc[, ..nm]))

  # UMAP coordinates (existing global embedding; Panel A subsets it)
  umap <- if ("umap" %in% Reductions(obj)) as.data.table(Embeddings(obj[["umap"]]))[strict_pos] else NULL
  if (!is.null(umap)) setnames(umap, c("umap1", "umap2"))

  out <- data.table(
    cohort = cohort,
    cell = strict_names,
    donor_id = md[[donor_col]][strict_pos],
    diagnosis = md$diagnosis[strict_pos],
    author_label = md$author_label[strict_pos],
    broad_class = md$broad_class[strict_pos],
    umap1 = if (!is.null(umap)) umap$umap1 else NA_real_,
    umap2 = if (!is.null(umap)) umap$umap2 else NA_real_,
    score_reference_TAL = sc[[grep("reference_TAL_UCell$", nm, value = TRUE)]],
    score_adaptive_TAL  = sc[[grep("adaptive_TAL_UCell$", nm, value = TRUE)]],
    score_degenerative_TAL = sc[[grep("degenerative_TAL_UCell$", nm, value = TRUE)]],
    score_adaptive_sensitivity = sc[[grep("adaptive_sensitivity_UCell$", nm, value = TRUE)]]
  )
  rm(obj, counts, sub_counts, sub_obj, scored, sc); gc()
  out
}

# ---------------------------------------------------------------------------
# Donor aggregation (frozen): mean on >=50 (primary), median on >=50,
# mean on >=20 (sensitivity). <20 excluded.
# ---------------------------------------------------------------------------
aggregate_donor <- function(pc) {
  pc[, n_cells := .N, by = .(cohort, donor_id)]
  pc[, tier := fifelse(n_cells >= 50, "ge50_primary",
                       fifelse(n_cells >= 20, "ge20_sensitivity", "lt20_low-confidence"))]

  long <- melt(pc,
               id.vars = c("cohort", "donor_id", "n_cells", "tier"),
               measure.vars = c("score_reference_TAL", "score_adaptive_TAL",
                                "score_degenerative_TAL", "score_adaptive_sensitivity"),
               variable.name = "score_var", value.name = "score")
  long[, state := sub("^score_", "", score_var)]
  long[, set_type := fifelse(state == "adaptive_sensitivity", "sensitivity", "primary")]
  long[, state := fifelse(state == "adaptive_sensitivity", "adaptive_TAL", state)]
  long[, score_var := NULL]

  agg_mean_ge50   <- long[tier == "ge50_primary",
                          .(aggregation = "mean_ge50", score = mean(score)),
                          by = .(cohort, donor_id, state, set_type)]
  agg_median_ge50 <- long[tier == "ge50_primary",
                          .(aggregation = "median_ge50", score = median(score)),
                          by = .(cohort, donor_id, state, set_type)]
  agg_mean_ge20   <- long[tier %in% c("ge50_primary", "ge20_sensitivity"),
                          .(aggregation = "mean_ge20", score = mean(score)),
                          by = .(cohort, donor_id, state, set_type)]

  donor_scores <- rbindlist(list(agg_mean_ge50, agg_median_ge50, agg_mean_ge20))
  donor_scores <- merge(donor_scores,
                        unique(pc[, .(cohort, donor_id, n_cells, tier)]),
                        by = c("cohort", "donor_id"))
  donor_scores <- merge(donor_scores,
                        unique(pc[, .(cohort, donor_id, diagnosis)]),
                        by = c("cohort", "donor_id"))
  donor_scores[, n_strict_TAL := n_cells][, n_cells := NULL]
  donor_scores
}

# ---------------------------------------------------------------------------
# Model helpers (donor-level; effect, 95% CI, raw P)
# ---------------------------------------------------------------------------
fit_donor_lm <- function(dt, ref_level) {
  dt <- copy(dt)
  dt[, AKI := as.integer(diagnosis == "AKI")]
  dt[, group := fifelse(diagnosis == "AKI", "AKI", ref_level)]
  dt[, group := factor(group, levels = c(ref_level, "AKI"))]
  fit <- lm(score ~ group, data = dt)
  s <- summary(fit)$coefficients
  ci <- confint(fit)
  term_aki <- grep("AKI$", rownames(s), value = TRUE)
  if (length(term_aki) != 1) stop("AKI term not found or ambiguous in lm fit")
  data.table(
    estimate = unname(s[term_aki, "Estimate"]),
    std_error = unname(s[term_aki, "Std. Error"]),
    df_residual = fit$df.residual,
    statistic = unname(s[term_aki, "t value"]),
    p_value = unname(s[term_aki, "Pr(>|t|)"]),
    ci_low = unname(ci[term_aki, 1]),
    ci_high = unname(ci[term_aki, 2])
  )
}

# leave-one-donor-out stability (sign concordance with full-data effect)
lodo_summary <- function(dt, ref_level) {
  donors <- unique(dt$donor_id)
  eff <- numeric(length(donors))
  for (i in seq_along(donors)) {
    sub <- dt[donor_id != donors[i]]
    e <- tryCatch(fit_donor_lm(sub, ref_level)$estimate, error = function(e) NA_real_)
    eff[i] <- e
  }
  eff <- eff[is.finite(eff)]
  list(min = if (length(eff)) min(eff) else NA_real_,
       max = if (length(eff)) max(eff) else NA_real_,
       vals = eff, n = length(eff))
}

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
ref <- fread(REF_CSV, na.strings = c("NA", ""))
sig <- build_signatures(ref)
stopifnot(sig$primary[state == "reference_TAL", .N] == 17)
stopifnot(sig$primary[state == "adaptive_TAL", .N] == 281)
stopifnot(sig$primary[state == "degenerative_TAL", .N] == 767)
stopifnot(sig$sensitivity[, .N] == 132)
stopifnot(!any(sig$primary$gene %in% FROZEN_SIX))
stopifnot(!any(sig$sensitivity$gene %in% FROZEN_SIX))
cat("Signature sizes verified: 17/281/767 primary + 132 sensitivity; six-gene exclusion OK\n")

# ---- per-cell scoring (per cohort, separately) ----
pc183 <- score_cohort(OBJ_183, "patient", "GSE183276", sig$primary, sig$sensitivity)
pc622 <- score_cohort(OBJ_622, "donor",    "GSE210622", sig$primary, sig$sensitivity)

# ---- donor aggregation ----
donor_scores <- rbind(aggregate_donor(pc183), aggregate_donor(pc622))

# ---- donor annotations: six-gene frozen score (formal_evidence, primary rows) ----
fe <- fread(FE_CSV, na.strings = c("NA", ""))
fe_primary <- fe[tier == "primary_ge50" & comparison == "all_primary"]
fe_primary[, key := paste(dataset, donor_id, sep = "|")]
stopifnot(!anyDuplicated(fe_primary$key))
cat("formal_evidence primary_ge50/all_primary rows:", nrow(fe_primary),
    "| unique cohort-donor:", uniqueN(fe_primary$key), "\n")
sixgene_map <- unique(fe_primary[, .(dataset, donor_id,
                                     six_gene_TAL_identity = TAL_identity,
                                     covid_status, sampling_context)])
donor_scores <- merge(donor_scores, sixgene_map,
                      by.x = c("cohort", "donor_id"), by.y = c("dataset", "donor_id"),
                      all.x = TRUE)
donor_scores[cohort == "GSE183276" & is.na(covid_status),
             `:=`(covid_status = "not_reported", sampling_context = "Clinical biopsy")]

# ---- state effects models (per cohort x state x aggregation x contrast) ----
all_eff <- list(); all_lodo <- list(); primary_tests <- list()

add_state_row <- function(cohort, st, setype, agg, contr, note, sub) {
  n_aki <- sum(sub$diagnosis == "AKI"); n_ref <- sum(sub$diagnosis != "AKI")
  eff <- fit_donor_lm(sub, if (contr == "AKI_vs_CKD") "CKD" else if (contr == "AKI_vs_Reference") "Ref" else "non_AKI")
  lo  <- lodo_summary(sub, if (contr == "AKI_vs_CKD") "CKD" else if (contr == "AKI_vs_Reference") "Ref" else "non_AKI")
  row <- data.table(cohort = cohort, state = st, set_type = setype, aggregation = agg,
                    contrast = contr, n_aki = n_aki, n_ref = n_ref,
                    estimate = eff$estimate, std_error = eff$std_error,
                    df_residual = eff$df_residual, statistic = eff$statistic,
                    p_value = eff$p_value, ci_low = eff$ci_low, ci_high = eff$ci_high,
                    note = note)
  all_eff[[length(all_eff) + 1]] <<- row
  lrow <- data.table(cohort = cohort, state = st, set_type = setype, aggregation = agg,
                     contrast = contr, full_effect = eff$estimate,
                     lodo_min = lo$min, lodo_max = lo$max, n_donors = lo$n,
                     lodo_concordant_frac = if (is.finite(eff$estimate) && eff$estimate != 0)
                       mean(sign(lo$vals) == sign(eff$estimate)) else NA_real_,
                     lodo_negative_frac = if (is.finite(eff$estimate) && eff$estimate < 0)
                       mean(lo$vals < 0) else NA_real_,
                     note = note)
  all_lodo[[length(all_lodo) + 1]] <<- lrow
  if (agg == "mean_ge50" & setype == "primary" &
      ((cohort == "GSE183276" & contr == "AKI_vs_Reference") |
       (cohort == "GSE210622" & contr == "AKI_vs_nonAKI"))) {
    primary_tests[[length(primary_tests) + 1]] <<- row
  }
  invisible(NULL)
}

for (coh in c("GSE183276", "GSE210622")) {
  don <- donor_scores[cohort == coh]
  for (st in c("reference_TAL", "adaptive_TAL", "degenerative_TAL")) {
    for (setype in c("primary", "sensitivity")) {
      for (agg in c("mean_ge50", "median_ge50", "mean_ge20")) {
        for (contr in if (coh == "GSE183276") c("AKI_vs_Reference", "AKI_vs_CKD") else "AKI_vs_nonAKI") {
          sub <- don[state == st & set_type == setype & aggregation == agg & !is.na(score)]
          if (contr == "AKI_vs_Reference") sub <- sub[diagnosis %in% c("AKI", "Ref")]
          if (contr == "AKI_vs_CKD")       sub <- sub[diagnosis %in% c("AKI", "CKD")]
          if (contr == "AKI_vs_nonAKI")    sub <- sub[diagnosis %in% c("AKI", "non_AKI")]
          if (nrow(sub) < 4) next
          note <- if (coh == "GSE183276")
            ifelse(contr == "AKI_vs_Reference", "primary_AKI_vs_Ref", "secondary_AKI_vs_CKD")
          else "primary_AKI_vs_nonAKI"
          add_state_row(coh, st, setype, agg, contr, note, sub)
        }
      }
    }
  }
}
eff_all  <- rbindlist(all_eff)
lodo_all <- rbindlist(all_lodo)

# BH-FDR across the six primary cohort x state tests (mean_ge50, primary sets)
pt <- rbindlist(primary_tests)
stopifnot(nrow(pt) == 6)
pt[, fdr_bh := p.adjust(p_value, method = "BH")]
eff_all[, fdr_bh := NA_real_]
eff_all[pt, on = c("cohort", "state", "set_type", "aggregation", "contrast"),
        fdr_bh := i.fdr_bh]

# ---- six-gene decomposition ----
dec <- list()
dec_add <- function(row) dec[[length(dec) + 1]] <<- row

reflev <- function(cohort) ifelse(cohort == "GSE183276", "Ref", "non_AKI")

# unadjusted six-gene model per cohort (primary donors)
for (coh in c("GSE183276", "GSE210622")) {
  d6 <- unique(donor_scores[cohort == coh & !is.na(six_gene_TAL_identity) &
                            aggregation == "mean_ge50" & state == "reference_TAL",
                            .(donor_id, diagnosis, six_gene_TAL_identity)])
  if (coh == "GSE183276") d6 <- d6[diagnosis %in% c("AKI", "Ref")]
  if (coh == "GSE210622") d6 <- d6[diagnosis %in% c("AKI", "non_AKI")]
  d6[, score := six_gene_TAL_identity]
  e <- fit_donor_lm(d6, reflev(coh))
  dec_add(data.table(cohort = coh, model_type = "six_gene_unadjusted",
                     state = NA_character_, set_type = NA_character_, aggregation = "primary_ge50",
                     gene = NA_character_, term = "AKI",
                     n_aki = sum(d6$diagnosis == "AKI"), n_ref = sum(d6$diagnosis != "AKI"),
                     e, fdr_bh = NA_real_,
                     note = "frozen six-gene TAL_identity ~ AKI, primary donors"))
}

# adjusted six-gene models (per state program x aggregation), on donors with
# six-gene scores (primary_ge50 donors)
for (coh in c("GSE183276", "GSE210622")) {
  for (st in c("reference_TAL", "adaptive_TAL", "degenerative_TAL")) {
    for (setype in c("primary", "sensitivity")) {
      for (agg in c("mean_ge50", "median_ge50", "mean_ge20")) {
        sub <- donor_scores[cohort == coh & state == st & set_type == setype &
                            aggregation == agg & !is.na(six_gene_TAL_identity) & !is.na(score)]
        sub <- unique(sub[, .(donor_id, diagnosis, six_gene_TAL_identity, score)])
        if (coh == "GSE183276") sub <- sub[diagnosis %in% c("AKI", "Ref")]
        if (coh == "GSE210622") sub <- sub[diagnosis %in% c("AKI", "non_AKI")]
        if (nrow(sub) < 5) next
        sub[, AKI := as.integer(diagnosis == "AKI")]
        sub[, group := factor(ifelse(diagnosis == "AKI", "AKI", reflev(coh)),
                              levels = c(reflev(coh), "AKI"))]
        fit <- lm(six_gene_TAL_identity ~ group + score, data = sub)
        s <- summary(fit)$coefficients; ci <- confint(fit)
        term_aki <- grep("AKI$", rownames(s), value = TRUE)
        if (length(term_aki) != 1) stop("AKI term not found in adjusted model")
        dec_add(data.table(cohort = coh, model_type = "six_gene_adjusted_state",
                           state = st, set_type = setype, aggregation = agg,
                           gene = NA_character_, term = "AKI",
                           n_aki = sum(sub$diagnosis == "AKI"), n_ref = sum(sub$diagnosis != "AKI"),
                           estimate = unname(s[term_aki, "Estimate"]),
                           std_error = unname(s[term_aki, "Std. Error"]),
                           df_residual = fit$df.residual,
                           statistic = unname(s[term_aki, "t value"]),
                           p_value = unname(s[term_aki, "Pr(>|t|)"]),
                           ci_low = unname(ci[term_aki, 1]), ci_high = unname(ci[term_aki, 2]),
                           fdr_bh = NA_real_,
                           note = "six-gene ~ AKI + state_score; AKI term (not mediation)"))
        dec_add(data.table(cohort = coh, model_type = "six_gene_adjusted_state",
                           state = st, set_type = setype, aggregation = agg,
                           gene = NA_character_, term = "state_score",
                           n_aki = sum(sub$diagnosis == "AKI"), n_ref = sum(sub$diagnosis != "AKI"),
                           estimate = unname(s["score", "Estimate"]),
                           std_error = unname(s["score", "Std. Error"]),
                           df_residual = fit$df.residual,
                           statistic = unname(s["score", "t value"]),
                           p_value = unname(s["score", "Pr(>|t|)"]),
                           ci_low = unname(ci["score", 1]), ci_high = unname(ci["score", 2]),
                           fdr_bh = NA_real_,
                           note = "six-gene ~ AKI + state_score; state_score term"))
      }
    }
  }
}

# GSE210622 COVID and sampling-source sensitivities
{
  coh <- "GSE210622"
  d6 <- unique(donor_scores[cohort == coh & !is.na(six_gene_TAL_identity) &
                            aggregation == "mean_ge50" & state == "reference_TAL",
                            .(donor_id, diagnosis, six_gene_TAL_identity,
                              covid_status, sampling_context)])
  d6 <- d6[diagnosis %in% c("AKI", "non_AKI")]
  d6[, AKI := as.integer(diagnosis == "AKI")]
  # COVID-adjusted (mirror of frozen 02 adjusted_COVID)
  d6[, covid_factor := factor(covid_status, levels = c("not_reported", "COVID", "Non-COVID"))]
  f1 <- lm(six_gene_TAL_identity ~ AKI + covid_factor, data = d6)
  s1 <- summary(f1)$coefficients; ci1 <- confint(f1)
  covid_mm <- model.matrix(f1)
  covid_rank <- qr(covid_mm)$rank
  stopifnot(covid_rank < ncol(covid_mm))
  dec_add(data.table(cohort = coh, model_type = "six_gene_covid_subgroup_contrast",
                     state = NA_character_, set_type = NA_character_, aggregation = "primary_ge50",
                     gene = NA_character_, term = "AKI",
                     n_aki = sum(d6$diagnosis == "AKI"), n_ref = sum(d6$diagnosis != "AKI"),
                     estimate = unname(s1["AKI", "Estimate"]),
                     std_error = unname(s1["AKI", "Std. Error"]),
                     df_residual = f1$df.residual, statistic = unname(s1["AKI", "t value"]),
                     p_value = unname(s1["AKI", "Pr(>|t|)"]),
                     ci_low = unname(ci1["AKI", 1]), ci_high = unname(ci1["AKI", 2]),
                     fdr_bh = NA_real_,
                     note = sprintf(paste0("rank-deficient COVID model (%d/%d rank); AKI coefficient is ",
                                           "Non-COVID AKI vs non-AKI not_reported subgroup contrast, ",
                                           "NOT an overall COVID-adjusted AKI effect"),
                                    covid_rank, ncol(covid_mm))))
  # sampling-source adjusted (full-rank; AKI identified among postmortem only)
  d6[, source_factor := factor(sampling_context, levels = c("nephrectomy", "postmortem"))]
  f2 <- lm(six_gene_TAL_identity ~ AKI + source_factor, data = d6)
  s2 <- summary(f2)$coefficients; ci2 <- confint(f2)
  dec_add(data.table(cohort = coh, model_type = "six_gene_source_adjusted",
                     state = NA_character_, set_type = NA_character_, aggregation = "primary_ge50",
                     gene = NA_character_, term = "AKI",
                     n_aki = sum(d6$diagnosis == "AKI"), n_ref = sum(d6$diagnosis != "AKI"),
                     estimate = unname(s2["AKI", "Estimate"]),
                     std_error = unname(s2["AKI", "Std. Error"]),
                     df_residual = f2$df.residual, statistic = unname(s2["AKI", "t value"]),
                     p_value = unname(s2["AKI", "Pr(>|t|)"]),
                     ci_low = unname(ci2["AKI", 1]), ci_high = unname(ci2["AKI", 2]),
                     fdr_bh = NA_real_,
                     note = "six-gene ~ AKI + sampling_source (full-rank; AKI identified among postmortem only; single postmortem control CPM -> limited)"))
}

# ---- six-gene detection / positive-expression decomposition ----
det_expression <- function(cohort) {
  obj <- readRDS(if (cohort == "GSE183276") OBJ_183 else OBJ_622)
  md <- as.data.table(obj@meta.data)
  strict_pos <- which(md$author_label == "TAL")
  strict_names <- colnames(obj)[strict_pos]
  X <- GetAssayData(obj, layer = "counts")[FROZEN_SIX, strict_names, drop = FALSE]
  D <- GetAssayData(obj, layer = "data")[FROZEN_SIX, strict_names, drop = FALSE]
  donor_of_cell <- md[[if (cohort == "GSE183276") "patient" else "donor"]][strict_pos]
  rm(obj, md); gc()

  out_det <- list(); out_expr <- list()
  for (d in unique(donor_of_cell)) {
    idx <- which(donor_of_cell == d)
    det <- rowMeans(X[, idx, drop = FALSE] > 0)
    for (g in FROZEN_SIX) {
      pos <- X[g, idx] > 0
      n_pos <- sum(pos)
      expr <- if (n_pos > 0) mean(D[g, idx][pos]) else NA_real_   # no imputation
      out_det[[length(out_det) + 1]] <- data.table(cohort = cohort, donor_id = d,
                                                   gene = g, n_strict_TAL = length(idx),
                                                   detection_fraction = unname(det[[g]]))
      out_expr[[length(out_expr) + 1]] <- data.table(cohort = cohort, donor_id = d,
                                                     gene = g, n_strict_TAL = length(idx),
                                                     n_positive = n_pos,
                                                     positive_expression = expr)
    }
  }
  rm(X, D); gc()
  list(det = rbindlist(out_det), expr = rbindlist(out_expr))
}

det_results <- lapply(c("GSE183276", "GSE210622"), det_expression)
det_dt  <- rbindlist(lapply(det_results, `[[`, "det"))
expr_dt <- rbindlist(lapply(det_results, `[[`, "expr"))
rm(det_results); gc()

don_diag <- unique(donor_scores[, .(cohort, donor_id, diagnosis)])
agg_det <- det_dt[, .(agg_detection = mean(detection_fraction)), by = .(cohort, donor_id)]
agg_det <- merge(agg_det, unique(det_dt[, .(cohort, donor_id, n_strict_TAL)]), by = c("cohort", "donor_id"))
agg_det <- merge(agg_det, don_diag, by = c("cohort", "donor_id"))

for (coh in c("GSE183276", "GSE210622")) {
  sub <- agg_det[cohort == coh & n_strict_TAL >= 50]
  if (coh == "GSE183276") sub <- sub[diagnosis %in% c("AKI", "Ref")]
  if (coh == "GSE210622") sub <- sub[diagnosis %in% c("AKI", "non_AKI")]
  sub[, score := agg_detection]
  e <- fit_donor_lm(sub, reflev(coh))
  dec_add(data.table(cohort = coh, model_type = "agg_detection_primary",
                     state = NA_character_, set_type = NA_character_, aggregation = "primary_ge50",
                     gene = NA_character_, term = "AKI",
                     n_aki = sum(sub$diagnosis == "AKI"), n_ref = sum(sub$diagnosis != "AKI"),
                     e, fdr_bh = NA_real_,
                     note = "aggregate six-gene detection fraction (PRIMARY decomposition endpoint)"))
}

# per-gene detection and positive-expression (secondary; BH-FDR over six genes)
dec_extra <- list()
for (coh in c("GSE183276", "GSE210622")) {
  dd <- merge(det_dt[cohort == coh], don_diag, by = c("cohort", "donor_id"))
  dd <- dd[n_strict_TAL >= 50]
  if (coh == "GSE183276") dd <- dd[diagnosis %in% c("AKI", "Ref")]
  if (coh == "GSE210622") dd <- dd[diagnosis %in% c("AKI", "non_AKI")]
  for (g in FROZEN_SIX) {
    sub <- dd[gene == g]; sub[, score := detection_fraction]
    e <- fit_donor_lm(sub, reflev(coh))
    dec_extra[[length(dec_extra) + 1]] <- data.table(
      cohort = coh, model_type = "pergene_detection",
      state = NA_character_, set_type = NA_character_, aggregation = "primary_ge50",
      gene = g, term = "AKI",
      n_aki = sum(sub$diagnosis == "AKI"), n_ref = sum(sub$diagnosis != "AKI"),
      e, fdr_bh = NA_real_, note = "per-gene detection fraction (secondary)")
  }
  ee <- merge(expr_dt[cohort == coh], don_diag, by = c("cohort", "donor_id"))
  ee <- ee[n_strict_TAL >= 50]
  if (coh == "GSE183276") ee <- ee[diagnosis %in% c("AKI", "Ref")]
  if (coh == "GSE210622") ee <- ee[diagnosis %in% c("AKI", "non_AKI")]
  for (g in FROZEN_SIX) {
    sub <- ee[gene == g & !is.na(positive_expression)]
    if (nrow(sub) < 4) {
      dec_extra[[length(dec_extra) + 1]] <- data.table(
        cohort = coh, model_type = "pergene_positive_expression",
        state = NA_character_, set_type = NA_character_, aggregation = "primary_ge50",
        gene = g, term = "AKI",
        n_aki = sum(sub$diagnosis == "AKI"), n_ref = sum(sub$diagnosis != "AKI"),
        estimate = NA_real_, std_error = NA_real_, df_residual = NA_integer_,
        statistic = NA_real_, p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        fdr_bh = NA_real_, note = "not fitted (valid donor n < 4)")
      next
    }
    sub[, score := positive_expression]
    e <- fit_donor_lm(sub, reflev(coh))
    dec_extra[[length(dec_extra) + 1]] <- data.table(
      cohort = coh, model_type = "pergene_positive_expression",
      state = NA_character_, set_type = NA_character_, aggregation = "primary_ge50",
      gene = g, term = "AKI",
      n_aki = sum(sub$diagnosis == "AKI"), n_ref = sum(sub$diagnosis != "AKI"),
      e, fdr_bh = NA_real_,
      note = "per-gene positive expression among expressing cells (secondary)")
  }
}
dec_dt <- rbindlist(c(dec, dec_extra))
dec_dt[model_type == "pergene_detection", fdr_bh := p.adjust(p_value, method = "BH"), by = cohort]
dec_dt[model_type == "pergene_positive_expression", fdr_bh := p.adjust(p_value, method = "BH"), by = cohort]

# Collinearity diagnostics for every six-gene + state-score model. These are
# essential for interpreting coefficient attenuation; a changed AKI point
# estimate is not evidence of mediation or of a biological state transition.
dec_dt[, `:=`(cor_diagnosis_state = NA_real_, cor_sixgene_state = NA_real_,
              vif_diagnosis = NA_real_)]
for (coh in c("GSE183276", "GSE210622")) {
  for (st in c("reference_TAL", "adaptive_TAL", "degenerative_TAL")) {
    for (setype in c("primary", "sensitivity")) {
      for (agg in c("mean_ge50", "median_ge50", "mean_ge20")) {
        sub <- donor_scores[cohort == coh & state == st & set_type == setype &
                            aggregation == agg & !is.na(six_gene_TAL_identity) & !is.na(score)]
        sub <- unique(sub[, .(donor_id, diagnosis, six_gene_TAL_identity, score)])
        if (coh == "GSE183276") sub <- sub[diagnosis %in% c("AKI", "Ref")]
        else sub <- sub[diagnosis %in% c("AKI", "non_AKI")]
        if (nrow(sub) < 5L) next
        aki <- as.integer(sub$diagnosis == "AKI")
        r_dx <- suppressWarnings(cor(aki, sub$score))
        r_y <- suppressWarnings(cor(sub$six_gene_TAL_identity, sub$score))
        vif <- if (is.finite(r_dx) && abs(r_dx) < 1) 1 / (1 - r_dx^2) else Inf
        dec_dt[cohort == coh & model_type == "six_gene_adjusted_state" &
                 state == st & set_type == setype & aggregation == agg,
               `:=`(cor_diagnosis_state = r_dx, cor_sixgene_state = r_y,
                    vif_diagnosis = vif)]
      }
    }
  }
}

# ---- write outputs ----
setcolorder(donor_scores, c("cohort", "donor_id", "diagnosis", "n_strict_TAL", "tier",
                            "state", "set_type", "aggregation", "score",
                            "six_gene_TAL_identity", "covid_status", "sampling_context"))
setcolorder(eff_all, c("cohort", "state", "set_type", "aggregation", "contrast",
                       "n_aki", "n_ref", "estimate", "std_error", "df_residual",
                       "statistic", "p_value", "ci_low", "ci_high", "fdr_bh", "note"))
setcolorder(lodo_all, c("cohort", "state", "set_type", "aggregation", "contrast",
                        "full_effect", "lodo_min", "lodo_max", "lodo_concordant_frac",
                        "lodo_negative_frac", "n_donors", "note"))
setcolorder(dec_dt, c("cohort", "model_type", "state", "set_type", "aggregation", "gene",
                      "term", "n_aki", "n_ref", "estimate", "std_error", "df_residual",
                      "statistic", "p_value", "ci_low", "ci_high", "fdr_bh",
                      "cor_diagnosis_state", "cor_sixgene_state", "vif_diagnosis", "note"))

fwrite(donor_scores, OUT_DON)
fwrite(eff_all, OUT_EFF)
fwrite(lodo_all, OUT_LODO)
fwrite(dec_dt, OUT_DEC)

# reproducible cell-coordinate + continuous-score plotting object (Panel A)
plot_data <- rbind(pc183, pc622)
saveRDS(plot_data, OUT_PLOT)

cat("\n===== OUTPUT SUMMARY =====\n")
cat("donor_scores rows:", nrow(donor_scores), "\n")
cat("effects rows:", nrow(eff_all), "\n")
cat("lodo rows:", nrow(lodo_all), "\n")
cat("decomposition rows:", nrow(dec_dt), "\n")
cat("plotting data rows:", nrow(plot_data), " (", uniqueN(plot_data$cohort), " cohorts)\n")
cat("primary FDR tests:", nrow(pt), "\n")
cat("DONE\n")
