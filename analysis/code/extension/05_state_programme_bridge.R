#!/usr/bin/env Rscript
# =============================================================================
# 05_state_programme_bridge.R — Figure 4f (extension): state-programme bridge
#                                between Figure 3 and Figure 4 (formal analysis)
# =============================================================================
# INTERNAL EXPLANATORY BRIDGE between the TAL-derived AKI consensus module
# (Figure 3, results/extension/consensus_module_definition.csv, 1,027 genes)
# and the continuous TAL state programmes (Figure 4,
# results/extension/tal_state_donor_scores.csv + metadata/extension/
# tal_state_reference.csv) at the DONOR level within the SAME two cell-resolved
# cohorts (GSE183276, GSE210622). NOT an independent validation.
#
# Frozen inputs:
#   ../AKI_DCT_analysis/results/stage3b_pseudobulk.rds  (strict TAL
#       donor-pseudobulk counts + per-segment TMM factors, frozen stage3b)
#   ../AKI_DCT_analysis/results/stage3b_qc.csv          (donor QC tiers)
#   results/extension/consensus_module_definition.csv   (frozen 1,027-gene
#       consensus module, direction = frozen logFC sign)
#   results/extension/tal_state_donor_scores.csv        (frozen Figure 4 donor
#       scores: primary mean_ge50, sensitivity mean_ge20)
#   metadata/extension/tal_state_reference.csv          (frozen state reference)
#
# Method (per cohort, cohorts never merged):
#   1. Strict TAL pseudobulk only. Donors: primary = AKI-Reference/non-AKI
#      donors with >=50 strict TAL cells; sensitivity = same diagnoses with
#      >=20 cells.  Two cohorts analysed independently.
#   2. TMM logCPM (reusing the frozen per-segment TMM factors from
#      stage3b_pseudobulk.rds) -> gene-wise z-score across the analysis donors
#      of that tier (per cohort).
#   3. Consensus signed score = mean(z of module-up genes) - mean(z of
#      module-down genes) using the FROZEN module directions.
#      ANTI-CIRCULARITY: (a) the six frozen endpoint genes are ALWAYS excluded
#      from the consensus score; (b) for each state test, all consensus genes
#      overlapping that state's primary signature are additionally excluded;
#      the remaining gene count is reported.
#   4. Primary endpoint: donor-level association of consensus score with the
#      reference_TAL score, adjusted for diagnosis. Adaptive_TAL and
#      degenerative_TAL are secondary. Reported per cohort x state:
#      raw Spearman rho + P; diagnosis-adjusted linear coefficient beta with
#      95% CI and P; partial correlation (consensus | state given diagnosis);
#      BH-FDR across the SIX cohort x state tests; LODO direction stability.
#   5. Collinearity: cor(diagnosis, state), diagnosis VIF, condition number
#      of the [1, AKI, state] design matrix. A highly collinear cohort/state
#      (e.g. GSE210622 reference_TAL) keeps its results but is not
#      interpretable as an independent effect.
#   6. Figure 4f: donor-level scatter and diagnosis-adjusted residual scatter
#      (both cohorts), state-effect forest, LODO stability. No cell-level P.
#
# Claim boundary: internal explanatory association between two extension
# analyses in the same cohorts. NO mediation, state-transition, mechanism,
# causality, TAL-specificity, PT-independence or biological-abundance claims.
# All null, opposite-direction and non-identifiable results are retained.
#
# Outputs (extension dirs only):
#   results/extension/state_programme_bridge_{donor_scores,alignment,
#       overlap_audit,effects,lodo,collinearity}.csv
#   results/extension/state_programme_bridge_report.md
#   figures/extension/figure4f_state_programme_bridge.{png,pdf}
# =============================================================================

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1 && nzchar(cli_args[1])) {
  proj_root <- normalizePath(cli_args[1])
} else {
  cand <- c(getwd(), file.path(getwd(), "..", ".."))
  hit <- cand[vapply(cand, function(p) file.exists(file.path(p, "results", "extension")), logical(1))]
  if (length(hit) == 0) stop("Could not locate project root (no results/extension). Pass project root as first CLI arg.")
  proj_root <- normalizePath(hit[1])
}
setwd(proj_root)
rl <- file.path(proj_root, "metadata", "extension", "r-lib")
.libPaths(c(rl, .libPaths()))

# UTF-8 locale: without this, R writes non-ASCII characters (em/en dashes)
# as "<U+2014>" escapes because the default LC_CTYPE is "C".
invisible(tryCatch(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"),
                   error = function(e) warning("could not set UTF-8 locale: ", conditionMessage(e))))

suppressMessages({
  for (pkg in c("edgeR", "data.table", "ggplot2", "gridExtra", "grid", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package missing: ", pkg)
  }
  library(edgeR); library(data.table); library(ggplot2)
  library(gridExtra); library(grid); library(scales)
})

RES   <- file.path(proj_root, "results", "extension")
FIG   <- file.path(proj_root, "figures", "extension")
DCT   <- file.path(proj_root, "..", "AKI_DCT_analysis", "results")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

FROZEN_SIX <- c("SLC12A1", "UMOD", "CLDN16", "KCNJ10", "PTH1R", "CLDN10")
SET_SEED <- 20260809L
set.seed(SET_SEED)

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
fmt3 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))
fmt4 <- function(x) ifelse(is.na(x), "NA", sprintf("%.4f", x))
fmt2 <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))

cat("Software versions:\n")
cat("  R:", R.version.string, "\n")
cat("  edgeR:", as.character(packageVersion("edgeR")), "\n")
cat("  data.table:", as.character(packageVersion("data.table")), "\n")
cat("  ggplot2:", as.character(packageVersion("ggplot2")), "\n")

# =============================================================================
# 1. Load frozen inputs
# =============================================================================
msg("--- 1. Load frozen inputs ---")
stopifnot(file.exists(file.path(DCT, "stage3b_pseudobulk.rds")),
          file.exists(file.path(DCT, "stage3b_qc.csv")),
          file.exists(file.path(RES, "consensus_module_definition.csv")),
          file.exists(file.path(RES, "tal_state_donor_scores.csv")),
          file.exists(file.path(proj_root, "metadata", "extension", "tal_state_reference.csv")))

pb  <- readRDS(file.path(DCT, "stage3b_pseudobulk.rds"))
qc  <- fread(file.path(DCT, "stage3b_qc.csv"))
mod <- fread(file.path(RES, "consensus_module_definition.csv"))
don <- fread(file.path(RES, "tal_state_donor_scores.csv"))
ref <- fread(file.path(proj_root, "metadata", "extension", "tal_state_reference.csv"),
             na.strings = c("NA", ""))

stopifnot(nrow(mod) == 1027L)
stopifnot(all(mod$direction %in% c("up", "down")))
stopifnot(sum(mod$direction == "up") == 723L, sum(mod$direction == "down") == 304L)

# Rebuild state primary signatures with the frozen build_signatures rule
# (identical to 04b/04d).
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
sig <- build_signatures(ref)
sig_prim <- sig$primary
stopifnot(sig_prim[state == "reference_TAL", .N] == 17L)
stopifnot(sig_prim[state == "adaptive_TAL", .N] == 281L)
stopifnot(sig_prim[state == "degenerative_TAL", .N] == 767L)
stopifnot(!any(sig_prim$gene %in% FROZEN_SIX))

STATES <- c("reference_TAL", "adaptive_TAL", "degenerative_TAL")

# =============================================================================
# 2. Donor sets (strict TAL pseudobulk; per cohort; never merged)
# =============================================================================
msg("--- 2. Donor sets ---")
tal_qc <- qc[segment == "TAL"]
donor_tier <- function(coh) {
  d <- tal_qc[dataset == coh]
  d[, tier := fifelse(n_cells >= 50, "ge50_primary",
                      fifelse(n_cells >= 20, "ge20_sensitivity", "lt20_excluded"))]
  d[, .(donor_id = donor, diagnosis, n_strict_TAL = n_cells, tier)]
}
tier183 <- donor_tier("GSE183276")
tier622 <- donor_tier("GSE210622")
tiers <- rbind(tier183[, cohort := "GSE183276"], tier622[, cohort := "GSE210622"])

# State donor scores (frozen Figure 4) for primary (mean_ge50) and
# sensitivity (mean_ge20) aggregations; primary sets only.
don_scores <- don[set_type == "primary" & aggregation %in% c("mean_ge50", "mean_ge20")]
don_scores <- don_scores[, .(cohort, donor_id, diagnosis, state, aggregation, score,
                             six_gene_TAL_identity, covid_status, sampling_context)]
don_wide <- dcast(don_scores, cohort + donor_id + diagnosis + aggregation +
                    six_gene_TAL_identity + covid_status + sampling_context ~ state,
                  value.var = "score")

primary_donors <- function(coh) {
  comp <- if (coh == "GSE183276") "Ref" else "non_AKI"
  d <- tiers[cohort == coh & tier == "ge50_primary" & diagnosis %in% c("AKI", comp)]
  d
}
sensitivity_donors <- function(coh) {
  comp <- if (coh == "GSE183276") "Ref" else "non_AKI"
  d <- tiers[cohort == coh & tier %in% c("ge50_primary", "ge20_sensitivity") &
               diagnosis %in% c("AKI", comp)]
  d
}
prim_don <- rbind(primary_donors("GSE183276")[, tier := "primary_ge50"],
                  primary_donors("GSE210622")[, tier := "primary_ge50"])
sens_don <- rbind(sensitivity_donors("GSE183276")[, tier := "sensitivity_ge20"],
                  sensitivity_donors("GSE210622")[, tier := "sensitivity_ge20"])
msg(sprintf("  GSE183276 primary: %d AKI + %d Ref = %d donors; sensitivity adds %d",
            sum(prim_don[cohort == "GSE183276"]$diagnosis == "AKI"),
            sum(prim_don[cohort == "GSE183276"]$diagnosis != "AKI"),
            nrow(prim_don[cohort == "GSE183276"]),
            nrow(sens_don[cohort == "GSE183276"]) - nrow(prim_don[cohort == "GSE183276"])))
msg(sprintf("  GSE210622 primary: %d AKI + %d non-AKI = %d donors; sensitivity adds %d",
            sum(prim_don[cohort == "GSE210622"]$diagnosis == "AKI"),
            sum(prim_don[cohort == "GSE210622"]$diagnosis != "AKI"),
            nrow(prim_don[cohort == "GSE210622"]),
            nrow(sens_don[cohort == "GSE210622"]) - nrow(prim_don[cohort == "GSE210622"])))

# =============================================================================
# 3. Consensus signed scores from strict TAL pseudobulk
# =============================================================================
msg("--- 3. Consensus signed scores (TMM logCPM, gene-wise z) ---")

# module genes excluding the six frozen genes ALWAYS
mod_no_frozen <- mod[!gene %in% FROZEN_SIX]
stopifnot(nrow(mod_no_frozen) == nrow(mod) - 3L)   # SLC12A1, CLDN16, KCNJ10 in module
up_all   <- mod_no_frozen[direction == "up", gene]
down_all <- mod_no_frozen[direction == "down", gene]

# per-state exclusions (anti-circularity): consensus minus frozen minus
# consensus genes overlapping the state's primary signature
overlap_audit <- data.table(state = character(), n_sig_genes = integer(),
                            n_overlap = integer(), n_overlap_up = integer(),
                            n_overlap_down = integer(),
                            n_remaining_no_frozen = integer())
for (st in STATES) {
  sig_genes <- sig_prim[state == st, gene]
  ov <- intersect(sig_genes, mod_no_frozen$gene)
  overlap_audit <- rbind(overlap_audit, data.table(
    state = st, n_sig_genes = length(sig_genes), n_overlap = length(ov),
    n_overlap_up = sum(mod_no_frozen[gene %in% ov]$direction == "up"),
    n_overlap_down = sum(mod_no_frozen[gene %in% ov]$direction == "down"),
    n_remaining_no_frozen = nrow(mod_no_frozen) - length(ov)))
}

score_donors <- function(coh, donor_dt, zspace_donors) {
  m <- pb$pseudobulk[[coh]][["TAL"]]
  md <- as.data.table(pb$metadata[[coh]])
  tmm <- pb$tmm_factors[[coh]][["TAL"]]
  cols <- intersect(colnames(m), zspace_donors)
  stopifnot(setequal(cols, zspace_donors))
  dge <- DGEList(counts = m[, cols, drop = FALSE],
                 lib.size = Matrix::colSums(m[, cols, drop = FALSE]))
  dge$samples$norm.factors <- tmm[cols]
  logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  # gene-wise z across the z-space donors of this cohort/tier
  mu <- rowMeans(logcpm); s <- apply(logcpm, 1, sd)
  s[!is.finite(s) | s < 1e-8] <- NA_real_
  z <- sweep(sweep(logcpm, 1, mu, "-"), 1, s, "/")
  z <- z[is.finite(rowMeans(z)), , drop = FALSE]
  signed_score_genes <- function(up_genes, down_genes) {
    up <- intersect(up_genes, rownames(z)); dn <- intersect(down_genes, rownames(z))
    colMeans(z[up, , drop = FALSE]) - colMeans(z[dn, , drop = FALSE])
  }
  sc_all <- signed_score_genes(up_all, down_all)
  sc_state <- lapply(STATES, function(st) {
    sig_genes <- sig_prim[state == st, gene]
    ov <- intersect(sig_genes, mod_no_frozen$gene)
    keep <- setdiff(mod_no_frozen$gene, ov)
    signed_score_genes(keep[mod_no_frozen[gene %in% keep]$direction == "up"],
                       keep[mod_no_frozen[gene %in% keep]$direction == "down"])
  })
  names(sc_state) <- STATES
  out <- data.table(
    cohort = coh, donor_id = cols,
    consensus = unname(sc_all[cols]),
    consensus_excl_reference_TAL    = unname(sc_state[["reference_TAL"]][cols]),
    consensus_excl_adaptive_TAL     = unname(sc_state[["adaptive_TAL"]][cols]),
    consensus_excl_degenerative_TAL = unname(sc_state[["degenerative_TAL"]][cols]),
    n_consensus_genes_no_frozen = length(up_all) + length(down_all),
    n_up_used_all = length(intersect(up_all, rownames(z))),
    n_down_used_all = length(intersect(down_all, rownames(z))))
  out
}

# primary: z-space = primary donors; sensitivity: z-space = sensitivity donors
score_prim <- rbind(score_donors("GSE183276", prim_don, prim_don[cohort == "GSE183276"]$donor_id),
                    score_donors("GSE210622", prim_don, prim_don[cohort == "GSE210622"]$donor_id))
score_sens <- rbind(score_donors("GSE183276", sens_don, sens_don[cohort == "GSE183276"]$donor_id),
                    score_donors("GSE210622", sens_don, sens_don[cohort == "GSE210622"]$donor_id))

score_prim[, tier := "primary_ge50"]
score_sens[, tier := "sensitivity_ge20"]
score_all <- rbind(score_prim, score_sens)

# merge state scores + frozen six-gene identity
dw_prim <- don_wide[aggregation == "mean_ge50"]
dw_sens <- don_wide[aggregation == "mean_ge20"]
dw_prim[, tier := "primary_ge50"]
dw_sens[, tier := "sensitivity_ge20"]
dw_all <- rbind(dw_prim, dw_sens)
setnames(dw_all,
         c("reference_TAL", "adaptive_TAL", "degenerative_TAL"),
         c("reference_TAL_state", "adaptive_TAL_state", "degenerative_TAL_state"))
donor_scores <- merge(score_all, dw_all, by = c("cohort", "donor_id", "tier"),
                      all.x = TRUE, suffixes = c("", "_state"))
stopifnot(!anyNA(donor_scores$reference_TAL_state))
stopifnot(!anyNA(donor_scores$adaptive_TAL_state))
stopifnot(!anyNA(donor_scores$degenerative_TAL_state))

# donor-level alignment (residuals on diagnosis) used by Figure 4f Panel B
# Both variables are residualized on AKI; the consensus side uses the
# per-state overlap-excluded consensus so the panel matches the analysis table.
alignment <- list()
for (tier_i in c("primary_ge50", "sensitivity_ge20")) {
  for (coh in c("GSE183276", "GSE210622")) {
    sub <- donor_scores[tier == tier_i & cohort == coh]
    sub[, AKI := as.integer(diagnosis == "AKI")]
    for (st in STATES) {
      ycol <- paste0("consensus_excl_", st)
      x <- sub[[paste0(st, "_state")]]
      ok <- is.finite(sub[[ycol]]) & is.finite(x)
      sub2 <- sub[ok]
      ry <- resid(lm(as.formula(paste0("`", ycol, "` ~ AKI")), data = sub2))
      rx <- resid(lm(as.formula(paste0("`", st, "_state` ~ AKI")), data = sub2))
      alignment[[length(alignment) + 1]] <- data.table(
        cohort = coh, donor_id = sub2$donor_id, diagnosis = sub2$diagnosis,
        state = st, tier = tier_i, consensus_resid = ry, state_resid = rx)
    }
  }
}
alignment_dt <- rbindlist(alignment)

fwrite(donor_scores, file.path(RES, "state_programme_bridge_donor_scores.csv"))
fwrite(alignment_dt, file.path(RES, "state_programme_bridge_alignment.csv"))
fwrite(overlap_audit, file.path(RES, "state_programme_bridge_overlap_audit.csv"))
msg(sprintf("  WROTE donor_scores (%d rows), alignment (%d rows), overlap_audit",
            nrow(donor_scores), nrow(alignment_dt)))

# =============================================================================
# 4. Donor-level association models (primary endpoint: reference_TAL;
#    secondary: adaptive_TAL, degenerative_TAL)
# =============================================================================
msg("--- 4. Donor-level association models ---")

fit_assoc <- function(coh, st, tier_i, donor_dt, scores_dt) {
  comp <- if (coh == "GSE183276") "Ref" else "non_AKI"
  d <- donor_dt[cohort == coh & tier == tier_i & diagnosis %in% c("AKI", comp)]
  sc <- scores_dt[tier == tier_i & cohort == coh]
  d <- merge(d[, .(donor_id)], sc, by = "donor_id")
  d[, AKI := as.integer(diagnosis == "AKI")]
  d[, state_score := get(paste0(st, "_state"))]
  y <- d$consensus; x <- d$state_score
  ok <- is.finite(y) & is.finite(x)
  d <- d[ok]
  n <- nrow(d)
  if (n < 5L) return(NULL)
  excl_col <- paste0("consensus_excl_", st)
  d[, consensus_for_state := get(excl_col)]
  rho <- suppressWarnings(cor.test(d$consensus_for_state, d$state_score, method = "spearman"))
  fit <- lm(consensus_for_state ~ AKI + state_score, data = d)
  s <- summary(fit)$coefficients; ci <- confint(fit)
  term <- "state_score"
  beta  <- unname(s[term, "Estimate"])
  beta_lo <- unname(ci[term, 1]); beta_hi <- unname(ci[term, 2])
  beta_p <- unname(s[term, "Pr(>|t|)"])
  # partial correlation: correlation of residuals given diagnosis
  ry <- resid(lm(consensus_for_state ~ AKI, data = d))
  rx <- resid(lm(state_score ~ AKI, data = d))
  r_partial <- suppressWarnings(cor(ry, rx))
  # collinearity
  r_dx_state <- suppressWarnings(cor(d$AKI, d$state_score))
  vif_dx <- if (is.finite(r_dx_state) && abs(r_dx_state) < 1) 1 / (1 - r_dx_state^2) else Inf
  # condition number of the CENTERED, unit-length-scaled two-predictor design
  # [AKI, state]: kappa = sqrt(lambda_max/lambda_min) of Xc'Xc after centering
  # AND scaling each column to unit length, so Xc'Xc is the correlation matrix
  # and kappa = sqrt((1+|r|)/(1-|r|)) (scale-free). The raw [1, AKI, state]
  # design's condition number is dominated by the intercept column and is not
  # reported as a collinearity diagnostic.
  Xc <- cbind(d$AKI - mean(d$AKI), d$state_score - mean(d$state_score))
  Xcs <- sweep(Xc, 2, sqrt(colSums(Xc^2)), "/")
  evc <- eigen(crossprod(Xcs), symmetric = TRUE, only.values = TRUE)$values
  cond_number <- sqrt(max(evc) / min(evc))
  n_cons <- nrow(mod_no_frozen) - length(intersect(sig_prim[state == st, gene], mod_no_frozen$gene))
  data.table(cohort = coh, state = st, tier = tier_i,
             n_donors = n, n_aki = sum(d$AKI), n_comp = sum(d$AKI == 0),
             n_consensus_genes = n_cons,
             rho = unname(rho$estimate), rho_p = rho$p.value,
             beta = beta, beta_ci_low = beta_lo, beta_ci_high = beta_hi,
             beta_p = beta_p, partial_r = r_partial,
             r_dx_state = r_dx_state, vif_dx = vif_dx, cond_number = cond_number,
             note = "")
}

assoc_prim <- rbindlist(lapply(STATES, function(st)
  fit_assoc("GSE183276", st, "primary_ge50", prim_don, donor_scores)))
assoc_prim <- rbind(assoc_prim, rbindlist(lapply(STATES, function(st)
  fit_assoc("GSE210622", st, "primary_ge50", prim_don, donor_scores))))
stopifnot(nrow(assoc_prim) == 6L)
assoc_prim[, fdr_bh := p.adjust(beta_p, method = "BH")]
msg("  Primary association table (consensus vs state, diagnosis-adjusted):")
print(assoc_prim[, .(cohort, state, n_donors, rho, rho_p, beta, beta_ci_low,
                     beta_ci_high, beta_p, partial_r, fdr_bh)])

assoc_sens <- rbindlist(lapply(STATES, function(st)
  fit_assoc("GSE183276", st, "sensitivity_ge20", sens_don, donor_scores)))
assoc_sens <- rbind(assoc_sens, rbindlist(lapply(STATES, function(st)
  fit_assoc("GSE210622", st, "sensitivity_ge20", sens_don, donor_scores))))
assoc_sens[, fdr_bh := NA_real_]

effects_all <- rbind(
  assoc_prim[, endpoint := "primary"],
  assoc_sens[, endpoint := "sensitivity"])
fwrite(effects_all, file.path(RES, "state_programme_bridge_effects.csv"))

# collinearity summary table (instability flag driven by VIF >= 10; the
# centred condition number is complementary: kappa >= ~6.2 is the VIF>=10
# equivalent for a two-predictor design)
collin <- effects_all[, .(cohort, state, tier = endpoint, n_donors,
                          r_dx_state, vif_dx, cond_number)]
collin[, unstable := vif_dx >= 10]
fwrite(collin, file.path(RES, "state_programme_bridge_collinearity.csv"))

# =============================================================================
# 5. LODO direction stability (leave-one-donor-out on the adjusted model)
# =============================================================================
msg("--- 5. LODO ---")
lodo_one <- function(coh, st, tier_i, donor_dt, scores_dt) {
  comp <- if (coh == "GSE183276") "Ref" else "non_AKI"
  d <- donor_dt[cohort == coh & tier == tier_i & diagnosis %in% c("AKI", comp)]
  sc <- scores_dt[tier == tier_i & cohort == coh]
  d <- merge(d[, .(donor_id)], sc, by = "donor_id")
  d[, AKI := as.integer(diagnosis == "AKI")]
  d[, state_score := get(paste0(st, "_state"))]
  d[, consensus_for_state := get(paste0("consensus_excl_", st))]
  d <- d[is.finite(consensus_for_state) & is.finite(state_score)]
  donors <- d$donor_id
  betas <- numeric(length(donors))
  for (i in seq_along(donors)) {
    sub <- d[donor_id != donors[i]]
    fit <- lm(consensus_for_state ~ AKI + state_score, data = sub)
    betas[i] <- unname(coef(fit)["state_score"])
  }
  full <- lm(consensus_for_state ~ AKI + state_score, data = d)
  full_beta <- unname(coef(full)["state_score"])
  data.table(cohort = coh, state = st, tier = tier_i,
             full_beta = full_beta, lodo_min = min(betas), lodo_max = max(betas),
             lodo_same_sign_frac = if (full_beta != 0)
               mean(sign(betas) == sign(full_beta)) else NA_real_,
             n_lodo = length(betas))
}
lodo_prim <- rbindlist(lapply(STATES, function(st)
  lodo_one("GSE183276", st, "primary_ge50", prim_don, donor_scores)))
lodo_prim <- rbind(lodo_prim, rbindlist(lapply(STATES, function(st)
  lodo_one("GSE210622", st, "primary_ge50", prim_don, donor_scores))))
lodo_sens <- rbindlist(lapply(STATES, function(st)
  lodo_one("GSE183276", st, "sensitivity_ge20", sens_don, donor_scores)))
lodo_sens <- rbind(lodo_sens, rbindlist(lapply(STATES, function(st)
  lodo_one("GSE210622", st, "sensitivity_ge20", sens_don, donor_scores))))
lodo_all <- rbind(lodo_prim[, endpoint := "primary"], lodo_sens[, endpoint := "sensitivity"])
fwrite(lodo_all, file.path(RES, "state_programme_bridge_lodo.csv"))

# per-donor LODO detail (long)
lodo_detail <- list()
for (tier_i in c("primary_ge50", "sensitivity_ge20")) {
  for (coh in c("GSE183276", "GSE210622")) {
    for (st in STATES) {
      comp <- if (coh == "GSE183276") "Ref" else "non_AKI"
      d <- (if (tier_i == "primary_ge50") prim_don else sens_don)[cohort == coh & diagnosis %in% c("AKI", comp)]
      sc <- donor_scores[tier == tier_i & cohort == coh]
      d <- merge(d[, .(donor_id)], sc, by = "donor_id")
      d[, AKI := as.integer(diagnosis == "AKI")]
      d[, state_score := get(paste0(st, "_state"))]
      d[, consensus_for_state := get(paste0("consensus_excl_", st))]
      d <- d[is.finite(consensus_for_state) & is.finite(state_score)]
      full <- lm(consensus_for_state ~ AKI + state_score, data = d)
      full_beta <- unname(coef(full)["state_score"])
      for (dd in d$donor_id) {
        sub <- d[donor_id != dd]
        fit <- lm(consensus_for_state ~ AKI + state_score, data = sub)
        lodo_detail[[length(lodo_detail) + 1]] <- data.table(
          cohort = coh, state = st, tier = tier_i, dropped_donor = dd,
          beta = unname(coef(fit)["state_score"]),
          same_sign_as_full = sign(unname(coef(fit)["state_score"])) == sign(full_beta))
      }
    }
  }
}
lodo_detail_dt <- rbindlist(lodo_detail)
fwrite(lodo_detail_dt, file.path(RES, "state_programme_bridge_lodo_detail.csv"))
msg("  LODO done")

# =============================================================================
# 6. Figure 4f (extension): donor-level scatter / residual scatter, forest,
#    stability. No cell-level P values anywhere.
# =============================================================================
msg("--- 6. Figure 4f ---")
col_aki <- "#C0392B"; col_comp <- "#7F8C8D"
col_ref <- "#2166AC"; col_adap <- "#B35806"; col_degen <- "#762A83"
state_cols <- c(reference_TAL = col_ref, adaptive_TAL = col_adap,
                degenerative_TAL = col_degen)
lab_diag <- function(x, coh) ifelse(x == "AKI", "AKI",
                              ifelse(coh == "GSE183276", "Reference", "non-AKI"))
state_lab <- function(x) c(reference_TAL = "Reference TAL",
                           adaptive_TAL = "Adaptive TAL",
                           degenerative_TAL = "Degenerative TAL")[x]

prim_plot <- donor_scores[tier == "primary_ge50"]
prim_plot[, diag_lab := factor(lab_diag(diagnosis, cohort),
                               levels = c("AKI", "Reference", "non-AKI"))]
prim_plot[, cohort_lab := factor(cohort, levels = c("GSE183276", "GSE210622"))]

# helper: annotation text for the reference-TAL test (panel A), using the
# SAME per-state overlap-excluded consensus as the analysis table
ref_test_text <- function(dat, coh) {
  dat <- copy(dat)
  dat[, AKI := as.integer(diagnosis == "AKI")]
  dat[, cons := consensus_excl_reference_TAL]
  fit <- lm(cons ~ AKI + reference_TAL_state, data = dat)
  b <- coef(fit)
  rho <- suppressWarnings(cor.test(dat$cons, dat$reference_TAL_state,
                                   method = "spearman"))
  rr <- suppressWarnings(cor(resid(lm(cons ~ AKI, data = dat)),
                             resid(lm(reference_TAL_state ~ AKI, data = dat))))
  list(txt = sprintf("rho = %.2f\nadj r = %.2f\nbeta = %.3f\nn = %d",
                     rho$estimate, rr, b[["reference_TAL_state"]], nrow(dat)),
       slope = b[["reference_TAL_state"]],
       intercept_mid = b[["(Intercept)"]] + b[["AKI"]] * 0.5)
}

# Panel A: donor-level scatter, per-state consensus (module minus frozen minus
# reference_TAL overlap) vs reference_TAL score, coloured by diagnosis;
# solid line = diagnosis-adjusted association, dashed = raw per-group lm.
pa_data  <- prim_plot[cohort == "GSE183276" & diagnosis %in% c("AKI", "Ref")]
pa_data2 <- prim_plot[cohort == "GSE210622" & diagnosis %in% c("AKI", "non_AKI")]
panelA_plot <- function(dat, coh, lab_letter, extra_note = NULL) {
  a <- ref_test_text(dat, coh)
  dat <- copy(dat); dat[, cons := consensus_excl_reference_TAL]
  sub <- paste0(coh, ": donor-level; solid = diagnosis-adjusted association, dashed = raw per-group lm")
  if (!is.null(extra_note)) sub <- paste0(sub, "\n", extra_note)
  ggplot(dat, aes(reference_TAL_state, cons, colour = diag_lab)) +
    geom_abline(slope = a$slope, intercept = a$intercept_mid,
                colour = "black", linewidth = 0.6) +
    geom_smooth(aes(y = cons, group = diag_lab), method = "lm", se = FALSE,
                linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
    geom_point(size = 2.2, alpha = 0.9) +
    annotate("text", x = -Inf, y = Inf, label = lab_letter,
             hjust = -0.4, vjust = 1.6, fontface = "bold", size = 5) +
    annotate("label", x = Inf, y = -Inf, label = a$txt, hjust = 1.05, vjust = -0.3,
             size = 3, label.size = 0.2, fill = "white", colour = "grey15") +
    scale_colour_manual(values = c(AKI = col_aki, Reference = col_comp,
                                   `non-AKI` = col_comp), name = "Diagnosis") +
    labs(x = "Reference-TAL score (donor mean UCell)",
         y = "Consensus score (z; frozen six genes and\nreference-TAL overlap excluded)",
         subtitle = sub) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8, colour = "grey25"))
}
pA1 <- panelA_plot(pa_data, "GSE183276", "A")
pA2 <- panelA_plot(pa_data2, "GSE210622", "",
                   extra_note = "GSE210622 reference: FDR=0.0775; VIF=12.17 (collinear)")

# Panel B: residual scatter (both variables residualized on diagnosis),
# 2 cohorts x 3 states, using per-state overlap-excluded consensus
pb_dat <- alignment_dt[tier == "primary_ge50"]
pb_dat[, cohort_lab := factor(cohort, levels = c("GSE183276", "GSE210622"))]
pb_dat[, state_lab := factor(state_lab(state),
                             levels = state_lab(STATES))]
pb_dat[, diag_lab := factor(lab_diag(diagnosis, cohort),
                            levels = c("AKI", "Reference", "non-AKI"))]
panelB_plot <- function(st, coh, show_y, show_x) {
  dat <- pb_dat[state == st & cohort == coh]
  fit <- lm(consensus_resid ~ state_resid + 0, data = dat)
  b <- unname(coef(fit)["state_resid"])
  rr <- suppressWarnings(cor(dat$consensus_resid, dat$state_resid))
  txt <- sprintf("slope = %.3f\npartial r = %.2f\nn = %d", b, rr, nrow(dat))
  primary <- (st == "reference_TAL")
  tag <- if (primary) "  [primary]" else "  [secondary]"
  p <- ggplot(dat, aes(state_resid, consensus_resid)) +
    geom_abline(slope = b, intercept = 0, colour = "black", linewidth = 0.5) +
    annotate("label", x = Inf, y = -Inf, label = txt, hjust = 1.05, vjust = -0.3,
             size = 2.6, label.size = 0.2, fill = "white", colour = "grey15")
  if (primary) {
    p <- p + geom_point(aes(colour = diag_lab), size = 1.8, alpha = 0.9) +
      scale_colour_manual(values = c(AKI = col_aki, Reference = col_comp,
                                     `non-AKI` = col_comp), guide = "none")
  } else {
    p <- p + geom_point(size = 1.6, alpha = 0.55, colour = "grey55")
  }
  p + labs(x = if (show_x) "State score residual (on dx)" else "",
           y = if (show_y) "Consensus residual" else "",
           title = paste0(state_lab(st), " \u2014 ", coh, tag)) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "plain"),
          panel.grid.minor = element_blank())
}
pB <- list()
for (st in STATES) {
  pB[[paste0("B_183_", st)]] <- panelB_plot(st, "GSE183276",
                                            show_y = TRUE, show_x = (st == "degenerative_TAL"))
  pB[[paste0("B_622_", st)]] <- panelB_plot(st, "GSE210622",
                                            show_y = FALSE, show_x = (st == "degenerative_TAL"))
}

# Panel C: forest of diagnosis-adjusted state coefficients (six primary tests)
eff_c <- effects_all[endpoint == "primary"]
eff_c[, cohort_lab := factor(cohort, levels = c("GSE183276", "GSE210622"))]
eff_c[, state_lab := factor(state_lab(state), levels = state_lab(STATES))]
eff_c[, sig := fifelse(!is.na(fdr_bh) & fdr_bh < 0.05, "FDR<0.05", "ns")]
eff_c[, primary := state == "reference_TAL"]
eff_c[, pt_col := ifelse(primary, col_ref, "grey60")]
eff_c[, rowlab := paste0(cohort_lab, " | ", state_lab)]
eff_c[cohort == "GSE210622" & state == "reference_TAL",
      rowlab := "GSE210622 | Reference TAL (FDR=0.078; VIF=12.17, collinear)"]
eff_c[, rowlab := factor(rowlab, levels = rev(unique(rowlab)))]
pC <- ggplot(eff_c, aes(beta, rowlab, colour = pt_col, shape = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_pointrange(aes(xmin = beta_ci_low, xmax = beta_ci_high),
                  size = 0.6, linewidth = 0.7) +
  scale_colour_identity() +
  scale_shape_manual(values = c("FDR<0.05" = 16, "ns" = 1), name = NULL) +
  labs(x = "Diagnosis-adjusted state coefficient (consensus ~ AKI + state)",
       y = NULL,
       title = "State effect on consensus score (primary, six tests)",
       subtitle = "95% CI; BH-FDR across the six tests. Reference TAL (blue) is the primary bridge; adaptive/degenerative (grey) are secondary.") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey25"),
        legend.position = "bottom")

# Panel D: LODO stability
lodo_c <- lodo_prim
lodo_c[, cohort_lab := factor(cohort, levels = c("GSE183276", "GSE210622"))]
lodo_c[, state_lab := factor(state_lab(state), levels = state_lab(STATES))]
lodo_c[, primary := state == "reference_TAL"]
lodo_c[, pt_col := ifelse(primary, col_ref, "grey60")]
lodo_c[, rowlab := paste0(cohort_lab, " | ", state_lab)]
lodo_c[, rowlab := factor(rowlab, levels = rev(unique(rowlab)))]
pD <- ggplot(lodo_c, aes(lodo_same_sign_frac, rowlab, colour = pt_col)) +
  geom_segment(aes(x = 0, xend = 1, y = rowlab, yend = rowlab), colour = "grey85") +
  geom_point(size = 2.8) +
  geom_text(aes(label = sprintf("%.2f", lodo_same_sign_frac)), hjust = -0.4,
            size = 2.8, colour = "grey20") +
  scale_x_continuous(limits = c(0, 1.35), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_colour_identity() +
  labs(x = "LODO same-sign fraction (diagnosis-adjusted beta)",
       y = NULL,
       title = "Leave-one-donor-out stability",
       subtitle = "Fraction of LODO refits with same sign as full-data beta; n = 27 (GSE183276) / 12 (GSE210622). Reference TAL (blue) is the primary bridge; adaptive/degenerative (grey) are secondary.") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8, colour = "grey25"))

# assemble with gridExtra (patchwork unavailable in this env)
gA <- arrangeGrob(pA1, pA2, ncol = 2, widths = c(1, 1))
gB <- arrangeGrob(grobs = list(pB[["B_183_reference_TAL"]], pB[["B_622_reference_TAL"]],
                               pB[["B_183_adaptive_TAL"]], pB[["B_622_adaptive_TAL"]],
                               pB[["B_183_degenerative_TAL"]], pB[["B_622_degenerative_TAL"]]),
                  ncol = 2, nrow = 3)
gB <- arrangeGrob(gB, top = textGrob("B  Diagnosis-adjusted residual scatter (consensus vs state, both variables residualized on AKI)",
                                     gp = gpar(fontsize = 10)))
gC <- arrangeGrob(pC, top = textGrob("C  Diagnosis-adjusted state effects (forest)",
                                     gp = gpar(fontsize = 10)))
gD <- arrangeGrob(pD, top = textGrob("D  LODO direction stability",
                                     gp = gpar(fontsize = 10)))

bottom_caption <- paste0(
  "Consensus score = TAL-derived AKI consensus module (z, frozen directions); six frozen genes always excluded;\n",
  "per-state overlap genes additionally excluded. Reference TAL is the primary bridge; adaptive/degenerative are secondary.\n",
  "Association only: no mediation, mechanism, causality, TAL-specificity, PT-independence or abundance claims.")

fig <- arrangeGrob(
  gA, gB, gC, gD, ncol = 1,
  heights = c(1.15, 2.1, 1.0, 0.85),
  top = textGrob("Figure 4f \u2014 State\u2013programme bridge (internal explanatory analysis; not independent validation)",
                 gp = gpar(fontface = "bold", fontsize = 12)),
  bottom = textGrob(bottom_caption, gp = gpar(fontsize = 8, col = "grey30", lineheight = 0.9))
)
ggsave(file.path(FIG, "figure4f_state_programme_bridge.png"), fig,
       width = 11, height = 13, dpi = 300, bg = "white")
ggsave(file.path(FIG, "figure4f_state_programme_bridge.pdf"), fig,
       width = 11, height = 13, bg = "white")
msg("  WROTE figures/extension/figure4f_state_programme_bridge.{png,pdf}")

# =============================================================================
# 7. Report
# =============================================================================
msg("--- 7. Report ---")
rep <- c(
  "# Figure 4f \u2014 State\u2013programme bridge (Figure 3 consensus module vs Figure 4 TAL states)",
  "",
  "**Status:** FINAL (extension analysis)",
  paste0("**Date:** ", format(Sys.time(), "%Y-%m-%d")),
  "**Analysis:** `code/extension/05_state_programme_bridge.R`",
  "**Figure:** `figures/extension/figure4f_state_programme_bridge.png/.pdf`",
  "",
  "## 1. Purpose and boundary",
  "",
  "This is an **internal explanatory bridge** connecting the Figure 3 TAL-derived AKI consensus module",
  "(1,027 genes; 723 up / 304 down; frozen directions) to the Figure 4 continuous TAL state programmes",
  "(reference-TAL 17 genes; adaptive-TAL 281 genes; degenerative-TAL 767 genes; frozen primary sets)",
  "**within the same two cell-resolved cohorts at the donor level**. It is **not** an independent",
  "validation. It supports no mediation, state-transition, mechanism, causality, TAL-specificity,",
  "PT-independence or biological-abundance claim. All null, opposite-direction and non-identifiable",
  "results are retained as-is.",
  "",
  "## 2. Data and donor sets",
  "",
  "- **Pseudobulk:** strict TAL donor-pseudobulk counts from `../AKI_DCT_analysis/results/stage3b_pseudobulk.rds`",
  "  (frozen stage3b), with the frozen per-segment TMM factors reused.",
  "- **Primary donors (>=50 strict TAL cells, AKI vs comparator):** GSE183276 10 AKI vs 17 Reference",
  "  (n=27); GSE210622 8 AKI vs 4 non-AKI (n=12).",
  "- **Sensitivity donors (>=20 cells):** GSE183276 adds 1153EO1 (Ref, 47 cells) and 30-10034 (AKI, 29 cells)",
  "  (n=29); GSE210622 unchanged (n=12).",
  "- **Cohorts are analysed independently; donors are never pooled across cohorts.**",
  "",
  "## 3. Scoring and anti-circularity",
  "",
  "- Per cohort: TMM logCPM on strict TAL pseudobulk (TMM factors reused from stage3b), then gene-wise",
  "  z-score across the analysis donors of that tier.",
  "- Consensus signed score = mean(z of module-up genes) - mean(z of module-down genes) using the frozen",
  "  module directions.",
  "- **Anti-circularity:** (a) the six frozen endpoint genes (SLC12A1, UMOD, CLDN16, KCNJ10, PTH1R, CLDN10)",
  "  are **always** excluded from the consensus score (3 are members of the 1,027-gene module); (b) for each",
  "  state test, all consensus genes overlapping that state's primary signature are additionally excluded.",
  "- State scores are the frozen Figure 4 donor scores (primary `mean_ge50`; sensitivity `mean_ge20`); the",
  "  state signatures exclude the six frozen genes by construction.",
  "",
  "### Overlap audit (consensus module vs state primary signatures)",
  "",
  "| State | State signature genes | Overlap with consensus | Overlap up | Overlap down | Consensus genes remaining after frozen + overlap exclusion |",
  "|---|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(overlap_audit))) {
  oa <- overlap_audit[i]
  rep <- c(rep, sprintf("| %s | %d | %d | %d | %d | %d |",
                        oa$state, oa$n_sig_genes, oa$n_overlap,
                        oa$n_overlap_up, oa$n_overlap_down,
                        oa$n_remaining_no_frozen))
}
rep <- c(rep,
  "",
  "## 4. Primary and secondary donor-level associations (>=50 tier)",
  "",
  "Primary endpoint: consensus score vs **reference_TAL** score, diagnosis-adjusted. Adaptive_TAL and",
  "degenerative_TAL are secondary. BH-FDR across the six cohort x state tests.",
  "",
  "| Cohort | State | n | rho | rho P | beta (95% CI) | beta P | partial r | FDR |",
  "|---|---|---:|---:|---:|---|---:|---:|---:|"
)
for (i in seq_len(nrow(assoc_prim))) {
  a <- assoc_prim[i]
  rep <- c(rep, sprintf("| %s | %s | %d | %.3f | %.3g | %.3f (%.3f, %.3f) | %.3g | %.3f | %.3g |",
                        a$cohort, a$state, a$n_donors, a$rho, a$rho_p,
                        a$beta, a$beta_ci_low, a$beta_ci_high, a$beta_p,
                        a$partial_r, a$fdr_bh))
}
rep <- c(rep,
  "",
  "`beta` is the state-score coefficient in `consensus ~ AKI + state_score`; `partial r` is the correlation",
  "of the two variables after residualizing on AKI. **Positive beta** means a higher state programme score",
  "accompanies a higher (more AKI-like) consensus score after diagnosis adjustment.",
  "",
  "### Bridge conclusions (primary tier)",
  "",
  "1. **reference-TAL is strongly supported in GSE183276:** rho = -0.87, adjusted beta = -13.20",
  "   (95% CI -16.23, -10.17), P = 3.8e-09, BH-FDR = 2.3e-08; LODO same-sign 1.00.",
  "2. **GSE210622 is directionally consistent but provides only collinear directional support:**",
  "   rho = -0.85, adjusted beta = -9.40 (95% CI -18.19, -0.60), but BH-FDR = 0.0775 and",
  "   VIF = 12.17 (Section 6). Because the adjusted test is not FDR-significant and the",
  "   diagnosis-reference-TAL pair is strongly collinear, this is **not** independent or complete",
  "   cross-cohort replication; it corroborates the direction of the GSE183276 finding only.",
  "3. **adaptive_TAL and degenerative_TAL do not replicate:** neither is FDR-significant in",
  "   GSE210622. GSE183276 degenerative_TAL reaches FDR = 0.043, but its raw and adjusted",
  "   directions are opposite (Section 7) and it does not reproduce in GSE210622, so it is",
  "   treated as an exploratory result only.",
  "",
  "## 5. Sensitivity (>=20 tier)",
  "",
  "| Cohort | State | n | rho | rho P | beta (95% CI) | beta P | partial r |",
  "|---|---|---:|---:|---:|---|---:|---:|"
)
for (i in seq_len(nrow(assoc_sens))) {
  a <- assoc_sens[i]
  rep <- c(rep, sprintf("| %s | %s | %d | %.3f | %.3g | %.3f (%.3f, %.3f) | %.3g | %.3f |",
                        a$cohort, a$state, a$n_donors, a$rho, a$rho_p,
                        a$beta, a$beta_ci_low, a$beta_ci_high, a$beta_p, a$partial_r))
}
rep <- c(rep,
  "",
  "## 6. Collinearity diagnostics (>=50 tier)",
  "",
  "| Cohort | State | cor(AKI, state) | VIF (diagnosis) | Condition number (centred [AKI, state]) | Unstable? |",
  "|---|---|---:|---:|---:|---|"
)
for (i in seq_len(nrow(collin[tier == "primary"]))) {
  c2 <- collin[tier == "primary"][i]
  rep <- c(rep, sprintf("| %s | %s | %.3f | %.2f | %.2f | %s |",
                        c2$cohort, c2$state, c2$r_dx_state, c2$vif_dx,
                        c2$cond_number, ifelse(c2$unstable, "yes", "no")))
}
rep <- c(rep,
  "",
  "`cor(AKI, state)` is the point-biserial correlation of the AKI indicator with the state score;",
  "`VIF` is the variance-inflation factor of the diagnosis term in `consensus ~ AKI + state`",
  "(= 1/(1-r^2) for this two-predictor design); the condition number is computed on the centred",
  "two-predictor design [AKI, state] (scale-free; equals sqrt((1+|r|)/(1-|r|))), not on the raw",
  "[1, AKI, state] matrix whose condition number is dominated by the intercept column.",
  "",
  "A VIF >= 10 (equivalently |cor(AKI, state)| >= ~0.95) marks a diagnosis-state pair too collinear",
  "for the adjusted coefficient to be interpreted as an independent effect. **GSE210622",
  "reference_TAL is in this class** (cor -0.958, VIF 12.17); its adjusted beta is retained in the",
  "table and figure but must not be read as an independent association or as evidence of",
  "reversal / mediation. Other pairs have moderate |r| (0.31-0.66) and VIF < 2: adjustment is",
  "better-conditioned there, though all state tests remain small-n.",
  "",
  "## 7. LODO direction stability (>=50 tier)",
  "",
  "| Cohort | State | Full beta | LODO min | LODO max | Same-sign fraction | n LODO |",
  "|---|---|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(lodo_prim))) {
  l <- lodo_prim[i]
  rep <- c(rep, sprintf("| %s | %s | %.3f | %.3f | %.3f | %.3f | %d |",
                        l$cohort, l$state, l$full_beta, l$lodo_min, l$lodo_max,
                        l$lodo_same_sign_frac, l$n_lodo))
}
rep <- c(rep,
  "",
  "### Raw vs diagnosis-adjusted direction (explicit disclosure)",
  "",
  "For 4 of 6 primary tests the raw Spearman rho and the diagnosis-adjusted linear coefficient differ",
  "in sign, and the sensitivity tier reproduces the same sign flips. These flips are retained and",
  "disclosed. Their most plausible reading is **suppression under diagnosis adjustment with small-n",
  "instability**: the state scores correlate with diagnosis (Section 6), so the partial association",
  "between consensus and state after conditioning on AKI can switch sign even when the marginal",
  "association is stable. GSE210622 reference_TAL is strongly collinear with diagnosis but does not",
  "show a raw-versus-adjusted sign flip. The four observed sign flips occur only for adaptive_TAL and",
  "degenerative_TAL in the two cohorts. No causal or mediation interpretation is drawn.",
  "",
  "| Cohort | State | rho (raw) | beta (adjusted) | Sign flip? |",
  "|---|---|---:|---:|---|"
)
for (i in seq_len(nrow(effects_all[endpoint == "primary"]))) {
  e <- effects_all[endpoint == "primary"][i]
  flip <- (e$rho > 0) != (e$beta > 0)
  rep <- c(rep, sprintf("| %s | %s | %.3f | %.3f | %s |",
                        e$cohort, e$state, e$rho, e$beta, ifelse(flip, "yes", "no")))
}
rep <- c(rep,
  "",
  "## 8. Figure 4f panels",
  "",
  "- **A** \u2014 donor-level scatter of consensus score (six frozen genes excluded) vs reference-TAL score,",
  "  both cohorts; solid line = diagnosis-adjusted association, dashed = raw per-group lm; the GSE210622",
  "  panel is annotated 'FDR=0.0775; VIF=12.17 (collinear)'.",
  "- **B** \u2014 diagnosis-adjusted residual scatter (both variables residualized on AKI) for the three",
  "  state programmes in both cohorts; slope = adjusted association; reference TAL panels are marked",
  "  [primary] with diagnosis-coloured points, adaptive/degenerative are grey [secondary].",
  "- **C** \u2014 forest of the six diagnosis-adjusted state coefficients with 95% CI and BH-FDR annotation;",
  "  reference TAL in blue (primary), adaptive/degenerative in grey (secondary); the GSE210622 reference",
  "  row is annotated with its FDR and VIF.",
  "- **D** \u2014 leave-one-donor-out same-sign fraction for the six tests; reference TAL in blue (primary).",
  "- No cell-level P values are used anywhere in the figure.",
  "",
  "## 9. Interpretation boundary (non-negotiable)",
  "",
  "Supported: a strong donor-level diagnosis-adjusted reference-TAL association in GSE183276, with",
  "directionally concordant but collinear support in GSE210622; adaptive and degenerative associations",
  "do not replicate across cohorts.",
  "",
  "**Not claimed:** independent validation, mediation, state transition, mechanism, causality, TAL",
  "specificity, PT independence, or biological abundance. The two cohorts are not combined into a",
  "single effect. The reference_TAL association in GSE210622 is collinear with diagnosis (VIF 12.17)",
  "and is not FDR-significant (BH-FDR 0.0775); it provides collinear directional support only, not",
  "independent or complete cross-cohort replication of the GSE183276 finding. adaptive_TAL and",
  "degenerative_TAL do not replicate across cohorts; the GSE183276 degenerative_TAL FDR = 0.043 is",
  "exploratory only (raw/adjusted direction flip, no GSE210622 replication).",
  "",
  "## 10. Software and reproduction",
  "",
  paste0("- R ", R.version.string, "; edgeR ", as.character(packageVersion("edgeR")),
         "; data.table ", as.character(packageVersion("data.table")),
         "; ggplot2 ", as.character(packageVersion("ggplot2")), "."),
  "- Reproduction: `Rscript code/extension/05_state_programme_bridge.R <project-root>`",
  "- All inputs are frozen files; no frozen core file, manuscript, or existing Figure 1-6 was modified."
)
writeLines(rep, file.path(RES, "state_programme_bridge_report.md"))
msg("  WROTE results/extension/state_programme_bridge_report.md")

cat("===== OUTPUT SUMMARY =====\n")
cat("donor_scores rows:", nrow(donor_scores), "\n")
cat("alignment rows:", nrow(alignment_dt), "\n")
cat("effects rows:", nrow(effects_all), "\n")
cat("lodo rows:", nrow(lodo_all), "| lodo_detail rows:", nrow(lodo_detail_dt), "\n")
cat("overlap_audit rows:", nrow(overlap_audit), "\n")
cat("DONE\n")
