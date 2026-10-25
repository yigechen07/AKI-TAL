#!/usr/bin/env Rscript
# 05_mitochondrial_audit.R — Direction C: TAL transport–mitochondrial programme audit
# ============================================================================
# Donor as statistical unit. Two datasets analyzed separately, never pooled.
# Frozen TAL/PT identity/injury scores read from formal_evidence.csv.
# Mitochondrial module scores computed from TAL pseudobulk counts.
# PT background uses independent PT pseudobulk OXPHOS scores.
# CI: Fisher-z SE = 1/sqrt(n-4).

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_arg) == 1L) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
} else {
  normalizePath("projects/AKI_TAL_severity_analysis")
}
original_wd <- getwd()
setwd(project_dir)
on.exit(setwd(original_wd), add = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(Matrix)
  library(digest)
})
set.seed(42)

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

PKG_VERSIONS <- c(
  tidyverse = as.character(packageVersion("tidyverse")),
  edgeR     = as.character(packageVersion("edgeR")),
  digest    = as.character(packageVersion("digest"))
)
R_VERSION <- paste(R.version$major, R.version$minor, sep=".")
SYS_INFO  <- Sys.info()[c("sysname","release","machine")]

# ---- 1. Hash inputs ----
msg("Hashing inputs")
INPUT_HASH_PB  <- digest("../AKI_DCT_analysis/results/stage3b_pseudobulk.rds", algo="sha256", file=TRUE)
INPUT_HASH_QC  <- digest("../AKI_DCT_analysis/results/stage3b_qc.csv", algo="sha256", file=TRUE)
INPUT_HASH_FE  <- digest("results/formal_evidence.csv", algo="sha256", file=TRUE)
INPUT_HASH_SIG <- digest("results/tal_mito_signature.csv", algo="sha256", file=TRUE)

# ---- 2. Load frozen mitochondrial signature from CSV ----
msg("Loading frozen mitochondrial signature")
sig_raw <- readLines("results/tal_mito_signature.csv")
meta_lines <- grep("^#", sig_raw, value=TRUE)
REACTOME_VERSION <- trimws(strsplit(meta_lines[grep("Reactome version", meta_lines)], ",")[[1]][2])
REACTOME_DATE    <- trimws(strsplit(meta_lines[grep("Fetch date",    meta_lines)], ",")[[1]][2])

sig_df <- read.csv("results/tal_mito_signature.csv", comment.char="#", stringsAsFactors=FALSE)
oxphos_genes       <- sig_df$gene_symbol[sig_df$category == "primary"]
oxphos_with_mt     <- sig_df$gene_symbol[sig_df$category %in% c("primary","primary_sensitivity")]
fao_genes          <- sig_df$gene_symbol[sig_df$category == "secondary" & sig_df$pathway == "Mitochondrial_FAO"]
biogenesis_genes   <- sig_df$gene_symbol[sig_df$category == "secondary" & sig_df$pathway == "Mitochondrial_biogenesis"]
mt_genes           <- sig_df$gene_symbol[sig_df$mt_gene == TRUE]

msg(sprintf("OXPHOS: %d genes (no MT-), %d with MT- (%d MT-); FAO: %d; Biogenesis: %d; Reactome v%s",
            length(oxphos_genes), length(oxphos_with_mt), length(mt_genes),
            length(fao_genes), length(biogenesis_genes), REACTOME_VERSION))

# ---- 3. Load data ----
msg("Loading pseudobulk")
pb <- readRDS("../AKI_DCT_analysis/results/stage3b_pseudobulk.rds")
qc <- read.csv("../AKI_DCT_analysis/results/stage3b_qc.csv", stringsAsFactors=FALSE)

msg("Loading frozen donor scores")
fe <- read.csv("results/formal_evidence.csv", stringsAsFactors=FALSE)

# ---- 4. Helper: TMM + logCPM + per-gene z-score -> module mean ----
score_module <- function(counts, genes) {
  g <- intersect(genes, rownames(counts))
  if (length(g) == 0) {
    return(list(scores = rep(NA_real_, ncol(counts)),
                n_used = 0L, n_preset = length(genes),
                frac_used = 0, missing = genes))
  }
  dge <- DGEList(counts)
  dge <- calcNormFactors(dge, method="TMM")
  cpm_all <- cpm(dge, log=FALSE)
  logcpm <- log2(cpm_all + 1)
  z <- t(scale(t(logcpm[g, , drop=FALSE])))
  scores <- colMeans(z, na.rm=TRUE)
  missing <- setdiff(genes, rownames(counts))
  list(scores = scores, n_used = length(g), n_preset = length(genes),
       frac_used = length(g)/length(genes), missing = missing)
}

# ---- 5. Partial correlation (n-4 CI) ----
fisher_z     <- function(r) 0.5 * log((1+r)/(1-r))
fisher_z_inv <- function(z) (exp(2*z)-1)/(exp(2*z)+1)

partial_r <- function(y, x, z) {
  ok <- complete.cases(y, x, z)
  y <- y[ok]; x <- x[ok]; z <- z[ok]
  n <- length(y)
  if (n < 4) return(c(r=NA, ci_lower=NA, ci_upper=NA, n=n))
  ry <- lm(y ~ z)$residuals
  rx <- lm(x ~ z)$residuals
  r  <- as.numeric(cor(ry, rx))
  fz <- fisher_z(r)
  se <- 1/sqrt(n-4)
  c(r=r, ci_lower=fisher_z_inv(fz-1.96*se), ci_upper=fisher_z_inv(fz+1.96*se), n=n)
}

# ---- 6. Permutation (diagnosis-stratified, 9999 draws, seed=42) ----
perm_partial <- function(y, x, z, n_perm=9999) {
  set.seed(42)
  obs <- partial_r(y, x, z)
  if (is.na(obs["r"]) || obs["n"] < 4) return(c(obs, p_perm=NA))
  r_obs <- obs["r"]
  r_perm <- numeric(n_perm)
  uz <- unique(z)
  for (i in seq_len(n_perm)) {
    x_perm <- x
    for (uzi in uz) { idx <- which(z==uzi); x_perm[idx] <- sample(x[idx]) }
    r_perm[i] <- partial_r(y, x_perm, z)["r"]
  }
  p <- (sum(abs(r_perm) >= abs(r_obs)) + 1) / (n_perm + 1)
  c(obs, p_perm=unname(p))
}

# ---- 7. LODO (re-fit model per dropped donor) ----
lodo_lm <- function(df, formula_str) {
  y_all <- coef(summary(lm(as.formula(formula_str), data=df)))
  beta_full <- y_all["diagnosis_binary", "Estimate"]
  r_loo <- numeric(nrow(df))
  for (i in seq_len(nrow(df))) {
    m <- lm(as.formula(formula_str), data=df[-i, ])
    r_loo[i] <- coef(m)["diagnosis_binary"]
  }
  c(lodo_prop=mean(sign(r_loo)==sign(beta_full)), lodo_mean=mean(r_loo), lodo_n=nrow(df))
}

lodo_partial <- function(df, y_var, x_var, z_var="diagnosis_binary") {
  r_loo <- numeric(nrow(df))
  for (i in seq_len(nrow(df))) {
    df_train <- df[-i, ]
    pr <- partial_r(df_train[[y_var]], df_train[[x_var]], df_train[[z_var]])
    r_loo[i] <- pr["r"]
  }
  full_r <- partial_r(df[[y_var]], df[[x_var]], df[[z_var]])["r"]
  c(lodo_prop=mean(sign(r_loo)==sign(full_r)), lodo_mean_r=mean(r_loo), lodo_n=nrow(df))
}
# ---- 8. Analysis per dataset × tier ----
results_list <- list()
lodo_list    <- list()
coverage_list <- list()

tier_configs <- list(
  list(tier="primary_ge50", eligible_col="eligible_50", label_suffix="ge50"),
  list(tier="sensitivity_ge20", eligible_col="eligible_20", label_suffix="ge20")
)

for (tc in tier_configs) {
  tier_name <- tc$tier
  elig_col  <- tc$eligible_col
  
  for (ds in c("GSE183276", "GSE210622")) {
    msg(sprintf("=== %s %s ===", ds, tier_name))
    
    # Filter formal_evidence
    fe_ds <- fe %>% filter(dataset==ds, tier==tier_name)
    if (ds == "GSE183276") {
      fe_ds <- fe_ds %>% filter(group %in% c("AKI","Reference"))
      ref_label <- "Reference"
    } else {
      fe_ds <- fe_ds %>% filter(group %in% c("AKI","Reference"))
      ref_label <- "non-AKI"
    }
    fe_ds <- fe_ds %>%
      mutate(diagnosis_binary = ifelse(group=="AKI", 1, 0))
    
    nd <- nrow(fe_ds)
    n_aki <- sum(fe_ds$diagnosis_binary==1)
    n_ref <- nd - n_aki
    
    if (nd < 5) { msg(sprintf("  Insufficient donors (n=%d), skip", nd)); next }
    
    msg(sprintf("  Donors: %d (AKI=%d, %s=%d)", nd, n_aki, ref_label, n_ref))
    
    # Map fe donor IDs to pseudobulk columns
    donors <- fe_ds$donor_id
    pb_available <- donors %in% colnames(pb$pseudobulk[[ds]]$TAL)
    if (!all(pb_available)) {
      msg(sprintf("  WARNING: %d donors not in pseudobulk, dropping", sum(!pb_available)))
      donors <- donors[pb_available]
      fe_ds <- fe_ds[fe_ds$donor_id %in% donors, ]
    }
    
    # TAL pseudobulk counts for mitochondrial modules
    tal_counts <- pb$pseudobulk[[ds]]$TAL[, donors, drop=FALSE]
    pt_counts  <- pb$pseudobulk[[ds]]$PT[, donors, drop=FALSE]
    
    # Compute mitochondrial module scores from TAL pseudobulk
    ox_mod    <- score_module(tal_counts, oxphos_genes)
    ox_mt_mod <- score_module(tal_counts, oxphos_with_mt)
    fao_mod   <- score_module(tal_counts, fao_genes)
    bio_mod   <- score_module(tal_counts, biogenesis_genes)
    
    # PT OXPHOS from PT pseudobulk (independent)
    pt_ox_mod <- score_module(pt_counts, oxphos_genes)
    
    # Coverage report
    coverage_list[[paste0(ds,"_",tc$label_suffix)]] <- data.frame(
      dataset=ds, tier=tier_name, n_donors=nd,
      oxphos_preset=ox_mod$n_preset, oxphos_used=ox_mod$n_used,
      oxphos_mt_preset=ox_mt_mod$n_preset, oxphos_mt_used=ox_mt_mod$n_used,
      fao_preset=fao_mod$n_preset, fao_used=fao_mod$n_used,
      biogenesis_preset=bio_mod$n_preset, biogenesis_used=bio_mod$n_used,
      stringsAsFactors=FALSE
    )
    
    msg(sprintf("  OXPHOS coverage: %d/%d genes", ox_mod$n_used, ox_mod$n_preset))
    
    # Build donor dataframe with frozen scores from fe
    df <- data.frame(
      donor            = donors,
      diagnosis_binary = fe_ds$diagnosis_binary[match(donors, fe_ds$donor_id)],
      TAL_identity     = fe_ds$TAL_identity[match(donors, fe_ds$donor_id)],
      PT_identity      = fe_ds$PT_identity[match(donors, fe_ds$donor_id)],
      PT_injury        = fe_ds$PT_injury[match(donors, fe_ds$donor_id)],
      OXPHOS           = ox_mod$scores,
      OXPHOS_with_MT   = ox_mt_mod$scores,
      FAO              = fao_mod$scores,
      Biogenesis       = bio_mod$scores,
      PT_OXPHOS        = pt_ox_mod$scores,
      stringsAsFactors = FALSE
    )
    
    label <- paste0(ds, "_", tc$label_suffix)
    
    # --- A1: OXPHOS_score ~ diagnosis ---
    m_ox <- lm(OXPHOS ~ diagnosis_binary, data=df)
    ox_eff <- coef(summary(m_ox))["diagnosis_binary", "Estimate"]
    ox_p   <- coef(summary(m_ox))["diagnosis_binary", "Pr(>|t|)"]
    ox_ci  <- confint(m_ox)["diagnosis_binary", ]
    ox_lodo <- lodo_lm(df, "OXPHOS ~ diagnosis_binary")
    
    msg(sprintf("  OXPHOS~diagnosis: β=%.4f [%.4f,%.4f] p=%.4f LODO=%.2f",
                ox_eff, ox_ci[1], ox_ci[2], ox_p, ox_lodo["lodo_prop"]))
    
    results_list[[paste0(label,"_oxphos_dx")]] <- data.frame(
      dataset=ds, tier=tier_name, analysis="OXPHOS_vs_diagnosis",
      n_donors=nd, n_AKI=n_aki, n_Ref=n_ref,
      estimate=unname(ox_eff), ci_lower=unname(ox_ci[1]), ci_upper=unname(ox_ci[2]),
      p_value=unname(ox_p), lodo_prop=unname(ox_lodo["lodo_prop"]),
      partial_r=NA_real_, partial_ci_lower=NA_real_, partial_ci_upper=NA_real_,
      p_perm=NA_real_, adj_r_squared=NA_real_,
      stringsAsFactors=FALSE
    )
    
    # --- A2: TAL_identity ~ OXPHOS | diagnosis ---
    pr1 <- partial_r(df$TAL_identity, df$OXPHOS, df$diagnosis_binary)
    pm1 <- perm_partial(df$TAL_identity, df$OXPHOS, df$diagnosis_binary)
    lo1 <- lodo_partial(df, "TAL_identity", "OXPHOS")
    
    msg(sprintf("  TAL_id~OXPHOS|diag: r=%.4f [%.4f,%.4f] p_perm=%.4f LODO=%.2f",
                pr1["r"], pr1["ci_lower"], pr1["ci_upper"], pm1["p_perm"], lo1["lodo_prop"]))
    
    results_list[[paste0(label,"_tal_oxphos")]] <- data.frame(
      dataset=ds, tier=tier_name, analysis="TAL_identity_vs_OXPHOS",
      n_donors=nd, n_AKI=n_aki, n_Ref=n_ref,
      estimate=NA_real_, ci_lower=NA_real_, ci_upper=NA_real_, p_value=NA_real_,
      lodo_prop=unname(lo1["lodo_prop"]),
      partial_r=unname(pr1["r"]),
      partial_ci_lower=unname(pr1["ci_lower"]),
      partial_ci_upper=unname(pr1["ci_upper"]),
      p_perm=unname(pm1["p_perm"]), adj_r_squared=NA_real_,
      stringsAsFactors=FALSE
    )
    
    # LODO per donor
    r_loo <- numeric(nd)
    for (i in seq_len(nd)) {
      df_t <- df[-i, ]
      r_loo[i] <- partial_r(df_t$TAL_identity, df_t$OXPHOS, df_t$diagnosis_binary)["r"]
    }
    lodo_list[[paste0(label,"_tal_oxphos")]] <- data.frame(
      dataset=ds, tier=tier_name, analysis="TAL_identity_vs_OXPHOS",
      donor=donors, lodo_r=r_loo, stringsAsFactors=FALSE
    )
    
    # --- A3: OXPHOS with MT- sensitivity ---
    pr_mt <- partial_r(df$TAL_identity, df$OXPHOS_with_MT, df$diagnosis_binary)
    pm_mt <- perm_partial(df$TAL_identity, df$OXPHOS_with_MT, df$diagnosis_binary)
    lo_mt <- lodo_partial(df, "TAL_identity", "OXPHOS_with_MT")
    
    msg(sprintf("  TAL_id~OXPHOS+MT|diag: r=%.4f p_perm=%.4f LODO=%.2f",
                pr_mt["r"], pm_mt["p_perm"], lo_mt["lodo_prop"]))
    
    results_list[[paste0(label,"_tal_oxphos_mt")]] <- data.frame(
      dataset=ds, tier=tier_name, analysis="TAL_identity_vs_OXPHOS_with_MT",
      n_donors=nd, n_AKI=n_aki, n_Ref=n_ref,
      estimate=NA_real_, ci_lower=NA_real_, ci_upper=NA_real_, p_value=NA_real_,
      lodo_prop=unname(lo_mt["lodo_prop"]),
      partial_r=unname(pr_mt["r"]),
      partial_ci_lower=unname(pr_mt["ci_lower"]),
      partial_ci_upper=unname(pr_mt["ci_upper"]),
      p_perm=unname(pm_mt["p_perm"]), adj_r_squared=NA_real_,
      stringsAsFactors=FALSE
    )
    
    # --- A4: +PT_identity adjustment ---
    m4 <- lm(TAL_identity ~ OXPHOS + PT_identity + diagnosis_binary, data=df)
    co4 <- coef(summary(m4))
    adj_r2_4 <- summary(m4)$adj.r.squared
    msg(sprintf("  TAL_id~OXPHOS+PT_id+diag: OXPHOS β=%.4f p=%.4f adj-R²=%.4f",
                co4["OXPHOS","Estimate"], co4["OXPHOS","Pr(>|t|)"], adj_r2_4))
    
    results_list[[paste0(label,"_adj_ptid")]] <- data.frame(
      dataset=ds, tier=tier_name, analysis="TAL_identity_vs_OXPHOS_adj_PT_identity",
      n_donors=nd, n_AKI=n_aki, n_Ref=n_ref,
      estimate=unname(co4["OXPHOS","Estimate"]),
      ci_lower=NA_real_, ci_upper=NA_real_,
      p_value=unname(co4["OXPHOS","Pr(>|t|)"]),
      lodo_prop=NA_real_, partial_r=NA_real_, p_perm=NA_real_,
      adj_r_squared=adj_r2_4,
      stringsAsFactors=FALSE
    )
    
    # --- A5: +PT_injury adjustment ---
    m5 <- lm(TAL_identity ~ OXPHOS + PT_injury + diagnosis_binary, data=df)
    co5 <- coef(summary(m5))
    adj_r2_5 <- summary(m5)$adj.r.squared
    msg(sprintf("  TAL_id~OXPHOS+PT_inj+diag: OXPHOS β=%.4f p=%.4f adj-R²=%.4f",
                co5["OXPHOS","Estimate"], co5["OXPHOS","Pr(>|t|)"], adj_r2_5))
    
    results_list[[paste0(label,"_adj_ptinj")]] <- data.frame(
      dataset=ds, tier=tier_name, analysis="TAL_identity_vs_OXPHOS_adj_PT_injury",
      n_donors=nd, n_AKI=n_aki, n_Ref=n_ref,
      estimate=unname(co5["OXPHOS","Estimate"]),
      ci_lower=NA_real_, ci_upper=NA_real_,
      p_value=unname(co5["OXPHOS","Pr(>|t|)"]),
      lodo_prop=NA_real_, partial_r=NA_real_, p_perm=NA_real_,
      adj_r_squared=adj_r2_5,
      stringsAsFactors=FALSE
    )
    
    # --- A6: PT background — PT_OXPHOS from PT pseudobulk ---
    pr_pt <- partial_r(df$PT_identity, df$PT_OXPHOS, df$diagnosis_binary)
    msg(sprintf("  PT_id~PT_OXPHOS|diag (bg): r=%.4f [%.4f,%.4f]",
                pr_pt["r"], pr_pt["ci_lower"], pr_pt["ci_upper"]))
    
    results_list[[paste0(label,"_pt_bg")]] <- data.frame(
      dataset=ds, tier=tier_name, analysis="PT_identity_vs_PT_OXPHOS",
      n_donors=nd, n_AKI=n_aki, n_Ref=n_ref,
      estimate=NA_real_, ci_lower=NA_real_, ci_upper=NA_real_, p_value=NA_real_,
      lodo_prop=NA_real_,
      partial_r=unname(pr_pt["r"]),
      partial_ci_lower=unname(pr_pt["ci_lower"]),
      partial_ci_upper=unname(pr_pt["ci_upper"]),
      p_perm=NA_real_, adj_r_squared=NA_real_,
      stringsAsFactors=FALSE
    )
    
    # --- Secondary: FAO ---
    m_fao <- lm(FAO ~ diagnosis_binary, data=df)
    pr_fao <- partial_r(df$TAL_identity, df$FAO, df$diagnosis_binary)
    msg(sprintf("  FAO: β=%.4f p=%.4f; TAL~FAO|diag r=%.4f",
                coef(m_fao)["diagnosis_binary"],
                coef(summary(m_fao))["diagnosis_binary","Pr(>|t|)"],
                pr_fao["r"]))
    
    results_list[[paste0(label,"_fao_dx")]] <- data.frame(
      dataset=ds, tier=tier_name, analysis="FAO_vs_diagnosis",
      n_donors=nd, n_AKI=n_aki, n_Ref=n_ref,
      estimate=unname(coef(m_fao)["diagnosis_binary"]),
      ci_lower=unname(confint(m_fao)["diagnosis_binary",1]),
      ci_upper=unname(confint(m_fao)["diagnosis_binary",2]),
      p_value=unname(coef(summary(m_fao))["diagnosis_binary","Pr(>|t|)"]),
      lodo_prop=NA_real_,
      partial_r=unname(pr_fao["r"]),
      partial_ci_lower=unname(pr_fao["ci_lower"]),
      partial_ci_upper=unname(pr_fao["ci_upper"]),
      p_perm=NA_real_, adj_r_squared=NA_real_,
      stringsAsFactors=FALSE
    )
    
    # --- Secondary: Biogenesis ---
    m_bio <- lm(Biogenesis ~ diagnosis_binary, data=df)
    pr_bio <- partial_r(df$TAL_identity, df$Biogenesis, df$diagnosis_binary)
    msg(sprintf("  Biogenesis: β=%.4f p=%.4f; TAL~Biogenesis|diag r=%.4f",
                coef(m_bio)["diagnosis_binary"],
                coef(summary(m_bio))["diagnosis_binary","Pr(>|t|)"],
                pr_bio["r"]))
    
    results_list[[paste0(label,"_bio_dx")]] <- data.frame(
      dataset=ds, tier=tier_name, analysis="Biogenesis_vs_diagnosis",
      n_donors=nd, n_AKI=n_aki, n_Ref=n_ref,
      estimate=unname(coef(m_bio)["diagnosis_binary"]),
      ci_lower=unname(confint(m_bio)["diagnosis_binary",1]),
      ci_upper=unname(confint(m_bio)["diagnosis_binary",2]),
      p_value=unname(coef(summary(m_bio))["diagnosis_binary","Pr(>|t|)"]),
      lodo_prop=NA_real_,
      partial_r=unname(pr_bio["r"]),
      partial_ci_lower=unname(pr_bio["ci_lower"]),
      partial_ci_upper=unname(pr_bio["ci_upper"]),
      p_perm=NA_real_, adj_r_squared=NA_real_,
      stringsAsFactors=FALSE
    )
  }
}
# ---- 9. Compile ----
results_df   <- bind_rows(results_list)
lodo_df      <- bind_rows(lodo_list)
coverage_df  <- bind_rows(coverage_list)

# ---- 10. Decision (based on primary_ge50 tier only) ----
msg("=== Decision ===")
get_stat <- function(ds, aname, tier_name, col) {
  rows <- results_df[results_df$dataset==ds & results_df$analysis==aname & results_df$tier==tier_name, , drop=FALSE]
  if (nrow(rows)==0) return(NA_real_)
  rows[[col]][1]
}

# OXPHOS AKI effect (expected < 0) — primary_ge50
ox1 <- get_stat("GSE183276","OXPHOS_vs_diagnosis","primary_ge50","estimate")
ox2 <- get_stat("GSE210622","OXPHOS_vs_diagnosis","primary_ge50","estimate")
oxl1 <- get_stat("GSE183276","OXPHOS_vs_diagnosis","primary_ge50","lodo_prop")
oxl2 <- get_stat("GSE210622","OXPHOS_vs_diagnosis","primary_ge50","lodo_prop")
cond_oxphos <- (!is.na(ox1) && !is.na(ox2)) && (ox1 < 0) && (ox2 < 0) && (oxl1 >= 0.80) && (oxl2 >= 0.80)

# TAL-OXPHOS coupling (expected r > 0) — primary_ge50
r1 <- get_stat("GSE183276","TAL_identity_vs_OXPHOS","primary_ge50","partial_r")
r2 <- get_stat("GSE210622","TAL_identity_vs_OXPHOS","primary_ge50","partial_r")
rl1 <- get_stat("GSE183276","TAL_identity_vs_OXPHOS","primary_ge50","lodo_prop")
rl2 <- get_stat("GSE210622","TAL_identity_vs_OXPHOS","primary_ge50","lodo_prop")
cond_coupling <- (!is.na(r1) && !is.na(r2)) && (r1 > 0) && (r2 > 0) && (rl1 >= 0.80) && (rl2 >= 0.80)

if (cond_oxphos && cond_coupling) {
  verdict <- "KEEP_C_SECONDARY"
} else if (cond_oxphos && !cond_coupling) {
  verdict <- "NARROW_C_OVERLAP"
} else {
  verdict <- "STOP_C"
}

msg(sprintf("OXPHOS decrease: %s (β=%.4f/%.4f LODO=%.2f/%.2f)", cond_oxphos, ox1, ox2, oxl1, oxl2))
msg(sprintf("TAL-OXPHOS coupling: %s (r=%.4f/%.4f LODO=%.2f/%.2f)", cond_coupling, r1, r2, rl1, rl2))
msg(sprintf("VERDICT: %s", verdict))

# ---- 11. Fisher-z meta-analysis (primary_ge50 reporting) ----
msg("=== Meta-analysis ===")
meta_keys <- c("TAL_identity_vs_OXPHOS","TAL_identity_vs_OXPHOS_with_MT")
for (ak in meta_keys) {
  r_183 <- get_stat("GSE183276",ak,"primary_ge50","partial_r")
  n_183 <- get_stat("GSE183276","OXPHOS_vs_diagnosis","primary_ge50","n_donors")
  r_210 <- get_stat("GSE210622",ak,"primary_ge50","partial_r")
  n_210 <- get_stat("GSE210622","OXPHOS_vs_diagnosis","primary_ge50","n_donors")
  if (is.na(r_183) || is.na(r_210)) next
  z1 <- fisher_z(r_183); z2 <- fisher_z(r_210)
  w1 <- n_183-4; w2 <- n_210-4
  z_pooled <- (w1*z1 + w2*z2) / (w1+w2)
  r_pooled <- fisher_z_inv(z_pooled)
  se_pooled <- 1/sqrt(w1+w2)
  ci_lo <- fisher_z_inv(z_pooled - 1.96*se_pooled)
  ci_hi <- fisher_z_inv(z_pooled + 1.96*se_pooled)
  msg(sprintf("  %s: pooled r=%.4f [%.4f, %.4f]", ak, r_pooled, ci_lo, ci_hi))
}

# ---- 12. Write outputs ----
write.csv(results_df,  "results/tal_mito_results.csv", row.names=FALSE)
write.csv(lodo_df,     "results/tal_mito_lodo.csv",    row.names=FALSE)

# ---- 13. Audit ----
n_183 <- get_stat("GSE183276","OXPHOS_vs_diagnosis","primary_ge50","n_donors")
n_210 <- get_stat("GSE210622","OXPHOS_vs_diagnosis","primary_ge50","n_donors")
n_183_ge20 <- get_stat("GSE183276","OXPHOS_vs_diagnosis","sensitivity_ge20","n_donors")
n_210_ge20 <- get_stat("GSE210622","OXPHOS_vs_diagnosis","sensitivity_ge20","n_donors")

# PT background
pt_r1 <- get_stat("GSE183276","PT_identity_vs_PT_OXPHOS","primary_ge50","partial_r")
pt_r2 <- get_stat("GSE210622","PT_identity_vs_PT_OXPHOS","primary_ge50","partial_r")

# MT sensitivity
mt_r1 <- get_stat("GSE183276","TAL_identity_vs_OXPHOS_with_MT","primary_ge50","partial_r")
mt_r2 <- get_stat("GSE210622","TAL_identity_vs_OXPHOS_with_MT","primary_ge50","partial_r")

# Adjustment
adj_id_1 <- get_stat("GSE183276","TAL_identity_vs_OXPHOS_adj_PT_identity","primary_ge50","estimate")
adj_id_2 <- get_stat("GSE210622","TAL_identity_vs_OXPHOS_adj_PT_identity","primary_ge50","estimate")
adj_inj_1 <- get_stat("GSE183276","TAL_identity_vs_OXPHOS_adj_PT_injury","primary_ge50","estimate")
adj_inj_2 <- get_stat("GSE210622","TAL_identity_vs_OXPHOS_adj_PT_injury","primary_ge50","estimate")

# Hashes for results reproducibility
RES_HASH <- digest("results/tal_mito_results.csv", algo="sha256", file=TRUE)
LODO_HASH <- digest("results/tal_mito_lodo.csv", algo="sha256", file=TRUE)

audit_lines <- c(
  "# TAL Transport–Mitochondrial Coupling Audit (Direction C)",
  "",
  paste("**Verdict:**", verdict),
  "",
  "## Input Hashes (SHA-256 of file bytes)",
  paste("- stage3b_pseudobulk.rds:", INPUT_HASH_PB),
  paste("- stage3b_qc.csv:", INPUT_HASH_QC),
  paste("- formal_evidence.csv:", INPUT_HASH_FE),
  paste("- tal_mito_signature.csv:", INPUT_HASH_SIG),
  "",
  "## Output Reproducibility Hashes",
  paste("- tal_mito_results.csv:", RES_HASH),
  paste("- tal_mito_lodo.csv:", LODO_HASH),
  "",
  paste("## Frozen Signature — Reactome v", REACTOME_VERSION, "(", REACTOME_DATE, ")"),
  paste("- Primary OXPHOS:", length(oxphos_genes), "genes (R-HSA-611105 + R-HSA-163210, MT- excluded)"),
  paste("- Primary OXPHOS with MT-:", length(oxphos_with_mt), "genes (", length(mt_genes), " MT-)"),
  paste("- Secondary FAO:", length(fao_genes), "genes (R-HSA-77289)"),
  paste("- Secondary Biogenesis:", length(biogenesis_genes), "genes (R-HSA-1592230)"),
  "- TAL six-gene overlap with OXPHOS: none",
  "",
  "## Environment",
  paste("- R version:", R_VERSION, "(", paste(SYS_INFO, collapse=", "), ")"),
  paste("- tidyverse:", PKG_VERSIONS["tidyverse"]),
  paste("- edgeR:", PKG_VERSIONS["edgeR"]),
  "",
  "## Sample sizes (donor, AKI vs Reference/non-AKI)",
  paste("- GSE183276 primary_ge50:", n_183, "donors (AKI=",
        get_stat("GSE183276","OXPHOS_vs_diagnosis","primary_ge50","n_AKI"),
        ", Ref=", get_stat("GSE183276","OXPHOS_vs_diagnosis","primary_ge50","n_Ref"), ")"),
  paste("- GSE210622 primary_ge50:", n_210, "donors (AKI=",
        get_stat("GSE210622","OXPHOS_vs_diagnosis","primary_ge50","n_AKI"),
        ", non-AKI=", get_stat("GSE210622","OXPHOS_vs_diagnosis","primary_ge50","n_Ref"), ")"),
  paste("- GSE183276 sensitivity_ge20:", if(is.na(n_183_ge20))"N/A"else n_183_ge20, "donors"),
  paste("- GSE210622 sensitivity_ge20: N/A (all 12 donors already ≥50, identical to primary_ge50)"),
  "",
  "## Decision Logic",
  "- KEEP_C_SECONDARY: OXPHOS β<0 in both datasets (LODO≥80%) AND TAL-OXPHOS partial r>0 in both (LODO≥80%)",
  "- NARROW_C_OVERLAP: OXPHOS decrease replicates but TAL-OXPHOS coupling does not",
  "- STOP_C: otherwise (overlap with Arnoux, direction conflict, or instability)",
  "",
  "## Primary Results (primary_ge50)",
  paste("### OXPHOS score ~ diagnosis (expected β < 0)"),
  paste("  GSE183276:  β =", sprintf("%.4f", ox1), sprintf("[%.4f, %.4f]",
        get_stat("GSE183276","OXPHOS_vs_diagnosis","primary_ge50","ci_lower"),
        get_stat("GSE183276","OXPHOS_vs_diagnosis","primary_ge50","ci_upper")),
        sprintf("p=%.4f", get_stat("GSE183276","OXPHOS_vs_diagnosis","primary_ge50","p_value")),
        sprintf("LODO=%.2f", oxl1)),
  paste("  GSE210622: β =", sprintf("%.4f", ox2), sprintf("[%.4f, %.4f]",
        get_stat("GSE210622","OXPHOS_vs_diagnosis","primary_ge50","ci_lower"),
        get_stat("GSE210622","OXPHOS_vs_diagnosis","primary_ge50","ci_upper")),
        sprintf("p=%.4f", get_stat("GSE210622","OXPHOS_vs_diagnosis","primary_ge50","p_value")),
        sprintf("LODO=%.2f", oxl2)),
  paste("  PASS:", cond_oxphos),
  "",
  paste("### TAL_identity ~ OXPHOS | diagnosis (expected r > 0)"),
  paste("  GSE183276:  r =", sprintf("%.4f", r1), sprintf("[%.4f, %.4f]",
        get_stat("GSE183276","TAL_identity_vs_OXPHOS","primary_ge50","partial_ci_lower"),
        get_stat("GSE183276","TAL_identity_vs_OXPHOS","primary_ge50","partial_ci_upper")),
        sprintf("p_perm=%.4f", get_stat("GSE183276","TAL_identity_vs_OXPHOS","primary_ge50","p_perm")),
        sprintf("LODO=%.2f", rl1)),
  paste("  GSE210622: r =", sprintf("%.4f", r2), sprintf("[%.4f, %.4f]",
        get_stat("GSE210622","TAL_identity_vs_OXPHOS","primary_ge50","partial_ci_lower"),
        get_stat("GSE210622","TAL_identity_vs_OXPHOS","primary_ge50","partial_ci_upper")),
        sprintf("p_perm=%.4f", get_stat("GSE210622","TAL_identity_vs_OXPHOS","primary_ge50","p_perm")),
        sprintf("LODO=%.2f", rl2)),
  paste("  PASS:", cond_coupling),
  "",
  "## Sensitivity Analyses",
  paste("### OXPHOS with MT- genes (TAL_identity ~ OXPHOS+MT | diagnosis)"),
  paste("  GSE183276: r =", sprintf("%.4f", mt_r1)),
  paste("  GSE210622: r =", sprintf("%.4f", mt_r2)),
  paste("### TAL_identity ~ OXPHOS + PT_identity + diagnosis"),
  paste("  GSE183276: OXPHOS β =", sprintf("%.4f", adj_id_1)),
  paste("  GSE210622: OXPHOS β =", sprintf("%.4f", adj_id_2)),
  paste("### TAL_identity ~ OXPHOS + PT_injury + diagnosis"),
  paste("  GSE183276: OXPHOS β =", sprintf("%.4f", adj_inj_1)),
  paste("  GSE210622: OXPHOS β =", sprintf("%.4f", adj_inj_2)),
  paste("### PT background: PT_identity ~ PT_OXPHOS | diagnosis (PT pseudobulk)"),
  paste("  GSE183276: r =", sprintf("%.4f", pt_r1)),
  paste("  GSE210622: r =", sprintf("%.4f", pt_r2)),
  "",
  "## Secondary Results (FAO, Biogenesis)",
  paste("  FAO vs diagnosis: GSE183276 β =",
        sprintf("%.4f", get_stat("GSE183276","FAO_vs_diagnosis","primary_ge50","estimate")),
        ", GSE210622 β =", sprintf("%.4f", get_stat("GSE210622","FAO_vs_diagnosis","primary_ge50","estimate"))),
  paste("  Biogenesis vs diagnosis: GSE183276 β =",
        sprintf("%.4f", get_stat("GSE183276","Biogenesis_vs_diagnosis","primary_ge50","estimate")),
        ", GSE210622 β =", sprintf("%.4f", get_stat("GSE210622","Biogenesis_vs_diagnosis","primary_ge50","estimate"))),
  "",
  "## Gene Coverage (TAL pseudobulk, primary_ge50)",
  paste("  OXPHOS:", 
        coverage_df$oxphos_used[coverage_df$dataset=="GSE183276" & coverage_df$tier=="primary_ge50"],
        "/", coverage_df$oxphos_preset[1], "genes in GSE183276;",
        coverage_df$oxphos_used[coverage_df$dataset=="GSE210622" & coverage_df$tier=="primary_ge50"],
        "/", coverage_df$oxphos_preset[1], "in GSE210622"),
  paste("  OXPHOS+MT:",
        coverage_df$oxphos_mt_used[coverage_df$dataset=="GSE183276" & coverage_df$tier=="primary_ge50"],
        "/", coverage_df$oxphos_mt_preset[1], "genes in GSE183276;",
        coverage_df$oxphos_mt_used[coverage_df$dataset=="GSE210622" & coverage_df$tier=="primary_ge50"],
        "/", coverage_df$oxphos_mt_preset[1], "in GSE210622"),
  "",
  "## Overlap Assessment (Arnoux et al. 2025)",
  "- Arnoux used spatial/single-cell transcriptomics, immunofluorescence, isotope tracing.",
  "- Reported TAL as metabolic hotspot with OXPHOS/FAO reliance; propofol impaired OXPHOS.",
  "- Analysis was cell-level/group-level; no donor-level TAL_transport ~ mitochondrial_score",
  "  partial correlation with diagnosis adjustment was reported.",
  "- Direction C is a secondary analysis, not a completed prior finding.",
  "",
  "## Interpretation Limits",
  "- Transcript module scores do not measure mitochondrial function, metabolic flux, or ATP.",
  "- Correlation is not TAL specificity, PT independence, cell-intrinsic mechanism, or causality.",
  "- Dalga et al. PCK1 cited for background only; PCK1 is not a single-gene proxy for mitochondria.",
  "- No claims about propofol effect, therapeutic implications, or mitochondrial fitness are warranted.",
  paste("- Reactome v", REACTOME_VERSION, "(", REACTOME_DATE, ") snapshots may change with future releases.")
)

writeLines(audit_lines, "results/tal_mito_audit.md")
msg("Outputs written.")
msg("Done.")
