#!/usr/bin/env Rscript
# =============================================================================
# 03d_consensus_module_external_validation.R
# External validation of the 1,027-gene TAL-derived AKI consensus
# transcriptional response module in three independent bulk cohorts.
#
# Frozen inputs (results/extension/):
#   tal_pseudobulk_concordance.csv  -> replicated==TRUE (1,027 genes);
#                                      723 up / 304 down by the common
#                                      direction of the two cohort logFCs.
#                                      0 discordant genes in the replicated set.
#
# Cohort pipelines are REUSED verbatim from the frozen scripts:
#   GSE139061, GSE30718 -> code/01_bulk_validation.R (matrix, annotation,
#                          sample selection, donor collapse)
#   GSE126805           -> code/03_longitudinal_validation.R (matrix,
#                          metadata, paired-allograft design)
# No re-downloading. No gene filtering or direction adjustment based on
# external results.
#
# Main score (per cohort): per-gene z-score within the cohort, then
#   signed score = mean(z of up genes) - mean(z of down genes)
#   (up/down averaged separately so gene-count imbalance is not an issue).
# Sensitivity: per-gene percentile (rank/(n+1)) signed score.
#
# Models (statistical unit = sample/donor / unique donor / paired allograft):
#   GSE139061: AKI vs Reference, sample/donor unit
#   GSE30718:  frozen unique-donor primary set, transplant-AKI vs Stable
#   GSE126805: paired baseline_post - baseline_pre, paired allograft unit;
#              sign-flip permutation + leave-one-kidney-out
# Three main raw P values are Holm-corrected across cohorts.
# Connectivity only: sample-level correlation of the consensus signed score
# with the frozen six-gene TAL score.
#
# Claim boundary: bulk tissue can only support tissue-level detectability
# of the TAL-derived AKI consensus transcriptional response. No
# TAL-specific, cell-intrinsic, mechanism, causality, or biomarker claims.
# Negative or inconsistent results are preserved as-is.
#
# Outputs (extension dirs only):
#   results/extension/consensus_module_{definition,coverage,scores,effects,lodo}.csv
#   results/extension/consensus_module_report.md
#   figures/extension/figure3e_consensus_external_validation.{png,pdf}
# =============================================================================

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_arg) == 1L) {
  dirname(dirname(dirname(normalizePath(sub("^--file=", "", script_arg)))))
} else {
  normalizePath("projects/AKI_TAL_severity_analysis")
}
setwd(project_dir)

res_dir <- file.path(project_dir, "results", "extension")
fig_dir <- file.path(project_dir, "figures", "extension")
data_bulk <- file.path(project_dir, "data", "bulk")
data_gse  <- file.path(project_dir, "data", "gse126805")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

SET_SEED <- 42L
N_PERM   <- 9999L
set.seed(SET_SEED)

AKI_COL  <- "#B2182B"
BLUE     <- "#2166AC"
ORANGE   <- "#D95F02"
GREY     <- "#BDBDBD"

FROZEN_TAL <- c("SLC12A1", "UMOD", "CLDN16", "KCNJ10", "PTH1R", "CLDN10")

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))
fmt3 <- function(x) if (is.na(x)) "NA" else sprintf("%.3f", x)
fmt4 <- function(x) if (is.na(x)) "NA" else sprintf("%.4f", x)

# =============================================================================
# 1. Consensus module definition (frozen)
# =============================================================================
msg("--- 1. Consensus module definition ---")
conc <- read.csv(file.path(res_dir, "tal_pseudobulk_concordance.csv"),
                 stringsAsFactors = FALSE)
mod <- conc[conc$replicated, c("gene", "logFC_GSE183276", "logFC_GSE210622",
                               "FDR_GSE183276", "FDR_GSE210622", "max_fdr")]
stopifnot(nrow(mod) == 1027L)
mod$direction <- ifelse(mod$logFC_GSE183276 > 0, "up", "down")
stopifnot(all(sign(mod$logFC_GSE183276) == sign(mod$logFC_GSE210622)))
stopifnot(sum(mod$direction == "up") == 723L, sum(mod$direction == "down") == 304L)
n_up <- sum(mod$direction == "up"); n_dn <- sum(mod$direction == "down")
msg(sprintf("  Consensus module: %d genes (%d up, %d down)", nrow(mod), n_up, n_dn))

definition <- data.frame(
  gene = mod$gene,
  direction = mod$direction,
  logFC_GSE183276 = mod$logFC_GSE183276,
  logFC_GSE210622 = mod$logFC_GSE210622,
  FDR_GSE183276 = mod$FDR_GSE183276,
  FDR_GSE210622 = mod$FDR_GSE210622,
  max_fdr = mod$max_fdr,
  stringsAsFactors = FALSE
)
write.csv(definition, file.path(res_dir, "consensus_module_definition.csv"),
          row.names = FALSE)
msg("  WROTE results/extension/consensus_module_definition.csv")

# =============================================================================
# 2. Shared scoring helpers (frozen conventions from 01/03)
# =============================================================================
zscore_genes <- function(mat) {
  # genes x samples; mean/sd across the columns provided (the cohort's
  # z-scoring space, identical to the frozen scripts)
  mu <- rowMeans(mat, na.rm = TRUE)
  s  <- apply(mat, 1, sd, na.rm = TRUE)
  s[!is.finite(s) | s < 1e-8] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, s, "/")
}

percentile_genes <- function(mat) {
  # per-gene percentile across the same columns: rank/(n+1)
  res <- t(apply(mat, 1, function(x) {
    ok <- !is.na(x)
    r <- rep(NA_real_, length(x))
    r[ok] <- rank(x[ok], ties.method = "average") / (sum(ok) + 1)
    r
  }))
  rownames(res) <- rownames(mat)
  colnames(res) <- colnames(mat)
  res
}

signed_score <- function(zmat, up_genes, down_genes) {
  up <- intersect(up_genes, rownames(zmat))
  dn <- intersect(down_genes, rownames(zmat))
  colMeans(zmat[up, , drop = FALSE], na.rm = TRUE) -
    colMeans(zmat[dn, , drop = FALSE], na.rm = TRUE)
}

# Welch two-sided t-test + 95% CI on two group score vectors (unit-level)
two_group_effect <- function(x_grp1, x_grp2, label1, label2) {
  tt <- t.test(x_grp1, x_grp2)
  data.frame(
    estimate = unname(diff(rev(tt$estimate))),  # grp1 - grp2
    ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
    p_value = tt$p.value, se = unname(tt$stderr),
    mean_grp1 = unname(tt$estimate[1]), mean_grp2 = unname(tt$estimate[2]),
    n_grp1 = length(x_grp1), n_grp2 = length(x_grp2),
    stringsAsFactors = FALSE
  )
}

coverage_row <- function(cohort, mat, analysis_cols, up_genes, down_genes) {
  present <- intersect(rownames(mat), c(up_genes, down_genes))
  sub <- mat[present, analysis_cols, drop = FALSE]
  sds <- apply(sub, 1, sd, na.rm = TRUE)
  usable <- present[is.finite(sds) & sds > 1e-8]
  data.frame(
    cohort = cohort,
    n_consensus = n_up + n_dn,
    n_present = length(present),
    n_up_present = length(intersect(up_genes, present)),
    n_down_present = length(intersect(down_genes, present)),
    n_usable = length(usable),
    frac_present = length(present) / (n_up + n_dn),
    frac_usable = length(usable) / (n_up + n_dn),
    n_zero_var_present = length(present) - length(usable),
    note = "usable = present and sd>1e-8 across the analysis samples",
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 3. GSE139061 — Eadon bulk biopsy: AKI vs nephrectomy reference
#    Reuse frozen matrix/sample selection from 01_bulk_validation.R
# =============================================================================
msg("--- 2. GSE139061 ---")
df139 <- read.csv(file.path(data_bulk, "GSE139061_Eadon_processed_QN_101419.csv.gz"),
                  stringsAsFactors = FALSE, check.names = FALSE)
rownames(df139) <- df139$Geneid
df139$Geneid <- NULL
aki_cols_139 <- grep("^AKI", colnames(df139), value = TRUE)
ref_cols_139 <- grep("^REF", colnames(df139), value = TRUE)
mat139 <- as.matrix(df139)
msg(sprintf("  Matrix %d genes x %d samples (AKI=%d, REF=%d)",
            nrow(mat139), ncol(mat139), length(aki_cols_139), length(ref_cols_139)))
stopifnot(length(aki_cols_139) == 39L, length(ref_cols_139) == 9L)

# Frozen six-gene score (same computation as 01_bulk_validation.R)
tal139 <- intersect(FROZEN_TAL, rownames(mat139))
six139 <- colMeans(zscore_genes(mat139[tal139, , drop = FALSE]))

# Consensus module scores; z-space = all 48 samples (frozen convention)
z139 <- zscore_genes(mat139)
pct139 <- percentile_genes(mat139)
sc139 <- signed_score(z139, mod$gene[mod$direction == "up"],
                      mod$gene[mod$direction == "down"])
sr139 <- signed_score(pct139, mod$gene[mod$direction == "up"],
                      mod$gene[mod$direction == "down"])

cov139 <- coverage_row("GSE139061", mat139, colnames(mat139),
                       mod$gene[mod$direction == "up"],
                       mod$gene[mod$direction == "down"])
msg(sprintf("  Coverage: %d/%d present, %d usable", cov139$n_present,
            cov139$n_consensus, cov139$n_usable))
msg(sprintf("  Frozen six genes present: %d/6 (%s)", length(tal139),
            paste(tal139, collapse = ",")))

# Main effect: AKI - Reference (sample/donor unit)
eff139 <- two_group_effect(sc139[aki_cols_139], sc139[ref_cols_139], "AKI", "REF")
eff139_rank <- two_group_effect(sr139[aki_cols_139], sr139[ref_cols_139], "AKI", "REF")
msg(sprintf("  MAIN signed score AKI-REF: %.3f (95%% CI %.3f to %.3f), P=%.4f",
            eff139$estimate, eff139$ci_low, eff139$ci_high, eff139$p_value))

# LODO: leave-one-sample-out (sample == donor for these bulk biopsies)
lodo139 <- sapply(seq_along(c(aki_cols_139, ref_cols_139)), function(i) {
  g <- c(rep("AKI", length(aki_cols_139)), rep("REF", length(ref_cols_139)))
  v <- c(sc139[aki_cols_139], sc139[ref_cols_139])
  keep <- seq_along(v) != i
  mean(v[keep & g == "AKI"]) - mean(v[keep & g == "REF"])
})
lodo139_frac <- mean(sign(lodo139) == sign(eff139$estimate), na.rm = TRUE)
msg(sprintf("  LODO same-direction fraction: %.3f (%d/%d)", lodo139_frac,
            sum(sign(lodo139) == sign(eff139$estimate), na.rm = TRUE), length(lodo139)))

# =============================================================================
# 4. GSE30718 — Famulski transplant AKI vs stable transplant
#    Reuse frozen parsing, primary unique-donor set from 01_bulk_validation.R
# =============================================================================
msg("--- 3. GSE30718 ---")
gz30718 <- gzfile(file.path(data_bulk, "GSE30718_series_matrix.txt.gz"), "rt")
lines30718 <- readLines(gz30718); close(gz30718)
title_line <- grep("^!Sample_title", lines30718, value = TRUE)[1]
id_line    <- grep("^!Sample_geo_accession", lines30718, value = TRUE)[1]
parse_tabs <- function(x) strsplit(x, "\t")[[1]][-1]
titles30718_raw <- parse_tabs(title_line)
ids30718_raw    <- parse_tabs(id_line)
titles30718 <- gsub('"', '', titles30718_raw)
ids30718    <- gsub('"', '', ids30718_raw)

find_grp <- function(t) {
  if (grepl("Acute kidney injury|^AKI", t, ignore.case = TRUE)) return("AKI")
  if (grepl("pristine protocol biopsy|^pPB", t, ignore.case = TRUE)) return("Stable")
  if (grepl("Nephrectomy|^NEPH", t, ignore.case = TRUE)) return("Nephrectomy")
  return("Other")
}
grps30718 <- sapply(titles30718, find_grp)

parse_donor <- function(t, grp) {
  if (grp == "Stable") {
    m <- regmatches(t, regexpr("pPB\\d+", t))
    if (length(m) > 0) return(m)
  }
  if (grp == "AKI") {
    parts <- strsplit(t, " ")[[1]]
    aki_token <- parts[length(parts)]
    m <- regmatches(aki_token, regexpr("^\\d+", aki_token))
    if (length(m) > 0) return(m)
  }
  return(t)
}
donor30718 <- mapply(parse_donor, titles30718, grps30718)
meta30718 <- data.frame(sample_id = ids30718, donor_id = donor30718,
                        group = grps30718, title = titles30718,
                        stringsAsFactors = FALSE)

data_start <- grep("^!series_matrix_table_begin", lines30718) + 1
data_lines <- lines30718[data_start:length(lines30718)]
data_lines <- data_lines[data_lines != "!series_matrix_table_end"]
probe_ids <- character(length(data_lines) - 1)
mat_rows <- matrix(NA_real_, nrow = length(data_lines) - 1, ncol = nrow(meta30718))
for (i in seq_along(data_lines)[-1]) {
  parts <- strsplit(data_lines[i], "\t")[[1]]
  probe_ids[i - 1] <- gsub('"', '', parts[1])
  mat_rows[i - 1, ] <- suppressWarnings(as.numeric(parts[-1]))
}
rownames(mat_rows) <- probe_ids
colnames(mat_rows) <- meta30718$sample_id
msg(sprintf("  Matrix %d probes x %d samples", nrow(mat_rows), ncol(mat_rows)))

annot_gz <- gzfile(file.path(data_bulk, "GPL570.annot.gz"), "rt")
annot_lines <- readLines(annot_gz); close(annot_gz)
gene_to_probes <- list()
for (line in annot_lines) {
  if (nchar(line) == 0 || substr(line, 1, 1) %in% c("^", "!", "#")) next
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 3) {
    pid <- parts[1]; gs <- parts[3]
    if (pid != "" && gs != "" && gs != "---" && gs != "Gene symbol") {
      for (g in strsplit(gs, " /// ")[[1]]) {
        gene_to_probes[[g]] <- c(gene_to_probes[[g]], pid)
      }
    }
  }
}

collapse_probes_fast <- function(mat, genes, gene_to_probes) {
  genes <- intersect(genes, names(gene_to_probes))
  out <- matrix(NA_real_, nrow = 0, ncol = ncol(mat))
  for (gene in genes) {
    probes <- intersect(gene_to_probes[[gene]], rownames(mat))
    if (length(probes) == 0) next
    if (length(probes) == 1) {
      row_vals <- mat[probes, , drop = FALSE]
    } else {
      row_vals <- matrix(apply(mat[probes, , drop = FALSE], 2, median, na.rm = TRUE), nrow = 1)
    }
    rownames(row_vals) <- gene
    out <- rbind(out, row_vals)
  }
  out
}

all_mod_genes <- mod$gene
mat307_gene <- collapse_probes_fast(mat_rows, all_mod_genes, gene_to_probes)
msg(sprintf("  Consensus genes collapsed from probes: %d", nrow(mat307_gene)))

# Frozen primary unique-donor set (earliest biopsy per donor)
keep30718 <- meta30718$group %in% c("AKI", "Stable")
meta_sub <- meta30718[keep30718, ]
meta_sub$biopsy_order <- "single"
for (i in seq_len(nrow(meta_sub))) {
  if (meta_sub$donor_id[i] %in% c("49", "56")) {
    t <- meta_sub$title[i]
    parts <- strsplit(t, " ")[[1]]
    aki_token <- parts[length(parts)]
    if (grepl("A1", aki_token)) meta_sub$biopsy_order[i] <- "A1"
    else if (grepl("A2", aki_token)) meta_sub$biopsy_order[i] <- "A2"
  }
}
primary_keep <- meta_sub$biopsy_order %in% c("A1", "single") | meta_sub$group == "Stable"
meta_primary <- meta_sub[primary_keep, ]
primary_cols <- meta_primary$sample_id
msg(sprintf("  Primary unique-donor set: AKI=%d Stable=%d (total %d)",
            sum(meta_primary$group == "AKI"), sum(meta_primary$group == "Stable"),
            nrow(meta_primary)))

mat307_pri <- mat307_gene[, primary_cols, drop = FALSE]
# Frozen six-gene score (same computation as 01_bulk_validation.R):
# collapse ALL six frozen TAL genes from probes (UMOD/PTH1R/CLDN10 are in the
# annotation but not in the 1027 consensus module, so they are absent from
# mat307_gene and must be collapsed separately for the connectivity score).
tal307_mat <- collapse_probes_fast(mat_rows, FROZEN_TAL, gene_to_probes)
tal307 <- rownames(tal307_mat)
stopifnot(length(tal307) == 6L)
six307 <- colMeans(zscore_genes(tal307_mat[, primary_cols, drop = FALSE]))

# Consensus module scores; z-space = primary cohort only (frozen convention)
z307 <- zscore_genes(mat307_pri)
pct307 <- percentile_genes(mat307_pri)
sc307 <- signed_score(z307, mod$gene[mod$direction == "up"],
                      mod$gene[mod$direction == "down"])
sr307 <- signed_score(pct307, mod$gene[mod$direction == "up"],
                      mod$gene[mod$direction == "down"])

cov307 <- coverage_row("GSE30718", mat307_pri, primary_cols,
                       mod$gene[mod$direction == "up"],
                       mod$gene[mod$direction == "down"])
msg(sprintf("  Coverage: %d/%d present, %d usable", cov307$n_present,
            cov307$n_consensus, cov307$n_usable))
msg(sprintf("  Frozen six genes present: %d/6 (%s)", length(tal307),
            paste(tal307, collapse = ",")))

aki_idx307 <- meta_primary$group == "AKI"
stab_idx307 <- meta_primary$group == "Stable"
eff307 <- two_group_effect(sc307[aki_idx307], sc307[stab_idx307], "AKI", "Stable")
eff307_rank <- two_group_effect(sr307[aki_idx307], sr307[stab_idx307], "AKI", "Stable")
msg(sprintf("  MAIN signed score AKI-Stable: %.3f (95%% CI %.3f to %.3f), P=%.4f",
            eff307$estimate, eff307$ci_low, eff307$ci_high, eff307$p_value))

# LODO: leave-one-donor-out
donors_aki <- unique(meta_primary$donor_id[aki_idx307])
donors_stab <- unique(meta_primary$donor_id[stab_idx307])
all_donors307 <- c(donors_aki, donors_stab)
lodo307 <- sapply(all_donors307, function(d) {
  ka <- aki_idx307 & meta_primary$donor_id != d
  ks <- stab_idx307 & meta_primary$donor_id != d
  if (sum(ka) < 2 || sum(ks) < 2) return(NA_real_)
  mean(sc307[aki_idx307][ka[aki_idx307]]) - mean(sc307[stab_idx307][ks[stab_idx307]])
})
lodo307_frac <- mean(sign(lodo307) == sign(eff307$estimate), na.rm = TRUE)
msg(sprintf("  LODO same-direction fraction: %.3f (%d/%d)", lodo307_frac,
            sum(sign(lodo307) == sign(eff307$estimate), na.rm = TRUE), length(lodo307)))

# =============================================================================
# 5. GSE126805 — paired transplant-IRI allografts
#    Reuse frozen matrix/metadata/pairing from 03_longitudinal_validation.R
# =============================================================================
msg("--- 4. GSE126805 ---")
cts <- read.delim(gzfile(file.path(data_gse, "GSE126805_all.gene_counts.txt.gz")),
                  stringsAsFactors = FALSE, check.names = FALSE)
sample_cols <- grep("^sample\\.kidney_", colnames(cts), value = TRUE)
cnt <- as.matrix(cts[, sample_cols]); mode(cnt) <- "numeric"
gene_names <- cts$external_gene_name
gene_mean <- rowMeans(cnt, na.rm = TRUE)
ord <- order(gene_names, -gene_mean)
cnt <- cnt[ord, ]; gene_names <- gene_names[ord]
keep <- !duplicated(gene_names)
cnt <- cnt[keep, ]; gene_names <- gene_names[keep]
rownames(cnt) <- gene_names; colnames(cnt) <- sample_cols

parse_sample <- function(s) {
  parts <- strsplit(sub("^sample\\.", "", s), "_")[[1]]
  c(kidney = parts[2], timepoint = paste(parts[3:length(parts)], collapse = "_"))
}
meta126 <- as.data.frame(do.call(rbind, lapply(sample_cols, function(s) {
  p <- parse_sample(s)
  c(sample = s, kidney = unname(p["kidney"]), timepoint = unname(p["timepoint"]))
})), stringsAsFactors = FALSE)
stopifnot(all(meta126$timepoint %in% c("baseline_pre", "baseline_post", "3months", "1year")))

dge <- DGEList(counts = cnt)
dge <- calcNormFactors(dge, method = "TMM")
lcpm <- cpm(dge, log = TRUE, prior.count = 1)

tal126 <- intersect(FROZEN_TAL, rownames(lcpm))
six126 <- colMeans(zscore_genes(lcpm[tal126, , drop = FALSE]))

# Consensus module scores; z-space = all 163 qualified samples (frozen convention)
z126 <- zscore_genes(lcpm)
pct126 <- percentile_genes(lcpm)
sc126 <- signed_score(z126, mod$gene[mod$direction == "up"],
                      mod$gene[mod$direction == "down"])
sr126 <- signed_score(pct126, mod$gene[mod$direction == "up"],
                      mod$gene[mod$direction == "down"])

cov126 <- coverage_row("GSE126805", lcpm, sample_cols,
                       mod$gene[mod$direction == "up"],
                       mod$gene[mod$direction == "down"])
msg(sprintf("  Coverage: %d/%d present, %d usable", cov126$n_present,
            cov126$n_consensus, cov126$n_usable))
msg(sprintf("  Frozen six genes present: %d/6 (%s)", length(tal126),
            paste(tal126, collapse = ",")))

has_pair <- sapply(unique(meta126$kidney), function(k)
  all(c("baseline_pre", "baseline_post") %in% meta126$timepoint[meta126$kidney == k]))
kidneys_pair <- names(has_pair)[has_pair]
n_pairs <- length(kidneys_pair)
msg(sprintf("  Paired allografts: %d/42 (kidney 29 lacks baseline_pre)", n_pairs))
stopifnot(n_pairs >= 25L)

get_sc <- function(kid, tp, score_vec) {
  r <- meta126[meta126$kidney == kid & meta126$timepoint == tp, ]
  if (nrow(r) == 0) NA_real_ else score_vec[r$sample[1]]
}
pre_sc126  <- sapply(kidneys_pair, get_sc, "baseline_pre", sc126)
post_sc126 <- sapply(kidneys_pair, get_sc, "baseline_post", sc126)
deltas126 <- post_sc126 - pre_sc126
obs_d126 <- mean(deltas126)
tt126 <- t.test(deltas126)
eff126 <- data.frame(
  estimate = unname(tt126$estimate), ci_low = tt126$conf.int[1],
  ci_high = tt126$conf.int[2], p_value = tt126$p.value, se = unname(tt126$stderr),
  mean_grp1 = unname(tt126$estimate), mean_grp2 = 0,
  n_grp1 = n_pairs, n_grp2 = n_pairs, stringsAsFactors = FALSE
)
msg(sprintf("  MAIN paired post-pre delta: %.3f (95%% CI %.3f to %.3f), P=%.4f, n=%d",
            eff126$estimate, eff126$ci_low, eff126$ci_high, eff126$p_value, n_pairs))

pre126_r  <- sapply(kidneys_pair, get_sc, "baseline_pre", sr126)
post126_r <- sapply(kidneys_pair, get_sc, "baseline_post", sr126)
deltas126_r <- post126_r - pre126_r
tt126_r <- t.test(deltas126_r)
eff126_rank <- data.frame(
  estimate = unname(tt126_r$estimate), ci_low = tt126_r$conf.int[1],
  ci_high = tt126_r$conf.int[2], p_value = tt126_r$p.value,
  se = unname(tt126_r$stderr), mean_grp1 = unname(tt126_r$estimate),
  mean_grp2 = 0, n_grp1 = n_pairs, n_grp2 = n_pairs, stringsAsFactors = FALSE
)

# Sign-flip permutation (9999, seed 42, two-sided)
set.seed(SET_SEED)
ext <- 0L
for (b in seq_len(N_PERM)) {
  if (abs(mean(deltas126 * sample(c(-1, 1), n_pairs, TRUE))) >= abs(obs_d126))
    ext <- ext + 1L
}
perm126 <- (ext + 1) / (N_PERM + 1)
msg(sprintf("  Sign-flip permutation P: %.4f", perm126))

# LOKO: leave-one-kidney-out
loko126 <- sapply(seq_len(n_pairs), function(i) mean(deltas126[-i]))
loko126_frac <- mean(sign(loko126) == sign(eff126$estimate), na.rm = TRUE)
msg(sprintf("  LOKO same-direction fraction: %.3f (%d/%d)", loko126_frac,
            sum(sign(loko126) == sign(eff126$estimate), na.rm = TRUE), length(loko126)))

# =============================================================================
# 6. Holm across the three main P values
# =============================================================================
msg("--- 5. Holm correction ---")
raw_p <- c(GSE139061 = eff139$p_value, GSE30718 = eff307$p_value,
           GSE126805 = eff126$p_value)
holm_p <- p.adjust(raw_p, method = "holm")
msg(sprintf("  raw P: %.4f / %.4f / %.4f; Holm P: %.4f / %.4f / %.4f",
            raw_p[1], raw_p[2], raw_p[3], holm_p[1], holm_p[2], holm_p[3]))

# =============================================================================
# 5b. Leave-six-genes-out sensitivity
#     Exclude the frozen six TAL genes (SLC12A1, UMOD, CLDN16, KCNJ10, PTH1R,
#     CLDN10) that are actually present in the consensus module and recompute
#     effects, 95% CI, raw P, three-cohort Holm P and LODO/LOKO with the exact
#     same scoring method and the same per-cohort z-spaces. Output is an
#     independent CSV; the main analysis results are not modified.
# =============================================================================
msg("--- 5b. Leave-six-genes-out sensitivity ---")
frozen_in_module <- intersect(FROZEN_TAL, mod$gene)
l6o_up <- mod$gene[mod$direction == "up" & !(mod$gene %in% frozen_in_module)]
l6o_dn <- mod$gene[mod$direction == "down" & !(mod$gene %in% frozen_in_module)]
msg(sprintf("  Excluded frozen six genes present in module: %s",
            paste(frozen_in_module, collapse = ", ")))
msg(sprintf("  L6O module: %d genes (%d up, %d down)",
            length(l6o_up) + length(l6o_dn), length(l6o_up), length(l6o_dn)))
stopifnot(length(l6o_up) == n_up)
stopifnot(length(l6o_dn) == n_dn - length(frozen_in_module))

# GSE139061
l6o_sc139 <- signed_score(z139, l6o_up, l6o_dn)
l6o_eff139 <- two_group_effect(l6o_sc139[aki_cols_139], l6o_sc139[ref_cols_139], "AKI", "REF")
l6o_lodo139 <- sapply(seq_along(c(aki_cols_139, ref_cols_139)), function(i) {
  g <- c(rep("AKI", length(aki_cols_139)), rep("REF", length(ref_cols_139)))
  v <- c(l6o_sc139[aki_cols_139], l6o_sc139[ref_cols_139])
  keep <- seq_along(v) != i
  mean(v[keep & g == "AKI"]) - mean(v[keep & g == "REF"])
})
l6o_lodo139_frac <- mean(sign(l6o_lodo139) == sign(l6o_eff139$estimate), na.rm = TRUE)
msg(sprintf("  GSE139061 L6O AKI-REF: %.3f (95%% CI %.3f to %.3f), P=%.4f, LODO=%.3f",
            l6o_eff139$estimate, l6o_eff139$ci_low, l6o_eff139$ci_high,
            l6o_eff139$p_value, l6o_lodo139_frac))

# GSE30718
l6o_sc307 <- signed_score(z307, l6o_up, l6o_dn)
l6o_eff307 <- two_group_effect(l6o_sc307[aki_idx307], l6o_sc307[stab_idx307], "AKI", "Stable")
l6o_lodo307 <- sapply(all_donors307, function(d) {
  ka <- aki_idx307 & meta_primary$donor_id != d
  ks <- stab_idx307 & meta_primary$donor_id != d
  if (sum(ka) < 2 || sum(ks) < 2) return(NA_real_)
  mean(l6o_sc307[aki_idx307][ka[aki_idx307]]) - mean(l6o_sc307[stab_idx307][ks[stab_idx307]])
})
l6o_lodo307_frac <- mean(sign(l6o_lodo307) == sign(l6o_eff307$estimate), na.rm = TRUE)
msg(sprintf("  GSE30718 L6O AKI-Stable: %.3f (95%% CI %.3f to %.3f), P=%.4f, LODO=%.3f",
            l6o_eff307$estimate, l6o_eff307$ci_low, l6o_eff307$ci_high,
            l6o_eff307$p_value, l6o_lodo307_frac))

# GSE126805
l6o_sc126 <- signed_score(z126, l6o_up, l6o_dn)
l6o_pre126  <- sapply(kidneys_pair, get_sc, "baseline_pre", l6o_sc126)
l6o_post126 <- sapply(kidneys_pair, get_sc, "baseline_post", l6o_sc126)
l6o_deltas126 <- l6o_post126 - l6o_pre126
l6o_obs_d126 <- mean(l6o_deltas126)
l6o_tt126 <- t.test(l6o_deltas126)
l6o_eff126 <- data.frame(
  estimate = unname(l6o_tt126$estimate), ci_low = l6o_tt126$conf.int[1],
  ci_high = l6o_tt126$conf.int[2], p_value = l6o_tt126$p.value,
  se = unname(l6o_tt126$stderr), mean_grp1 = unname(l6o_tt126$estimate),
  mean_grp2 = 0, n_grp1 = n_pairs, n_grp2 = n_pairs, stringsAsFactors = FALSE
)
msg(sprintf("  GSE126805 L6O post-pre: %.3f (95%% CI %.3f to %.3f), P=%.4f, n=%d",
            l6o_eff126$estimate, l6o_eff126$ci_low, l6o_eff126$ci_high,
            l6o_eff126$p_value, n_pairs))
set.seed(SET_SEED)
l6o_ext <- 0L
for (b in seq_len(N_PERM)) {
  if (abs(mean(l6o_deltas126 * sample(c(-1, 1), n_pairs, TRUE))) >= abs(l6o_obs_d126))
    l6o_ext <- l6o_ext + 1L
}
l6o_perm126 <- (l6o_ext + 1) / (N_PERM + 1)
l6o_loko126 <- sapply(seq_len(n_pairs), function(i) mean(l6o_deltas126[-i]))
l6o_loko126_frac <- mean(sign(l6o_loko126) == sign(l6o_eff126$estimate), na.rm = TRUE)
msg(sprintf("  GSE126805 L6O perm P=%.4f, LOKO=%.3f (%d/%d)", l6o_perm126,
            l6o_loko126_frac,
            sum(sign(l6o_loko126) == sign(l6o_eff126$estimate), na.rm = TRUE),
            length(l6o_loko126)))

# Holm across the three L6O P values
l6o_raw_p <- c(GSE139061 = l6o_eff139$p_value, GSE30718 = l6o_eff307$p_value,
               GSE126805 = l6o_eff126$p_value)
l6o_holm_p <- p.adjust(l6o_raw_p, method = "holm")
msg(sprintf("  L6O raw P: %.4f / %.4f / %.4f; Holm P: %.4f / %.4f / %.4f",
            l6o_raw_p[1], l6o_raw_p[2], l6o_raw_p[3],
            l6o_holm_p[1], l6o_holm_p[2], l6o_holm_p[3]))

# Independent L6O CSV (same schema as consensus_module_effects.csv)
excl_note <- sprintf("leave-six-genes-out; excluded module-present frozen six genes: %s",
                     paste(frozen_in_module, collapse = ", "))
l6o_effects <- rbind(
  data.frame(cohort = "META", model = "excluded_genes",
             estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
             p_value = NA_real_, p_holm = NA_real_,
             n_group1 = NA_real_, n_group2 = NA_real_,
             group1 = paste(frozen_in_module, collapse = ", "), group2 = "frozen six genes",
             unit = "module-present only", sensitivity = "meta",
             note = excl_note, stringsAsFactors = FALSE),
  data.frame(cohort = "GSE139061", model = "l6o_signed_z",
             estimate = l6o_eff139$estimate, ci_low = l6o_eff139$ci_low,
             ci_high = l6o_eff139$ci_high, p_value = l6o_eff139$p_value,
             p_holm = l6o_holm_p["GSE139061"],
             n_group1 = l6o_eff139$n_grp1, n_group2 = l6o_eff139$n_grp2,
             group1 = "AKI", group2 = "Reference",
             unit = "sample (biopsy; sample == donor)",
             sensitivity = "sensitivity_leave6", note = excl_note,
             stringsAsFactors = FALSE),
  data.frame(cohort = "GSE30718", model = "l6o_signed_z",
             estimate = l6o_eff307$estimate, ci_low = l6o_eff307$ci_low,
             ci_high = l6o_eff307$ci_high, p_value = l6o_eff307$p_value,
             p_holm = l6o_holm_p["GSE30718"],
             n_group1 = l6o_eff307$n_grp1, n_group2 = l6o_eff307$n_grp2,
             group1 = "AKI", group2 = "Stable",
             unit = "unique donor (frozen primary set)",
             sensitivity = "sensitivity_leave6", note = excl_note,
             stringsAsFactors = FALSE),
  data.frame(cohort = "GSE126805", model = "l6o_signed_z",
             estimate = l6o_eff126$estimate, ci_low = l6o_eff126$ci_low,
             ci_high = l6o_eff126$ci_high, p_value = l6o_eff126$p_value,
             p_holm = l6o_holm_p["GSE126805"],
             n_group1 = l6o_eff126$n_grp1, n_group2 = l6o_eff126$n_grp2,
             group1 = "baseline_post", group2 = "baseline_pre",
             unit = "paired kidney allograft",
             sensitivity = "sensitivity_leave6", note = excl_note,
             stringsAsFactors = FALSE),
  data.frame(cohort = "GSE126805", model = "l6o_signflip_permutation",
             estimate = l6o_obs_d126, ci_low = NA_real_, ci_high = NA_real_,
             p_value = l6o_perm126, p_holm = NA_real_,
             n_group1 = n_pairs, n_group2 = n_pairs,
             group1 = "baseline_post", group2 = "baseline_pre",
             unit = "paired kidney allograft",
             sensitivity = "sensitivity_leave6",
             note = sprintf("9999 sign flips, seed %d; %s", SET_SEED, excl_note),
             stringsAsFactors = FALSE),
  data.frame(cohort = "GSE139061", model = "l6o_lodo_direction_fraction",
             estimate = l6o_lodo139_frac, ci_low = NA_real_, ci_high = NA_real_,
             p_value = NA_real_, p_holm = NA_real_,
             n_group1 = sum(sign(l6o_lodo139) == sign(l6o_eff139$estimate), na.rm = TRUE),
             n_group2 = length(l6o_lodo139),
             group1 = "same_direction", group2 = "omitted_units",
             unit = "sample (sample == donor)",
             sensitivity = "sensitivity_leave6", note = excl_note,
             stringsAsFactors = FALSE),
  data.frame(cohort = "GSE30718", model = "l6o_lodo_direction_fraction",
             estimate = l6o_lodo307_frac, ci_low = NA_real_, ci_high = NA_real_,
             p_value = NA_real_, p_holm = NA_real_,
             n_group1 = sum(sign(l6o_lodo307) == sign(l6o_eff307$estimate), na.rm = TRUE),
             n_group2 = length(l6o_lodo307),
             group1 = "same_direction", group2 = "omitted_donors",
             unit = "unique donor",
             sensitivity = "sensitivity_leave6", note = excl_note,
             stringsAsFactors = FALSE),
  data.frame(cohort = "GSE126805", model = "l6o_loko_direction_fraction",
             estimate = l6o_loko126_frac, ci_low = NA_real_, ci_high = NA_real_,
             p_value = NA_real_, p_holm = NA_real_,
             n_group1 = sum(sign(l6o_loko126) == sign(l6o_eff126$estimate), na.rm = TRUE),
             n_group2 = length(l6o_loko126),
             group1 = "same_direction", group2 = "omitted_allografts",
             unit = "paired kidney allograft",
             sensitivity = "sensitivity_leave6", note = excl_note,
             stringsAsFactors = FALSE)
)
write.csv(l6o_effects, file.path(res_dir, "consensus_module_l6o.csv"), row.names = FALSE)
msg("  WROTE results/extension/consensus_module_l6o.csv")

# =============================================================================
# 7. Connectivity: consensus signed score vs frozen six-gene score (sample-level)
# =============================================================================
msg("--- 6. Connectivity (consensus vs frozen six-gene score) ---")
conn <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  pr <- cor.test(x[ok], y[ok], method = "pearson")
  sp <- cor.test(x[ok], y[ok], method = "spearman")
  data.frame(n = sum(ok), pearson_r = pr$estimate, pearson_p = pr$p.value,
             spearman_rho = sp$estimate, spearman_p = sp$p.value,
             stringsAsFactors = FALSE)
}
conn139 <- conn(sc139, six139)
conn307 <- conn(sc307[primary_cols], six307)
conn126 <- conn(sc126, six126)
msg(sprintf("  GSE139061 Pearson r=%.3f (P=%.3g) | GSE30718 r=%.3f (P=%.3g) | GSE126805 r=%.3f (P=%.3g)",
            conn139$pearson_r, conn139$pearson_p, conn307$pearson_r,
            conn307$pearson_p, conn126$pearson_r, conn126$pearson_p))

# =============================================================================
# 8. Assemble scores.csv (long format)
# =============================================================================
msg("--- 7. Writing consensus_module_scores.csv ---")
scores_long <- rbind(
  data.frame(
    cohort = "GSE139061", sample = colnames(mat139),
    unit_id = colnames(mat139),
    group = ifelse(colnames(mat139) %in% aki_cols_139, "AKI", "Reference"),
    timepoint = NA_character_,
    signed_score = as.numeric(sc139), signed_score_percentile = as.numeric(sr139),
    six_gene_score = as.numeric(six139),
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE30718", sample = primary_cols,
    unit_id = meta_primary$donor_id,
    group = meta_primary$group,
    timepoint = NA_character_,
    signed_score = as.numeric(sc307[primary_cols]),
    signed_score_percentile = as.numeric(sr307[primary_cols]),
    six_gene_score = as.numeric(six307[primary_cols]),
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE126805", sample = meta126$sample,
    unit_id = meta126$kidney,
    group = meta126$timepoint,
    timepoint = meta126$timepoint,
    signed_score = as.numeric(sc126[meta126$sample]),
    signed_score_percentile = as.numeric(sr126[meta126$sample]),
    six_gene_score = as.numeric(six126[meta126$sample]),
    stringsAsFactors = FALSE
  )
)
write.csv(scores_long, file.path(res_dir, "consensus_module_scores.csv"),
          row.names = FALSE)

# =============================================================================
# 9. Effects table (tidy)
# =============================================================================
msg("--- 8. Writing consensus_module_effects.csv ---")
effects <- rbind(
  data.frame(
    cohort = "GSE139061", model = "main_signed_z",
    estimate = eff139$estimate, ci_low = eff139$ci_low, ci_high = eff139$ci_high,
    p_value = eff139$p_value, p_holm = holm_p["GSE139061"],
    n_group1 = eff139$n_grp1, n_group2 = eff139$n_grp2,
    group1 = "AKI", group2 = "Reference",
    unit = "sample (biopsy; sample == donor)",
    sensitivity = "primary", note = "expected direction: AKI higher",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE139061", model = "sensitivity_percentile",
    estimate = eff139_rank$estimate, ci_low = eff139_rank$ci_low,
    ci_high = eff139_rank$ci_high, p_value = eff139_rank$p_value, p_holm = NA_real_,
    n_group1 = eff139_rank$n_grp1, n_group2 = eff139_rank$n_grp2,
    group1 = "AKI", group2 = "Reference",
    unit = "sample (biopsy; sample == donor)",
    sensitivity = "sensitivity", note = "per-gene percentile signed score",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE30718", model = "main_signed_z",
    estimate = eff307$estimate, ci_low = eff307$ci_low, ci_high = eff307$ci_high,
    p_value = eff307$p_value, p_holm = holm_p["GSE30718"],
    n_group1 = eff307$n_grp1, n_group2 = eff307$n_grp2,
    group1 = "AKI", group2 = "Stable",
    unit = "unique donor (frozen primary set)",
    sensitivity = "primary", note = "expected direction: AKI higher",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE30718", model = "sensitivity_percentile",
    estimate = eff307_rank$estimate, ci_low = eff307_rank$ci_low,
    ci_high = eff307_rank$ci_high, p_value = eff307_rank$p_value, p_holm = NA_real_,
    n_group1 = eff307_rank$n_grp1, n_group2 = eff307_rank$n_grp2,
    group1 = "AKI", group2 = "Stable",
    unit = "unique donor (frozen primary set)",
    sensitivity = "sensitivity", note = "per-gene percentile signed score",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE126805", model = "main_signed_z",
    estimate = eff126$estimate, ci_low = eff126$ci_low, ci_high = eff126$ci_high,
    p_value = eff126$p_value, p_holm = holm_p["GSE126805"],
    n_group1 = eff126$n_grp1, n_group2 = eff126$n_grp2,
    group1 = "baseline_post", group2 = "baseline_pre",
    unit = "paired kidney allograft",
    sensitivity = "primary", note = "paired post-pre; expected direction: post higher",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE126805", model = "sensitivity_percentile",
    estimate = eff126_rank$estimate, ci_low = eff126_rank$ci_low,
    ci_high = eff126_rank$ci_high, p_value = eff126_rank$p_value, p_holm = NA_real_,
    n_group1 = eff126_rank$n_grp1, n_group2 = eff126_rank$n_grp2,
    group1 = "baseline_post", group2 = "baseline_pre",
    unit = "paired kidney allograft",
    sensitivity = "sensitivity", note = "per-gene percentile signed score",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE126805", model = "signflip_permutation",
    estimate = obs_d126, ci_low = NA_real_, ci_high = NA_real_,
    p_value = perm126, p_holm = NA_real_,
    n_group1 = n_pairs, n_group2 = n_pairs,
    group1 = "baseline_post", group2 = "baseline_pre",
    unit = "paired kidney allograft",
    sensitivity = "primary", note = sprintf("9999 sign flips, seed %d", SET_SEED),
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE139061", model = "lodo_direction_fraction",
    estimate = lodo139_frac, ci_low = NA_real_, ci_high = NA_real_,
    p_value = NA_real_, p_holm = NA_real_,
    n_group1 = sum(sign(lodo139) == sign(eff139$estimate), na.rm = TRUE),
    n_group2 = length(lodo139), group1 = "same_direction", group2 = "omitted_units",
    unit = "sample (sample == donor)",
    sensitivity = "stability", note = "leave-one-sample-out",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE30718", model = "lodo_direction_fraction",
    estimate = lodo307_frac, ci_low = NA_real_, ci_high = NA_real_,
    p_value = NA_real_, p_holm = NA_real_,
    n_group1 = sum(sign(lodo307) == sign(eff307$estimate), na.rm = TRUE),
    n_group2 = length(lodo307), group1 = "same_direction", group2 = "omitted_donors",
    unit = "unique donor",
    sensitivity = "stability", note = "leave-one-donor-out",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE126805", model = "loko_direction_fraction",
    estimate = loko126_frac, ci_low = NA_real_, ci_high = NA_real_,
    p_value = NA_real_, p_holm = NA_real_,
    n_group1 = sum(sign(loko126) == sign(eff126$estimate), na.rm = TRUE),
    n_group2 = length(loko126), group1 = "same_direction", group2 = "omitted_allografts",
    unit = "paired kidney allograft",
    sensitivity = "stability", note = "leave-one-kidney-out",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE139061", model = "connectivity_pearson",
    estimate = conn139$pearson_r, ci_low = NA_real_, ci_high = NA_real_,
    p_value = conn139$pearson_p, p_holm = NA_real_,
    n_group1 = conn139$n, n_group2 = NA_real_,
    group1 = "samples", group2 = NA_character_,
    unit = "sample", sensitivity = "connectivity",
    note = "consensus signed score vs frozen six-gene score",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE30718", model = "connectivity_pearson",
    estimate = conn307$pearson_r, ci_low = NA_real_, ci_high = NA_real_,
    p_value = conn307$pearson_p, p_holm = NA_real_,
    n_group1 = conn307$n, n_group2 = NA_real_,
    group1 = "samples", group2 = NA_character_,
    unit = "unique donor (primary set)", sensitivity = "connectivity",
    note = "consensus signed score vs frozen six-gene score",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE126805", model = "connectivity_pearson",
    estimate = conn126$pearson_r, ci_low = NA_real_, ci_high = NA_real_,
    p_value = conn126$pearson_p, p_holm = NA_real_,
    n_group1 = conn126$n, n_group2 = NA_real_,
    group1 = "samples", group2 = NA_character_,
    unit = "sample (all timepoints)", sensitivity = "connectivity",
    note = "consensus signed score vs frozen six-gene score",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE139061", model = "connectivity_spearman",
    estimate = conn139$spearman_rho, ci_low = NA_real_, ci_high = NA_real_,
    p_value = conn139$spearman_p, p_holm = NA_real_,
    n_group1 = conn139$n, n_group2 = NA_real_,
    group1 = "samples", group2 = NA_character_,
    unit = "sample", sensitivity = "connectivity",
    note = "consensus signed score vs frozen six-gene score",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE30718", model = "connectivity_spearman",
    estimate = conn307$spearman_rho, ci_low = NA_real_, ci_high = NA_real_,
    p_value = conn307$spearman_p, p_holm = NA_real_,
    n_group1 = conn307$n, n_group2 = NA_real_,
    group1 = "samples", group2 = NA_character_,
    unit = "unique donor (primary set)", sensitivity = "connectivity",
    note = "consensus signed score vs frozen six-gene score",
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "GSE126805", model = "connectivity_spearman",
    estimate = conn126$spearman_rho, ci_low = NA_real_, ci_high = NA_real_,
    p_value = conn126$spearman_p, p_holm = NA_real_,
    n_group1 = conn126$n, n_group2 = NA_real_,
    group1 = "samples", group2 = NA_character_,
    unit = "sample (all timepoints)", sensitivity = "connectivity",
    note = "consensus signed score vs frozen six-gene score",
    stringsAsFactors = FALSE
  )
)
write.csv(effects, file.path(res_dir, "consensus_module_effects.csv"), row.names = FALSE)

# =============================================================================
# 10. LODO/LOKO detail table
# =============================================================================
msg("--- 9. Writing consensus_module_lodo.csv ---")
lodo_long <- rbind(
  data.frame(cohort = "GSE139061", analysis = "lodo",
             omitted_unit = c(aki_cols_139, ref_cols_139),
             estimate = as.numeric(lodo139), stringsAsFactors = FALSE),
  data.frame(cohort = "GSE30718", analysis = "lodo",
             omitted_unit = all_donors307,
             estimate = as.numeric(lodo307), stringsAsFactors = FALSE),
  data.frame(cohort = "GSE126805", analysis = "loko",
             omitted_unit = kidneys_pair,
             estimate = as.numeric(loko126), stringsAsFactors = FALSE)
)
write.csv(lodo_long, file.path(res_dir, "consensus_module_lodo.csv"), row.names = FALSE)

# =============================================================================
# 11. Coverage table
# =============================================================================
msg("--- 10. Writing consensus_module_coverage.csv ---")
coverage <- rbind(cov139, cov307, cov126)
write.csv(coverage, file.path(res_dir, "consensus_module_coverage.csv"),
          row.names = FALSE)

# =============================================================================
# 12. Figure 3E: consensus module external validation
# =============================================================================
msg("--- 11. Figure 3E ---")
forest_rows <- data.frame(
  cohort = c("GSE139061", "GSE30718", "GSE126805"),
  label = c("AKI - Reference", "AKI - Stable", "post - pre (paired)"),
  estimate = c(eff139$estimate, eff307$estimate, eff126$estimate),
  ci_low = c(eff139$ci_low, eff307$ci_low, eff126$ci_low),
  ci_high = c(eff139$ci_high, eff307$ci_high, eff126$ci_high),
  unit = c("statistical unit: sample (biopsy; sample == donor)",
           "statistical unit: unique donor (frozen primary set)",
           "statistical unit: paired kidney allograft"),
  n_text = c(sprintf("n = 39 AKI, 9 Reference"), sprintf("n = 26 AKI, 11 Stable"),
             sprintf("n = 41 allografts (paired)")),
  p_text = c(sprintf("raw P = %.4f; Holm P = %.4f", raw_p["GSE139061"], holm_p["GSE139061"]),
             sprintf("raw P = %.4f; Holm P = %.4f", raw_p["GSE30718"], holm_p["GSE30718"]),
             sprintf("raw P = %.4f; Holm P = %.4f; perm P = %.4f",
                     raw_p["GSE126805"], holm_p["GSE126805"], perm126)),
  stringsAsFactors = FALSE
)
forest_rows$cohort <- factor(forest_rows$cohort,
                             levels = c("GSE139061", "GSE30718", "GSE126805"))

forest_panel <- function(dat) {
  xr <- range(c(dat$ci_low, dat$ci_high, 0))
  pad <- max(abs(xr)) * 0.12
  xlim <- c(xr[1] - pad, xr[2] + pad)
  ggplot(dat, aes(x = estimate, y = 1)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18,
                   colour = "grey25", linewidth = 0.7) +
    geom_point(size = 4.2, colour = AKI_COL) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.02, vjust = 2.0,
             label = paste0(dat$n_text, "\n", dat$p_text),
             size = 3.1, colour = "grey25") +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = c(0.35, 1.65)) +
    labs(title = dat$cohort, subtitle = dat$label,
         x = "Signed consensus score effect (AKI/post minus reference/pre)",
         y = NULL, caption = dat$unit) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
}

pA <- forest_panel(forest_rows[1, ])
pB <- forest_panel(forest_rows[2, ])
pC <- forest_panel(forest_rows[3, ])

# Paired trajectories for GSE126805 (per allograft)
traj_df <- data.frame(
  kidney = rep(kidneys_pair, 2),
  timepoint = factor(rep(c("baseline_pre", "baseline_post"), each = n_pairs),
                     levels = c("baseline_pre", "baseline_post")),
  score = c(pre_sc126, post_sc126),
  stringsAsFactors = FALSE
)
mean_pre <- mean(pre_sc126); mean_post <- mean(post_sc126)

pD <- ggplot(traj_df, aes(x = timepoint, y = score, group = kidney)) +
  geom_line(colour = "grey72", linewidth = 0.35, alpha = 0.7) +
  geom_point(size = 1.1, colour = "grey35") +
  annotate("segment", x = 1, xend = 2, y = mean_pre, yend = mean_post,
           colour = AKI_COL, linewidth = 1.1) +
  annotate("point", x = 1, y = mean_pre, colour = AKI_COL, size = 3) +
  annotate("point", x = 2, y = mean_post, colour = AKI_COL, size = 3) +
  scale_x_discrete(labels = c("baseline_pre" = "Before reperfusion\n(baseline_pre)",
                              "baseline_post" = "After reperfusion\n(baseline_post)")) +
  labs(title = "GSE126805", subtitle = "Paired score trajectories (per allograft)",
       y = "Signed consensus score", x = NULL) +
  annotate("text", x = 1.5, y = max(traj_df$score), hjust = 0.5, vjust = 1.2,
           label = sprintf("Mean paired change = %.3f (95%% CI %.3f to %.3f)",
                           eff126$estimate, eff126$ci_low, eff126$ci_high),
           size = 3.1, colour = "grey25") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

fig3e <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(
    title = "Figure 3E: 1,027-gene TAL-derived AKI consensus transcriptional response - external validation",
    subtitle = paste0(
      "Signed score = mean(z of 723 up genes) - mean(z of 304 down genes), per-gene z within cohort. ",
      "Expected direction: AKI/post higher. No patient merging; genes/cells are not replicates.\n",
      "Connectivity (sample-level, descriptive only): the consensus score and the frozen six-gene TAL score move in opposite directions,\n",
      "so the negative correlations in GSE30718/GSE126805 are direction-consistent; GSE139061 shows no correlation."
    ),
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.5, colour = "grey25"),
                  plot.tag = element_text(face = "bold", size = 11))
  )

ggsave(file.path(fig_dir, "figure3e_consensus_external_validation.png"), fig3e,
       width = 12, height = 7.2, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figure3e_consensus_external_validation.pdf"), fig3e,
       width = 12, height = 7.2, bg = "white")
msg("  WROTE figures/extension/figure3e_consensus_external_validation.png/.pdf")

# =============================================================================
# 13. Report
# =============================================================================
msg("--- 12. Writing consensus_module_report.md ---")
direction_status <- function(est) if (est > 0) "expected (positive)" else "negative (opposite of expected)"
dir139 <- direction_status(eff139$estimate)
dir307 <- direction_status(eff307$estimate)
dir126 <- direction_status(eff126$estimate)

rpt <- c(
  "# External validation of the 1,027-gene TAL-derived AKI consensus transcriptional module (GSE139061 / GSE30718 / GSE126805)",
  "",
  sprintf("**Definition:** %d-gene TAL-derived AKI consensus transcriptional response module, frozen from `results/extension/tal_pseudobulk_concordance.csv` (`replicated == TRUE`; %d up, %d down by the common direction of the two cohort logFCs; 0 discordant). No gene was filtered or re-directed based on external results.", nrow(mod), n_up, n_dn),
  "",
  "**Main score:** per-gene z-score within each cohort (frozen z-scoring space), then",
  "`signed score = mean(z of up genes) - mean(z of down genes)`; up/down are averaged",
  "separately so gene-count imbalance does not bias the score. Expected direction:",
  "signed score higher in AKI / after reperfusion (positive effect).",
  "**Sensitivity score:** same formula with per-gene percentile (rank/(n+1)).",
  sprintf("**Leave-six-genes-out (L6O):** the frozen six TAL genes present in the module (%s) are excluded and effects, 95%% CI, raw P, three-cohort Holm P and LODO/LOKO are recomputed with the same method and the same per-cohort z-spaces (independent `consensus_module_l6o.csv`).", paste(frozen_in_module, collapse = ", ")),
  "",
  "**Statistical units:** GSE139061 sample (biopsy; sample == donor), GSE30718 unique donor",
  "(frozen primary set), GSE126805 paired kidney allograft. Patients are not merged and",
  "genes/cells are not treated as replicates.",
  "",
  "## 1. Gene coverage per cohort",
  "",
  "| Cohort | Consensus | Present | Up present | Down present | Usable (sd>0) | Fraction present |",
  "|---|---|---|---|---|---|---|",
  sprintf("| GSE139061 | %d | %d | %d | %d | %d | %.3f |", cov139$n_consensus, cov139$n_present, cov139$n_up_present, cov139$n_down_present, cov139$n_usable, cov139$frac_present),
  sprintf("| GSE30718 | %d | %d | %d | %d | %d | %.3f |", cov307$n_consensus, cov307$n_present, cov307$n_up_present, cov307$n_down_present, cov307$n_usable, cov307$frac_present),
  sprintf("| GSE126805 | %d | %d | %d | %d | %d | %.3f |", cov126$n_consensus, cov126$n_present, cov126$n_up_present, cov126$n_down_present, cov126$n_usable, cov126$frac_present),
  "",
  "## 2. Main effects (unit-level, two-sided)",
  "",
  "| Cohort | Model | Effect (95% CI) | Raw P | Holm P | Direction |",
  "|---|---|---|---|---|---|",
  sprintf("| GSE139061 | AKI - Reference | %.3f (%.3f to %.3f) | %.4f | %.4f | %s |", eff139$estimate, eff139$ci_low, eff139$ci_high, eff139$p_value, holm_p["GSE139061"], dir139),
  sprintf("| GSE30718 | AKI - Stable | %.3f (%.3f to %.3f) | %.4f | %.4f | %s |", eff307$estimate, eff307$ci_low, eff307$ci_high, eff307$p_value, holm_p["GSE30718"], dir307),
  sprintf("| GSE126805 | post - pre (paired, n=%d) | %.3f (%.3f to %.3f) | %.4f | %.4f | %s |", n_pairs, eff126$estimate, eff126$ci_low, eff126$ci_high, eff126$p_value, holm_p["GSE126805"], dir126),
  "",
  "Holm correction applied across the three main P values.",
  "",
  "## 3. Sensitivity (per-gene percentile signed score)",
  "",
  sprintf("- GSE139061: %.3f (%.3f to %.3f), P = %.4f", eff139_rank$estimate, eff139_rank$ci_low, eff139_rank$ci_high, eff139_rank$p_value),
  sprintf("- GSE30718: %.3f (%.3f to %.3f), P = %.4f", eff307_rank$estimate, eff307_rank$ci_low, eff307_rank$ci_high, eff307_rank$p_value),
  sprintf("- GSE126805: %.3f (%.3f to %.3f), P = %.4f", eff126_rank$estimate, eff126_rank$ci_low, eff126_rank$ci_high, eff126_rank$p_value),
  "",
  "## 4. GSE126805 paired permutation and LOKO",
  "",
  sprintf("- Paired sign-flip permutation (%d flips, seed %d, two-sided): P = %.4f", N_PERM, SET_SEED, perm126),
  sprintf("- Leave-one-kidney-out: same-direction fraction %.3f (%d/%d)", loko126_frac, sum(sign(loko126) == sign(eff126$estimate), na.rm = TRUE), length(loko126)),
  "",
  "## 5. Leave-six-genes-out sensitivity (excludes module-present frozen six genes)",
  "",
  sprintf("- Excluded genes (present in the module): %s (all down).", paste(frozen_in_module, collapse = ", ")),
  sprintf("- GSE139061 AKI - Reference: %.3f (%.3f to %.3f), raw P = %.4f, Holm P = %.4f", l6o_eff139$estimate, l6o_eff139$ci_low, l6o_eff139$ci_high, l6o_eff139$p_value, l6o_holm_p["GSE139061"]),
  sprintf("- GSE30718 AKI - Stable: %.3f (%.3f to %.3f), raw P = %.4f, Holm P = %.4f", l6o_eff307$estimate, l6o_eff307$ci_low, l6o_eff307$ci_high, l6o_eff307$p_value, l6o_holm_p["GSE30718"]),
  sprintf("- GSE126805 post - pre (paired, n=%d): %.3f (%.3f to %.3f), raw P = %.4f, Holm P = %.4f; sign-flip permutation P = %.4f", n_pairs, l6o_eff126$estimate, l6o_eff126$ci_low, l6o_eff126$ci_high, l6o_eff126$p_value, l6o_holm_p["GSE126805"], l6o_perm126),
  sprintf("- LODO/LOKO same-direction fractions: GSE139061 %.3f (%d/%d), GSE30718 %.3f (%d/%d), GSE126805 %.3f (%d/%d)", l6o_lodo139_frac, sum(sign(l6o_lodo139) == sign(l6o_eff139$estimate), na.rm = TRUE), length(l6o_lodo139), l6o_lodo307_frac, sum(sign(l6o_lodo307) == sign(l6o_eff307$estimate), na.rm = TRUE), length(l6o_lodo307), l6o_loko126_frac, sum(sign(l6o_loko126) == sign(l6o_eff126$estimate), na.rm = TRUE), length(l6o_loko126)),
  "Holm correction applied across the three L6O P values.",
  "",
  "## 6. LODO/LOKO direction stability (main score)",
  "",
  "| Cohort | Analysis | Same-direction fraction | Omitted units |",
  "|---|---|---|---|",
  sprintf("| GSE139061 | leave-one-sample-out (sample == donor) | %.3f (%d/%d) | %d |", lodo139_frac, sum(sign(lodo139) == sign(eff139$estimate), na.rm = TRUE), length(lodo139), length(lodo139)),
  sprintf("| GSE30718 | leave-one-donor-out | %.3f (%d/%d) | %d |", lodo307_frac, sum(sign(lodo307) == sign(eff307$estimate), na.rm = TRUE), length(lodo307), length(lodo307)),
  sprintf("| GSE126805 | leave-one-kidney-out | %.3f (%d/%d) | %d |", loko126_frac, sum(sign(loko126) == sign(eff126$estimate), na.rm = TRUE), length(loko126), length(loko126)),
  "",
  "## 7. Connectivity with the frozen six-gene TAL score (sample-level, descriptive only)",
  "",
  "| Cohort | n samples | Pearson r (P) | Spearman rho (P) |",
  "|---|---|---|---|",
  sprintf("| GSE139061 | %d | %.3f (%.3g) | %.3f (%.3g) |", conn139$n, conn139$pearson_r, conn139$pearson_p, conn139$spearman_rho, conn139$spearman_p),
  sprintf("| GSE30718 | %d | %.3f (%.3g) | %.3f (%.3g) |", conn307$n, conn307$pearson_r, conn307$pearson_p, conn307$spearman_rho, conn307$spearman_p),
  sprintf("| GSE126805 | %d | %.3f (%.3g) | %.3f (%.3g) |", conn126$n, conn126$pearson_r, conn126$pearson_p, conn126$spearman_rho, conn126$spearman_p),
  "",
  "These correlations are connectivity checks only; they are not validation tests.",
  "The consensus signed score and the frozen six-gene TAL identity score are constructed in",
  "opposite directions (the consensus score rises in AKI/post while the frozen identity score",
  "falls), so the negative sample-level correlations in GSE30718 and GSE126805 are",
  sprintf("direction-consistent with the two scores moving oppositely (GSE30718 r = %.3f, GSE126805 r = %.3f).", conn307$pearson_r, conn126$pearson_r),
  sprintf("GSE139061 shows no sample-level correlation (r = %.3f).", conn139$pearson_r),
  "",
  "## 8. Conclusion boundaries",
  "",
  "- Bulk tissue analyses measure **tissue-level detectability** of the TAL-derived AKI",
  "  consensus transcriptional response; they do not establish TAL-specific, cell-intrinsic,",
  "  mechanistic, causal, or biomarker claims.",
  "- Negative or inconsistent results are preserved as-is and not promoted.",
  "- No patient merging; the statistical unit is sample/donor (GSE139061), unique",
  "  donor (GSE30718), or paired allograft (GSE126805); genes/cells are never",
  "  treated as replicates.",
  "- All three cohorts come from the frozen pipelines (`01_bulk_validation.R`,",
  "  `03_longitudinal_validation.R`); no data were re-downloaded.",
  "",
  "## 9. Files",
  "",
  "- `results/extension/consensus_module_definition.csv` (1,027 genes, frozen direction)",
  "- `results/extension/consensus_module_coverage.csv`",
  "- `results/extension/consensus_module_scores.csv` (sample-level scores, all cohorts)",
  "- `results/extension/consensus_module_effects.csv` (effects, CI, raw/Holm P, stability, connectivity)",
  "- `results/extension/consensus_module_l6o.csv` (leave-six-genes-out sensitivity)",
  "- `results/extension/consensus_module_lodo.csv` (LODO/LOKO detail)",
  "- `figures/extension/figure3e_consensus_external_validation.png/.pdf`",
  ""
)
writeLines(rpt, file.path(res_dir, "consensus_module_report.md"))
cat(paste(rpt, collapse = "\n"), "\n")
msg("DONE")
