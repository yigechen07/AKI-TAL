#!/usr/bin/env Rscript
# Submission sensitivity analyses for the frozen AKI-TAL evidence package.
#
# This script does not overwrite the frozen primary outputs. It addresses two
# submission-facing technical questions:
#   1. Figure 3: GSE183276 pseudobulk values are SoupX-corrected and fractional.
#      Compare the frozen edgeR result with rounded-pseudobulk edgeR and a
#      limma-trend analysis of TMM logCPM.
#   2. GSE126805: recompute the frozen six-gene and 1,027-gene signed scores
#      after z-standardising only the 41 paired pre/post allografts; replace
#      repeated-sample connectivity P values with allograft-level summaries.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) == 1L) {
  project_dir <- dirname(dirname(dirname(normalizePath(sub("^--file=", "", script_arg)))))
} else {
  project_dir <- getwd()
}
setwd(project_dir)

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
})

set.seed(20260818)
res_dir <- file.path(project_dir, "results", "submission")
fig_dir <- file.path(project_dir, "figures", "supplement")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pb_path <- file.path(project_dir, "..", "AKI_DCT_analysis", "results", "stage3b_pseudobulk.rds")
qc_path <- file.path(project_dir, "..", "AKI_DCT_analysis", "results", "stage3b_qc.csv")
gse126_dir <- file.path(project_dir, "data", "external", "GSE126805")
module_path <- file.path(project_dir, "results", "extension", "consensus_module_definition.csv")
frozen_conc_path <- file.path(project_dir, "results", "extension", "tal_pseudobulk_concordance.csv")

stopifnot(file.exists(pb_path), file.exists(qc_path), file.exists(module_path),
          file.exists(frozen_conc_path))

pb <- readRDS(pb_path)
qc <- read.csv(qc_path, stringsAsFactors = FALSE)
module <- read.csv(module_path, stringsAsFactors = FALSE)
frozen_conc <- read.csv(frozen_conc_path, stringsAsFactors = FALSE)
frozen_1027 <- frozen_conc$gene[frozen_conc$replicated]
stopifnot(length(frozen_1027) == 1027L)

primary_data <- function(dataset) {
  meta <- pb$metadata[[dataset]]
  q <- qc[qc$dataset == dataset & qc$segment == "TAL", ]
  q <- q[match(meta$donor, q$donor), ]
  stopifnot(all(q$donor == meta$donor))
  if (dataset == "GSE183276") {
    keep <- q$eligible_50 & meta$diagnosis %in% c("AKI", "Ref")
    control <- "Ref"
  } else {
    keep <- q$eligible_50 & meta$diagnosis %in% c("AKI", "non_AKI")
    control <- "non_AKI"
  }
  donors <- meta$donor[keep]
  info <- data.frame(
    donor = donors,
    group = factor(ifelse(meta$diagnosis[keep] == "AKI", "AKI", control),
                   levels = c(control, "AKI")),
    stringsAsFactors = FALSE
  )
  counts <- as.matrix(pb$pseudobulk[[dataset]]$TAL[, donors, drop = FALSE])
  list(counts = counts, info = info)
}

run_edger <- function(counts, info) {
  y <- DGEList(counts = counts, samples = info)
  keep <- filterByExpr(y, group = info$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~ group, data = info)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  tab <- topTags(glmQLFTest(fit, coef = "groupAKI"), n = Inf,
                 sort.by = "none")$table
  data.frame(gene = rownames(tab), logFC = tab$logFC, PValue = tab$PValue,
             FDR = tab$FDR, stringsAsFactors = FALSE)
}

run_limma_trend <- function(counts, info) {
  y <- DGEList(counts = counts, samples = info)
  keep <- filterByExpr(y, group = info$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  lcpm <- cpm(y, log = TRUE, prior.count = 2)
  design <- model.matrix(~ group, data = info)
  fit <- eBayes(lmFit(lcpm, design), trend = TRUE, robust = TRUE)
  tab <- topTable(fit, coef = "groupAKI", n = Inf, sort.by = "none")
  data.frame(gene = rownames(tab), logFC = tab$logFC, PValue = tab$P.Value,
             FDR = tab$adj.P.Val, stringsAsFactors = FALSE)
}

merge_method <- function(base, extra, prefix) {
  names(extra)[names(extra) != "gene"] <- paste0(prefix, names(extra)[names(extra) != "gene"])
  merge(base, extra, by = "gene", all = FALSE, sort = FALSE)
}

# -------------------------------------------------------------------------
# Figure 3 count-assumption sensitivity
# -------------------------------------------------------------------------
d183 <- primary_data("GSE183276")
d210 <- primary_data("GSE210622")
nonzero <- d183$counts[d183$counts != 0]
fractional_n <- sum(abs(nonzero - round(nonzero)) > 1e-8)
fractional_fraction <- fractional_n / length(nonzero)

edger_fractional <- run_edger(d183$counts, d183$info)
edger_rounded <- run_edger(round(d183$counts), d183$info)
limma_fractional <- run_limma_trend(d183$counts, d183$info)
edger_210 <- run_edger(d210$counts, d210$info)

gene_sens <- edger_fractional
names(gene_sens)[-1] <- paste0("fractional_edger_", names(gene_sens)[-1])
gene_sens <- merge_method(gene_sens, edger_rounded, "rounded_edger_")
gene_sens <- merge_method(gene_sens, limma_fractional, "fractional_limma_")
gene_sens <- merge_method(gene_sens, edger_210, "GSE210622_edger_")
gene_sens$in_frozen_1027 <- gene_sens$gene %in% frozen_1027
gene_sens$rounded_same_direction <-
  sign(gene_sens$fractional_edger_logFC) == sign(gene_sens$rounded_edger_logFC)
gene_sens$limma_same_direction <-
  sign(gene_sens$fractional_edger_logFC) == sign(gene_sens$fractional_limma_logFC)
gene_sens$rounded_strict_cross_cohort <- gene_sens$rounded_same_direction &
  sign(gene_sens$rounded_edger_logFC) == sign(gene_sens$GSE210622_edger_logFC) &
  gene_sens$rounded_edger_FDR < 0.05 & gene_sens$GSE210622_edger_FDR < 0.05
gene_sens$limma_strict_cross_cohort <- gene_sens$limma_same_direction &
  sign(gene_sens$fractional_limma_logFC) == sign(gene_sens$GSE210622_edger_logFC) &
  gene_sens$fractional_limma_FDR < 0.05 & gene_sens$GSE210622_edger_FDR < 0.05
write.csv(gene_sens, file.path(res_dir, "figure3_count_assumption_gene_sensitivity.csv"),
          row.names = FALSE)

in1027 <- gene_sens$in_frozen_1027
fig3_summary <- data.frame(
  metric = c(
    "GSE183276_nonzero_pseudobulk_values",
    "GSE183276_fractional_nonzero_values",
    "GSE183276_fractional_nonzero_fraction",
    "all_tested_rounded_logFC_sign_concordance",
    "all_tested_limma_logFC_sign_concordance",
    "frozen_1027_rounded_logFC_sign_retention",
    "frozen_1027_limma_logFC_sign_retention",
    "rounded_edgeR_strict_cross_cohort_genes",
    "fractional_limma_strict_cross_cohort_genes"
  ),
  value = c(
    length(nonzero), fractional_n, fractional_fraction,
    mean(gene_sens$rounded_same_direction),
    mean(gene_sens$limma_same_direction),
    mean(gene_sens$rounded_same_direction[in1027]),
    mean(gene_sens$limma_same_direction[in1027]),
    sum(gene_sens$rounded_strict_cross_cohort),
    sum(gene_sens$limma_strict_cross_cohort)
  ),
  stringsAsFactors = FALSE
)
write.csv(fig3_summary, file.path(res_dir, "figure3_count_assumption_summary.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# GSE126805 paired-only z-standardisation and allograft-level connectivity
# -------------------------------------------------------------------------
count_candidates <- c(
  file.path(project_dir, "data", "gse126805", "GSE126805_all.gene_counts.txt.gz"),
  file.path(gse126_dir, "GSE126805_all.gene_counts.txt.gz"),
  file.path(project_dir, "data", "GSE126805_all.gene_counts.txt.gz")
)
count_path <- count_candidates[file.exists(count_candidates)][1]
if (is.na(count_path)) stop("GSE126805 count matrix not found")

cts <- read.delim(gzfile(count_path), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("external_gene_name" %in% names(cts))
sample_cols <- grep("^sample\\.kidney_", names(cts), value = TRUE)
cnt126 <- as.matrix(cts[, sample_cols, drop = FALSE])
storage.mode(cnt126) <- "numeric"
gene_names <- as.character(cts$external_gene_name)
gene_mean <- rowMeans(cnt126, na.rm = TRUE)
ord <- order(gene_names, -gene_mean)
cnt126 <- cnt126[ord, , drop = FALSE]
gene_names <- gene_names[ord]
keep_gene <- !duplicated(gene_names) & !is.na(gene_names) & nzchar(gene_names)
cnt126 <- cnt126[keep_gene, , drop = FALSE]
rownames(cnt126) <- gene_names[keep_gene]

parse_sample <- function(s) {
  parts <- strsplit(sub("^sample\\.", "", s), "_")[[1]]
  c(kidney = parts[2], timepoint = paste(parts[3:length(parts)], collapse = "_"))
}
meta126 <- as.data.frame(do.call(rbind, lapply(sample_cols, function(s) {
  p <- parse_sample(s)
  c(sample = s, kidney = unname(p["kidney"]), timepoint = unname(p["timepoint"]))
})), stringsAsFactors = FALSE)

has_pair <- sapply(unique(meta126$kidney), function(k)
  all(c("baseline_pre", "baseline_post") %in% meta126$timepoint[meta126$kidney == k]))
paired_kidneys <- names(has_pair)[has_pair]
paired_rows <- meta126$kidney %in% paired_kidneys &
  meta126$timepoint %in% c("baseline_pre", "baseline_post")
paired_samples <- meta126$sample[paired_rows]
stopifnot(length(paired_kidneys) == 41L, length(paired_samples) == 82L)

y126 <- calcNormFactors(DGEList(counts = cnt126), method = "TMM")
lcpm126 <- cpm(y126, log = TRUE, prior.count = 1)

zscore_rows <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  s <- apply(m, 1, sd, na.rm = TRUE)
  s[!is.finite(s) | s < 1e-8] <- 1
  sweep(sweep(m, 1, mu, "-"), 1, s, "/")
}

score_set <- function(reference_samples) {
  z <- zscore_rows(lcpm126[, reference_samples, drop = FALSE])
  frozen <- intersect(c("SLC12A1", "UMOD", "CLDN16", "KCNJ10", "PTH1R", "CLDN10"), rownames(z))
  up <- intersect(module$gene[module$direction == "up"], rownames(z))
  down <- intersect(module$gene[module$direction == "down"], rownames(z))
  list(
    six = colMeans(z[frozen, , drop = FALSE]),
    consensus = colMeans(z[up, , drop = FALSE]) - colMeans(z[down, , drop = FALSE]),
    n_frozen = length(frozen), n_up = length(up), n_down = length(down)
  )
}

scores_all <- score_set(sample_cols)
scores_pair <- score_set(paired_samples)

paired_vectors <- function(scores, samples_available) {
  m <- meta126[meta126$sample %in% samples_available & meta126$kidney %in% paired_kidneys &
                 meta126$timepoint %in% c("baseline_pre", "baseline_post"), ]
  get <- function(k, tp, v) {
    s <- m$sample[m$kidney == k & m$timepoint == tp][1]
    unname(v[s])
  }
  six_pre <- sapply(paired_kidneys, get, "baseline_pre", scores$six)
  six_post <- sapply(paired_kidneys, get, "baseline_post", scores$six)
  con_pre <- sapply(paired_kidneys, get, "baseline_pre", scores$consensus)
  con_post <- sapply(paired_kidneys, get, "baseline_post", scores$consensus)
  data.frame(
    kidney = paired_kidneys,
    six_pre = six_pre, six_post = six_post, six_delta = six_post - six_pre,
    consensus_pre = con_pre, consensus_post = con_post,
    consensus_delta = con_post - con_pre,
    stringsAsFactors = FALSE
  )
}

pair_all <- paired_vectors(scores_all, sample_cols)
pair_only <- paired_vectors(scores_pair, paired_samples)

effect_row <- function(delta, endpoint, standardisation) {
  tt <- t.test(delta)
  nperm <- 99999L
  signs <- matrix(sample(c(-1, 1), length(delta) * nperm, replace = TRUE),
                  nrow = nperm)
  perm <- rowMeans(signs * matrix(delta, nrow = nperm, ncol = length(delta), byrow = TRUE))
  pperm <- (1 + sum(abs(perm) >= abs(mean(delta)))) / (nperm + 1)
  loko <- sapply(seq_along(delta), function(i) mean(delta[-i]))
  data.frame(
    endpoint = endpoint, standardisation = standardisation,
    n_allografts = length(delta), estimate = mean(delta),
    ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], p_t = tt$p.value,
    p_signflip_monte_carlo = pperm,
    loko_same_direction_fraction = mean(sign(loko) == sign(mean(delta))),
    stringsAsFactors = FALSE
  )
}

paired_effects <- rbind(
  effect_row(pair_all$six_delta, "frozen_six_gene", "all_163_samples"),
  effect_row(pair_only$six_delta, "frozen_six_gene", "paired_pre_post_82_samples"),
  effect_row(pair_all$consensus_delta, "signed_1027_gene", "all_163_samples"),
  effect_row(pair_only$consensus_delta, "signed_1027_gene", "paired_pre_post_82_samples")
)
write.csv(paired_effects, file.path(res_dir, "GSE126805_prepost_only_effects.csv"),
          row.names = FALSE)
write.csv(pair_only, file.path(res_dir, "GSE126805_prepost_only_allograft_scores.csv"),
          row.names = FALSE)

connectivity <- data.frame(
  analysis = c("paired_pre_post_mean", "acute_post_minus_pre_delta"),
  n_allografts = c(nrow(pair_only), nrow(pair_only)),
  spearman_rho = c(
    cor((pair_only$six_pre + pair_only$six_post) / 2,
        (pair_only$consensus_pre + pair_only$consensus_post) / 2,
        method = "spearman"),
    cor(pair_only$six_delta, pair_only$consensus_delta, method = "spearman")
  ),
  interpretation = c(
    "descriptive between-allograft association; one value per allograft",
    "descriptive association between acute paired changes; one value per allograft"
  ),
  stringsAsFactors = FALSE
)
write.csv(connectivity, file.path(res_dir, "GSE126805_allograft_level_connectivity.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# Compact supplementary figure
# -------------------------------------------------------------------------
png(file.path(fig_dir, "figureS1_submission_sensitivity.png"), width = 2400,
    height = 1200, res = 220)
par(mfrow = c(1, 2), mar = c(4.6, 4.8, 2.8, 1.2), las = 1)
idx <- gene_sens$in_frozen_1027
plot(gene_sens$fractional_edger_logFC[idx], gene_sens$rounded_edger_logFC[idx],
     pch = 16, cex = 0.45, col = grDevices::adjustcolor("#2166AC", alpha.f = 0.45),
     xlab = "Fractional SoupX pseudobulk: edgeR logFC",
     ylab = "Rounded donor pseudobulk: edgeR logFC",
     main = "A  Figure 3 count-assumption sensitivity")
abline(0, 1, col = "#B2182B", lwd = 1.5)
legend("topleft", legend = sprintf("Frozen 1,027 genes\nSign retained: %.1f%%",
       100 * mean(gene_sens$rounded_same_direction[idx])), bty = "n", cex = 0.85)

pe <- paired_effects
xx <- c(0.9, 1.1, 1.9, 2.1)
cols <- c("#777777", "#2166AC", "#777777", "#2166AC")
plot(NA, xlim = c(0.5, 2.5), ylim = range(c(pe$ci_low, pe$ci_high)), xaxt = "n",
     xlab = "", ylab = "Post minus pre paired effect (95% CI)",
     main = "B  GSE126805 standardisation sensitivity")
abline(h = 0, lty = 2, col = "#999999")
segments(xx, pe$ci_low, xx, pe$ci_high, col = cols, lwd = 2)
points(xx, pe$estimate, pch = 16, col = cols, cex = 1.1)
axis(1, at = c(1, 2), labels = c("Frozen six-gene", "Signed 1,027-gene"))
legend("topleft", legend = c("All 163 samples", "Paired pre/post only"),
       col = c("#777777", "#2166AC"), pch = 16, bty = "n", cex = 0.85)
dev.off()

pdf(file.path(fig_dir, "figureS1_submission_sensitivity.pdf"), width = 10.8, height = 5.3)
par(mfrow = c(1, 2), mar = c(4.6, 4.8, 2.8, 1.2), las = 1)
plot(gene_sens$fractional_edger_logFC[idx], gene_sens$rounded_edger_logFC[idx],
     pch = 16, cex = 0.45, col = grDevices::adjustcolor("#2166AC", alpha.f = 0.45),
     xlab = "Fractional SoupX pseudobulk: edgeR logFC",
     ylab = "Rounded donor pseudobulk: edgeR logFC",
     main = "A  Figure 3 count-assumption sensitivity")
abline(0, 1, col = "#B2182B", lwd = 1.5)
legend("topleft", legend = sprintf("Frozen 1,027 genes\nSign retained: %.1f%%",
       100 * mean(gene_sens$rounded_same_direction[idx])), bty = "n", cex = 0.85)
plot(NA, xlim = c(0.5, 2.5), ylim = range(c(pe$ci_low, pe$ci_high)), xaxt = "n",
     xlab = "", ylab = "Post minus pre paired effect (95% CI)",
     main = "B  GSE126805 standardisation sensitivity")
abline(h = 0, lty = 2, col = "#999999")
segments(xx, pe$ci_low, xx, pe$ci_high, col = cols, lwd = 2)
points(xx, pe$estimate, pch = 16, col = cols, cex = 1.1)
axis(1, at = c(1, 2), labels = c("Frozen six-gene", "Signed 1,027-gene"))
legend("topleft", legend = c("All 163 samples", "Paired pre/post only"),
       col = c("#777777", "#2166AC"), pch = 16, bty = "n", cex = 0.85)
dev.off()

report <- c(
  "# Submission sensitivity audit",
  "",
  "## Figure 3 count-assumption sensitivity",
  "",
  sprintf("- GSE183276 donor pseudobulk contains %d fractional non-zero values (%.1f%% of non-zero values), reflecting author-provided SoupX correction.",
          fractional_n, 100 * fractional_fraction),
  sprintf("- Rounded-pseudobulk edgeR retained the logFC direction for %.1f%% of the frozen 1,027 genes and yielded %d strict cross-cohort genes.",
          100 * mean(gene_sens$rounded_same_direction[in1027]),
          sum(gene_sens$rounded_strict_cross_cohort)),
  sprintf("- Fractional-logCPM limma-trend retained the logFC direction for %.1f%% of the frozen 1,027 genes and yielded %d strict cross-cohort genes.",
          100 * mean(gene_sens$limma_same_direction[in1027]),
          sum(gene_sens$limma_strict_cross_cohort)),
  "- These are method sensitivities for exploratory Figure 3; they do not redefine the frozen six-gene primary endpoint.",
  "",
  "## GSE126805 paired-only standardisation",
  "",
  sprintf("- Frozen six-gene score, all-sample z reference: %.3f (95%% CI %.3f to %.3f), sign-flip P=%.4f.",
          pe$estimate[1], pe$ci_low[1], pe$ci_high[1], pe$p_signflip_monte_carlo[1]),
  sprintf("- Frozen six-gene score, paired-pre/post-only z reference: %.3f (95%% CI %.3f to %.3f), sign-flip P=%.4f.",
          pe$estimate[2], pe$ci_low[2], pe$ci_high[2], pe$p_signflip_monte_carlo[2]),
  sprintf("- Signed 1,027-gene score, all-sample z reference: %.3f (95%% CI %.3f to %.3f), sign-flip P=%.4f.",
          pe$estimate[3], pe$ci_low[3], pe$ci_high[3], pe$p_signflip_monte_carlo[3]),
  sprintf("- Signed 1,027-gene score, paired-pre/post-only z reference: %.3f (95%% CI %.3f to %.3f), sign-flip P=%.4f.",
          pe$estimate[4], pe$ci_low[4], pe$ci_high[4], pe$p_signflip_monte_carlo[4]),
  "- Connectivity is now summarised once per allograft; no P value from 163 repeated samples is used.",
  "",
  "## Interpretation boundary",
  "",
  "These analyses assess technical robustness only. They do not establish TAL selectivity, cell-intrinsic regulation, PT independence, causality, or a mitochondrial mechanism."
)
writeLines(report, file.path(res_dir, "submission_sensitivity_report.md"))

writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo_submission_sensitivity.txt"))
cat("Submission sensitivity analysis complete.\n")
