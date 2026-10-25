#!/usr/bin/env Rscript
# =============================================================================
# figure3_tal_pseudobulk.R — Figure 3 (extension): TAL donor-pseudobulk DEG
# =============================================================================
# Scope (analysis_extension_plan.md WP3):
#   * Sum raw TAL counts within donor (stage3b_pseudobulk.rds) and model each
#     cohort SEPARATELY with edgeR quasi-likelihood (template:
#     ../AKI_DCT_analysis/code/phase1_analysis.R and
#     stage2_gse210622_validation.R adapted to TAL).
#   * Primary tier = donors with >=50 TAL cells (eligible_50) and diagnosis in
#     AKI vs Reference (GSE183276) / AKI vs non-AKI (GSE210622); sensitivity
#     tier = >=20 TAL cells.
#   * Filter genes by donor-level expression support (filterByExpr), not by
#     final P value. Report effect size, FDR, direction agreement and
#     leave-one-donor-out stability.
#   * Strict replication requires the same direction and FDR < 0.05 in both
#     cohorts. Same-direction significance in only one cohort is reported as
#     single-cohort support, not independent replication.
#   * Pathway enrichment panels are DEFERRED (no local gene-set resource).
#
# Outputs (results/extension/):
#   tal_pseudobulk_deg_GSE183276.csv / tal_pseudobulk_deg_GSE210622.csv
#   tal_pseudobulk_concordance.csv, tal_pseudobulk_summary.csv
#   tal_pseudobulk_lodo.csv
#   figures/extension/figure3_tal_pseudobulk.png/.pdf
# =============================================================================

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
  library(ggplot2)
  library(patchwork)
})

pb_path <- file.path(project_dir, "..", "AKI_DCT_analysis", "results", "stage3b_pseudobulk.rds")
qc_path <- file.path(project_dir, "..", "AKI_DCT_analysis", "results", "stage3b_qc.csv")
res_dir <- file.path(project_dir, "results", "extension")
fig_dir <- file.path(project_dir, "figures", "extension")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

pb <- readRDS(pb_path)
qc <- read.csv(qc_path, stringsAsFactors = FALSE)

FROZEN <- c("SLC12A1","UMOD","CLDN16","KCNJ10","PTH1R","CLDN10")

# ---------------------------------------------------------------------------
# edgeR wrapper: cohort-specific AKI vs control on TAL donor pseudobulk
# ---------------------------------------------------------------------------
run_tal_edger <- function(dataset, tier) {
  meta <- pb$metadata[[dataset]]
  qc_ds <- qc[qc$dataset == dataset & qc$segment == "TAL", ]
  qc_ds <- qc_ds[match(meta$donor, qc_ds$donor), ]
  stopifnot(all(qc_ds$donor == meta$donor))

  if (tier == "primary_ge50") keep_tier <- qc_ds$eligible_50
  else keep_tier <- qc_ds$eligible_20

  if (dataset == "GSE183276") {
    keep_dx <- meta$diagnosis %in% c("AKI", "Ref")
    ctrl_name <- "Ref"
  } else {
    keep_dx <- meta$diagnosis %in% c("AKI", "non_AKI")
    ctrl_name <- "non_AKI"
  }
  keep <- keep_tier & keep_dx
  donors <- meta$donor[keep]
  info <- data.frame(
    donor = donors,
    group = factor(ifelse(meta$diagnosis[keep] == "AKI", "AKI", ctrl_name),
                   levels = c(ctrl_name, "AKI")),
    stringsAsFactors = FALSE
  )
  counts <- as.matrix(pb$pseudobulk[[dataset]]$TAL[, donors, drop = FALSE])
  rownames(counts) <- rownames(pb$pseudobulk[[dataset]]$TAL)

  y <- DGEList(counts = counts, samples = info)
  keep_genes <- filterByExpr(y, group = info$group)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~ group, data = info)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "groupAKI")
  tt <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tt$gene <- rownames(tt)
  tt$dataset <- dataset
  tt$tier <- tier
  tt$n_AKI <- sum(info$group == "AKI")
  tt$n_ctrl <- sum(info$group != "AKI")
  tt$direction <- ifelse(tt$logFC < 0, "down", "up")
  tt <- tt[, c("gene","logFC","logCPM","F","PValue","FDR","direction",
               "dataset","tier","n_AKI","n_ctrl")]
  rownames(tt) <- NULL

  # LODO: refit with each donor removed, record sign(logFC) per gene
  lodo_sign <- matrix(NA_real_, nrow = length(donors), ncol = nrow(tt),
                      dimnames = list(donors, tt$gene))
  for (d in donors) {
    info_loo <- info[info$donor != d, ]
    counts_loo <- counts[, info_loo$donor, drop = FALSE]
    yl <- DGEList(counts = counts_loo, samples = info_loo)
    yl <- yl[keep_genes, , keep.lib.sizes = FALSE]
    yl <- calcNormFactors(yl, method = "TMM")
    design_loo <- model.matrix(~ group, data = info_loo)
    yl <- estimateDisp(yl, design_loo, robust = TRUE)
    fit_loo <- glmQLFit(yl, design_loo, robust = TRUE)
    qlf_loo <- glmQLFTest(fit_loo, coef = "groupAKI")
    tt_loo <- topTags(qlf_loo, n = Inf, sort.by = "PValue")$table
    lodo_sign[d, intersect(colnames(lodo_sign), rownames(tt_loo))] <-
      sign(tt_loo[intersect(colnames(lodo_sign), rownames(tt_loo)), "logFC"])
  }
  full_sign <- sign(tt$logFC); names(full_sign) <- tt$gene
  lodo_consist <- colMeans(lodo_sign == matrix(full_sign[colnames(lodo_sign)],
                                               nrow = nrow(lodo_sign),
                                               ncol = ncol(lodo_sign), byrow = TRUE),
                           na.rm = TRUE)
  tt$lodo_consistency <- lodo_consist[tt$gene]

  list(table = tt, lodo = lodo_sign, n_donors = length(donors))
}

# ---------------------------------------------------------------------------
# Run both cohorts, both tiers
# ---------------------------------------------------------------------------
res <- list()
for (ds in c("GSE183276", "GSE210622")) {
  for (tier in c("primary_ge50", "sensitivity_ge20")) {
    msg(sprintf("edgeR %s %s", ds, tier))
    res[[paste(ds, tier)]] <- run_tal_edger(ds, tier)
  }
}

deg183 <- res[["GSE183276 primary_ge50"]]$table
deg210 <- res[["GSE210622 primary_ge50"]]$table
write.csv(deg183, file.path(res_dir, "tal_pseudobulk_deg_GSE183276.csv"), row.names = FALSE)
write.csv(deg210, file.path(res_dir, "tal_pseudobulk_deg_GSE210622.csv"), row.names = FALSE)
write.csv(res[["GSE183276 sensitivity_ge20"]]$table,
          file.path(res_dir, "tal_pseudobulk_deg_GSE183276_ge20.csv"), row.names = FALSE)
write.csv(res[["GSE210622 sensitivity_ge20"]]$table,
          file.path(res_dir, "tal_pseudobulk_deg_GSE210622_ge20.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Concordance (shared genes, primary tier)
# ---------------------------------------------------------------------------
common <- intersect(deg183$gene, deg210$gene)
conc <- data.frame(
  gene = common,
  logFC_GSE183276 = deg183$logFC[match(common, deg183$gene)],
  FDR_GSE183276 = deg183$FDR[match(common, deg183$gene)],
  logFC_GSE210622 = deg210$logFC[match(common, deg210$gene)],
  FDR_GSE210622 = deg210$FDR[match(common, deg210$gene)],
  stringsAsFactors = FALSE
)
conc$same_direction <- sign(conc$logFC_GSE183276) == sign(conc$logFC_GSE210622)
conc$fdr_any <- conc$FDR_GSE183276 < 0.05 | conc$FDR_GSE210622 < 0.05
conc$fdr_both <- conc$FDR_GSE183276 < 0.05 & conc$FDR_GSE210622 < 0.05
conc$replicated <- conc$same_direction & conc$fdr_both
conc$consensus <- conc$same_direction & conc$fdr_any  # backward-compatible exploratory tier
conc$support_tier <- ifelse(conc$replicated, "FDR<0.05 both",
                            ifelse(conc$consensus, "FDR<0.05 one cohort", "not supported"))
conc$max_fdr <- pmax(conc$FDR_GSE183276, conc$FDR_GSE210622)
tier_rank <- match(conc$support_tier,
                   c("FDR<0.05 both", "FDR<0.05 one cohort", "not supported"))
conc <- conc[order(tier_rank, conc$max_fdr,
                   pmin(conc$FDR_GSE183276, conc$FDR_GSE210622)), ]
write.csv(conc, file.path(res_dir, "tal_pseudobulk_concordance.csv"), row.names = FALSE)

replicated_genes <- conc$gene[conc$replicated]
supported_genes <- conc$gene[conc$consensus]
msg(sprintf("shared=%d, same direction=%d, strict replicated=%d, single-cohort supported=%d",
            nrow(conc), sum(conc$same_direction), length(replicated_genes),
            sum(conc$support_tier == "FDR<0.05 one cohort")))

# ---------------------------------------------------------------------------
# Donor-level heatmap matrix (TMM logCPM, z-scored per gene within cohort)
# ---------------------------------------------------------------------------
donor_logcpm <- function(dataset, genes) {
  meta <- pb$metadata[[dataset]]
  qc_ds <- qc[qc$dataset == dataset & qc$segment == "TAL", ]
  qc_ds <- qc_ds[match(meta$donor, qc_ds$donor), ]
  keep <- qc_ds$eligible_50
  if (dataset == "GSE183276") keep <- keep & meta$diagnosis %in% c("AKI","Ref")
  else keep <- keep & meta$diagnosis %in% c("AKI","non_AKI")
  donors <- meta$donor[keep]
  counts <- pb$pseudobulk[[dataset]]$TAL[genes, donors, drop = FALSE]
  counts <- as.matrix(counts)
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  lcpm <- cpm(dge, log = TRUE, prior.count = 2)
  # z-score per gene across donors
  z <- t(scale(t(lcpm)))
  group_label <- as.character(meta$diagnosis[keep])
  colnames(z) <- paste(group_label, donors, sep = "|")
  z
}
heat_genes <- if (length(replicated_genes) > 0) replicated_genes[1:min(30, length(replicated_genes))] else FROZEN
z183 <- donor_logcpm("GSE183276", heat_genes)
z210 <- donor_logcpm("GSE210622", heat_genes)

# LODO export (primary tier, all genes -> keep memory small: only consensus + frozen)
lodo_rows <- list()
for (ds in c("GSE183276","GSE210622")) {
  r <- res[[paste(ds, "primary_ge50")]]
  genes_sel <- intersect(colnames(r$lodo), unique(c(replicated_genes, FROZEN)))
  for (g in genes_sel) {
    lodo_rows[[length(lodo_rows) + 1]] <- data.frame(
      dataset = ds, gene = g, donor = rownames(r$lodo),
      lodo_logFC_sign = r$lodo[, g], stringsAsFactors = FALSE)
  }
}
lodo_df <- do.call(rbind, lodo_rows)
write.csv(lodo_df, file.path(res_dir, "tal_pseudobulk_lodo.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
frozen_summary <- conc[conc$gene %in% FROZEN, ]
summary_df <- data.frame(
  metric = c(
    "cohort_GSE183276_n_donors_primary", "cohort_GSE210622_n_donors_primary",
    "GSE183276_n_genes_tested", "GSE210622_n_genes_tested",
    "GSE183276_n_FDR005", "GSE210622_n_FDR005",
    "shared_genes", "same_direction_shared", "fdr_any_shared", "supported_genes",
    "strict_replicated_genes",
    "frozen_6_same_direction",
    "pathway_analysis", "pathway_missing_requirement"
  ),
  value = c(
    res[["GSE183276 primary_ge50"]]$n_donors,
    res[["GSE210622 primary_ge50"]]$n_donors,
    nrow(deg183), nrow(deg210),
    sum(deg183$FDR < 0.05), sum(deg210$FDR < 0.05),
    nrow(conc), sum(conc$same_direction), sum(conc$fdr_any), length(supported_genes),
    length(replicated_genes),
    paste(frozen_summary$gene[frozen_summary$same_direction], collapse = ";"),
    "COMPLETE_PRIMARY_NEGATIVE",
    "none; see tal_pathway_report.md (0 replicated/supportive pathways under primary estimated-correlation CAMERA)"
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(res_dir, "tal_pseudobulk_summary.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Figure 3
# ---------------------------------------------------------------------------
AKI_COL <- "#B2182B"; CTRL_COL <- "#2166AC"; CONS_COL <- "#D95F02"

# Panel B: effect scatter
scatter_df <- conc
pB <- ggplot(scatter_df, aes(logFC_GSE183276, logFC_GSE210622)) +
  geom_point(aes(colour = support_tier), size = 1.0, alpha = 0.65) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  scale_colour_manual(values = c("FDR<0.05 both" = "#B2182B",
                                 "FDR<0.05 one cohort" = CONS_COL,
                                 "not supported" = "grey70"),
                      breaks = c("FDR<0.05 both", "FDR<0.05 one cohort", "not supported"),
                      labels = c("Both", "One", "None"),
                      name = "FDR support") +
  labs(x = "GSE183276 logFC (discovery)", y = "GSE210622 logFC (validation)",
       title = sprintf("Gene-effect concordance (strict n=%d)",
                       length(replicated_genes))) +
  coord_fixed() +
  theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# Panel A: MA summaries (one per cohort, combined small multiples)
vol <- rbind(
  transform(deg183, cohort = "GSE183276"),
  transform(deg210, cohort = "GSE210622")
)
vol$sig <- ifelse(vol$FDR < 0.05, "FDR<0.05", "ns")
pA <- ggplot(vol, aes(logCPM, logFC, colour = sig)) +
  geom_point(size = 0.35, alpha = 0.6) +
  facet_wrap(~ cohort, ncol = 2) +
  scale_colour_manual(values = c("FDR<0.05" = AKI_COL, "ns" = "grey70")) +
  labs(x = "Mean logCPM", y = "logFC (AKI vs control)",
       title = "TAL donor-pseudobulk MA summaries (>=50 exact author-TAL cells)") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# Panel C: consensus-gene donor heatmap (z-scored logCPM)
heat_df <- rbind(
  data.frame(cohort = "GSE183276", donor = rep(colnames(z183), each = nrow(z183)),
             gene = rep(rownames(z183), ncol(z183)), z = as.numeric(z183)),
  data.frame(cohort = "GSE210622", donor = rep(colnames(z210), each = nrow(z210)),
             gene = rep(rownames(z210), ncol(z210)), z = as.numeric(z210))
)
heat_df$gene <- factor(heat_df$gene, levels = rev(unique(heat_df$gene)))
heat_df$donor <- factor(heat_df$donor, levels = unique(heat_df$donor))
pC <- ggplot(heat_df, aes(donor, gene, fill = z)) +
  geom_tile() +
  facet_grid(~ cohort, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "z(logCPM)") +
  labs(x = "Diagnosis | donor", y = NULL,
       title = "Top strict-replication genes, donor-level expression") +
  theme_minimal(base_size = 7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.2),
        axis.text.y = element_text(size = 3.6),
        strip.text = element_text(face = "bold", size = 7),
        legend.key.size = unit(0.28, "cm"), legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 6))

# Panel D: frozen six-gene logFC by cohort (real data; replaces deferred pathways)
fz <- conc[conc$gene %in% FROZEN, ]
fz_long <- rbind(
  data.frame(gene = fz$gene, cohort = "GSE183276", logFC = fz$logFC_GSE183276,
             FDR = fz$FDR_GSE183276),
  data.frame(gene = fz$gene, cohort = "GSE210622", logFC = fz$logFC_GSE210622,
             FDR = fz$FDR_GSE210622)
)
fz_long$gene <- factor(fz_long$gene, levels = FROZEN)
fz_long$lab <- ifelse(fz_long$FDR < 0.05, "*", "")
pD <- ggplot(fz_long, aes(gene, logFC, fill = cohort)) +
  geom_col(position = position_dodge(0.75), width = 0.65) +
  geom_text(aes(label = lab), position = position_dodge(0.75), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("GSE183276" = "#8C6BB1", "GSE210622" = "#6AAED6")) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  labs(x = NULL, y = "logFC (AKI vs control)",
       title = "Frozen six-gene core in TAL pseudobulk (* FDR<0.05)") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        panel.grid.minor = element_blank(), legend.position = "bottom")

p_all <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.55)) +
  plot_annotation(
    title = "Figure 3  |  Exploratory TAL donor-pseudobulk extension across two cohorts",
    subtitle = sprintf(paste0("Cohorts analysed separately: GSE183276 n=%d; GSE210622 n=%d. ",
                             "Strict replication = same direction and FDR<0.05 in both. ",
                             "GSE210622 remains limited by mixed source/COVID context; no pathway result is shown."),
                       res[["GSE183276 primary_ge50"]]$n_donors,
                       res[["GSE210622 primary_ge50"]]$n_donors),
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 7, colour = "grey30"),
                  plot.tag = element_text(face = "bold", size = 10)))

ggsave(file.path(fig_dir, "figure3_tal_pseudobulk.png"), p_all, width = 11.7, height = 7.8, dpi = 300)
ggsave(file.path(fig_dir, "figure3_tal_pseudobulk.pdf"), p_all, width = 11.7, height = 7.8)
msg("wrote figures/extension/figure3_tal_pseudobulk.png/.pdf")
msg("DONE")
