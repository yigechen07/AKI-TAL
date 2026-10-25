#!/usr/bin/env Rscript
# Deterministic Figure 3 renderer: gene-level pseudobulk replication plus the
# pathway-robustness audit. Reads results only; performs no new inference.

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
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

res_dir <- file.path(project_dir, "results", "extension")
fig_dir <- file.path(project_dir, "figures", "extension")
pb_path <- file.path(project_dir, "..", "AKI_DCT_analysis", "results", "stage3b_pseudobulk.rds")
qc_path <- file.path(project_dir, "..", "AKI_DCT_analysis", "results", "stage3b_qc.csv")

deg183 <- fread(file.path(res_dir, "tal_pseudobulk_deg_GSE183276.csv"))
deg210 <- fread(file.path(res_dir, "tal_pseudobulk_deg_GSE210622.csv"))
conc <- fread(file.path(res_dir, "tal_pseudobulk_concordance.csv"))
path_primary <- fread(file.path(res_dir, "tal_pathway_summary.csv"))
path_sens <- fread(file.path(res_dir, "tal_pathway_summary_fixed_cor_sensitivity.csv"))
pb <- readRDS(pb_path)
qc <- fread(qc_path)

AKI_COL <- "#B2182B"
ORANGE <- "#D95F02"
BLUE <- "#2166AC"
GREY <- "#BDBDBD"

# Panel A: cohort-specific MA summaries.
ma <- rbindlist(list(
  copy(deg183)[, cohort := "GSE183276"],
  copy(deg210)[, cohort := "GSE210622"]
))
ma[, significance := fifelse(FDR < 0.05, "FDR<0.05", "Not FDR-significant")]
pA <- ggplot(ma, aes(logCPM, logFC, colour = significance)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.3) +
  geom_point(size = 0.35, alpha = 0.58) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_colour_manual(values = c("FDR<0.05" = AKI_COL,
                                 "Not FDR-significant" = "grey75"), name = NULL) +
  labs(x = "Mean logCPM", y = "logFC (AKI vs comparator)",
       subtitle = "Cohort-specific edgeR QL; donor is the statistical unit") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold"))

# Panel B: gene-effect concordance.
conc[, support_label := factor(
  fifelse(replicated, "Dual-FDR replication",
    fifelse(support_tier == "FDR<0.05 one cohort", "One-cohort FDR support", "Other")),
  levels = c("Dual-FDR replication", "One-cohort FDR support", "Other")
)]
rho <- cor(conc$logFC_GSE183276, conc$logFC_GSE210622, method = "spearman")
pB <- ggplot(conc, aes(logFC_GSE183276, logFC_GSE210622, colour = support_label)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(size = 0.65, alpha = 0.62) +
  scale_colour_manual(values = c("Dual-FDR replication" = AKI_COL,
                                 "One-cohort FDR support" = ORANGE,
                                 "Other" = "grey78"), name = NULL) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("Spearman rho = %.2f\nstrict replicated genes = %d",
                           rho, sum(conc$replicated)),
           hjust = 1.05, vjust = -0.35, size = 3, colour = "grey25") +
  labs(x = "GSE183276 logFC", y = "GSE210622 logFC",
       subtitle = "Gene effects compared across separately fitted cohorts") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# Panel C: top 24 strict-replication genes, selected deterministically by the
# existing concordance table order (max cross-cohort FDR, then evidence).
heat_genes <- conc[replicated == TRUE, head(gene, 24L)]
stopifnot(length(heat_genes) == 24L)

donor_logcpm <- function(dataset, genes) {
  meta <- pb$metadata[[dataset]]
  qc_ds <- qc[qc$dataset == dataset & qc$segment == "TAL"]
  qc_ds <- qc_ds[match(meta$donor, qc_ds$donor)]
  stopifnot(all(qc_ds$donor == meta$donor))
  keep <- qc_ds$eligible_50
  if (dataset == "GSE183276") keep <- keep & meta$diagnosis %in% c("AKI", "Ref")
  else keep <- keep & meta$diagnosis %in% c("AKI", "non_AKI")
  donors <- meta$donor[keep]
  counts <- as.matrix(pb$pseudobulk[[dataset]]$TAL[genes, donors, drop = FALSE])
  y <- calcNormFactors(DGEList(counts = counts), method = "TMM")
  z <- t(scale(t(cpm(y, log = TRUE, prior.count = 2))))
  diagnosis <- ifelse(meta$diagnosis[keep] == "AKI", "AKI", "Comparator")
  colnames(z) <- paste(diagnosis, donors, sep = " | ")
  z
}

z183 <- donor_logcpm("GSE183276", heat_genes)
z210 <- donor_logcpm("GSE210622", heat_genes)
heat <- rbindlist(list(
  data.table(cohort = "GSE183276", donor = rep(colnames(z183), each = nrow(z183)),
             gene = rep(rownames(z183), ncol(z183)), z = as.numeric(z183)),
  data.table(cohort = "GSE210622", donor = rep(colnames(z210), each = nrow(z210)),
             gene = rep(rownames(z210), ncol(z210)), z = as.numeric(z210))
))
heat[, gene := factor(gene, levels = rev(heat_genes))]
heat[, donor := factor(donor, levels = unique(donor))]
pC <- ggplot(heat, aes(donor, gene, fill = z)) +
  geom_tile() +
  facet_grid(~ cohort, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = BLUE, mid = "white", high = AKI_COL,
                       limits = c(-2.5, 2.5), oob = scales::squish,
                       name = "z(logCPM)") +
  labs(x = "Diagnosis | donor", y = NULL,
       subtitle = "Top strict-replication genes; z-scored within cohort") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.8),
        axis.text.y = element_text(size = 5.6), panel.grid = element_blank(),
        strip.text = element_text(face = "bold"), legend.position = "right")

# Panel D: pathway robustness audit. The primary estimated-correlation CAMERA
# is deliberately displayed even when every count is zero. Fixed 0.01 is a
# sensitivity and cannot override the primary negative result.
path_primary[, method := "Estimated correlation\n(primary)"]
path_sens[, method := "Fixed correlation 0.01\n(sensitivity)"]
path <- rbindlist(list(path_primary, path_sens), use.names = TRUE)
path_long <- melt(path,
                  id.vars = c("collection", "method", "n_tested"),
                  measure.vars = c("n_replicated", "n_supportive"),
                  variable.name = "evidence", value.name = "n_pathways")
path_long[, evidence := factor(evidence,
                               levels = c("n_replicated", "n_supportive"),
                               labels = c("Replicated (dual FDR)", "Supportive"))]
path_long[, method := factor(method,
                             levels = c("Estimated correlation\n(primary)",
                                        "Fixed correlation 0.01\n(sensitivity)"))]
path_long[, collection := factor(collection, levels = c("GO_BP", "Hallmark", "Reactome"),
                                 labels = c("GO BP", "Hallmark", "Reactome"))]
path_long[, percent_tested := 100 * n_pathways / n_tested]
pD <- ggplot(path_long, aes(collection, percent_tested, fill = evidence)) +
  geom_col(position = "stack", width = 0.66) +
  geom_text(aes(label = n_pathways), position = position_stack(vjust = 0.5),
            colour = "white", size = 3, show.legend = FALSE) +
  geom_text(data = path_long[n_pathways == 0], aes(x = collection, y = 0.7, label = "0"),
            inherit.aes = FALSE, colour = "grey35", size = 3) +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_manual(values = c("Replicated (dual FDR)" = AKI_COL,
                               "Supportive" = ORANGE), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "% of tested pathways",
       subtitle = paste0("Primary CAMERA: 0 replicated/supportive pathways.\n",
                         "Fixed-correlation candidates are method-sensitive and are not promoted.")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom")

fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.35), widths = c(1.18, 1)) +
  plot_annotation(
    title = "Figure 3: Reproducible TAL pseudobulk genes and pathway robustness audit",
    subtitle = paste0(
      "GSE183276: 10 AKI vs 17 Reference; GSE210622: 8 AKI vs 4 non-AKI. ",
      "Strict gene replication requires same direction and FDR<0.05 in both cohorts."
    ),
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8, colour = "grey25"),
                  plot.tag = element_text(face = "bold", size = 13))
  )

ggsave(file.path(fig_dir, "figure3_tal_pseudobulk.png"), fig,
       width = 12, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "figure3_tal_pseudobulk.pdf"), fig,
       width = 12, height = 8.5, bg = "white")
message("WROTE figures/extension/figure3_tal_pseudobulk.png/.pdf")
message("DONE")
