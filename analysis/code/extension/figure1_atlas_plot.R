#!/usr/bin/env Rscript
# =============================================================================
# figure1_atlas_plot.R — Figure 1 (extension): study design + cell atlas
# =============================================================================
# Panels (analysis_extension_plan.md WP1 / manuscript_v4 storyboard):
#   A  study and evidence-flow schematic (donor/sample counts; descriptive)
#   B  GSE183276 global UMAP, broad kidney classes (author labels harmonised)
#   C  GSE210622 global UMAP, same palette
#   D  annotation marker dot plot (rows=prespecified markers, cols=broad classes)
#   E  donor-by-cell-type contribution (proportional dots, no P values)
#   F  recovered TAL/LOH cells per donor with >=20 / >=50 thresholds
# UMAPs are descriptive; no cell-level P values; no abundance claims.
#
# Inputs: results/extension/atlas_*.csv (produced by figure1_atlas_prep.R)
# Outputs: figures/extension/figure1_atlas.png/.pdf
# =============================================================================

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) == 1L) {
  # script lives at code/extension/<name>.R -> project root is 3 levels up
  project_dir <- dirname(dirname(dirname(normalizePath(sub("^--file=", "", script_arg)))))
} else {
  project_dir <- getwd()
}
setwd(project_dir)

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(dplyr)
})

res_dir <- file.path(project_dir, "results", "extension")
fig_dir <- file.path(project_dir, "figures", "extension")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Load atlas outputs
# ---------------------------------------------------------------------------
u183 <- read.csv(file.path(res_dir, "atlas_GSE183276_umap.csv"), stringsAsFactors = FALSE)
u210 <- read.csv(file.path(res_dir, "atlas_GSE210622_umap.csv"), stringsAsFactors = FALSE)
mk183 <- read.csv(file.path(res_dir, "atlas_GSE183276_markers.csv"), stringsAsFactors = FALSE)
mk210 <- read.csv(file.path(res_dir, "atlas_GSE210622_markers.csv"), stringsAsFactors = FALSE)
cnt183 <- read.csv(file.path(res_dir, "atlas_GSE183276_counts.csv"), stringsAsFactors = FALSE)
cnt210 <- read.csv(file.path(res_dir, "atlas_GSE210622_counts.csv"), stringsAsFactors = FALSE)
sum183 <- read.csv(file.path(res_dir, "atlas_GSE183276_summary.csv"), stringsAsFactors = FALSE)
sum210 <- read.csv(file.path(res_dir, "atlas_GSE210622_summary.csv"), stringsAsFactors = FALSE)
pal <- read.csv(file.path(res_dir, "atlas_broad_palette.csv"), stringsAsFactors = FALSE)
palette_vec <- setNames(pal$colour, pal$broad_class)

# Class order (fixed)
class_order <- c("PT","TAL/LOH","DCT","CNT","CD","POD","EC","Stroma","IMM","Other")
class_order <- class_order[class_order %in% unique(c(u183$broad_class, u210$broad_class))]

# Condition palette (visual system): AKI muted red, reference blue
cond_col <- c("AKI" = "#B2182B", "Ref" = "#2166AC", "non_AKI" = "#2166AC", "CKD" = "#E0A800")
AKI_COL <- "#B2182B"; REF_COL <- "#2166AC"

# ---------------------------------------------------------------------------
# Panel A — study design schematic (drawn with ggplot annotation)
# ---------------------------------------------------------------------------
ov <- read.csv(file.path(res_dir, "atlas_cohort_overview.csv"), stringsAsFactors = FALSE)
pa <- ggplot() +
  annotate("rect", xmin = 0.2, xmax = 4.6, ymin = 6.6, ymax = 9.6, fill = "#F7F7F7", colour = "grey40") +
  annotate("rect", xmin = 5.4, xmax = 9.8, ymin = 6.6, ymax = 9.6, fill = "#F7F7F7", colour = "grey40") +
  annotate("text", x = 2.4, y = 9.2, label = "GSE183276 \u2014 scRNA-seq", fontface = "bold", size = 3.0) +
  annotate("text", x = 7.6, y = 9.2, label = "GSE210622 \u2014 snRNA-seq", fontface = "bold", size = 3.0) +
  annotate("text", x = 2.4, y = 8.15, label = sprintf("AKI n=%d  \u00b7  Ref n=%d  \u00b7  CKD n=%d", sum183$n_AKI, sum183$n_Ref, sum183$n_CKD), size = 2.2) +
  annotate("text", x = 7.6, y = 8.15, label = sprintf("AKI n=%d  \u00b7  non-AKI n=%d", sum210$n_AKI, sum210$n_nonAKI), size = 2.2) +
  annotate("text", x = 2.4, y = 7.1, label = sprintf("cells n=%s  \u00b7  donors n=%d  \u00b7  TAL/LOH n=%s",
                                                     comma(sum183$n_cells), sum183$n_donors, comma(sum183$n_TAL_LOH_cells)), size = 2.2) +
  annotate("text", x = 7.6, y = 7.1, label = sprintf("nuclei n=%s  \u00b7  donors n=%d  \u00b7  TAL/LOH n=%s",
                                                     comma(sum210$n_cells), sum210$n_donors, comma(sum210$n_TAL_LOH_cells)), size = 2.2) +
  annotate("segment", x = 4.7, xend = 5.3, y = 7.9, yend = 7.9, arrow = arrow(length = unit(0.08, "in")), linewidth = 0.4) +
  annotate("rect", xmin = 1.2, xmax = 4.4, ymin = 3.5, ymax = 5.8, fill = "#FFFFFF", colour = "grey60") +
  annotate("rect", xmin = 5.6, xmax = 8.8, ymin = 3.5, ymax = 5.8, fill = "#FFFFFF", colour = "grey60") +
  annotate("text", x = 2.8, y = 5.5, label = "Discovery (donor-pseudobulk)", size = 2.4) +
  annotate("text", x = 7.2, y = 5.5, label = "Independent validation (donor-pseudobulk)", size = 2.4) +
  annotate("text", x = 2.8, y = 4.7, label = "6-gene TAL score, effect + CI", size = 2.2) +
  annotate("text", x = 7.2, y = 4.7, label = "6-gene TAL score, effect + CI", size = 2.2) +
  annotate("text", x = 2.8, y = 3.9, label = "TAL pseudobulk DEG (Fig 3)", size = 2.2) +
  annotate("text", x = 7.2, y = 3.9, label = "TAL pseudobulk DEG (Fig 3)", size = 2.2) +
  annotate("rect", xmin = 0.8, xmax = 5.0, ymin = 0.7, ymax = 2.9, fill = "#FBF3E4", colour = "#E0A800") +
  annotate("text", x = 2.9, y = 2.35, label = "Cross-sectional tissue support", fontface = "bold", size = 2.4) +
  annotate("text", x = 2.9, y = 1.35, label = "GSE139061 bulk RNA-seq\nGSE30718 microarray", size = 2.1) +
  annotate("rect", xmin = 5.2, xmax = 9.4, ymin = 0.7, ymax = 2.9, fill = "#EDF4FA", colour = "#2166AC") +
  annotate("text", x = 7.3, y = 2.35, label = "Paired transplant-IRI support", fontface = "bold", size = 2.4) +
  annotate("text", x = 7.3, y = 1.35, label = "GSE126805: 41 pre/post allograft pairs\n42 grafts \u00b7 163 protocol biopsies", size = 2.1) +
  annotate("segment", x = 3.5, xend = 3.5, y = 5.8, yend = 6.6, arrow = arrow(length = unit(0.06, "in")), colour = "grey50", linewidth = 0.3) +
  annotate("segment", x = 6.5, xend = 6.5, y = 5.8, yend = 6.6, arrow = arrow(length = unit(0.06, "in")), colour = "grey50", linewidth = 0.3) +
  annotate("segment", x = 3.0, xend = 3.0, y = 2.9, yend = 3.5, arrow = arrow(length = unit(0.06, "in")), colour = "grey50", linewidth = 0.3) +
  annotate("segment", x = 7.0, xend = 7.0, y = 2.9, yend = 3.5, arrow = arrow(length = unit(0.06, "in")), colour = "grey50", linewidth = 0.3) +
  theme_void() +
  annotate("text", x = 5, y = 10.05, label = "Separate cell-resolved cohorts; donors are the statistical unit", size = 2.0, fontface = "italic") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10), clip = "off")

# ---------------------------------------------------------------------------
# Panels B/C — global UMAPs
# ---------------------------------------------------------------------------
umap_plot <- function(u, cohort_label, sub) {
  u$broad_class <- factor(u$broad_class, levels = class_order)
  ggplot(u, aes(umap1, umap2, colour = broad_class)) +
    geom_point(size = 0.12, alpha = 0.55) +
    scale_colour_manual(values = palette_vec, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1),
                                 nrow = 2, byrow = TRUE)) +
    labs(title = cohort_label, subtitle = sub, colour = "Broad class") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 7),
          legend.position = "bottom", legend.key.size = unit(0.28, "cm"),
          legend.text = element_text(size = 6.4), legend.title = element_text(size = 7),
          legend.margin = margin(0, 0, 0, 0))
}

pB <- umap_plot(u183, "GSE183276", sprintf("scRNA-seq; %s cells, %d donors", comma(nrow(u183)), length(unique(u183$donor))))
pC <- umap_plot(u210, "GSE210622", sprintf("snRNA-seq; %s nuclei, %d donors", comma(nrow(u210)), length(unique(u210$donor))))

# ---------------------------------------------------------------------------
# Panel D — marker dot plot (both cohorts side by side)
# ---------------------------------------------------------------------------
mk <- rbind(mk183, mk210)
mk$broad_class <- factor(mk$broad_class, levels = class_order)
safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
mk$mean_expr_z <- ave(mk$mean_expr, interaction(mk$cohort, mk$gene), FUN = safe_z)
mk$gene <- factor(mk$gene, levels = rev(unique(marker_genes <- c(
  "LRP2","SLC34A1","SLC5A2","SLC12A1","UMOD","CLDN16","KCNJ10","PTH1R","CLDN10",
  "SLC12A3","PVALB","CALB1","SCNN1G","ATP6V1B1","SLC4A1","NPHS1","NPHS2",
  "PECAM1","EMCN","DCN","COL1A1","PTPRC"))))
pD <- ggplot(mk, aes(broad_class, gene)) +
  geom_point(aes(size = pct_expr, colour = mean_expr_z)) +
  facet_grid(~ cohort, scales = "free_x", space = "free_x") +
  scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, limits = c(-2.5, 2.5), oob = squish,
                         name = "mean expression\n(within-gene z)") +
  scale_size_continuous(range = c(0.3, 3.4), name = "pct expr") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.6),
        axis.text.y = element_text(size = 6.4),
        strip.text = element_text(face = "bold", size = 8),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5, margin = margin(b = 2)),
        panel.grid.minor = element_blank())

# ---------------------------------------------------------------------------
# Panel E — donor-by-cell-type contribution
# ---------------------------------------------------------------------------
cnt <- rbind(cnt183, cnt210)
cnt$broad_class <- factor(cnt$broad_class, levels = class_order)
tot <- aggregate(n_cells ~ donor + cohort, data = cnt, FUN = sum)
cnt <- merge(cnt, tot, by = c("donor","cohort"), suffixes = c("", "_total"))
cnt$frac <- cnt$n_cells / cnt$n_cells_total
cnt$diagnosis <- factor(cnt$diagnosis, levels = c("Ref","non_AKI","AKI","CKD"))
cnt <- cnt[order(cnt$diagnosis, cnt$donor), ]
# donor id alone is unique within each cohort facet; drop the cohort prefix so
# the rotated labels stay short and clear the bottom edge.
cnt$donor_f <- factor(cnt$donor, levels = unique(cnt$donor))
pE <- ggplot(cnt, aes(donor_f, broad_class)) +
  geom_point(aes(size = frac, colour = diagnosis)) +
  facet_grid(~ cohort, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = cond_col, name = "Diagnosis") +
  scale_size_continuous(range = c(0.4, 4.2), name = "Fraction of donor cells") +
  labs(x = "Donor", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.2),
        axis.text.y = element_text(size = 6.4),
        strip.text = element_text(face = "bold", size = 8),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5),
        plot.margin = margin(0, 0, 0.25, 0, "in"),
        panel.grid.minor = element_blank())

# ---------------------------------------------------------------------------
# Panel F — TAL/LOH recovery per donor with >=20 / >=50 thresholds
# ---------------------------------------------------------------------------
tal_cnt <- cnt[cnt$broad_class == "TAL/LOH", ]
tal_cnt <- tal_cnt[order(tal_cnt$diagnosis, tal_cnt$donor), ]
tal_cnt$donor_f <- factor(tal_cnt$donor, levels = unique(tal_cnt$donor))
pF <- ggplot(tal_cnt, aes(donor_f, n_cells, colour = diagnosis)) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_hline(yintercept = 20, linetype = "dotted", colour = "grey65", linewidth = 0.4) +
  facet_grid(~ cohort, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = cond_col, name = "Diagnosis") +
  scale_y_continuous(trans = log1p_trans(), labels = comma) +
  labs(x = "Donor", y = "Recovered broad TAL/LOH cells per donor (log1p)",
       subtitle = "Atlas recovery/QC only; not biological abundance or the exact-author-TAL analysis tier") +
  annotate("text", x = -Inf, y = 50, label = "  \u226550", hjust = 0, vjust = -0.3, size = 2.3, colour = "grey45") +
  annotate("text", x = -Inf, y = 20, label = "  \u226520", hjust = 0, vjust = -0.3, size = 2.3, colour = "grey65") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.2),
        strip.text = element_text(face = "bold", size = 8),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5),
        plot.margin = margin(0, 0, 0.25, 0, "in"),
        panel.grid.minor = element_blank())

# ---------------------------------------------------------------------------
# Compose
# ---------------------------------------------------------------------------
layout <- "
AAAAAAA
BBBBCCC
BBBBCCC
DDDDDDD
DDDDDDD
EEEEFFF
"
p_all <- pa + pB + pC + pD + pE + pF +
  plot_layout(design = layout) +
  plot_annotation(
    title = "Figure 1  |  Author-annotated, marker-audited human kidney cell atlases (descriptive)",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.tag = element_text(face = "bold", size = 10))
  )

ggsave(file.path(fig_dir, "figure1_atlas.png"), p_all, width = 11.7, height = 12.0, dpi = 300)
ggsave(file.path(fig_dir, "figure1_atlas.pdf"), p_all, width = 11.7, height = 12.0)
msg <- function(...) cat(paste0(..., "\n"))
msg("wrote figures/extension/figure1_atlas.png/.pdf")
