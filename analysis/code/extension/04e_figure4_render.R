#!/usr/bin/env Rscript
# =============================================================================
# 04e_figure4_render.R — Figure 4 (extension): render panels A-E (deterministic)
# =============================================================================
# Reads ONLY the saved analysis outputs:
#   results/extension/figure4_plotting_data.rds   (cell-level; Panel A)
#   results/extension/tal_state_donor_scores.csv  (Panel B)
#   results/extension/tal_state_effects.csv       (Panel C)
#   results/extension/tal_state_decomposition.csv (Panels D, E)
# Writes figures/extension/figure4_tal_states.png and .pdf
#
# Panel A : TAL-only UMAP, SUBSET of the existing global UMAP (no recomputation),
#           coloured by continuous reference-TAL UCell score (no hard classification)
# Panel B : donor-level state scores by diagnosis (primary mean_ge50 donors)
# Panel C : six primary cohort x state effects, forest, facet per state,
#           independent axes per state — raw UCell units, magnitudes NOT
#           comparable across programmes (explicit label)
# Panel D : six-gene AKI effect, unadjusted vs adjusted for state score
# Panel E : per-gene detection-fraction and positive-expression effects
#           (secondary endpoints, BH-FDR annotated)
# =============================================================================

suppressMessages({
  for (pkg in c("data.table", "ggplot2", "patchwork", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package missing: ", pkg)
  }
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1 && nzchar(cli_args[1])) {
  proj_root <- normalizePath(cli_args[1])
} else {
  cand <- c(getwd(), file.path(getwd(), "..", ".."))
  hit <- cand[vapply(cand, function(p) file.exists(file.path(p, "results", "extension")), logical(1))]
  if (length(hit) == 0) stop("Could not locate project root (no results/extension). Pass project root as first CLI arg.")
  proj_root <- normalizePath(hit[1])
}

set.seed(20260809)   # deterministic jitter in Panel B
RES  <- file.path(proj_root, "results", "extension")
FIG  <- file.path(proj_root, "figures", "extension")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

plot_data <- readRDS(file.path(RES, "figure4_plotting_data.rds"))
don <- fread(file.path(RES, "tal_state_donor_scores.csv"))
eff <- fread(file.path(RES, "tal_state_effects.csv"))
dec <- fread(file.path(RES, "tal_state_decomposition.csv"))

# ---- colour threads (CVD-safe, diverging-family; see figure-style skill) ----
col_aki   <- "#C0392B"
col_ref   <- "#7F8C8D"   # Reference / non_AKI (grey)
col_ckd   <- "#2E86C1"
col_ref_TAL    <- "#2166AC"
col_adaptive   <- "#B35806"
col_degenerate <- "#762A83"

# diagnosis display labels
lab <- function(x) ifelse(x == "AKI", "AKI",
                   ifelse(x == "Ref", "Reference",
                   ifelse(x == "CKD", "CKD", "non-AKI")))
don[, diagnosis_lab := factor(lab(diagnosis), levels = c("AKI", "Reference", "CKD", "non-AKI"))]
plot_data[, cohort_lab := factor(cohort, levels = c("GSE183276", "GSE210622"))]

state_lab <- function(x) c(
  reference_TAL = "Reference TAL (17 genes)",
  adaptive_TAL  = "Adaptive TAL (281 genes)",
  degenerative_TAL = "Degenerative TAL (767 genes)"
)[x]

# ============================ PANEL A =======================================
# TAL-only subset of the existing global UMAP; continuous reference-TAL score
pa_dat <- plot_data[!is.na(umap1) & !is.na(umap2)]
pa <- ggplot(pa_dat, aes(umap1, umap2)) +
  geom_point(aes(colour = score_reference_TAL), size = 0.25, stroke = 0) +
  scale_colour_viridis_c(option = "C", name = "Reference TAL\nUCell score") +
  facet_wrap(~ cohort_lab, nrow = 1) +
  labs(subtitle = "Subset of global UMAP (coordinates reused, not recomputed); strict author TAL cells") +
  theme_minimal(base_size = 10) +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.background = element_rect(fill = "grey97", colour = NA),
        strip.text = element_text(face = "bold"),
        legend.position = "right", legend.key.height = unit(1.4, "cm"))

# ============================ PANEL B =======================================
# donor-level state scores by diagnosis (primary: mean_ge50, primary sets)
pb_dat <- don[aggregation == "mean_ge50" & set_type == "primary" &
                ((cohort == "GSE183276" & diagnosis %in% c("AKI", "Ref")) |
                 (cohort == "GSE210622" & diagnosis %in% c("AKI", "non_AKI")))]
pb_dat[, comparison_group := factor(ifelse(diagnosis == "AKI", "AKI", "Comparator"),
                                     levels = c("AKI", "Comparator"))]
pb_dat[, state_strip := factor(
  c(reference_TAL = "Reference TAL\n17 genes",
    adaptive_TAL = "Adaptive TAL\n281 genes",
    degenerative_TAL = "Degenerative TAL\n767 genes")[state],
  levels = c("Reference TAL\n17 genes", "Adaptive TAL\n281 genes",
             "Degenerative TAL\n767 genes"))]
pb_dat[, cohort_lab := factor(cohort, levels = c("GSE183276","GSE210622"))]
# jitter for small n
set.seed(20260809)
pb_dat[, jx := as.numeric(comparison_group) + runif(.N, -0.16, 0.16)]
pb <- ggplot(pb_dat, aes(comparison_group, score)) +
  geom_boxplot(aes(fill = comparison_group), outlier.shape = NA, width = 0.5, alpha = 0.35) +
  geom_point(aes(jx, score, colour = comparison_group), size = 1.6, alpha = 0.85) +
  facet_grid(state_strip ~ cohort_lab, scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(AKI = col_aki, Comparator = col_ref),
                    guide = "none") +
  scale_colour_manual(values = c(AKI = col_aki, Comparator = col_ref), name = NULL) +
  labs(y = "Donor mean UCell score (raw units)",
       subtitle = paste0("Donor is the unit (>=50 TAL cells); Comparator = Reference ",
                         "in GSE183276 and non-AKI in GSE210622")) +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7.5),
        strip.placement = "outside",
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# ============================ PANEL C =======================================
# six primary effects, forest, facet per state with independent axes
pc_dat <- eff[aggregation == "mean_ge50" & set_type == "primary" &
              contrast %in% c("AKI_vs_Reference", "AKI_vs_nonAKI")]
pc_dat[, state_lab := factor(state_lab(state), levels = state_lab(c("reference_TAL","adaptive_TAL","degenerative_TAL")))]
pc_dat[, cohort_lab := factor(cohort, levels = c("GSE183276","GSE210622"))]
pc_dat[, sig := fifelse(!is.na(fdr_bh) & fdr_bh < 0.05, "FDR<0.05", "ns")]
pc <- ggplot(pc_dat, aes(estimate, cohort_lab, colour = state)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high, shape = sig), size = 0.55, linewidth = 0.7) +
  facet_grid(. ~ state_lab, scales = "free_x") +
  scale_colour_manual(values = c(reference_TAL = col_ref_TAL, adaptive_TAL = col_adaptive,
                                 degenerative_TAL = col_degenerate), guide = "none") +
  scale_shape_manual(values = c("FDR<0.05" = 16, "ns" = 1), name = NULL) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
  labs(y = NULL, x = "AKI effect on donor mean UCell score (raw units)",
       subtitle = "Primary analysis; per-state axis; raw magnitudes are NOT comparable across programmes") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom", panel.grid.minor = element_blank())

# ============================ PANEL D =======================================
# six-gene AKI effect unadjusted vs adjusted for state score (primary mean_ge50)
pd_dat <- rbind(
  dec[model_type == "six_gene_unadjusted",
      .(cohort, model = "Unadjusted", state = NA_character_, estimate, ci_low, ci_high, p_value)],
  dec[model_type == "six_gene_adjusted_state" & set_type == "primary" & aggregation == "mean_ge50" & term == "AKI",
      .(cohort, model = paste0("Adjusted for ", state_lab(state)), state, estimate, ci_low, ci_high, p_value)]
)
pd_dat[, cohort_lab := factor(cohort, levels = c("GSE183276","GSE210622"))]
pd_dat[, model := factor(model, levels = c("Unadjusted",
                        "Adjusted for Reference TAL (17 genes)",
                        "Adjusted for Adaptive TAL (281 genes)",
                        "Adjusted for Degenerative TAL (767 genes)"))]
pd <- ggplot(pd_dat, aes(estimate, model, colour = model)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high), size = 0.6, linewidth = 0.7) +
  facet_wrap(~ cohort_lab, nrow = 1) +
  scale_colour_manual(values = c("Unadjusted" = "grey40",
                                 "Adjusted for Reference TAL (17 genes)" = col_ref_TAL,
                                 "Adjusted for Adaptive TAL (281 genes)" = col_adaptive,
                                 "Adjusted for Degenerative TAL (767 genes)" = col_degenerate),
                      guide = "none") +
  labs(y = NULL, x = "AKI effect on six-gene TAL identity (frozen units)",
       subtitle = "Change in AKI effect after adding state score; association only, not mediation") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())

# ============================ PANEL E =======================================
# Aggregate detection is the primary decomposition endpoint; per-gene
# detection and positive expression are secondary.
pe_dat <- rbind(
  dec[model_type == "agg_detection_primary",
      .(cohort, model_type, gene = "6-gene aggregate", estimate, ci_low, ci_high,
        p_value, fdr_bh)],
  dec[model_type %in% c("pergene_detection", "pergene_positive_expression"),
      .(cohort, model_type, gene, estimate, ci_low, ci_high, p_value, fdr_bh)],
  fill = TRUE
)
pe_dat[, endpoint := factor(ifelse(model_type == "pergene_detection",
                                   "Detection fraction",
                            ifelse(model_type == "agg_detection_primary",
                                   "Detection fraction", "Positive expression")),
                            levels = c("Detection fraction", "Positive expression"))]
pe_dat[, cohort_lab := factor(cohort, levels = c("GSE183276","GSE210622"))]
pe_dat[, evidence := fifelse(model_type == "agg_detection_primary" & p_value < 0.05,
                             "Primary aggregate P<0.05",
                       fifelse(!is.na(fdr_bh) & fdr_bh < 0.05,
                               "Secondary FDR<0.05", "Not significant"))]
pe_dat[, gene := factor(gene, levels = rev(c("6-gene aggregate", "SLC12A1", "UMOD",
                                             "CLDN16", "KCNJ10", "PTH1R", "CLDN10")))]
pe <- ggplot(pe_dat, aes(estimate, gene, colour = evidence, shape = evidence)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high), size = 0.55, linewidth = 0.6) +
  facet_grid(endpoint ~ cohort_lab, scales = "free", space = "free_y") +
  scale_colour_manual(values = c("Primary aggregate P<0.05" = col_ref_TAL,
                                 "Secondary FDR<0.05" = "#C0392B",
                                 "Not significant" = "grey50"), name = NULL) +
  scale_shape_manual(values = c("Primary aggregate P<0.05" = 18,
                                "Secondary FDR<0.05" = 16,
                                "Not significant" = 16), name = NULL) +
  labs(y = NULL, x = "AKI effect (detection fraction / mean expression, raw units)",
       subtitle = paste0("Within recovered strict TAL; aggregate detection is primary; ",
                         "per-gene endpoints are secondary with BH-FDR")) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom",
        panel.grid.minor = element_blank())

# ============================ COMPOSE =======================================
fig <- (pa / pb / pc / pd / pe) +
  plot_annotation(
    title = "Figure 4: TAL state context and six-gene score decomposition (extension analysis)",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.tag = element_text(face = "bold", size = 14))
  ) &
  theme(plot.subtitle = element_text(size = 8, colour = "grey25"))

ggsave(file.path(FIG, "figure4_tal_states.png"), fig, width = 11, height = 19, dpi = 300, bg = "white")
ggsave(file.path(FIG, "figure4_tal_states.pdf"), fig, width = 11, height = 19, bg = "white")
cat("WROTE", file.path(FIG, "figure4_tal_states.png"), "\n")
cat("WROTE", file.path(FIG, "figure4_tal_states.pdf"), "\n")
cat("panels: A cells", nrow(pa_dat), "| B donors", uniqueN(pb_dat[, paste(cohort, donor_id)]),
    "| C effects", nrow(pc_dat), "| D models", nrow(pd_dat), "| E rows", nrow(pe_dat), "\n")
cat("DONE\n")
