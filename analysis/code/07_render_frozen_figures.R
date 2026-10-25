#!/usr/bin/env Rscript
# 07_render_frozen_figures.R — render Figure 4 from frozen allograft-level outputs.
# This script does not reanalyse expression data; run 03_longitudinal_validation.R first
# when the edgeR/limma analysis environment is available.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_arg) == 1L) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
} else {
  normalizePath("projects/AKI_TAL_severity_analysis")
}
result_dir <- file.path(project_dir, "results")
meta_dir <- file.path(project_dir, "metadata")
figure_dir <- file.path(project_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.csv(file.path(meta_dir, "gse126805_samples.csv"), stringsAsFactors = FALSE)
results <- read.csv(file.path(result_dir, "gse126805_results.csv"), stringsAsFactors = FALSE)
samples$timepoint <- as.character(samples$timepoint)

paired <- samples[samples$in_primary_pair %in% TRUE &
                  samples$timepoint %in% c("baseline_pre", "baseline_post"), ]
pre <- paired[paired$timepoint == "baseline_pre", c("kidney", "TAL_score")]
post <- paired[paired$timepoint == "baseline_post", c("kidney", "TAL_score")]
names(pre)[2] <- "pre"; names(post)[2] <- "post"
pairs <- merge(pre, post, by = "kidney", sort = TRUE)
stopifnot(nrow(pairs) == 41L)

value <- function(label, field = "estimate") {
  x <- results[results$label == label, field]
  stopifnot(length(x) == 1L)
  x
}

draw_longitudinal_figure <- function(file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") pdf(file, width = 11, height = 5.4, useDingbats = FALSE)
  else png(file, width = 2200, height = 1080, res = 200)
  on.exit(dev.off(), add = TRUE)
  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.8, 4.8, 3.0, 1.2), las = 1)

  score_range <- range(c(pairs$pre, pairs$post), finite = TRUE)
  pad <- diff(score_range) * 0.10
  if (!is.finite(pad) || pad == 0) pad <- 0.2
  plot(c(1, 2), score_range + c(-pad, pad), type = "n", xaxt = "n",
       xlab = "Protocol-biopsy timepoint", ylab = "Frozen TAL programme score",
       main = "Paired acute change (n = 41 allografts)")
  axis(1, at = c(1, 2), labels = c("Before reperfusion", "After reperfusion"))
  for (i in seq_len(nrow(pairs))) segments(1, pairs$pre[i], 2, pairs$post[i], col = "grey70", lwd = 0.8)
  points(rep(1, nrow(pairs)), pairs$pre, pch = 16, cex = 0.65, col = "grey25")
  points(rep(2, nrow(pairs)), pairs$post, pch = 16, cex = 0.65, col = "grey25")
  segments(0.88, mean(pairs$pre), 1.12, mean(pairs$pre), lwd = 3, col = "#1B6CA8")
  segments(1.88, mean(pairs$post), 2.12, mean(pairs$post), lwd = 3, col = "#1B6CA8")
  mtext(sprintf("Mean paired change = %.3f (95%% CI %.3f to %.3f); P = %.4f",
                value("baseline_post_vs_pre_delta"),
                value("baseline_post_vs_pre_delta", "ci_low"),
                value("baseline_post_vs_pre_delta", "ci_high"),
                value("baseline_post_vs_pre_delta", "p_value")),
        side = 3, line = 0.25, cex = 0.75)

  forest <- data.frame(
    label = c("Post - pre (primary)", "3 months - post (exploratory)",
              "1 year - post (exploratory)", "1 year - pre (exploratory)"),
    estimate = c(value("baseline_post_vs_pre_delta"), value("3mo_vs_post"),
                 value("1y_vs_post"), value("1y_vs_pre")),
    lower = c(value("baseline_post_vs_pre_delta", "ci_low"), value("3mo_vs_post", "ci_low"),
              value("1y_vs_post", "ci_low"), value("1y_vs_pre", "ci_low")),
    upper = c(value("baseline_post_vs_pre_delta", "ci_high"), value("3mo_vs_post", "ci_high"),
              value("1y_vs_post", "ci_high"), value("1y_vs_pre", "ci_high")),
    p_text = c(sprintf("P = %.4f", value("baseline_post_vs_pre_delta", "p_value")),
               sprintf("Holm P = %.4f", value("3mo_vs_post_holm_P")),
               sprintf("Holm P = %.4f", value("1y_vs_post_holm_P")),
               sprintf("Holm P = %.4f", value("1y_vs_pre_holm_P"))),
    stringsAsFactors = FALSE
  )
  yy <- rev(seq_len(nrow(forest)))
  xlim <- range(c(forest$lower, forest$upper, 0)) + c(-0.11, 0.11)
  plot(forest$estimate, yy, xlim = xlim, ylim = c(0.5, 4.5), pch = 16,
       yaxt = "n", xlab = "Mean paired score change", ylab = "",
       main = "Paired allograft contrasts")
  abline(v = 0, lty = 2, col = "grey50")
  segments(forest$lower, yy, forest$upper, yy, lwd = 1.5, col = "grey25")
  points(forest$estimate, yy, pch = 16, cex = 1, col = "#1B6CA8")
  axis(2, at = yy, labels = forest$label, las = 1, cex.axis = 0.78)
  text(x = xlim[2], y = yy, labels = forest$p_text, pos = 2, cex = 0.72)
  mtext("Only post - pre is primary; later contrasts are exploratory.", side = 3, line = 0.25, cex = 0.75)
}

draw_longitudinal_figure(file.path(figure_dir, "longitudinal_gse126805.pdf"), "pdf")
draw_longitudinal_figure(file.path(figure_dir, "longitudinal_gse126805.png"), "png")
message("Rendered Figure 4: figures/longitudinal_gse126805.pdf/png")
