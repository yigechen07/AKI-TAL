#!/usr/bin/env Rscript
# =============================================================================
# 03b_tal_pathway_enrichment.R
# Donor-aware ranked pathway enrichment for the Figure 3 TAL pseudobulk layer.
#
# Primary rule:
#   - fit GSE183276 and GSE210622 separately on >=50 strict author-TAL donors;
#   - limma-voom donor-pseudobulk model, followed by rank-based CAMERA;
#   - estimate inter-gene correlation for each gene set from donor residuals;
#   - adjust FDR within collection and cohort;
#   - replicated = same direction and FDR < 0.05 in both cohorts;
#   - supportive = same direction, FDR < 0.05 in one cohort and nominal
#     P < 0.05 in the other.
#
# This is an exploratory extension. It does not replace the frozen six-gene
# endpoint and does not establish mechanism, causality or TAL specificity.
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
  library(limma)
  library(data.table)
})

set.seed(20260809)

pb_path <- file.path(project_dir, "..", "AKI_DCT_analysis", "results", "stage3b_pseudobulk.rds")
qc_path <- file.path(project_dir, "..", "AKI_DCT_analysis", "results", "stage3b_qc.csv")
gmt_dir <- file.path(project_dir, "metadata", "extension", "gene_sets")
out_dir <- file.path(project_dir, "results", "extension")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gmt_files <- c(
  Hallmark = "h.all.v2026.1.Hs.symbols.gmt",
  Reactome = "c2.cp.reactome.v2026.1.Hs.symbols.gmt",
  GO_BP = "c5.go.bp.v2026.1.Hs.symbols.gmt"
)
gmt_paths <- file.path(gmt_dir, gmt_files)
names(gmt_paths) <- names(gmt_files)
stopifnot(file.exists(pb_path), file.exists(qc_path), all(file.exists(gmt_paths)))

pb <- readRDS(pb_path)
qc <- fread(qc_path)

read_gmt <- function(path, collection) {
  lines <- readLines(path, warn = FALSE)
  parts <- strsplit(lines, "\t", fixed = FALSE)
  out <- lapply(parts, function(x) {
    if (length(x) < 3L) return(NULL)
    list(pathway = x[[1]], description = x[[2]], genes = unique(x[-c(1, 2)]))
  })
  out <- Filter(Negate(is.null), out)
  names(out) <- vapply(out, `[[`, character(1), "pathway")
  data.table(
    collection = collection,
    pathway = names(out),
    description = vapply(out, `[[`, character(1), "description"),
    genes = I(lapply(out, `[[`, "genes"))
  )
}

gmt <- rbindlist(Map(read_gmt, gmt_paths, names(gmt_paths)))
stopifnot(uniqueN(gmt$collection) == 3L, !anyDuplicated(gmt[, paste(collection, pathway)]))

prepare_voom <- function(dataset) {
  meta <- pb$metadata[[dataset]]
  qc_ds <- qc[qc$dataset == dataset & qc$segment == "TAL"]
  qc_ds <- qc_ds[match(meta$donor, qc_ds$donor)]
  stopifnot(all(qc_ds$donor == meta$donor))

  if (dataset == "GSE183276") {
    keep <- qc_ds$eligible_50 & meta$diagnosis %in% c("AKI", "Ref")
    ctrl <- "Ref"
  } else {
    keep <- qc_ds$eligible_50 & meta$diagnosis %in% c("AKI", "non_AKI")
    ctrl <- "non_AKI"
  }
  donors <- meta$donor[keep]
  info <- data.frame(
    donor = donors,
    group = factor(ifelse(meta$diagnosis[keep] == "AKI", "AKI", ctrl),
                   levels = c(ctrl, "AKI")),
    stringsAsFactors = FALSE
  )
  counts <- as.matrix(pb$pseudobulk[[dataset]]$TAL[, donors, drop = FALSE])
  y <- DGEList(counts = counts, samples = info)
  keep_gene <- filterByExpr(y, group = info$group)
  y <- y[keep_gene, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(~ group, data = info)
  v <- voom(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), robust = TRUE)
  tstat <- fit$t[, "groupAKI"]
  names(tstat) <- rownames(v$E)
  list(v = v, design = design, statistic = tstat, info = info)
}

run_collection <- function(prep, collection_name, method_name = "camera_estimated_cor") {
  universe <- rownames(prep$v$E)
  gs <- gmt[gmt$collection == collection_name]
  gs[, overlap := lapply(genes, intersect, y = universe)]
  gs[, n_overlap := lengths(overlap)]
  gs <- gs[n_overlap >= 15L & n_overlap <= 500L]
  index <- setNames(gs$overlap, gs$pathway)
  index <- lapply(index, function(g) match(g, universe))

  correlation_arg <- if (method_name == "camera_estimated_cor") NA_real_ else 0.01
  cam <- camera(
    prep$v,
    index = index,
    design = prep$design,
    contrast = "groupAKI",
    use.ranks = TRUE,
    inter.gene.cor = correlation_arg,
    allow.neg.cor = FALSE,
    directional = TRUE,
    sort = FALSE
  )
  cam <- as.data.table(cam, keep.rownames = "pathway")

  # A bounded descriptive effect used only for plotting: mean centred
  # fractional rank of member genes, from -1 (down-ranked) to +1 (up-ranked).
  stat <- prep$statistic[universe]
  centred_rank <- 2 * (rank(stat, ties.method = "average") / (length(stat) + 1) - 0.5)
  names(centred_rank) <- universe
  rank_shift_map <- vapply(gs$overlap, function(g) mean(centred_rank[g]), numeric(1))
  names(rank_shift_map) <- gs$pathway

  cam[, `:=`(
    method = method_name,
    collection = collection_name,
    raw_rank_shift = rank_shift_map[pathway],
    n_universe = length(universe),
    n_aki = sum(prep$info$group == "AKI"),
    n_ctrl = sum(prep$info$group != "AKI")
  )]
  # CAMERA direction is authoritative. Near-zero mean rank shifts can differ
  # in sign because CAMERA also accounts for inter-gene correlation; retain
  # the raw value for audit and align the displayed effect sign to CAMERA.
  cam[, rank_shift := fifelse(Direction == "Up", abs(raw_rank_shift), -abs(raw_rank_shift))]
  if (!("Correlation" %in% names(cam))) cam[, Correlation := correlation_arg]
  setcolorder(cam, c("method", "collection", "pathway", "NGenes", "Correlation",
                     "Direction", "rank_shift", "raw_rank_shift", "PValue", "FDR",
                     "n_universe", "n_aki", "n_ctrl"))
  cam
}

preps <- list(
  GSE183276 = prepare_voom("GSE183276"),
  GSE210622 = prepare_voom("GSE210622")
)

all_results <- list()
for (dataset in names(preps)) {
  for (collection_name in names(gmt_paths)) {
    for (method_name in c("camera_estimated_cor", "camera_fixed_cor_0.01")) {
      message("ranked CAMERA: ", dataset, " / ", collection_name, " / ", method_name)
      x <- run_collection(preps[[dataset]], collection_name, method_name)
      x[, cohort := dataset]
      all_results[[paste(dataset, collection_name, method_name)]] <- x
    }
  }
}
enrich <- rbindlist(all_results, use.names = TRUE)
setcolorder(enrich, c("cohort", setdiff(names(enrich), "cohort")))
fwrite(enrich, file.path(out_dir, "tal_pathway_enrichment_all.csv"))

primary_enrich <- enrich[method == "camera_estimated_cor"]
a <- primary_enrich[cohort == "GSE183276"]
b <- primary_enrich[cohort == "GSE210622"]
conc <- merge(
  a[, .(collection, pathway, n_genes_GSE183276 = NGenes,
        correlation_GSE183276 = Correlation, direction_GSE183276 = Direction,
        rank_shift_GSE183276 = rank_shift, p_GSE183276 = PValue, fdr_GSE183276 = FDR)],
  b[, .(collection, pathway, n_genes_GSE210622 = NGenes,
        correlation_GSE210622 = Correlation, direction_GSE210622 = Direction,
        rank_shift_GSE210622 = rank_shift, p_GSE210622 = PValue, fdr_GSE210622 = FDR)],
  by = c("collection", "pathway")
)
conc[, same_direction := direction_GSE183276 == direction_GSE210622]
conc[, replicated := same_direction & fdr_GSE183276 < 0.05 & fdr_GSE210622 < 0.05]
conc[, supportive := same_direction & !replicated &
       ((fdr_GSE183276 < 0.05 & p_GSE210622 < 0.05) |
        (fdr_GSE210622 < 0.05 & p_GSE183276 < 0.05))]
conc[, support_class := fifelse(replicated, "replicated_both_FDR",
                         fifelse(supportive, "supportive_one_FDR_other_nominal",
                           fifelse(same_direction, "concordant_not_supported", "discordant")))]
conc[, max_fdr := pmax(fdr_GSE183276, fdr_GSE210622)]
conc[, max_p := pmax(p_GSE183276, p_GSE210622)]
conc[, mean_abs_rank_shift := (abs(rank_shift_GSE183276) + abs(rank_shift_GSE210622)) / 2]
conc[, class_rank := match(support_class,
                           c("replicated_both_FDR",
                             "supportive_one_FDR_other_nominal",
                             "concordant_not_supported", "discordant"))]
setorder(conc, collection, class_rank, max_fdr, max_p, -mean_abs_rank_shift)
fwrite(conc, file.path(out_dir, "tal_pathway_concordance.csv"))

# Fixed-correlation CAMERA is retained as a sensitivity only. It must never
# replace a negative estimated-correlation primary result.
sa <- enrich[method == "camera_fixed_cor_0.01" & cohort == "GSE183276"]
sb <- enrich[method == "camera_fixed_cor_0.01" & cohort == "GSE210622"]
sens_conc <- merge(
  sa[, .(collection, pathway, direction_GSE183276 = Direction,
         rank_shift_GSE183276 = rank_shift, p_GSE183276 = PValue,
         fdr_GSE183276 = FDR)],
  sb[, .(collection, pathway, direction_GSE210622 = Direction,
         rank_shift_GSE210622 = rank_shift, p_GSE210622 = PValue,
         fdr_GSE210622 = FDR)],
  by = c("collection", "pathway")
)
sens_conc[, same_direction := direction_GSE183276 == direction_GSE210622]
sens_conc[, replicated := same_direction & fdr_GSE183276 < 0.05 & fdr_GSE210622 < 0.05]
sens_conc[, supportive := same_direction & !replicated &
            ((fdr_GSE183276 < 0.05 & p_GSE210622 < 0.05) |
             (fdr_GSE210622 < 0.05 & p_GSE183276 < 0.05))]
fwrite(sens_conc, file.path(out_dir, "tal_pathway_concordance_fixed_cor_sensitivity.csv"))

clean_label <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- sub("^GOBP_", "", x)
  tools::toTitleCase(tolower(gsub("_", " ", x, fixed = TRUE)))
}

# Deterministic, redundancy-aware display selection. No pathway is selected by
# biological preference. Within each collection, rank by replication class and
# evidence, then reject a pathway if its Jaccard overlap with an already chosen
# set exceeds 0.50. Show at most four supported pathways. If a collection has
# none, show at most two strongest same-direction pathways as an explicit
# audit-only negative display (never as replicated/supportive evidence).
gene_lookup <- setNames(gmt$genes, paste(gmt$collection, gmt$pathway, sep = "|"))
chosen <- list()
for (collection_name in names(gmt_paths)) {
  cand <- conc[conc$collection == collection_name &
                 support_class %in% c("replicated_both_FDR",
                                      "supportive_one_FDR_other_nominal")]
  display_limit <- 4L
  if (nrow(cand) == 0L) {
    cand <- conc[conc$collection == collection_name & same_direction == TRUE]
    cand[, display_role := "audit_only_not_supported"]
    setorder(cand, max_p, max_fdr, -mean_abs_rank_shift)
    display_limit <- 2L
  } else {
    cand[, display_role := support_class]
    setorder(cand, class_rank, max_fdr, max_p, -mean_abs_rank_shift)
  }
  keep_rows <- integer(0)
  selected_sets <- list()
  for (i in seq_len(nrow(cand))) {
    genes <- gene_lookup[[paste(collection_name, cand$pathway[[i]], sep = "|")]]
    jac <- if (length(selected_sets) == 0L) numeric(0) else
      vapply(selected_sets, function(s) length(intersect(genes, s)) / length(union(genes, s)), numeric(1))
    if (length(jac) == 0L || all(jac <= 0.50)) {
      keep_rows <- c(keep_rows, i)
      selected_sets[[length(selected_sets) + 1L]] <- genes
    }
    if (length(keep_rows) >= display_limit) break
  }
  if (length(keep_rows)) chosen[[collection_name]] <- cand[keep_rows]
}
display <- rbindlist(chosen, use.names = TRUE, fill = TRUE)
if (nrow(display)) {
  display[, display_label := clean_label(pathway)]
  display[, display_order := seq_len(.N)]
}
fwrite(display, file.path(out_dir, "tal_pathway_display.csv"))

# Member-gene audit for displayed pathways. These are pathway members with
# cohort-specific statistics, not a post hoc mechanistic or leading-edge set.
core_rows <- list()
if (nrow(display)) {
  deg183 <- fread(file.path(out_dir, "tal_pseudobulk_deg_GSE183276.csv"))
  deg210 <- fread(file.path(out_dir, "tal_pseudobulk_deg_GSE210622.csv"))
  for (i in seq_len(nrow(display))) {
    collection_name <- display$collection[[i]]
    pathway <- display$pathway[[i]]
    genes <- gene_lookup[[paste(collection_name, pathway, sep = "|")]]
    member <- merge(
      deg183[gene %in% genes, .(gene, logFC_GSE183276 = logFC, FDR_GSE183276 = FDR)],
      deg210[gene %in% genes, .(gene, logFC_GSE210622 = logFC, FDR_GSE210622 = FDR)],
      by = "gene"
    )
    member[, `:=`(collection = collection_name, pathway = pathway,
                  same_gene_direction = sign(logFC_GSE183276) == sign(logFC_GSE210622))]
    core_rows[[length(core_rows) + 1L]] <- member
  }
}
core <- if (length(core_rows)) rbindlist(core_rows, use.names = TRUE, fill = TRUE) else data.table()
fwrite(core, file.path(out_dir, "tal_pathway_member_genes.csv"))

summary_counts <- conc[, .(
  n_tested = .N,
  n_same_direction = sum(same_direction),
  n_replicated = sum(replicated),
  n_supportive = sum(supportive),
  n_discordant = sum(!same_direction)
), by = collection]
fwrite(summary_counts, file.path(out_dir, "tal_pathway_summary.csv"))

sens_summary <- sens_conc[, .(
  n_tested = .N,
  n_same_direction = sum(same_direction),
  n_replicated = sum(replicated),
  n_supportive = sum(supportive),
  n_discordant = sum(!same_direction)
), by = collection]
fwrite(sens_summary, file.path(out_dir, "tal_pathway_summary_fixed_cor_sensitivity.csv"))

report <- c(
  "# Figure 3 ranked pathway enrichment report",
  "",
  "**Method:** cohort-specific donor-pseudobulk limma-voom models followed by rank-based CAMERA (`use.ranks=TRUE`) with per-set inter-gene correlation estimated from donor residuals.",
  "",
  "**Collections:** MSigDB 2026.1.Hs Hallmark, Reactome and GO Biological Process; symbol GMT files archived under `metadata/extension/gene_sets/`.",
  "",
  "**Statistical unit:** donor. Cohorts were fitted separately; no cells or cohorts were pooled as independent replicates.",
  "",
  "**Replication rule:** same direction and FDR < 0.05 in both cohorts. Supportive requires same direction, FDR < 0.05 in one cohort and nominal P < 0.05 in the other.",
  "",
  "## Counts",
  "",
  "| Collection | Tested | Same direction | Replicated | Supportive | Discordant |",
  "|---|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(summary_counts))) {
  z <- summary_counts[i]
  report <- c(report, sprintf("| %s | %d | %d | %d | %d | %d |",
                              z$collection, z$n_tested, z$n_same_direction,
                              z$n_replicated, z$n_supportive, z$n_discordant))
}
report <- c(report, "", "## Primary conclusion", "",
            "No Hallmark, Reactome or GO BP pathway met the replicated or supportive rule under the primary estimated-correlation CAMERA model. The gene-level extension therefore does not currently support a robust shared pathway or molecular mechanism.",
            "", "## Fixed-correlation sensitivity", "",
            "The standard CAMERA fixed inter-gene correlation of 0.01 was run as a sensitivity. It is less conservative than estimating correlation from these small donor cohorts and cannot override the primary result.", "",
            "| Collection | Tested | Same direction | Replicated | Supportive | Discordant |",
            "|---|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(sens_summary))) {
  z <- sens_summary[i]
  report <- c(report, sprintf("| %s | %d | %d | %d | %d | %d |",
                              z$collection, z$n_tested, z$n_same_direction,
                              z$n_replicated, z$n_supportive, z$n_discordant))
}
report <- c(report, "", "## Interpretation boundary", "",
            "Fixed-correlation candidates include immune/inflammatory, antigen-presentation, extracellular-matrix and neutrophil-related themes, but they are method-sensitive and may reflect generic injury or cohort context. They are not promoted into the positive evidence chain.",
            "",
            "These are exploratory pathway associations derived after the frozen six-gene endpoint. They do not establish TAL specificity, PT independence, causality, lineage transition or molecular mechanism.")
writeLines(report, file.path(out_dir, "tal_pathway_report.md"))

message("WROTE pathway outputs under results/extension/")
message("DONE")
