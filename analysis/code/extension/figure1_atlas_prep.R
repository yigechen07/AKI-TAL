#!/usr/bin/env Rscript
# =============================================================================
# figure1_atlas_prep.R — Figure 1 (extension): cohort-specific cell atlas build
# =============================================================================
# Scope (analysis_extension_plan.md WP1):
#   * Reconstruct GSE183276 (scRNA) and GSE210622 (snRNA) UMAP/annotation panels
#     separately. No joint scRNA-snRNA embedding.
#   * Author labels are the primary annotation source; they are harmonised to
#     broad kidney classes ONLY for cross-cohort presentation, and the exact
#     author-label -> broad-class mapping is exported.
#   * UMAPs are descriptive. No cell-level P values, no abundance claims.
#
# Statistical unit for all downstream inference remains the DONOR; this script
# only produces descriptive atlas summaries.
#
# Inputs (from ../AKI_DCT_analysis/data):
#   GSE183276_Counts_decomp.RDS, GSE183276_Metadata.txt.gz
#   GSE210622_RAW.tar, GSE210622_metadata.txt.gz
#   ../AKI_DCT_analysis/results/stage2_gse210622_donor_audit.csv (donor map)
# Barcode logic mirrors ../AKI_DCT_analysis/code/stage3b_build_pseudobulk.R
#   (last-18-character barcode matching for GSE210622).
#
# Outputs (results/extension/):
#   atlas_GSE183276_umap.csv, atlas_GSE183276_counts.csv,
#   atlas_GSE183276_markers.csv, atlas_GSE183276_summary.csv
#   atlas_GSE210622_umap.csv, atlas_GSE210622_counts.csv,
#   atlas_GSE210622_markers.csv, atlas_GSE210622_summary.csv
#   atlas_broad_mapping.csv, atlas_cohort_overview.csv
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
  library(Matrix)
  library(Seurat)
})

data_dir  <- file.path(project_dir, "..", "AKI_DCT_analysis", "data")
aki_res   <- file.path(project_dir, "..", "AKI_DCT_analysis", "results")
out_dir   <- file.path(project_dir, "results", "extension")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"))

# ---------------------------------------------------------------------------
# 0. Broad-class harmonisation (author label -> broad class), fixed palette
# ---------------------------------------------------------------------------
broad_map <- function(cohort) {
  if (cohort == "GSE183276") {
    data.frame(
      author_label = c("PT","TAL","DTL","ATL","DCT","CNT","PC","IC","POD",
                       "EC","VSM/P","FIB","IMM","PEC","PapE","NEU"),
      broad_class  = c("PT","TAL/LOH","TAL/LOH","TAL/LOH","DCT","CNT","CD","CD","POD",
                       "EC","Stroma","Stroma","IMM","Other","Other","Other"),
      stringsAsFactors = FALSE
    )
  } else if (cohort == "GSE210622") {
    data.frame(
      author_label = c("PT","TAL","tL","DCT","CNT","CD-PC","CD-IC-A","CD-IC-B",
                       "Podo","EC","IntC","Leuk"),
      broad_class  = c("PT","TAL/LOH","TAL/LOH","DCT","CNT","CD","CD","CD",
                       "POD","EC","Stroma","IMM"),
      stringsAsFactors = FALSE
    )
  }
}

# Fixed broad-class palette used in every extension figure (visual system).
broad_palette <- c(
  "PT"      = "#1B9E77",
  "TAL/LOH" = "#D95F02",
  "DCT"     = "#7570B3",
  "CNT"     = "#E7298A",
  "CD"      = "#66A61E",
  "POD"     = "#E6AB02",
  "EC"      = "#1F78B4",
  "Stroma"  = "#A6761D",
  "IMM"     = "#666666",
  "Other"   = "#CCCCCC"
)

# Prespecified kidney marker panel (rows of the annotation dot plot).
marker_panel <- data.frame(
  gene = c("LRP2","SLC34A1","SLC5A2",
           "SLC12A1","UMOD","CLDN16","KCNJ10","PTH1R","CLDN10",
           "SLC12A3","PVALB",
           "CALB1","SCNN1G",
           "ATP6V1B1","SLC4A1",
           "NPHS1","NPHS2",
           "PECAM1","EMCN",
           "DCN","COL1A1",
           "PTPRC"),
  group = c(rep("PT",3), rep("TAL",6), rep("DCT",2), rep("CNT",2),
            rep("CD/IC",2), rep("POD",2), rep("EC",2), rep("Stroma",2),
            "IMM"),
  stringsAsFactors = FALSE
)
marker_genes <- marker_panel$gene

# ---------------------------------------------------------------------------
# 1. Seurat atlas function (common pipeline)
# ---------------------------------------------------------------------------
run_seurat <- function(counts, meta_df, seed = 42L) {
  obj <- CreateSeuratObject(counts = counts, meta.data = meta_df,
                            min.cells = 0, min.features = 0)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 2000, verbose = FALSE)
  # Scale ONLY the 2000 HVGs (Seurat v5 would otherwise materialise a dense
  # 33k x 107k scale.data layer, which blows the 16 GB R vector limit).
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, features = VariableFeatures(obj), verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, k.param = 20, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, n.neighbors = 20, min.dist = 0.3,
                 seed.use = seed, verbose = FALSE)
  obj
}

# Marker summary (mean log-norm expression and detection fraction per class)
marker_summary <- function(obj, class_col = "broad_class") {
  classes <- sort(unique(obj@meta.data[[class_col]]))
  expr <- GetAssayData(obj, layer = "data")
  found <- intersect(marker_genes, rownames(expr))
  if (length(found) < length(marker_genes)) {
    warning("Markers missing from this cohort: ",
            paste(setdiff(marker_genes, found), collapse = ", "))
  }
  rows <- list()
  for (cl in classes) {
    cells <- which(obj@meta.data[[class_col]] == cl)
    for (g in found) {
      v <- expr[g, cells]
      rows[[length(rows) + 1]] <- data.frame(
        broad_class = cl, gene = g, n_cells = length(cells),
        mean_expr = as.numeric(mean(v)), pct_expr = as.numeric(mean(v > 0)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# 2. GSE183276 (scRNA, 45 donors)
# ---------------------------------------------------------------------------
force_atlas_rebuild <- identical(Sys.getenv("FORCE_ATLAS_REBUILD", unset = "0"), "1")
resume_183 <- !force_atlas_rebuild && all(file.exists(file.path(out_dir, c(
  "atlas_GSE183276_umap.csv", "atlas_GSE183276_counts.csv",
  "atlas_GSE183276_markers.csv", "atlas_GSE183276_summary.csv",
  "atlas_GSE183276_seurat.rds"))))
if (resume_183) {
  msg("=== GSE183276 outputs present; skipping recomputation ===")
  sum183 <- read.csv(file.path(out_dir, "atlas_GSE183276_summary.csv"),
                     stringsAsFactors = FALSE)
} else {
  msg("=== GSE183276 ===")
meta183 <- read.delim(gzfile(file.path(data_dir, "GSE183276_Metadata.txt.gz")),
                      row.names = 1, check.names = FALSE, quote = "\"",
                      stringsAsFactors = FALSE)
counts183 <- readRDS(file.path(data_dir, "GSE183276_Counts_decomp.RDS"))
if (!inherits(counts183, "sparseMatrix")) counts183 <- as(counts183, "dgCMatrix")
counts183 <- drop0(counts183)
stopifnot(all(colnames(counts183) %in% rownames(meta183)))
meta183 <- meta183[colnames(counts183), , drop = FALSE]
msg(sprintf("cells=%d donors=%d genes=%d", ncol(counts183),
            length(unique(meta183$patient)), nrow(counts183)))

map183 <- broad_map("GSE183276")
meta183$broad_class <- map183$broad_class[match(meta183$subclass.l1, map183$author_label)]
meta183$author_label <- meta183$subclass.l1
meta183$cohort <- "GSE183276"
meta183$diagnosis <- meta183$condition.l1
stopifnot(!any(is.na(meta183$broad_class)))

meta_sub <- meta183[, c("cohort","patient","diagnosis","author_label","broad_class",
                        "condition.l2","condition.l3","state","tissue_type","assay",
                        "percent.mt","nCount_RNA","nFeature_RNA")]
obj183 <- run_seurat(counts183, meta_sub)
umap183 <- data.frame(
  cell = rownames(obj183@meta.data),
  umap1 = obj183@reductions$umap@cell.embeddings[, 1],
  umap2 = obj183@reductions$umap@cell.embeddings[, 2],
  obj183@meta.data[, c("cohort","patient","diagnosis","author_label","broad_class")],
  stringsAsFactors = FALSE, check.names = FALSE
)
# column order: 1 cell, 2 umap1, 3 umap2, 4 cohort, 5 patient, 6 diagnosis, ...
colnames(umap183)[5] <- "donor"
write.csv(umap183, file.path(out_dir, "atlas_GSE183276_umap.csv"), row.names = FALSE)

# Save the Seurat object (drop scale.data layer to keep the file lean; counts,
# log-normalised data, metadata and UMAP coordinates are retained)
if ("scale.data" %in% names(obj183[["RNA"]]@layers)) {
  obj183[["RNA"]]@layers[["scale.data"]] <- NULL
}
saveRDS(obj183, file.path(out_dir, "atlas_GSE183276_seurat.rds"))
msg(sprintf("saved atlas_GSE183276_seurat.rds (%.1f MB)",
            file.info(file.path(out_dir, "atlas_GSE183276_seurat.rds"))$size / 1e6))

counts183_df <- as.data.frame(table(umap183$donor, umap183$broad_class))
colnames(counts183_df) <- c("donor","broad_class","n_cells")
counts183_df$cohort <- "GSE183276"
diag183 <- unique(umap183[, c("donor","diagnosis")])
counts183_df <- merge(counts183_df, diag183, by = "donor", all.x = TRUE)
write.csv(counts183_df, file.path(out_dir, "atlas_GSE183276_counts.csv"), row.names = FALSE)

mk183 <- marker_summary(obj183)
mk183$cohort <- "GSE183276"
write.csv(mk183, file.path(out_dir, "atlas_GSE183276_markers.csv"), row.names = FALSE)

tal183 <- umap183[umap183$broad_class == "TAL/LOH", ]
tal_donor183 <- aggregate(cell ~ donor, data = tal183, FUN = length)
sum183 <- data.frame(
  cohort = "GSE183276", assay = "scRNA-seq", n_cells = nrow(umap183),
  n_donors = length(unique(umap183$donor)),
  n_TAL_LOH_cells = nrow(tal183),
  n_donors_TAL_ge20 = sum(tal_donor183$cell >= 20),
  n_donors_TAL_ge50 = sum(tal_donor183$cell >= 50),
  n_AKI = sum(diag183$diagnosis == "AKI"),
  n_Ref = sum(diag183$diagnosis == "Ref"),
  n_CKD = sum(diag183$diagnosis == "CKD"),
  stringsAsFactors = FALSE
)
write.csv(sum183, file.path(out_dir, "atlas_GSE183276_summary.csv"), row.names = FALSE)

rm(counts183, obj183, meta183); gc()
}  # end GSE183276 (or resume)

# ---------------------------------------------------------------------------
# 3. GSE210622 (snRNA, 12 effective donors / 14 samples)
# ---------------------------------------------------------------------------
msg("=== GSE210622 ===")

# 3a. Metadata (space-separated; orig.ident may contain spaces, e.g. "AKI 1")
con <- gzfile(file.path(data_dir, "GSE210622_metadata.txt.gz"), "r")
hdr <- readLines(con, n = 1); ml_raw <- readLines(con); close(con)
parse_line <- function(line) {
  parts <- strsplit(line, "\\s+")[[1]]
  n <- length(parts)
  if (n < 7) return(NULL)
  # token layout (R 1-indexed): 1 barcode, 2..(n-5) orig.ident (may span
  # multiple tokens, e.g. "AKI 1"), then exactly 3 QC fields, then
  # celltype and group as the last two tokens. Anchor on the tail so
  # multi-token idents do not shift the QC columns.
  data.frame(
    barcode    = parts[1],
    orig_ident = paste(parts[2:(n - 5)], collapse = " "),
    nCount_RNA = as.numeric(parts[n - 4]),
    nFeature   = as.numeric(parts[n - 3]),
    percent.mt = as.numeric(parts[n - 2]),
    celltype   = parts[n - 1],
    author_group = parts[n],
    stringsAsFactors = FALSE
  )
}
meta210 <- do.call(rbind, lapply(ml_raw, parse_line))
stopifnot(nrow(meta210) == 106971)
msg(sprintf("metadata cells=%d", nrow(meta210)))

# 3b. Donor map from the stage-2 audit
audit <- read.csv(file.path(aki_res, "stage2_gse210622_donor_audit.csv"),
                  stringsAsFactors = FALSE)
donor_map <- list()
for (i in seq_len(nrow(audit))) {
  samples <- strsplit(audit$sample_ids[i], ";")[[1]]
  gsms    <- strsplit(audit$GSM_accessions[i], ";")[[1]]
  for (j in seq_along(samples)) {
    donor_map[[samples[j]]] <- data.frame(
      sample = samples[j], gsm = gsms[j], donor = audit$donor_id[i],
      AKI_status = audit$AKI_status[i], tissue_source = audit$tissue_source[i],
      postmortem_status = audit$postmortem_status[i],
      COVID_status = audit$COVID_status[i], stringsAsFactors = FALSE)
  }
}
donor_df <- do.call(rbind, donor_map)
stopifnot(nrow(donor_df) == 14)
meta210 <- merge(meta210, donor_df, by.x = "orig_ident", by.y = "sample", all.x = TRUE)
stopifnot(!any(is.na(meta210$donor)))

map210 <- broad_map("GSE210622")
meta210$broad_class <- map210$broad_class[match(meta210$celltype, map210$author_label)]
meta210$author_label <- meta210$celltype
meta210$cohort <- "GSE210622"
meta210$diagnosis <- ifelse(meta210$AKI_status == "AKI", "AKI", "non_AKI")
stopifnot(!any(is.na(meta210$broad_class)))

# 3c. Extract RAW tar and read per-sample matrices (stage3b barcode logic)
tmp_dir <- tempfile("fig1_gse210622_")
dir.create(tmp_dir)
untar(file.path(data_dir, "GSE210622_RAW.tar"), exdir = tmp_dir)
gz_files <- list.files(tmp_dir, pattern = "\\.gz$")
stopifnot(length(gz_files) == 42)

read_sample <- function(gsm, sample_name) {
  stub <- gsub(" ", "", sample_name)
  ft  <- read.table(gzfile(file.path(tmp_dir, sprintf("%s_%s.features.tsv.gz", gsm, stub))),
                    sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
  bc  <- readLines(gzfile(file.path(tmp_dir, sprintf("%s_%s.barcodes.tsv.gz", gsm, stub))))
  mat <- readMM(gzfile(file.path(tmp_dir, sprintf("%s_%s.matrix.mtx.gz", gsm, stub))))
  mat <- as(mat, "dgCMatrix")
  list(genes = ft[, 2], barcodes = bc, mat = mat)
}

# Gene universe from first sample; verify identical across samples
s1 <- read_sample(donor_df$gsm[1], donor_df$sample[1])
gene_symbols <- s1$genes
for (i in seq_len(nrow(donor_df))) {
  s <- read_sample(donor_df$gsm[i], donor_df$sample[i])
  stopifnot(identical(s$genes, gene_symbols))
  sm_meta <- meta210[meta210$orig_ident == donor_df$sample[i], ]
  sm_bc <- substr(sm_meta$barcode, nchar(sm_meta$barcode) - 17, nchar(sm_meta$barcode))
  stopifnot(sum(!is.na(match(sm_bc, s$barcodes))) == nrow(sm_meta))
}

# Collapse duplicate gene symbols by sum (deterministic; matches stage3b)
sym_to_rows <- split(seq_along(gene_symbols), gene_symbols)
unique_syms <- names(sym_to_rows)
# rowsum(..., reorder=TRUE) orders rows by sorted unique group levels, which
# must equal the sorted unique-symbol list used for rownames.
expected_genes <- sort(unique(gene_symbols))
stopifnot(identical(expected_genes, sort(unique_syms)))
msg(sprintf("unique symbols=%d (from %d features)", length(unique_syms), length(gene_symbols)))

cell_list <- list()
# One sparse collapse matrix C (unique_genes x features) built once: C[g, f] = 1
# iff feature f maps to unique gene g. Then per sample: C %*% sub gives the
# collapsed gene x cell dgCMatrix directly (handles any duplicate multiplicity,
# no sapply list drift, no rowsum export dependency).
n_feat <- length(gene_symbols)
g_levels <- unique(gene_symbols)             # order of appearance
g_idx <- match(gene_symbols, g_levels)       # 1..n_unique per feature
collapse_C <- sparseMatrix(i = g_idx, j = seq_len(n_feat), x = rep(1, n_feat),
                           dims = c(length(g_levels), n_feat))
rownames(collapse_C) <- g_levels
stopifnot(identical(sort(g_levels), sort(unique_syms)))
for (i in seq_len(nrow(donor_df))) {
  gsm <- donor_df$gsm[i]; sample_name <- donor_df$sample[i]
  s <- read_sample(gsm, sample_name)
  sm_meta <- meta210[meta210$orig_ident == sample_name, ]
  sm_bc <- substr(sm_meta$barcode, nchar(sm_meta$barcode) - 17, nchar(sm_meta$barcode))
  col_idx <- match(sm_bc, s$barcodes)
  sub <- s$mat[, col_idx, drop = FALSE]
  colnames(sub) <- sm_meta$barcode
  sub2 <- collapse_C %*% sub                 # unique_genes x cells, dgCMatrix
  cell_list[[i]] <- sub2
  msg(sprintf("  %-22s -> %-4s cells=%d", sample_name, donor_df$donor[i], ncol(sub2)))
}
counts210 <- do.call(cbind, cell_list)
stopifnot(ncol(counts210) == nrow(meta210))
unlink(tmp_dir, recursive = TRUE)
msg(sprintf("combined cells=%d donors=%d genes=%d", ncol(counts210),
            length(unique(meta210$donor)), nrow(counts210)))

# Release extraction/collapse intermediates before Seurat duplicates assay
# bookkeeping. This materially lowers the peak memory needed for a full rebuild.
rm(cell_list, collapse_C, s1, s, sub, sub2, sm_meta, sm_bc, col_idx)
gc()

meta_sub210 <- meta210[, c("cohort","donor","diagnosis","author_label","broad_class",
                           "orig_ident","tissue_source","postmortem_status",
                           "COVID_status","AKI_status","percent.mt","nCount_RNA","nFeature")]
colnames(meta_sub210)[colnames(meta_sub210) == "nFeature"] <- "nFeature_RNA"
# Make the metadata-to-count correspondence explicit. Named/reordered rows
# prevent silent drift if merge() or the donor audit order changes.
rownames(meta_sub210) <- meta210$barcode
stopifnot(setequal(rownames(meta_sub210), colnames(counts210)))
meta_sub210 <- meta_sub210[colnames(counts210), , drop = FALSE]
stopifnot(identical(rownames(meta_sub210), colnames(counts210)))
obj210 <- run_seurat(counts210, meta_sub210)
umap210 <- data.frame(
  cell = rownames(obj210@meta.data),
  umap1 = obj210@reductions$umap@cell.embeddings[, 1],
  umap2 = obj210@reductions$umap@cell.embeddings[, 2],
  obj210@meta.data[, c("cohort","donor","diagnosis","author_label","broad_class")],
  stringsAsFactors = FALSE, check.names = FALSE
)
write.csv(umap210, file.path(out_dir, "atlas_GSE210622_umap.csv"), row.names = FALSE)

if ("scale.data" %in% names(obj210[["RNA"]]@layers)) {
  obj210[["RNA"]]@layers[["scale.data"]] <- NULL
}
saveRDS(obj210, file.path(out_dir, "atlas_GSE210622_seurat.rds"))
msg(sprintf("saved atlas_GSE210622_seurat.rds (%.1f MB)",
            file.info(file.path(out_dir, "atlas_GSE210622_seurat.rds"))$size / 1e6))

counts210_df <- as.data.frame(table(umap210$donor, umap210$broad_class))
colnames(counts210_df) <- c("donor","broad_class","n_cells")
counts210_df$cohort <- "GSE210622"
diag210 <- unique(umap210[, c("donor","diagnosis")])
counts210_df <- merge(counts210_df, diag210, by = "donor", all.x = TRUE)
write.csv(counts210_df, file.path(out_dir, "atlas_GSE210622_counts.csv"), row.names = FALSE)

mk210 <- marker_summary(obj210)
mk210$cohort <- "GSE210622"
write.csv(mk210, file.path(out_dir, "atlas_GSE210622_markers.csv"), row.names = FALSE)

tal210 <- umap210[umap210$broad_class == "TAL/LOH", ]
tal_donor210 <- aggregate(cell ~ donor, data = tal210, FUN = length)
sum210 <- data.frame(
  cohort = "GSE210622", assay = "snRNA-seq", n_cells = nrow(umap210),
  n_donors = length(unique(umap210$donor)),
  n_TAL_LOH_cells = nrow(tal210),
  n_donors_TAL_ge20 = sum(tal_donor210$cell >= 20),
  n_donors_TAL_ge50 = sum(tal_donor210$cell >= 50),
  n_AKI = sum(diag210$diagnosis == "AKI"),
  n_nonAKI = sum(diag210$diagnosis == "non_AKI"),
  stringsAsFactors = FALSE
)
write.csv(sum210, file.path(out_dir, "atlas_GSE210622_summary.csv"), row.names = FALSE)

rm(counts210, obj210, meta210); gc()

# ---------------------------------------------------------------------------
# 4. Mapping + cohort overview exports
# ---------------------------------------------------------------------------
mapping_df <- rbind(
  transform(broad_map("GSE183276"), cohort = "GSE183276"),
  transform(broad_map("GSE210622"), cohort = "GSE210622")
)
write.csv(mapping_df, file.path(out_dir, "atlas_broad_mapping.csv"), row.names = FALSE)

overview <- rbind(
  data.frame(cohort = "GSE183276", assay = "scRNA-seq", diagnosis = c("AKI","Ref","CKD"),
             n_donors = c(sum183$n_AKI, sum183$n_Ref, sum183$n_CKD),
             stringsAsFactors = FALSE),
  data.frame(cohort = "GSE210622", assay = "snRNA-seq", diagnosis = c("AKI","non_AKI"),
             n_donors = c(sum210$n_AKI, sum210$n_nonAKI), stringsAsFactors = FALSE)
)
write.csv(overview, file.path(out_dir, "atlas_cohort_overview.csv"), row.names = FALSE)

# Palette CSV for the plot script
write.csv(data.frame(broad_class = names(broad_palette), colour = unname(broad_palette)),
          file.path(out_dir, "atlas_broad_palette.csv"), row.names = FALSE)

msg("=== DONE ===")
msg(paste(capture.output(sessionInfo()), collapse = "\n"))
