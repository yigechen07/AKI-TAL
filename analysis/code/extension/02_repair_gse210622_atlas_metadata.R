#!/usr/bin/env Rscript
# Repair the original GSE210622 atlas metadata row-order bug without recomputing
# the expression embedding. Author labels are rematched by exact cell barcode.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_arg) == 1L) {
  dirname(dirname(dirname(normalizePath(sub("^--file=", "", script_arg)))))
} else getwd()
setwd(project_dir)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

data_dir <- file.path(project_dir, "..", "AKI_DCT_analysis", "data")
dct_res <- file.path(project_dir, "..", "AKI_DCT_analysis", "results")
out_dir <- file.path(project_dir, "results", "extension")
rds_path <- file.path(out_dir, "atlas_GSE210622_seurat.rds")

map210 <- data.frame(
  author_label = c("PT","TAL","tL","DCT","CNT","CD-PC","CD-IC-A","CD-IC-B",
                   "Podo","EC","IntC","Leuk"),
  broad_class = c("PT","TAL/LOH","TAL/LOH","DCT","CNT","CD","CD","CD",
                  "POD","EC","Stroma","IMM"),
  stringsAsFactors = FALSE
)
marker_genes <- c("LRP2","SLC34A1","SLC5A2","SLC12A1","UMOD","CLDN16","KCNJ10",
                  "PTH1R","CLDN10","SLC12A3","PVALB","CALB1","SCNN1G","ATP6V1B1",
                  "SLC4A1","NPHS1","NPHS2","PECAM1","EMCN","DCN","COL1A1","PTPRC")

con <- gzfile(file.path(data_dir, "GSE210622_metadata.txt.gz"), "r")
invisible(readLines(con, n = 1)); raw_lines <- readLines(con); close(con)
parse_line <- function(line) {
  parts <- strsplit(line, "[[:space:]]+")[[1]]
  n <- length(parts)
  data.frame(
    barcode = parts[1],
    orig_ident = paste(parts[2:(n - 5)], collapse = " "),
    author_label = parts[n - 1],
    stringsAsFactors = FALSE
  )
}
raw_meta <- do.call(rbind, lapply(raw_lines, parse_line))
stopifnot(nrow(raw_meta) == 106971L, !anyDuplicated(raw_meta$barcode))

audit <- read.csv(file.path(dct_res, "stage2_gse210622_donor_audit.csv"),
                  stringsAsFactors = FALSE)
donor_rows <- list()
for (i in seq_len(nrow(audit))) {
  samples <- strsplit(audit$sample_ids[i], ";")[[1]]
  for (sample in samples) {
    donor_rows[[length(donor_rows) + 1L]] <- data.frame(
      orig_ident = sample, donor = audit$donor_id[i],
      diagnosis = ifelse(audit$AKI_status[i] == "AKI", "AKI", "non_AKI"),
      stringsAsFactors = FALSE
    )
  }
}
donor_map <- do.call(rbind, donor_rows)
raw_meta$donor <- donor_map$donor[match(raw_meta$orig_ident, donor_map$orig_ident)]
raw_meta$diagnosis <- donor_map$diagnosis[match(raw_meta$orig_ident, donor_map$orig_ident)]
raw_meta$broad_class <- map210$broad_class[match(raw_meta$author_label, map210$author_label)]
stopifnot(!anyNA(raw_meta$donor), !anyNA(raw_meta$broad_class))

message("Loading existing Seurat object (UMAP coordinates retained)...")
obj <- readRDS(rds_path)
idx <- match(colnames(obj), raw_meta$barcode)
stopifnot(!anyNA(idx))
fixed <- raw_meta[idx, , drop = FALSE]
stopifnot(identical(colnames(obj), fixed$barcode))
obj$cohort <- "GSE210622"
obj$donor <- fixed$donor
obj$diagnosis <- fixed$diagnosis
obj$author_label <- fixed$author_label
obj$broad_class <- fixed$broad_class
obj$orig_ident <- fixed$orig_ident

umap <- data.frame(
  cell = colnames(obj),
  umap1 = Embeddings(obj, "umap")[, 1],
  umap2 = Embeddings(obj, "umap")[, 2],
  cohort = "GSE210622", donor = obj$donor, diagnosis = obj$diagnosis,
  author_label = obj$author_label, broad_class = obj$broad_class,
  stringsAsFactors = FALSE
)
write.csv(umap, file.path(out_dir, "atlas_GSE210622_umap.csv"), row.names = FALSE)

counts_df <- as.data.frame(table(umap$donor, umap$broad_class))
colnames(counts_df) <- c("donor", "broad_class", "n_cells")
counts_df$cohort <- "GSE210622"
diag_df <- unique(umap[, c("donor", "diagnosis")])
counts_df <- merge(counts_df, diag_df, by = "donor", all.x = TRUE)
write.csv(counts_df, file.path(out_dir, "atlas_GSE210622_counts.csv"), row.names = FALSE)

expr <- GetAssayData(obj, layer = "data")
found <- intersect(marker_genes, rownames(expr))
marker_rows <- list()
for (cl in sort(unique(obj$broad_class))) {
  cells <- which(obj$broad_class == cl)
  for (gene in found) {
    values <- expr[gene, cells]
    marker_rows[[length(marker_rows) + 1L]] <- data.frame(
      broad_class = cl, gene = gene, n_cells = length(cells),
      mean_expr = as.numeric(mean(values)), pct_expr = as.numeric(mean(values > 0)),
      cohort = "GSE210622", stringsAsFactors = FALSE
    )
  }
}
write.csv(do.call(rbind, marker_rows),
          file.path(out_dir, "atlas_GSE210622_markers.csv"), row.names = FALSE)

tal <- umap[umap$broad_class == "TAL/LOH", ]
tal_donor <- aggregate(cell ~ donor, data = tal, FUN = length)
summary_df <- data.frame(
  cohort = "GSE210622", assay = "snRNA-seq", n_cells = nrow(umap),
  n_donors = length(unique(umap$donor)), n_TAL_LOH_cells = nrow(tal),
  n_donors_TAL_ge20 = sum(tal_donor$cell >= 20),
  n_donors_TAL_ge50 = sum(tal_donor$cell >= 50),
  n_AKI = sum(diag_df$diagnosis == "AKI"),
  n_nonAKI = sum(diag_df$diagnosis == "non_AKI"),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(out_dir, "atlas_GSE210622_summary.csv"), row.names = FALSE)

message("Saving barcode-corrected Seurat metadata...")
saveRDS(obj, rds_path)
message("Repair complete; UMAP coordinates unchanged, labels and summaries rematched by barcode.")
