#!/usr/bin/env Rscript
# Validate that atlas cell identifiers retain the original author labels.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_arg) == 1L) {
  dirname(dirname(dirname(normalizePath(sub("^--file=", "", script_arg)))))
} else getwd()
setwd(project_dir)

data_dir <- file.path(project_dir, "..", "AKI_DCT_analysis", "data")
out_dir <- file.path(project_dir, "results", "extension")

u183 <- read.csv(file.path(out_dir, "atlas_GSE183276_umap.csv"), stringsAsFactors = FALSE)
m183 <- read.delim(gzfile(file.path(data_dir, "GSE183276_Metadata.txt.gz")),
                   row.names = 1, check.names = FALSE, quote = "\"",
                   stringsAsFactors = FALSE)
i183 <- match(u183$cell, rownames(m183))
ok183 <- !is.na(i183) & u183$author_label == m183$subclass.l1[i183]

con <- gzfile(file.path(data_dir, "GSE210622_metadata.txt.gz"), "r")
invisible(readLines(con, n = 1)); raw210 <- readLines(con); close(con)
parse_label <- function(line) {
  parts <- strsplit(line, "[[:space:]]+")[[1]]
  n <- length(parts)
  c(barcode = parts[1], author_label = parts[n - 1])
}
m210 <- do.call(rbind, lapply(raw210, parse_label))
u210 <- read.csv(file.path(out_dir, "atlas_GSE210622_umap.csv"), stringsAsFactors = FALSE)
i210 <- match(u210$cell, m210[, "barcode"])
ok210 <- !is.na(i210) & u210$author_label == unname(m210[i210, "author_label"])

receipt <- data.frame(
  cohort = c("GSE183276", "GSE210622"),
  n_atlas_cells = c(nrow(u183), nrow(u210)),
  n_barcode_matched = c(sum(!is.na(i183)), sum(!is.na(i210))),
  n_author_label_matched = c(sum(ok183), sum(ok210)),
  match_fraction = c(mean(ok183), mean(ok210)),
  status = ifelse(c(all(ok183), all(ok210)), "PASS", "FAIL")
)
write.csv(receipt, file.path(out_dir, "atlas_metadata_validation.csv"), row.names = FALSE)
print(receipt)
if (!all(receipt$status == "PASS")) quit(status = 1L)
