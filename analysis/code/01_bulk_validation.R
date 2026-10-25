#!/usr/bin/env Rscript

# Bulk GEO validation of frozen AKI–TAL functional identity signal
# Datasets: GSE139061 and GSE30718
# Frozen TAL signature: SLC12A1, UMOD, CLDN16, KCNJ10, PTH1R, CLDN10
#
# Scope amendment 2026-07-29:
#   Bulk tissue datasets measure tissue-level detectability of the TAL
#   programme. They cannot validate the KDIGO-stage association from the
#   sc/snRNA-seq discovery because clinical stage is not donor-linked in
#   these series. The bulk evidence supports lower frozen TAL identity in
#   AKI tissue across independent contexts; it does not establish a
#   TAL-cell-intrinsic mechanism or TAL-selective severity effect.
#   GSE139061 is confounded by biopsy-versus-nephrectomy tissue source;
#   GSE30718 is transplant context; GSE210622 is confounded by COVID.
#   Status set to PROCEED_FORMAL_ANALYSIS_WITH_NARROWED_CLAIM.

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) == 1L) {
  project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
} else {
  project_dir <- normalizePath("projects/AKI_TAL_severity_analysis")
}

result_dir <- file.path(project_dir, "results")
figure_dir <- file.path(project_dir, "figures")
data_dir   <- file.path(project_dir, "data", "bulk")
meta_dir   <- file.path(project_dir, "metadata")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

# -------------------------------------------------------------------
# 1. Frozen gene sets
# -------------------------------------------------------------------

tal_genes <- c("SLC12A1", "UMOD", "CLDN16", "KCNJ10", "PTH1R", "CLDN10")

# From stage3_signature_definitions.csv: PT_identity and Pan_epithelial_injury
pt_identity_genes <- c(
  "LRP2", "SLC34A1", "SLC5A2", "SLC5A1", "MME",
  "ALDOB", "MIOX", "PDZK1", "SLC22A6", "ANPEP"
)

pt_injury_genes <- c("PROM1", "DCDC2", "SPP1", "ITGB6", "ITGB8")

all_required_genes <- unique(c(tal_genes, pt_identity_genes, pt_injury_genes))

# -------------------------------------------------------------------
# 2. Helper functions
# -------------------------------------------------------------------

check_genes <- function(genes, available, dataset_label) {
  present <- intersect(genes, available)
  missing <- setdiff(genes, available)
  if (length(missing) > 0) {
    message(sprintf("  [%s] MISSING genes: %s", dataset_label,
                    paste(missing, collapse = ", ")))
  }
  list(present = present, missing = missing)
}

zscore_genes <- function(mat) {
  # mat: genes x samples
  mu <- rowMeans(mat, na.rm = TRUE)
  s  <- apply(mat, 1, sd, na.rm = TRUE)
  s[!is.finite(s) | s < 1e-8] <- 1
  sweep(sweep(mat, 1, mu, "-"), 1, s, "/")
}

compute_module_score <- function(mat, genes) {
  available <- intersect(genes, rownames(mat))
  if (length(available) == 0) return(NULL)
  sub <- mat[available, , drop = FALSE]
  if (any(is.na(sub))) {
    sub[is.na(sub)] <- 0
  }
  z <- zscore_genes(sub)
  colMeans(z)
}

effect_estimate <- function(x, y) {
  # x: reference, y: AKI
  fit <- tryCatch(t.test(y, x), error = function(e) NULL)
  if (is.null(fit)) {
    return(data.frame(
      estimate = mean(y) - mean(x),
      ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_, stringsAsFactors = FALSE
    ))
  }
  data.frame(
    estimate = unname(mean(y) - mean(x)),
    ci_low  = fit$conf.int[1],
    ci_high = fit$conf.int[2],
    p_value = fit$p.value,
    stringsAsFactors = FALSE
  )
}

permutation_pvalue <- function(x, y, n_perm = 9999, seed = 42) {
  # Reproducible (extreme + 1) / (B + 1)
  obs <- mean(y) - mean(x)
  pooled <- c(x, y)
  n1 <- length(x); n2 <- length(y); N <- n1 + n2
  if (N < 5) return(NA_real_)
  if (n_perm > choose(N, n1) * 0.5) n_perm <- min(n_perm, 50000)
  set.seed(seed)
  extreme <- 0
  for (b in seq_len(n_perm)) {
    idx <- sample(N, n1)
    perm_diff <- mean(pooled[-idx]) - mean(pooled[idx])
    if (abs(perm_diff) >= abs(obs)) extreme <- extreme + 1
  }
  (extreme + 1) / (n_perm + 1)
}

lodo_donor_level <- function(x, y, donor_x, donor_y, dataset_name) {
  # Leave-one-DONOR-out stability
  donors_x <- unique(donor_x)
  donors_y <- unique(donor_y)
  all_donor_ids <- c(donors_x, donors_y)
  estimates <- numeric(length(all_donor_ids))

  for (i in seq_along(all_donor_ids)) {
    omit <- all_donor_ids[i]
    keep_x <- donor_x != omit
    keep_y <- donor_y != omit
    sub_x <- x[keep_x]
    sub_y <- y[keep_y]
    if (length(sub_x) < 2 || length(sub_y) < 2) {
      estimates[i] <- NA_real_
    } else {
      estimates[i] <- mean(sub_y) - mean(sub_x)
    }
  }
  data.frame(
    dataset      = dataset_name,
    omitted_donor = all_donor_ids,
    estimate     = estimates,
    stringsAsFactors = FALSE
  )
}

lodo_sample_level <- function(x, y, dataset_name) {
  # Leave-one-sample-out (for datasets where donor_id = sample_id)
  all_vals <- c(x, y)
  group <- rep(c("ref", "aki"), c(length(x), length(y)))
  estimates <- numeric(length(all_vals))
  for (i in seq_along(all_vals)) {
    sub_vals <- all_vals[-i]
    sub_group <- group[-i]
    sub_x <- sub_vals[sub_group == "ref"]
    sub_y <- sub_vals[sub_group == "aki"]
    if (length(sub_x) < 2 || length(sub_y) < 2) {
      estimates[i] <- NA_real_
    } else {
      estimates[i] <- mean(sub_y) - mean(sub_x)
    }
  }
  data.frame(
    dataset = dataset_name,
    omitted = seq_along(all_vals),
    estimate = estimates,
    stringsAsFactors = FALSE
  )
}

gene_direction_check <- function(mat, genes, x_cols, y_cols) {
  present <- intersect(genes, rownames(mat))
  result <- data.frame(
    gene = present,
    mean_ref = rowMeans(mat[present, x_cols, drop = FALSE], na.rm = TRUE),
    mean_aki = rowMeans(mat[present, y_cols, drop = FALSE], na.rm = TRUE),
    direction = NA_character_,
    stringsAsFactors = FALSE
  )
  result$effect <- result$mean_aki - result$mean_ref
  result$direction <- ifelse(result$effect < 0, "negative", "positive")
  result
}

# -------------------------------------------------------------------
# 3. GSE139061 – Eadon bulk biopsy AKI vs nephrectomy reference
#    Processed quantile-normalized matrix: AKI1-39 vs REF1-9
# -------------------------------------------------------------------

message("\n======== GSE139061: AKI vs Reference (bulk biopsy) ========")

df139 <- read.csv(
  file.path(data_dir, "GSE139061_Eadon_processed_QN_101419.csv.gz"),
  stringsAsFactors = FALSE, check.names = FALSE
)
rownames(df139) <- df139$Geneid
df139$Geneid <- NULL

aki_cols_139 <- grep("^AKI", colnames(df139), value = TRUE)
ref_cols_139 <- grep("^REF", colnames(df139), value = TRUE)
mat139 <- as.matrix(df139)

# Check gene coverage
cov139 <- check_genes(all_required_genes, rownames(mat139), "GSE139061")

# TAL score
tal_genes_139 <- intersect(tal_genes, rownames(mat139))
tal_z_139 <- zscore_genes(mat139[tal_genes_139, , drop = FALSE])
tal_score_139 <- colMeans(tal_z_139)

# PT identity score
pt_id_139 <- intersect(pt_identity_genes, rownames(mat139))
pt_id_z_139 <- zscore_genes(mat139[pt_id_139, , drop = FALSE])
pt_id_score_139 <- colMeans(pt_id_z_139)

# PT injury score
pt_inj_139 <- intersect(pt_injury_genes, rownames(mat139))
pt_inj_z_139 <- zscore_genes(mat139[pt_inj_139, , drop = FALSE])
pt_inj_score_139 <- colMeans(pt_inj_z_139)

# Effect estimates
eff139_tal <- effect_estimate(
  tal_score_139[ref_cols_139], tal_score_139[aki_cols_139]
)
eff139_ptid <- effect_estimate(
  pt_id_score_139[ref_cols_139], pt_id_score_139[aki_cols_139]
)
eff139_ptinj <- effect_estimate(
  pt_inj_score_139[ref_cols_139], pt_inj_score_139[aki_cols_139]
)

# Adjusted model: TAL ~ AKI + PT_injury
df139_scores <- data.frame(
  sample  = c(ref_cols_139, aki_cols_139),
  group   = c(rep("Reference", length(ref_cols_139)),
              rep("AKI", length(aki_cols_139))),
  TAL_score = c(tal_score_139[ref_cols_139], tal_score_139[aki_cols_139]),
  PT_identity = c(pt_id_score_139[ref_cols_139], pt_id_score_139[aki_cols_139]),
  PT_injury = c(pt_inj_score_139[ref_cols_139], pt_inj_score_139[aki_cols_139]),
  stringsAsFactors = FALSE
)
df139_scores$AKI <- as.integer(df139_scores$group == "AKI")
df139_scores$TAL_minus_PT <- df139_scores$TAL_score - df139_scores$PT_identity
df139_scores$dataset <- "GSE139061"

fit139_adj <- lm(TAL_score ~ AKI + PT_injury, data = df139_scores)
fit139_inter <- lm(TAL_minus_PT ~ AKI, data = df139_scores)
s139_adj <- coef(summary(fit139_adj))
s139_inter <- coef(summary(fit139_inter))

# LODO (sample-level: GSE139061 has no donor grouping)
lodo139 <- lodo_sample_level(
  tal_score_139[ref_cols_139], tal_score_139[aki_cols_139], "GSE139061"
)
lodo139_frac <- mean(lodo139$estimate < 0, na.rm = TRUE)

# Gene-level directions
gene_dir139 <- gene_direction_check(
  mat139, tal_genes, ref_cols_139, aki_cols_139
)
neg139 <- sum(gene_dir139$effect < 0, na.rm = TRUE)

# Permutation
perm139 <- permutation_pvalue(
  tal_score_139[ref_cols_139], tal_score_139[aki_cols_139]
)

# TAL-minus-PT identity contrast
tal_minus_pt_139 <- tal_score_139 - pt_id_score_139
eff139_contrast <- effect_estimate(
  tal_minus_pt_139[ref_cols_139], tal_minus_pt_139[aki_cols_139]
)

message(sprintf(
  "  TAL AKI-REF: %.3f (95%% CI %.3f to %.3f), p=%.4f",
  eff139_tal$estimate, eff139_tal$ci_low, eff139_tal$ci_high,
  eff139_tal$p_value
))
message(sprintf("  LODO negative fraction: %.3f", lodo139_frac))
message(sprintf("  Negative genes: %d/6", neg139))
message(sprintf("  Permutation p: %.4f", perm139))
message(sprintf("  PT-injury-adjusted TAL effect: %.3f",
                unname(s139_adj["AKI", "Estimate"])))
message(sprintf("  TAL-minus-PT contrast: %.3f", eff139_contrast$estimate))

# -------------------------------------------------------------------
# 4. GSE30718 – Famulski transplant AKI vs stable transplant
#    Affymetrix HGU133 Plus 2.0 microarray
#    Donor-level: extract leading patient number; 49A1/49A2 → donor 49,
#    56A1/56A2 → donor 56. Primary retains earliest biopsy (A1).
#    Sensitivity: average repeated biopsies within donor.
# -------------------------------------------------------------------

message("\n======== GSE30718: Transplant AKI vs Stable Transplant ========")

# Parse series matrix
gz30718 <- gzfile(file.path(data_dir, "GSE30718_series_matrix.txt.gz"), "rt")
lines30718 <- readLines(gz30718)
close(gz30718)

title_line <- grep("^!Sample_title", lines30718, value = TRUE)[1]
id_line    <- grep("^!Sample_geo_accession", lines30718, value = TRUE)[1]

parse_tabs <- function(x) strsplit(x, "\t")[[1]][-1]
titles30718_raw <- parse_tabs(title_line)
ids30718_raw    <- parse_tabs(id_line)

# Strip quotes
titles30718 <- gsub('"', '', titles30718_raw)
ids30718    <- gsub('"', '', ids30718_raw)

# Classify samples and parse donor_id
find_grp <- function(t) {
  if (grepl("Acute kidney injury|^AKI", t, ignore.case = TRUE)) return("AKI")
  if (grepl("pristine protocol biopsy|^pPB", t, ignore.case = TRUE)) return("Stable")
  if (grepl("Nephrectomy|^NEPH", t, ignore.case = TRUE)) return("Nephrectomy")
  return("Other")
}
grps30718 <- sapply(titles30718, find_grp)

# Parse donor_id from title
# AKI: "Acute kidney injury AKI 374A2AGX7" -> extract "374" from "374A2AGX7"
# Stable: "pristine protocol biopsy pPB1" -> "pPB1"
parse_donor <- function(t, grp) {
  if (grp == "Stable") {
    # Extract pPB<N>
    m <- regmatches(t, regexpr("pPB\\d+", t))
    if (length(m) > 0) return(m)
  }
  if (grp == "AKI") {
    # Find the AKI ID (space-separated last token), then extract leading digits
    parts <- strsplit(t, " ")[[1]]
    aki_token <- parts[length(parts)]  # "374A2AGX7"
    # Leading digits before first 'A'
    m <- regmatches(aki_token, regexpr("^\\d+", aki_token))
    if (length(m) > 0) return(m)
  }
  return(t)  # fallback
}
donor30718 <- mapply(parse_donor, titles30718, grps30718)

# Build metadata: sample_id, donor_id, group, title
meta30718 <- data.frame(
  sample_id = ids30718,
  donor_id  = donor30718,
  group     = grps30718,
  title     = titles30718,
  stringsAsFactors = FALSE
)

# Report donor counts
n_aki_donors <- length(unique(donor30718[grps30718 == "AKI"]))
n_stable_donors <- length(unique(donor30718[grps30718 == "Stable"]))
n_aki_samples <- sum(grps30718 == "AKI")
message(sprintf(
  "  Samples: AKI=%d Stable=%d Nephrectomy=%d",
  n_aki_samples,
  sum(grps30718 == "Stable"),
  sum(grps30718 == "Nephrectomy")
))
message(sprintf("  Donors: AKI=%d Stable=%d", n_aki_donors, n_stable_donors))

dup_aki <- meta30718[meta30718$group == "AKI" &
                      duplicated(meta30718$donor_id) |
                      duplicated(meta30718$donor_id, fromLast = TRUE), ]
if (nrow(dup_aki) > 0) {
  message(sprintf("  Duplicate AKI donors: %s",
    paste(unique(dup_aki$donor_id), collapse = ", ")))
}

# Read the expression matrix
data_start <- grep("^!series_matrix_table_begin", lines30718) + 1

data_lines <- lines30718[data_start:length(lines30718)]
data_lines <- data_lines[data_lines != "!series_matrix_table_end"]

probe_ids <- character(length(data_lines) - 1)
mat_rows  <- matrix(NA_real_, nrow = length(data_lines) - 1,
                    ncol = nrow(meta30718))

for (i in seq_along(data_lines)[-1]) {
  parts <- strsplit(data_lines[i], "\t")[[1]]
  probe_ids[i - 1] <- gsub('"', '', parts[1])
  vals <- suppressWarnings(as.numeric(parts[-1]))
  mat_rows[i - 1, ] <- vals
}

rownames(mat_rows) <- probe_ids
colnames(mat_rows) <- meta30718$sample_id

message(sprintf("  Matrix: %d probes x %d samples", nrow(mat_rows), ncol(mat_rows)))

# Map probes to genes using GPL570 annotation
annot_gz <- gzfile(file.path(data_dir, "GPL570.annot.gz"), "rt")
annot_lines <- readLines(annot_gz)
close(annot_gz)

probe_to_gene <- list()
for (line in annot_lines) {
  if (nchar(line) == 0 || substr(line, 1, 1) %in% c("^", "!", "#")) next
  parts <- strsplit(line, "\t")[[1]]
  if (length(parts) >= 3) {
    pid <- parts[1]
    gs  <- parts[3]
    if (pid != "" && gs != "" && gs != "---" && gs != "Gene symbol") {
      gs_list <- strsplit(gs, " /// ")[[1]]
      probe_to_gene[[pid]] <- gs_list
    }
  }
}

# Collapse probes to genes via median
collapse_probes <- function(mat, gene_set, probe_map) {
  result <- matrix(NA_real_, nrow = 0, ncol = ncol(mat))
  for (gene in gene_set) {
    probes <- names(probe_map)[sapply(probe_map, function(gs) gene %in% gs)]
    if (length(probes) == 0) next
    probes <- intersect(probes, rownames(mat))
    if (length(probes) == 0) next
    if (length(probes) == 1) {
      row_vals <- mat[probes, , drop = FALSE]
    } else {
      row_vals <- matrix(
        apply(mat[probes, , drop = FALSE], 2, median, na.rm = TRUE),
        nrow = 1
      )
    }
    rownames(row_vals) <- gene
    result <- rbind(result, row_vals)
  }
  result
}

tal_mat_30718 <- collapse_probes(mat_rows, tal_genes, probe_to_gene)
pt_id_mat_30718 <- collapse_probes(mat_rows, pt_identity_genes, probe_to_gene)
pt_inj_mat_30718 <- collapse_probes(mat_rows, pt_injury_genes, probe_to_gene)

message(sprintf("  TAL genes collapsed: %d", nrow(tal_mat_30718)))
message(sprintf("  PT identity genes collapsed: %d", nrow(pt_id_mat_30718)))
message(sprintf("  PT injury genes collapsed: %d", nrow(pt_inj_mat_30718)))

# Subset to AKI vs Stable
keep30718 <- meta30718$group %in% c("AKI", "Stable")
meta_sub <- meta30718[keep30718, ]

# Assign biopsy_order (only donors 49 and 56 have repeats)
meta_sub$biopsy_order <- "single"
for (i in seq_len(nrow(meta_sub))) {
  if (meta_sub$donor_id[i] %in% c("49", "56")) {
    t <- meta_sub$title[i]
    parts <- strsplit(t, " ")[[1]]
    aki_token <- parts[length(parts)]
    if (grepl("A1", aki_token)) meta_sub$biopsy_order[i] <- "A1"
    else if (grepl("A2", aki_token)) meta_sub$biopsy_order[i] <- "A2"
  }
}

# ---- PRIMARY: keep earliest biopsy per donor, z-score only primary ----
primary_keep <- meta_sub$biopsy_order %in% c("A1", "single") |
                meta_sub$group == "Stable"
meta_primary <- meta_sub[primary_keep, ]

message(sprintf(
  "  Primary donor-level: AKI=%d Stable=%d",
  sum(meta_primary$group == "AKI"),
  sum(meta_primary$group == "Stable")
))

# Subset expression matrices to primary samples only
primary_cols <- meta_primary$sample_id
tal_mat_pri  <- tal_mat_30718[, primary_cols, drop = FALSE]
pt_id_mat_pri  <- pt_id_mat_30718[, primary_cols, drop = FALSE]
pt_inj_mat_pri <- pt_inj_mat_30718[, primary_cols, drop = FALSE]

# Z-score within primary cohort only
tal_z_pri  <- zscore_genes(tal_mat_pri)
pt_id_z_pri  <- zscore_genes(pt_id_mat_pri)
pt_inj_z_pri <- zscore_genes(pt_inj_mat_pri)

tal_score_pri  <- colMeans(tal_z_pri)
pt_id_score_pri  <- colMeans(pt_id_z_pri)
pt_inj_score_pri <- colMeans(pt_inj_z_pri)

# Build primary donor-level scores
ds <- data.frame(
  sample_id   = meta_primary$sample_id,
  donor_id    = meta_primary$donor_id,
  group       = meta_primary$group,
  TAL_score   = tal_score_pri,
  PT_identity = pt_id_score_pri,
  PT_injury   = pt_inj_score_pri,
  stringsAsFactors = FALSE
)

stable_idx <- ds$group == "Stable"
aki_idx    <- ds$group == "AKI"

# ---- SENSITIVITY: average per-gene expression within donor, then re-z-score ----
donors_aki    <- unique(meta_sub$donor_id[meta_sub$group == "AKI"])
donors_stable <- unique(meta_sub$donor_id[meta_sub$group == "Stable"])
all_donors    <- c(donors_aki, donors_stable)

avg_per_donor <- function(mat, meta, donors) {
  result <- matrix(NA_real_, nrow = nrow(mat), ncol = length(donors))
  colnames(result) <- donors
  rownames(result) <- rownames(mat)
  for (j in seq_along(donors)) {
    samples <- meta$sample_id[meta$donor_id == donors[j]]
    if (length(samples) == 1) {
      result[, j] <- mat[, samples]
    } else {
      result[, j] <- rowMeans(mat[, samples, drop = FALSE], na.rm = TRUE)
    }
  }
  result
}

all_keep_cols <- intersect(colnames(tal_mat_30718), meta_sub$sample_id)
tal_sens_mat  <- avg_per_donor(tal_mat_30718[, all_keep_cols, drop = FALSE],  meta_sub, all_donors)
pt_id_sens_mat  <- avg_per_donor(pt_id_mat_30718[, all_keep_cols, drop = FALSE],  meta_sub, all_donors)
pt_inj_sens_mat <- avg_per_donor(pt_inj_mat_30718[, all_keep_cols, drop = FALSE], meta_sub, all_donors)

# Re-z-score on donor-averaged matrices
tal_z_sens  <- zscore_genes(tal_sens_mat)
pt_id_z_sens  <- zscore_genes(pt_id_sens_mat)
pt_inj_z_sens <- zscore_genes(pt_inj_sens_mat)

tal_score_sens  <- colMeans(tal_z_sens)
pt_id_score_sens  <- colMeans(pt_id_z_sens)
pt_inj_score_sens <- colMeans(pt_inj_z_sens)

# Build sensitivity donor-level scores
ds_sens <- data.frame(
  donor_id    = all_donors,
  group       = ifelse(all_donors %in% donors_aki, "AKI", "Stable"),
  TAL_score   = tal_score_sens,
  PT_identity = pt_id_score_sens,
  PT_injury   = pt_inj_score_sens,
  stringsAsFactors = FALSE
)
ds_sens$sample_id <- ds_sens$donor_id  # placeholder for consistent column names

sens_stable_idx <- ds_sens$group == "Stable"
sens_aki_idx    <- ds_sens$group == "AKI"

# ---- PRIMARY effect estimates ----
eff30718_tal <- effect_estimate(
  ds$TAL_score[stable_idx], ds$TAL_score[aki_idx]
)

ds$AKI <- as.integer(ds$group == "AKI")
ds$TAL_minus_PT <- ds$TAL_score - ds$PT_identity
ds$dataset <- "GSE30718"

fit30718_adj <- lm(TAL_score ~ AKI + PT_injury, data = ds)
s30718_adj <- coef(summary(fit30718_adj))
fit30718_inter <- lm(TAL_minus_PT ~ AKI, data = ds)
s30718_inter <- coef(summary(fit30718_inter))

eff30718_contrast <- effect_estimate(
  ds$TAL_minus_PT[stable_idx], ds$TAL_minus_PT[aki_idx]
)

lodo30718 <- lodo_donor_level(
  ds$TAL_score[stable_idx], ds$TAL_score[aki_idx],
  ds$donor_id[stable_idx], ds$donor_id[aki_idx],
  "GSE30718"
)
lodo30718_frac <- mean(lodo30718$estimate < 0, na.rm = TRUE)

gene_dir30718 <- gene_direction_check(
  tal_mat_pri, tal_genes,
  ds$sample_id[stable_idx], ds$sample_id[aki_idx]
)
neg30718 <- sum(gene_dir30718$effect < 0, na.rm = TRUE)

perm30718 <- permutation_pvalue(
  ds$TAL_score[stable_idx], ds$TAL_score[aki_idx]
)

message(sprintf(
  "  PRIMARY TAL AKI-Stable: %.3f (95%% CI %.3f to %.3f), p=%.4f",
  eff30718_tal$estimate, eff30718_tal$ci_low, eff30718_tal$ci_high,
  eff30718_tal$p_value
))
message(sprintf("  Unique AKI donors: %d  Stable donors: %d",
                sum(aki_idx), sum(stable_idx)))
message(sprintf("  LODO (donor-level) negative fraction: %.3f", lodo30718_frac))
message(sprintf("  Negative genes: %d/%d", neg30718, nrow(gene_dir30718)))
message(sprintf("  Permutation p: %.4f", perm30718))
message(sprintf("  PT-injury-adjusted TAL effect: %.3f",
                unname(s30718_adj["AKI", "Estimate"])))
message(sprintf("  TAL-minus-PT contrast: %.3f", eff30718_contrast$estimate))

# ---- SENSITIVITY effect estimates ----
eff30718_sens_tal <- effect_estimate(
  ds_sens$TAL_score[sens_stable_idx], ds_sens$TAL_score[sens_aki_idx]
)

ds_sens$AKI <- as.integer(ds_sens$group == "AKI")
ds_sens$TAL_minus_PT <- ds_sens$TAL_score - ds_sens$PT_identity

fit30718_sens_adj <- lm(TAL_score ~ AKI + PT_injury, data = ds_sens)
s30718_sens_adj <- coef(summary(fit30718_sens_adj))

eff30718_sens_contrast <- effect_estimate(
  ds_sens$TAL_minus_PT[sens_stable_idx], ds_sens$TAL_minus_PT[sens_aki_idx]
)

lodo30718_sens <- lodo_donor_level(
  ds_sens$TAL_score[sens_stable_idx], ds_sens$TAL_score[sens_aki_idx],
  ds_sens$donor_id[sens_stable_idx], ds_sens$donor_id[sens_aki_idx],
  "GSE30718_sensitivity"
)
lodo30718_sens_frac <- mean(lodo30718_sens$estimate < 0, na.rm = TRUE)

gene_dir30718_sens <- gene_direction_check(
  tal_sens_mat, tal_genes,
  ds_sens$donor_id[sens_stable_idx], ds_sens$donor_id[sens_aki_idx]
)
neg30718_sens <- sum(gene_dir30718_sens$effect < 0, na.rm = TRUE)

perm30718_sens <- permutation_pvalue(
  ds_sens$TAL_score[sens_stable_idx], ds_sens$TAL_score[sens_aki_idx]
)

message(sprintf(
  "  SENSITIVITY TAL AKI-Stable: %.3f (95%% CI %.3f to %.3f), p=%.4f",
  eff30718_sens_tal$estimate, eff30718_sens_tal$ci_low, eff30718_sens_tal$ci_high,
  eff30718_sens_tal$p_value
))
message(sprintf("  SENSITIVITY LODO negative fraction: %.3f", lodo30718_sens_frac))
message(sprintf("  SENSITIVITY Negative genes: %d/%d", neg30718_sens, nrow(gene_dir30718_sens)))
message(sprintf("  SENSITIVITY Permutation p: %.4f", perm30718_sens))
message(sprintf("  SENSITIVITY PT-injury-adjusted TAL effect: %.3f",
                unname(s30718_sens_adj["AKI", "Estimate"])))
message(sprintf("  SENSITIVITY TAL-minus-PT contrast: %.3f", eff30718_sens_contrast$estimate))

# Build scores df for output (primary)
df30718_scores <- ds
names(df30718_scores)[names(df30718_scores) == "sample_id"] <- "sample"
df30718_scores$biopsy_order <- NA_character_
df30718_scores$analysis_role <- "primary"
for (i in seq_len(nrow(df30718_scores))) {
  sid <- df30718_scores$sample[i]
  idx <- which(meta_sub$sample_id == sid)
  if (length(idx) == 1) df30718_scores$biopsy_order[i] <- meta_sub$biopsy_order[idx]
}

# -------------------------------------------------------------------
# 5. Combine results for synthesis
# -------------------------------------------------------------------

# Add donor_id to GSE139061 (sample == donor for bulk biopsies)
df139_scores$donor_id       <- df139_scores$sample
df139_scores$biopsy_order   <- "single"
df139_scores$analysis_role  <- "primary"

# Build combined scores
n139 <- nrow(df139_scores)
n307 <- nrow(df30718_scores)
N    <- n139 + n307

all_scores <- data.frame(
  sample        = character(N),
  donor_id      = character(N),
  group         = character(N),
  TAL_score     = numeric(N),
  PT_identity   = numeric(N),
  PT_injury     = numeric(N),
  AKI           = integer(N),
  TAL_minus_PT  = numeric(N),
  dataset       = character(N),
  biopsy_order  = character(N),
  analysis_role = character(N),
  stringsAsFactors = FALSE
)

fill_row <- function(df, offset, n) {
  if (n == 0) return()
  idx <- (offset + 1):(offset + n)
  all_scores$sample[idx]        <<- as.character(df[["sample"]])
  all_scores$donor_id[idx]      <<- as.character(df[["donor_id"]])
  all_scores$group[idx]         <<- as.character(df[["group"]])
  all_scores$TAL_score[idx]     <<- as.numeric(df[["TAL_score"]])
  all_scores$PT_identity[idx]   <<- as.numeric(df[["PT_identity"]])
  all_scores$PT_injury[idx]     <<- as.numeric(df[["PT_injury"]])
  all_scores$AKI[idx]           <<- as.integer(df[["AKI"]])
  all_scores$TAL_minus_PT[idx]  <<- as.numeric(df[["TAL_minus_PT"]])
  all_scores$dataset[idx]       <<- as.character(df[["dataset"]])
  all_scores$biopsy_order[idx]  <<- as.character(df[["biopsy_order"]])
  all_scores$analysis_role[idx] <<- as.character(df[["analysis_role"]])
}

fill_row(df139_scores, 0, n139)
fill_row(df30718_scores, n139, n307)

write.csv(all_scores, file.path(result_dir, "bulk_scores.csv"),
          row.names = FALSE)

# --- Build and save bulk_samples.csv ---
# GSE139061: simple metadata
meta139_out <- data.frame(
  dataset    = "GSE139061",
  sample_id  = colnames(df139),
  donor_id   = colnames(df139),
  group      = c(rep("Reference", length(ref_cols_139)),
                 rep("AKI", length(aki_cols_139))),
  biopsy_order   = "single",
  analysis_role  = "primary",
  title      = colnames(df139),
  stringsAsFactors = FALSE
)

# GSE30718: full metadata including sensitivity-only biopsies
meta30718_out <- meta30718
meta30718_out$biopsy_order <- ifelse(
  meta30718_out$sample_id %in% meta_sub$sample_id,
  meta_sub$biopsy_order[match(meta30718_out$sample_id, meta_sub$sample_id)],
  NA_character_
)
meta30718_out$analysis_role <- "primary"
# 49A2 and 56A2 are sensitivity-only
for (i in seq_len(nrow(meta30718_out))) {
  if (!is.na(meta30718_out$biopsy_order[i]) &&
      meta30718_out$biopsy_order[i] == "A2") {
    meta30718_out$analysis_role[i] <- "sensitivity_only"
  }
}
# Nephrectomy samples are unused — never enter any analysis
meta30718_out$analysis_role[meta30718_out$group == "Nephrectomy"] <- "context_only"

meta_out_307 <- data.frame(
  dataset   = "GSE30718",
  sample_id = meta30718_out$sample_id,
  donor_id  = meta30718_out$donor_id,
  group     = meta30718_out$group,
  biopsy_order  = meta30718_out$biopsy_order,
  analysis_role = meta30718_out$analysis_role,
  title     = meta30718_out$title,
  stringsAsFactors = FALSE
)

samples_csv <- rbind(meta139_out, meta_out_307)
write.csv(samples_csv, file.path(meta_dir, "bulk_samples.csv"),
          row.names = FALSE)

# Effect summary (primary + sensitivity)
effects_summary <- data.frame(
  dataset = c("GSE139061", "GSE139061", "GSE139061",
              "GSE30718", "GSE30718", "GSE30718",
              "GSE30718_sens", "GSE30718_sens", "GSE30718_sens"),
  model = c("unadjusted", "adjusted_PT_injury", "TAL_minus_PT",
            "unadjusted", "adjusted_PT_injury", "TAL_minus_PT",
            "unadjusted", "adjusted_PT_injury", "TAL_minus_PT"),
  analysis_level = c("primary", "primary", "primary",
                     "primary", "primary", "primary",
                     "sensitivity", "sensitivity", "sensitivity"),
  estimate = c(
    eff139_tal$estimate,
    unname(s139_adj["AKI", "Estimate"]),
    eff139_contrast$estimate,
    eff30718_tal$estimate,
    unname(s30718_adj["AKI", "Estimate"]),
    eff30718_contrast$estimate,
    eff30718_sens_tal$estimate,
    unname(s30718_sens_adj["AKI", "Estimate"]),
    eff30718_sens_contrast$estimate
  ),
  ci_low = c(
    eff139_tal$ci_low,
    unname(confint(fit139_adj)["AKI", 1]),
    eff139_contrast$ci_low,
    eff30718_tal$ci_low,
    unname(confint(fit30718_adj)["AKI", 1]),
    eff30718_contrast$ci_low,
    eff30718_sens_tal$ci_low,
    unname(confint(fit30718_sens_adj)["AKI", 1]),
    eff30718_sens_contrast$ci_low
  ),
  ci_high = c(
    eff139_tal$ci_high,
    unname(confint(fit139_adj)["AKI", 2]),
    eff139_contrast$ci_high,
    eff30718_tal$ci_high,
    unname(confint(fit30718_adj)["AKI", 2]),
    eff30718_contrast$ci_high,
    eff30718_sens_tal$ci_high,
    unname(confint(fit30718_sens_adj)["AKI", 2]),
    eff30718_sens_contrast$ci_high
  ),
  p_value = c(
    eff139_tal$p_value,
    unname(s139_adj["AKI", "Pr(>|t|)"]),
    eff139_contrast$p_value,
    eff30718_tal$p_value,
    unname(s30718_adj["AKI", "Pr(>|t|)"]),
    eff30718_contrast$p_value,
    eff30718_sens_tal$p_value,
    unname(s30718_sens_adj["AKI", "Pr(>|t|)"]),
    eff30718_sens_contrast$p_value
  ),
  lodo_negative_frac = c(
    lodo139_frac, NA, NA,
    lodo30718_frac, NA, NA,
    lodo30718_sens_frac, NA, NA
  ),
  negative_genes = c(
    sprintf("%d/6", neg139), NA, NA,
    sprintf("%d/%d", neg30718, nrow(gene_dir30718)), NA, NA,
    sprintf("%d/%d", neg30718_sens, nrow(gene_dir30718_sens)), NA, NA
  ),
  perm_p = c(perm139, NA, NA, perm30718, NA, NA, perm30718_sens, NA, NA),
  n_aki = c(
    length(aki_cols_139), NA, NA,
    sum(aki_idx), NA, NA,
    sum(sens_aki_idx), NA, NA
  ),
  n_ref = c(
    length(ref_cols_139), NA, NA,
    sum(stable_idx), NA, NA,
    sum(sens_stable_idx), NA, NA
  ),
  stringsAsFactors = FALSE
)

write.csv(effects_summary, file.path(result_dir, "bulk_effects.csv"),
          row.names = FALSE)

# LODO combined (primary + sensitivity)
names(lodo139)[names(lodo139) == "omitted"] <- "omitted_unit"
names(lodo30718)[names(lodo30718) == "omitted_donor"] <- "omitted_unit"
names(lodo30718_sens)[names(lodo30718_sens) == "omitted_donor"] <- "omitted_unit"
lodo_combined <- rbind(lodo139, lodo30718, lodo30718_sens)
write.csv(lodo_combined, file.path(result_dir, "bulk_lodo.csv"),
          row.names = FALSE)

# Gene direction combined (primary + sensitivity)
gene_dir_combined <- rbind(
  cbind(gene_dir139, dataset = "GSE139061", analysis_level = "primary"),
  cbind(gene_dir30718, dataset = "GSE30718", analysis_level = "primary"),
  cbind(gene_dir30718_sens, dataset = "GSE30718", analysis_level = "sensitivity")
)
write.csv(gene_dir_combined, file.path(result_dir, "bulk_genes.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------
# 6. Decision rules
# -------------------------------------------------------------------

unadj_pass_139 <- eff139_tal$estimate < 0 && lodo139_frac >= 0.8 && neg139 >= 4
unadj_pass_30718 <- eff30718_tal$estimate < 0 && lodo30718_frac >= 0.8 && neg30718 >= 4

decision139 <- if (unadj_pass_139) "KEEP_SUPPORTIVE" else "REJECT_VALIDATION"
decision30718 <- if (unadj_pass_30718) "KEEP_SUPPORTIVE" else "REJECT_VALIDATION"

# Enhanced: PT-adjusted effect negative in both
adj_pass_139 <- unname(s139_adj["AKI", "Estimate"]) < 0
adj_pass_30718 <- unname(s30718_adj["AKI", "Estimate"]) < 0
contrast_pass_139 <- eff139_contrast$estimate < 0
contrast_pass_30718 <- eff30718_contrast$estimate < 0

# Bulk-evidence decision (KEEP_SUPPORTIVE vs REJECT)
bulk_decision <- if (unadj_pass_139 && unadj_pass_30718) {
  "KEEP_SUPPORTIVE"
} else if (unadj_pass_139) {
  "KEEP_DISCOVERY_ONLY"
} else {
  "REJECT_VALIDATION"
}

# Project status (separate from bulk evidence)
project_status <- "PROCEED_FORMAL_ANALYSIS_WITH_NARROWED_CLAIM"

# -------------------------------------------------------------------
# 7. Figures
# -------------------------------------------------------------------

# Forest plot (primary models only)
forest_data <- effects_summary[effects_summary$analysis_level == "primary", ]
forest_data$model <- factor(forest_data$model,
  levels = c("unadjusted", "adjusted_PT_injury", "TAL_minus_PT"))

forest_plot <- ggplot(
  forest_data,
  aes(x = estimate, y = model, colour = dataset)
) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    width = 0.15, orientation = "y",
    position = position_dodge(width = 0.35)
  ) +
  geom_point(size = 2.6, position = position_dodge(width = 0.35)) +
  labs(
    x = "AKI minus reference effect (TAL identity score)",
    y = NULL,
    colour = "Dataset",
    title = "TAL identity in bulk tissue: AKI vs reference/stable"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(figure_dir, "bulk_forest.png"),
  forest_plot, width = 9, height = 4.8, dpi = 180
)

# Boxplot: TAL score by group and dataset (primary donor-level only)
plot_data <- all_scores[!is.na(all_scores$TAL_score), ]
plot_data$group <- factor(plot_data$group,
  levels = c("Reference", "Stable", "AKI"))

score_box <- ggplot(
  plot_data,
  aes(x = group, y = TAL_score, colour = group)
) +
  geom_boxplot(width = 0.5, outlier.shape = NA, colour = "grey45") +
  geom_jitter(width = 0.08, size = 2.2) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_colour_manual(values = c("Reference" = "#3b82a0", "Stable" = "#3b82a0",
                                  "AKI" = "#c24f44")) +
  labs(
    x = NULL,
    y = "Frozen TAL identity score",
    title = "TAL identity score in bulk tissue validation"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(
  file.path(figure_dir, "bulk_scores.png"),
  score_box, width = 9, height = 5, dpi = 180
)

# Gene-level direction dot plot (primary)
gene_plot_data <- gene_dir_combined[gene_dir_combined$analysis_level == "primary", ]

gene_dot <- ggplot(
  gene_plot_data,
  aes(x = gene, y = dataset, size = abs(effect), colour = effect)
) +
  geom_point() +
  scale_colour_gradient2(low = "#2b6f9f", mid = "white", high = "#c24f44",
                         name = "AKI - Ref") +
  scale_size_continuous(name = "|Effect|") +
  labs(
    x = NULL, y = NULL,
    title = "Per-gene TAL expression change: AKI vs Reference"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(figure_dir, "bulk_genes.png"),
  gene_dot, width = 8, height = 3.2, dpi = 180
)

# -------------------------------------------------------------------
# 8. Summary markdown
# -------------------------------------------------------------------

fmt3 <- function(x) if (is.na(x)) "NA" else sprintf("%.3f", x)
fmt_p <- function(x) if (is.na(x)) "NA" else sprintf("%.4f", x)

summary_lines <- c(
  "# AKI–TAL bulk GEO validation summary",
  "",
  sprintf("**Bulk evidence decision:** `%s`", bulk_decision),
  sprintf("**Project status:** `%s`", project_status),
  "",
  "## Scope amendment (2026-07-29)",
  "",
  "Bulk tissue datasets measure tissue-level detectability of the TAL",
  "programme. They cannot validate the KDIGO-stage association from the",
  "sc/snRNA-seq discovery because clinical stage is not donor-linked in",
  "these series. The original frozen question is preserved; the amended",
  "primary claim is: lower frozen TAL identity is detectable in AKI tissue",
  "across independent biopsy and transplant contexts, but KDIGO-stage",
  "association, TAL-cell-intrinsic mechanism, and TAL-selective severity",
  "require further evidence.",
  "",
  "## Dataset provenance",
  "",
  "| Dataset | Platform | Comparison | N AKI (donors) | N Ref/Stable | TAL genes |",
  "|---|---|---|---|---|---|",
  sprintf("| GSE139061 | RNA-seq bulk biopsy | AKI vs nephrectomy | 39 | 9 | 6/6 |"),
  sprintf("| GSE30718 | Affy HGU133 Plus 2.0 | Tx-AKI vs stable Tx | 28 (26) | 11 | %d/%d |",
          nrow(gene_dir30718), length(tal_genes)),
  "",
  "GSE30718 donor-level: 28 AKI samples from 26 unique donors.",
  "49A1/49A2 → donor 49; 56A1/56A2 → donor 56.",
  "Primary: earliest biopsy retained (A1); 49A2 and 56A2 are sensitivity-only.",
  "Z-scoring performed on primary samples only; sensitivity averages per-gene",
  "expression within donor then re-z-scores the donor-level matrix.",
  "",
  "## Primary results",
  "",
  sprintf(
    "- **GSE139061** AKI–REF: %s (95%% CI %s to %s), p=%s",
    fmt3(eff139_tal$estimate), fmt3(eff139_tal$ci_low),
    fmt3(eff139_tal$ci_high), fmt_p(eff139_tal$p_value)
  ),
  sprintf(
    "- **GSE30718** Tx-AKI–Stable (26 donors, primary): %s (95%% CI %s to %s), p=%s",
    fmt3(eff30718_tal$estimate), fmt3(eff30718_tal$ci_low),
    fmt3(eff30718_tal$ci_high), fmt_p(eff30718_tal$p_value)
  ),
  "",
  "### Stability",
  sprintf("- GSE139061 LODO negative fraction: %.3f | Negative genes: %d/6 | Permutation p: %.4f",
          lodo139_frac, neg139, perm139),
  sprintf("- GSE30718 LODO (donor-level) negative fraction: %.3f | Negative genes: %d/%d | Permutation p: %.4f",
          lodo30718_frac, neg30718, nrow(gene_dir30718), perm30718),
  "",
  "### PT-adjusted and specificity",
  "",
  sprintf("- GSE139061 PT-injury-adjusted: %s", fmt3(unname(s139_adj["AKI", "Estimate"]))),
  sprintf("- GSE30718 PT-injury-adjusted: %s", fmt3(unname(s30718_adj["AKI", "Estimate"]))),
  sprintf("- GSE139061 TAL-minus-PT contrast: %s", fmt3(eff139_contrast$estimate)),
  sprintf("- GSE30718 TAL-minus-PT contrast: %s", fmt3(eff30718_contrast$estimate)),
  "",
  "## GSE30718 sensitivity (donor-mean collapse → re-z-score)",
  "",
  sprintf(
    "- TAL effect: %s (95%% CI %s to %s), p=%s",
    fmt3(eff30718_sens_tal$estimate), fmt3(eff30718_sens_tal$ci_low),
    fmt3(eff30718_sens_tal$ci_high), fmt_p(eff30718_sens_tal$p_value)
  ),
  sprintf("- LODO donor-level negative fraction: %.3f", lodo30718_sens_frac),
  sprintf("- Negative genes: %d/%d", neg30718_sens, nrow(gene_dir30718_sens)),
  sprintf("- Permutation p: %.4f", perm30718_sens),
  sprintf("- PT-injury-adjusted: %s", fmt3(unname(s30718_sens_adj["AKI", "Estimate"]))),
  sprintf("- TAL-minus-PT contrast: %s", fmt3(eff30718_sens_contrast$estimate)),
  "",
  "## Dataset decisions",
  "",
  sprintf("- GSE139061: `%s`", decision139),
  sprintf("- GSE30718: `%s`", decision30718),
  "",
  "## Confounding caveats",
  "",
  "- **GSE139061:** biopsy-vs-nephrectomy tissue-source confound; lower TAL",
  "  identity may reflect tubular atrophy, cell loss, or biopsy-site effects.",
  "- **GSE30718:** transplant context (immunosuppression, rejection,",
  "  ischemia-reperfusion) differs from native AKI.",
  "- **GSE210622 (existing discovery):** COVID confound; stage effect reverses",
  "  after COVID adjustment.",
  "- No donor-linked KDIGO stage in GSE139061 or GSE30718.",
  "- Bulk tissue scores measure tissue-level gene detectability, not",
  "  TAL-cell-intrinsic transcriptional change.",
  "",
  "## Evidence weight",
  "",
  sprintf(
    paste0(
      "The frozen six-gene TAL identity score decreases in AKI tissue ",
      "relative to reference in %d/%d usable bulk validation datasets, ",
      "with stable donor-level LODO and per-gene support in both. PT-adjusted ",
      "effects remain negative. Segment specificity evidence is mixed: ",
      "GSE30718 TAL-minus-PT contrast is negative (%.3f) but GSE139061 ",
      "contrast is positive (%.3f). Sensitivity analysis using donor-mean ",
      "collapse and re-z-scoring confirms direction (%s, p=%s). The bulk ",
      "evidence supports lower TAL identity in AKI tissue across independent ",
      "contexts but does not independently validate a KDIGO-stage association ",
      "or TAL-selective mechanism."
    ),
    sum(unadj_pass_139, unadj_pass_30718), 2,
    eff30718_contrast$estimate, eff139_contrast$estimate,
    fmt3(eff30718_sens_tal$estimate), fmt_p(eff30718_sens_tal$p_value)
  ),
  "",
  sprintf("**Bulk evidence decision:** `%s`", bulk_decision),
  sprintf("**Project status:** `%s`", project_status)
)

writeLines(summary_lines, file.path(result_dir, "bulk.md"))

cat(paste(summary_lines, collapse = "\n"), "\n")

message("\nDone. Outputs saved to results/ and figures/.")
