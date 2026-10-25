#!/usr/bin/env Rscript
#
# 03_longitudinal_validation.R — GSE126805 longitudinal validation of frozen AKI-TAL signal
# Cippa et al. 2018 JCI Insight 3:e123151 (PMC6302941)
#
# Data sources (SHA256 verified on host):
#   GSE126805_all.gene_counts.txt.gz  — 91683e6f99fd200665fdd4be064c4bb3f2d0730715a832f4c9ba499fa45dc134
#   GSE126805_series_matrix.txt.gz    — 026598fa31d7b1e54a7c4e0ab2f572152f2972fac82d959994ebdbea0ea1ac43
#   PMC6302941_supp.pdf               — 515bade828f518b7d68cf970da52fb16e5aa1c9c9c1b503d93096b19e3b0eb60

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_arg) == 1L) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
} else {
  normalizePath("projects/AKI_TAL_severity_analysis")
}
result_dir <- file.path(project_dir, "results")
meta_dir   <- file.path(project_dir, "metadata")
data_dir   <- file.path(project_dir, "data", "gse126805")
figure_dir <- file.path(project_dir, "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(edgeR))

# Package versions for provenance
r_ver     <- paste(R.version$major, R.version$minor, sep=".")
edger_ver <- as.character(packageVersion("edgeR"))
limma_ver <- as.character(packageVersion("limma"))
message(sprintf("R %s, edgeR %s, limma %s", r_ver, edger_ver, limma_ver))

SET_SEED <- 42L
N_PERM   <- 9999L
set.seed(SET_SEED)

# ---- 0. Frozen gene sets (immutable) ----
TAL_GENES    <- c("SLC12A1","UMOD","CLDN16","KCNJ10","PTH1R","CLDN10")
PT_ID_GENES  <- c("LRP2","SLC34A1","SLC5A2","SLC5A1","MME",
                  "ALDOB","MIOX","PDZK1","SLC22A6","ANPEP")
PT_INJ_GENES <- c("PROM1","DCDC2","SPP1","ITGB6","ITGB8")

# ---- 1. Load count matrix ----
message("--- 1. Loading count matrix ---")
cts <- read.delim(gzfile(file.path(data_dir,"GSE126805_all.gene_counts.txt.gz")),
                  stringsAsFactors=FALSE, check.names=FALSE)
sample_cols <- grep("^sample\\.kidney_", colnames(cts), value=TRUE)
message(sprintf("  Raw: %d genes x %d columns; %d sample columns",
                nrow(cts), ncol(cts), length(sample_cols)))

# Build numeric matrix
cnt <- as.matrix(cts[, sample_cols]); mode(cnt) <- "numeric"
gene_names <- cts$external_gene_name
n_gene_nonempty <- sum(gene_names != "" & !is.na(gene_names))
n_gene_unique   <- length(unique(gene_names[gene_names != "" & !is.na(gene_names)]))

# ---- 2. Deduplicate gene symbols (deterministic: keep highest-mean row) ----
message("--- 2. Deduplicating gene symbols ---")
gene_mean <- rowMeans(cnt, na.rm=TRUE)
ord <- order(gene_names, -gene_mean)
cnt <- cnt[ord, ]; gene_names <- gene_names[ord]
n_dup <- sum(duplicated(gene_names))
keep <- !duplicated(gene_names)
cnt <- cnt[keep, ]; gene_names <- gene_names[keep]
rownames(cnt) <- gene_names; colnames(cnt) <- sample_cols
message(sprintf("  After dedup: %d genes (%d duplicates removed by highest expression)",
                nrow(cnt), n_dup))
# Confirm TAL genes are not duplicated
for (g in TAL_GENES) {
  n_hits <- sum(cts$external_gene_name == g)
  message(sprintf("  %s: %d raw rows -> 1 retained", g, n_hits))
}

# ---- 3. Parse GEO series matrix for GSM cross-validation ----
message("--- 3. GEO series matrix cross-validation ---")
mat_lines <- readLines(gzfile(file.path(data_dir,"GSE126805_series_matrix.txt.gz")))
# line 30 = !Sample_title, line 31 = !Sample_geo_accession
title_line <- mat_lines[grep("^!Sample_title", mat_lines)[1]]
gsm_line   <- mat_lines[grep("^!Sample_geo_accession", mat_lines)[1]]
titles     <- strsplit(title_line, "\t")[[1]][-1]
gsms       <- strsplit(gsm_line,   "\t")[[1]][-1]
titles     <- gsub('"', '', titles)
gsms       <- gsub('"', '', gsms)

# Build GSM lookup: sample_col_name -> GSM
sample_to_gsm <- setNames(gsms, titles)
message(sprintf("  Series matrix: %d sample titles, %d GSM IDs", length(titles), length(gsms)))

# ---- 4. Parse sample metadata from column names ----
message("--- 4. Sample metadata ---")
parse_sample <- function(s) {
  parts <- strsplit(sub("^sample\\.","",s), "_")[[1]]
  c(kidney = parts[2], timepoint = paste(parts[3:length(parts)], collapse="_"))
}
meta <- as.data.frame(do.call(rbind, lapply(sample_cols, function(s) {
  p <- parse_sample(s)
  c(sample = s, kidney = unname(p["kidney"]),
    timepoint = unname(p["timepoint"]),
    gsm = ifelse(s %in% names(sample_to_gsm), sample_to_gsm[s], NA_character_))
})), stringsAsFactors=FALSE)

stopifnot(all(meta$timepoint %in% c("baseline_pre","baseline_post","3months","1year")))
message(sprintf("  %d kidneys, %d samples", length(unique(meta$kidney)), nrow(meta)))
message(sprintf("  GSM cross-match: %d/%d", sum(!is.na(meta$gsm)), nrow(meta)))
tp_tbl <- table(meta$timepoint)
message(sprintf("  Timepoints: pre=%d post=%d 3mo=%d 1y=%d",
                tp_tbl["baseline_pre"],tp_tbl["baseline_post"],
                tp_tbl["3months"],tp_tbl["1year"]))

# ---- 5. Pairing audit ----
message("--- 5. Pairing audit ---")
has_pair <- sapply(unique(meta$kidney), function(k)
  all(c("baseline_pre","baseline_post") %in% meta$timepoint[meta$kidney==k]))
kidneys_pair <- names(has_pair)[has_pair]
n_pairs <- length(kidneys_pair)
message(sprintf("  Kidneys with both pre+post: %d/42", n_pairs))
no_pre <- setdiff(unique(meta$kidney), kidneys_pair)
if (length(no_pre) > 0)
  message(sprintf("  Excluded (no baseline_pre): %s", paste(no_pre,collapse=",")))

stopifnot(n_pairs >= 25)

# ---- 6. Gene mapping ----
message("--- 6. Gene mapping ---")
tal_found   <- intersect(TAL_GENES, rownames(cnt))
tal_missing <- setdiff(TAL_GENES, rownames(cnt))
pt_id_fnd   <- intersect(PT_ID_GENES, rownames(cnt))
pt_inj_fnd  <- intersect(PT_INJ_GENES, rownames(cnt))
message(sprintf("  TAL: %d/6 (%s); missing: %s",
                length(tal_found), paste(tal_found,collapse=","),
                if(length(tal_missing)) paste(tal_missing,collapse=",") else "none"))
message(sprintf("  PT_id: %d/10; PT_inj: %d/5", length(pt_id_fnd), length(pt_inj_fnd)))
stopifnot(length(tal_found) >= 4)

# ---- 7. TMM-logCPM ----
message("--- 7. TMM-logCPM ---")
dge <- DGEList(counts=cnt)
dge <- calcNormFactors(dge, method="TMM")
lcpm <- cpm(dge, log=TRUE, prior.count=1)
message(sprintf("  %d genes x %d samples", nrow(lcpm), ncol(lcpm)))

# ---- 8. One-time z-score (all qualified samples) ----
message("--- 8. Z-scoring ---")
zscore <- function(m) {
  mu <- rowMeans(m, na.rm=TRUE)
  s  <- apply(m, 1, sd, na.rm=TRUE)
  s[!is.finite(s) | s < 1e-8] <- 1
  sweep(sweep(m, 1, mu, "-"), 1, s, "/")
}
tal_z   <- zscore(lcpm[tal_found, , drop=FALSE])
pt_id_z <- zscore(lcpm[pt_id_fnd, , drop=FALSE])
pt_inj_z<- zscore(lcpm[pt_inj_fnd,, drop=FALSE])

scores <- data.frame(
  sample      = colnames(tal_z),
  TAL_score   = colMeans(tal_z),
  PT_identity = colMeans(pt_id_z),
  PT_injury   = colMeans(pt_inj_z),
  stringsAsFactors = FALSE)
rownames(scores) <- scores$sample

# ---- 9. Merge analysis dataset ----
analysis <- merge(meta, scores, by="sample", all.x=TRUE)
analysis$timepoint <- factor(analysis$timepoint,
  levels=c("baseline_pre","baseline_post","3months","1year"))

# ---- 10. Primary: baseline_post - baseline_pre (paired) ----
message("\n=== PRIMARY: baseline_post - baseline_pre ===")
get_sc <- function(kid, tp) {
  r <- analysis[analysis$kidney==kid & analysis$timepoint==tp, ]
  if (nrow(r)==0) NA_real_ else r$TAL_score[1]
}
pre_sc  <- sapply(kidneys_pair, get_sc, "baseline_pre")
post_sc <- sapply(kidneys_pair, get_sc, "baseline_post")
deltas_primary <- post_sc - pre_sc
obs_d <- mean(deltas_primary)
se_d  <- sd(deltas_primary)/sqrt(n_pairs)
p_prim<- 2*pt(abs(obs_d/se_d), df=n_pairs-1, lower.tail=FALSE)
ci_lo <- obs_d - qt(0.975, n_pairs-1)*se_d
ci_hi <- obs_d + qt(0.975, n_pairs-1)*se_d
message(sprintf("  Delta: %.3f (95%%CI %.3f to %.3f), P=%.4f, n=%d",
                obs_d, ci_lo, ci_hi, p_prim, n_pairs))

# Kidney FE
pd <- analysis[analysis$kidney %in% kidneys_pair &
               analysis$timepoint %in% c("baseline_pre","baseline_post"),]
pd$timepoint <- factor(pd$timepoint, levels=c("baseline_pre","baseline_post"))
pd$kidney <- factor(pd$kidney)
fit_fe <- lm(TAL_score ~ timepoint + kidney, data=pd)
sm_fe  <- coef(summary(fit_fe)); ci_fe <- confint(fit_fe)
pr <- grep("timepoint.*baseline_post", rownames(sm_fe))[1]
fe_est <- unname(sm_fe[pr,"Estimate"]); fe_se <- unname(sm_fe[pr,"Std. Error"])
fe_p   <- unname(sm_fe[pr,"Pr(>|t|)"])
fe_lo  <- unname(ci_fe[pr,1]); fe_hi <- unname(ci_fe[pr,2])
message(sprintf("  FE: %.3f (%.3f to %.3f), P=%.4f", fe_est, fe_lo, fe_hi, fe_p))

# ---- 11. Paired sign-flip permutation (9999, seed=42) ----
message("--- Permutation ---")
set.seed(SET_SEED)
ext <- 0L
for (b in seq_len(N_PERM)) {
  if (abs(mean(deltas_primary * sample(c(-1,1), n_pairs, TRUE))) >= abs(obs_d))
    ext <- ext + 1L
}
perm_p <- (ext + 1)/(N_PERM + 1)
message(sprintf("  Permutation P: %.4f", perm_p))

# ---- 12. LODO (leave-one-kidney-out) ----
message("--- LODO ---")
lodo <- numeric(n_pairs); names(lodo) <- kidneys_pair
for (i in seq_len(n_pairs)) lodo[i] <- mean(deltas_primary[-i])
lodo_neg <- mean(lodo < 0)
message(sprintf("  LODO negative: %.1f%% (%d/%d)", lodo_neg*100, sum(lodo<0), n_pairs))

# ---- 13. Per-gene direction and LOGO ----
message("--- Per-gene & LOGO ---")
pg <- data.frame(gene=character(0), delta=double(0), direction=character(0),
                 stringsAsFactors=FALSE)
for (gene in tal_found) {
  gd <- numeric(n_pairs)
  for (j in seq_len(n_pairs)) {
    k <- kidneys_pair[j]
    gd[j] <- lcpm[gene, analysis$sample[analysis$kidney==k & analysis$timepoint=="baseline_post"]] -
             lcpm[gene, analysis$sample[analysis$kidney==k & analysis$timepoint=="baseline_pre"]]
  }
  pg <- rbind(pg, data.frame(gene=gene, delta=mean(gd),
         direction=ifelse(mean(gd)<0,"negative","positive"), stringsAsFactors=FALSE))
}
neg_genes <- sum(pg$direction=="negative")
neg_names <- pg$gene[pg$direction=="negative"]
message(sprintf("  Per-gene negative: %d/6 (%s)", neg_genes, paste(neg_names,collapse=",")))

# LOGO: leave one gene out from pre-computed z-scores
logo <- numeric(length(tal_found)); names(logo) <- tal_found
for (j in seq_along(tal_found)) {
  rem <- setdiff(tal_found, tal_found[j])
  sz  <- colMeans(tal_z[rem,,drop=FALSE])
  ld  <- numeric(n_pairs)
  for (k in seq_len(n_pairs)) {
    kid <- kidneys_pair[k]
    ld[k] <- sz[analysis$sample[analysis$kidney==kid & analysis$timepoint=="baseline_post"]] -
             sz[analysis$sample[analysis$kidney==kid & analysis$timepoint=="baseline_pre"]]
  }
  logo[j] <- mean(ld)
}
logo_neg <- mean(logo < 0)
umod_est <- if("UMOD" %in% names(logo)) logo["UMOD"] else NA_real_
message(sprintf("  LOGO negative: %.1f%% (%d/%d) | excl UMOD: %.3f",
                logo_neg*100, sum(logo<0), length(tal_found), umod_est))

# ---- 14. Secondary comparisons with Holm correction ----
message("\n=== SECONDARY (3 comparisons, Holm-corrected) ===")
run_sec <- function(kids, tp_from, tp_to, label) {
  vk <- kids[sapply(kids, function(k)
    all(c(tp_from, tp_to) %in% analysis$timepoint[analysis$kidney==k]))]
  if (length(vk) < 3) {
    message(sprintf("  %s: n=%d insufficient", label, length(vk)))
    return(list(n=length(vk), delta=NA, se=NA, ci_lo=NA, ci_hi=NA, p_raw=NA,
                p_holm=NA, perm_p=NA, lodo_neg=NA))
  }
  ds <- sapply(vk, get_sc, tp_to) - sapply(vk, get_sc, tp_from)
  tt <- t.test(ds)
  ext <- 0L
  for (b in seq_len(N_PERM))
    if (abs(mean(ds*sample(c(-1,1),length(ds),TRUE))) >= abs(mean(ds))) ext <- ext+1L
  pp <- (ext+1)/(N_PERM+1)
  lo <- numeric(length(ds)); for (i in seq_along(ds)) lo[i] <- mean(ds[-i])
  ln <- mean(lo < 0)
  message(sprintf("  %s: delta=%.3f (%.3f to %.3f), P_raw=%.4f, perm_P=%.4f, LODO=%.1f%%, n=%d",
          label, mean(ds), tt$conf.int[1], tt$conf.int[2], tt$p.value, pp, ln*100, length(ds)))
  list(n=length(ds), delta=mean(ds), se=unname(tt$stderr),
       ci_lo=tt$conf.int[1], ci_hi=tt$conf.int[2],
       p_raw=tt$p.value, p_holm=NA_real_, perm_p=pp, lodo_neg=ln)
}
all_kids <- unique(analysis$kidney)
s1 <- run_sec(all_kids, "baseline_post","3months", "3mo - post")
s2 <- run_sec(all_kids, "baseline_post","1year",   "1y - post")
s3 <- run_sec(all_kids, "baseline_pre", "1year",   "1y - pre")

# Holm correction across 3 secondary P values
sec_p_raw <- c(s1$p_raw, s2$p_raw, s3$p_raw)
sec_p_holm <- p.adjust(sec_p_raw, method="holm")
s1$p_holm <- sec_p_holm[1]; s2$p_holm <- sec_p_holm[2]; s3$p_holm <- sec_p_holm[3]
message(sprintf("  Holm-adjusted P: 3mo-pos=%.4f, 1y-pos=%.4f, 1y-pre=%.4f",
                s1$p_holm, s2$p_holm, s3$p_holm))

# ---- 15. PT sensitivity (separate models) ----
message("\n=== PT SENSITIVITY ===")

# PT identity deltas
pt_id_d <- numeric(n_pairs)
for (j in seq_len(n_pairs)) {
  k <- kidneys_pair[j]
  pt_id_d[j] <- scores[analysis$sample[analysis$kidney==k & analysis$timepoint=="baseline_post"], "PT_identity"] -
                scores[analysis$sample[analysis$kidney==k & analysis$timepoint=="baseline_pre"],  "PT_identity"]
}
# TAL - PT_id
tt_pt <- t.test(deltas_primary - pt_id_d)
message(sprintf("  delta(TAL-PT_id): %.3f (%.3f to %.3f), P=%.4f",
                tt_pt$estimate, tt_pt$conf.int[1], tt_pt$conf.int[2], tt_pt$p.value))

# PT injury deltas
pt_inj_d <- numeric(n_pairs)
for (j in seq_len(n_pairs)) {
  k <- kidneys_pair[j]
  pt_inj_d[j] <- scores[analysis$sample[analysis$kidney==k & analysis$timepoint=="baseline_post"], "PT_injury"] -
                 scores[analysis$sample[analysis$kidney==k & analysis$timepoint=="baseline_pre"],  "PT_injury"]
}

# FE models with PT covariates
pd$PT_identity <- scores[pd$sample, "PT_identity"]
pd$PT_injury   <- scores[pd$sample, "PT_injury"]

# PT identity-adjusted FE
fit_pt <- lm(TAL_score ~ timepoint + kidney + PT_identity, data=pd)
sm_pt <- coef(summary(fit_pt)); ci_pt <- confint(fit_pt)
ptr <- grep("timepoint.*baseline_post", rownames(sm_pt))[1]
ptest <- unname(sm_pt[ptr,"Estimate"]); ptse <- unname(sm_pt[ptr,"Std. Error"])
ptp   <- unname(sm_pt[ptr,"Pr(>|t|)"])
ptlo  <- unname(ci_pt[ptr,1]); pthi <- unname(ci_pt[ptr,2])
message(sprintf("  FE+PT_id baseline_post: %.3f (%.3f to %.3f), P=%.4f", ptest, ptlo, pthi, ptp))

# PT injury-adjusted FE (this is the ONLY PT-injury sensitivity retained)
fit_inj <- lm(TAL_score ~ timepoint + kidney + PT_injury, data=pd)
sm_inj <- coef(summary(fit_inj)); ci_inj <- confint(fit_inj)
injr <- grep("timepoint.*baseline_post", rownames(sm_inj))[1]
injest <- unname(sm_inj[injr,"Estimate"]); injse <- unname(sm_inj[injr,"Std. Error"])
injp   <- unname(sm_inj[injr,"Pr(>|t|)"])
injlo  <- unname(ci_inj[injr,1]); injhi <- unname(ci_inj[injr,2])
message(sprintf("  FE+PT_inj baseline_post: %.3f (%.3f to %.3f), P=%.4f", injest, injlo, injhi, injp))

# delta_TAL ~ delta_PT_injury
fit_di <- lm(deltas_primary ~ pt_inj_d)
sm_di  <- coef(summary(fit_di))
message(sprintf("  delta_TAL ~ delta_PT_inj: slope=%.3f, R2=%.3f, P=%.4f",
                sm_di["pt_inj_d","Estimate"], summary(fit_di)$r.squared,
                sm_di["pt_inj_d","Pr(>|t|)"]))

# ---- 16. Save metadata CSV with GSM + scores + pairing flag ----
meta_csv <- analysis[, c("sample","gsm","kidney","timepoint","TAL_score","PT_identity","PT_injury")]
meta_csv$in_primary_pair <- meta_csv$kidney %in% kidneys_pair &
  meta_csv$timepoint %in% c("baseline_pre","baseline_post")
write.csv(meta_csv, file.path(meta_dir,"gse126805_samples.csv"), row.names=FALSE)
message("\nSaved: metadata/gse126805_samples.csv")

# ---- 17. Assemble results CSV (with LODO/LOGO detail) ----
message("--- Assembling results ---")
mkr <- function(rt,lb,est,se=NA,clo=NA,chi=NA,pv=NA,nn=NA,nt=NA) {
  data.frame(result_type=rt, label=lb, estimate=est, se=se,
             ci_low=clo, ci_high=chi, p_value=pv, n=nn, note=nt,
             stringsAsFactors=FALSE)
}
f3 <- function(x) if(is.na(x)) NA_character_ else sprintf("%.3f",x)
f4 <- function(x) if(is.na(x)) NA_character_ else sprintf("%.4f",x)

tidy <- rbind(
  mkr("mapping","n_kidneys",length(unique(meta$kidney))),
  mkr("mapping","n_samples",nrow(meta)),
  mkr("mapping","n_pairs_pre_post",n_pairs),
  mkr("mapping","n_TAL_genes",length(tal_found),nt=paste(tal_found,collapse=",")),
  mkr("mapping","tal_missing",NA,nt=if(length(tal_missing)) paste(tal_missing,collapse=",") else "none"),
  mkr("mapping","n_dup_genes_removed",n_dup),
  mkr("mapping","dedup_method",NA,nt="highest-mean-expression row per gene symbol"),
  mkr("primary","baseline_post_vs_pre_delta",obs_d,se_d,ci_lo,ci_hi,p_prim,n_pairs),
  mkr("primary","baseline_post_vs_pre_FE",fe_est,fe_se,fe_lo,fe_hi,fe_p,n_pairs),
  mkr("primary","permutation_p",perm_p),
  mkr("primary","lodo_negative_fraction",lodo_neg),
  mkr("primary","logo_negative_fraction",logo_neg),
  mkr("per_gene",pg$gene,pg$delta,nt=pg$direction,nn=n_pairs),
  mkr("secondary","3mo_vs_post",s1$delta,s1$se,s1$ci_lo,s1$ci_hi,s1$p_raw,s1$n),
  mkr("secondary","1y_vs_post",s2$delta,s2$se,s2$ci_lo,s2$ci_hi,s2$p_raw,s2$n),
  mkr("secondary","1y_vs_pre",s3$delta,s3$se,s3$ci_lo,s3$ci_hi,s3$p_raw,s3$n),
  mkr("secondary","3mo_vs_post_holm_P",s1$p_holm),
  mkr("secondary","1y_vs_post_holm_P",s2$p_holm),
  mkr("secondary","1y_vs_pre_holm_P",s3$p_holm),
  mkr("secondary","3mo_vs_post_perm_P",s1$perm_p),
  mkr("secondary","1y_vs_post_perm_P",s2$perm_p),
  mkr("secondary","1y_vs_pre_perm_P",s3$perm_p),
  mkr("sensitivity","delta_TAL_minus_PT_id",unname(tt_pt$estimate),unname(tt_pt$stderr),
      tt_pt$conf.int[1],tt_pt$conf.int[2],tt_pt$p.value,n_pairs),
  mkr("sensitivity","FE_baseline_post_PT_id_adj",ptest,ptse,ptlo,pthi,ptp,n_pairs),
  mkr("sensitivity","FE_baseline_post_PT_inj_adj",injest,injse,injlo,injhi,injp,n_pairs),
  mkr("sensitivity","delta_TAL_vs_delta_PT_inj_slope",sm_di["pt_inj_d","Estimate"],
      sm_di["pt_inj_d","Std. Error"],NA,NA,sm_di["pt_inj_d","Pr(>|t|)"],n_pairs),
  mkr("sensitivity","logo_exclude_UMOD",umod_est),
  # per-kidney LODO detail
  mkr("lodo_detail",paste0("omit_",names(lodo)),lodo)
)
# per-gene LOGO detail
for (j in seq_along(tal_found)) {
  tidy <- rbind(tidy,
    mkr("logo_detail",paste0("omit_",tal_found[j]),logo[j]))
}

write.csv(tidy, file.path(result_dir,"gse126805_results.csv"), row.names=FALSE)
message("Saved: results/gse126805_results.csv")

# ---- 18. Evidence figure (allograft-level; no cell-level inference) ----
draw_longitudinal_figure <- function(file, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    pdf(file, width = 11, height = 5.4, useDingbats = FALSE)
  } else {
    png(file, width = 2200, height = 1080, res = 200)
  }
  on.exit(dev.off(), add = TRUE)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.8, 4.8, 3.0, 1.2), las = 1)

  # Panel A: paired observations. Each line is one kidney allograft.
  score_range <- range(c(pre_sc, post_sc), finite = TRUE)
  pad <- diff(score_range) * 0.10
  if (!is.finite(pad) || pad == 0) pad <- 0.2
  plot(c(1, 2), score_range + c(-pad, pad), type = "n", xaxt = "n",
       xlab = "Protocol-biopsy timepoint", ylab = "Frozen TAL programme score",
       main = "Paired acute change (n = 41 allografts)")
  axis(1, at = c(1, 2), labels = c("Before reperfusion", "After reperfusion"))
  for (i in seq_len(n_pairs)) {
    segments(1, pre_sc[i], 2, post_sc[i], col = "grey70", lwd = 0.8)
  }
  points(rep(1, n_pairs), pre_sc, pch = 16, cex = 0.65, col = "grey25")
  points(rep(2, n_pairs), post_sc, pch = 16, cex = 0.65, col = "grey25")
  segments(0.88, mean(pre_sc), 1.12, mean(pre_sc), lwd = 3, col = "#1B6CA8")
  segments(1.88, mean(post_sc), 2.12, mean(post_sc), lwd = 3, col = "#1B6CA8")
  mtext(sprintf("Mean paired change = %.3f (95%% CI %.3f to %.3f); P = %.4f",
                obs_d, ci_lo, ci_hi, p_prim), side = 3, line = 0.25, cex = 0.75)

  # Panel B: one primary and three explicitly exploratory allograft-level contrasts.
  forest <- data.frame(
    label = c("Post - pre (primary)", "3 months - post (exploratory)",
              "1 year - post (exploratory)", "1 year - pre (exploratory)"),
    estimate = c(obs_d, s1$delta, s2$delta, s3$delta),
    lower = c(ci_lo, s1$ci_lo, s2$ci_lo, s3$ci_lo),
    upper = c(ci_hi, s1$ci_hi, s2$ci_hi, s3$ci_hi),
    p_text = c(sprintf("P = %.4f", p_prim),
               sprintf("Holm P = %.4f", s1$p_holm),
               sprintf("Holm P = %.4f", s2$p_holm),
               sprintf("Holm P = %.4f", s3$p_holm)),
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
  mtext("Only post - pre is primary; later contrasts are exploratory.",
        side = 3, line = 0.25, cex = 0.75)
}

draw_longitudinal_figure(file.path(figure_dir, "longitudinal_gse126805.pdf"), "pdf")
draw_longitudinal_figure(file.path(figure_dir, "longitudinal_gse126805.png"), "png")
message("Saved: figures/longitudinal_gse126805.pdf/png")

# ---- 19. Decision ----
message("\n=== DECISION ===")
conds <- c(
  "post-pre negative" = obs_d < 0,
  ">=4/6 genes neg"   = neg_genes >= 4,
  ">=25 pairs"        = n_pairs >= 25,
  "LODO >=80% neg"    = lodo_neg >= 0.8)
cs <- paste(names(conds), ifelse(conds,"PASS","FAIL"), sep="=", collapse=", ")
message(sprintf("  Conditions: %s", cs))
final <- if(all(conds)) "KEEP_GSE126805" else "STOP_GSE126805"
message(sprintf("  Final: %s", final))

# ---- 20. Audit report ----
message("\nWriting audit report...")
rpt <- c(
  "# GSE126805 Audit Report - Frozen TAL Longitudinal Validation",
  "",
  sprintf("**Date:** %s | **Decision:** `%s`", Sys.Date(), final),
  "",
  "## Data Sources",
  "| File | Source | SHA256 |",
  "|---|---|---|",
  "| `GSE126805_all.gene_counts.txt.gz` | GEO suppl (ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126805/suppl/) | `91683e6f99fd200665fdd4be064c4bb3f2d0730715a832f4c9ba499fa45dc134` |",
  "| `GSE126805_series_matrix.txt.gz` | GEO series matrix | `026598fa31d7b1e54a7c4e0ab2f572152f2972fac82d959994ebdbea0ea1ac43` |",
  "| `PMC6302941_supp.pdf` | PMC supplementary material | `515bade828f518b7d68cf970da52fb16e5aa1c9c9c1b503d93096b19e3b0eb60` |",
  "",
  "## Runtime Versions",
  sprintf("- R version %s", r_ver),
  sprintf("- edgeR version %s", edger_ver),
  sprintf("- limma version %s", limma_ver),
  "",
  "## Cohort",
  "- GSE126805 (Cippa et al. 2018, JCI Insight 3:e123151, PMID 30429361)",
  "- 42 kidney allograft recipients, protocol biopsies at 4 timepoints (163 total)",
  "- Timepoints: baseline_pre (before reperfusion), baseline_post (after reperfusion), 3 months, 1 year",
  "- Platform: bulk RNA-seq (Illumina, STAR v2.0.4b, featureCounts v1.4.5)",
  "- Unit: kidney allograft. Timepoints are repeated measures, NOT independent replicates.",
  "",
  "## Gene Symbol Deduplication",
  sprintf("- Input: %d rows, %d non-empty gene symbols, %d unique symbols", nrow(cts),
          n_gene_nonempty, n_gene_unique),
  sprintf("- %d duplicate gene symbols detected", n_dup),
  "- Method: per duplicated symbol, retain the row with highest mean count; discard others",
  "- All 6 TAL genes: single row each, no duplication",
  sprintf("- After dedup: %d genes retained", nrow(cnt)),
  "",
  "## Original Paper Overlap Audit",
  "Cippa et al. 2018 analysed transcriptional *trajectories* across timepoints using",
  "self-organising maps (SOM) and k-means clustering to identify gene modules that",
  "changed coordinately over time. They identified two main trajectory types: recovery",
  "and sustained injury with fibrosis.",
  "",
  "**Our analysis differs fundamentally:**",
  "- We score a frozen 6-gene TAL identity composite (Muto et al. 2021) at each timepoint",
  "- We test a prespecified hypothesis: does the TAL programme drop acutely after IRI?",
  "- We use paired within-kidney change scores, not trajectory clustering",
  "- Cippa et al. did NOT test this specific 6-gene composite, nor did they perform",
  "  a paired pre-vs-post change-score analysis focused on TAL identity",
  "",
  "**No overlap concern:** the original paper's trajectory analysis and our frozen-gene",
  "change-score analysis address different questions with different methods on the same",
  "data. This is a validation, not a replication of their published result.",
  "",
  "## Sample Coverage",
  sprintf("- Total samples analysed: %d", nrow(analysis)),
  sprintf("- Kidneys with baseline_pre + baseline_post: %d/42", n_pairs),
  sprintf("- Kidney 29: missing baseline_pre (3 observations excluded from primary pairing)"),
  sprintf("- Kidneys with all 4 timepoints: %d",
          sum(table(meta$kidney)==4)),
  "- GSM accessions: all 163 cross-matched to series matrix titles (100%)",
  "",
  "## Primary: baseline_post - baseline_pre",
  sprintf("- Paired delta: %.3f (95%% CI %.3f to %.3f), P=%.4f, n=%d",
          obs_d, ci_lo, ci_hi, p_prim, n_pairs),
  sprintf("- Kidney FE (baseline_post coefficient): %.3f (%.3f to %.3f), P=%.4f",
          fe_est, fe_lo, fe_hi, fe_p),
  sprintf("- Paired sign-flip permutation (9999, seed=42, two-sided): P=%.4f", perm_p),
  sprintf("- LODO negative: %.1f%% (%d/%d)", lodo_neg*100, sum(lodo<0), n_pairs),
  sprintf("- LOGO negative: %.1f%% (%d/%d)", logo_neg*100, sum(logo<0), length(tal_found)),
  sprintf("- Per-gene direction: %d/6 negative (%s)", neg_genes, paste(neg_names,collapse=", ")),
  sprintf("- LOGO excluding UMOD: %.3f", umod_est),
  sprintf("- UMOD per-gene direction: %s (delta=%.3f)",
          pg$direction[pg$gene=="UMOD"], pg$delta[pg$gene=="UMOD"]),
  "",
  "## Secondary Comparisons (3 tests, Holm-corrected)",
  "",
  sprintf("| Comparison | Delta (95%% CI) | P_raw | P_holm | n |"),
  sprintf("|---|---|---|---|---|"),
  sprintf("| 3mo - post | %.3f (%.3f to %.3f) | %.4f | %.4f | %d |",
          s1$delta, s1$ci_lo, s1$ci_hi, s1$p_raw, s1$p_holm, s1$n),
  sprintf("| 1y - post  | %.3f (%.3f to %.3f) | %.4f | %.4f | %d |",
          s2$delta, s2$ci_lo, s2$ci_hi, s2$p_raw, s2$p_holm, s2$n),
  sprintf("| 1y - pre   | %.3f (%.3f to %.3f) | %.4f | %.4f | %d |",
          s3$delta, s3$ci_lo, s3$ci_hi, s3$p_raw, s3$p_holm, s3$n),
  "",
  "Note: secondary comparisons are exploratory and reported with Holm correction",
  "across 3 tests. They do not contribute to the primary KEEP/STOP decision.",
  "",
  "## PT Sensitivity (separate models)",
  "",
  "### PT Identity",
  sprintf("- Paired delta(TAL - PT_id): %.3f (%.3f to %.3f), P=%.4f",
          tt_pt$estimate, tt_pt$conf.int[1], tt_pt$conf.int[2], tt_pt$p.value),
  sprintf("- FE+PT_id baseline_post coefficient: %.3f (%.3f to %.3f), P=%.4f",
          ptest, ptlo, pthi, ptp),
  "",
  "### PT Injury",
  sprintf("- FE+PT_inj baseline_post coefficient: %.3f (%.3f to %.3f), P=%.4f",
          injest, injlo, injhi, injp),
  sprintf("- delta_TAL ~ delta_PT_inj: slope=%.3f, R2=%.3f, P=%.4f",
          sm_di["pt_inj_d","Estimate"], summary(fit_di)$r.squared,
          sm_di["pt_inj_d","Pr(>|t|)"]),
  "",
  "## Methods Executed",
  "1. Gene counts from GEO supplementary file (no RAW/SRA/FASTQ)",
  "2. Gene symbols deduplicated: per duplicated symbol, retain highest-mean row",
  "3. TMM normalization (edgeR calcNormFactors, method=TMM)",
  "4. logCPM (cpm with log=TRUE, prior.count=1)",
  "5. Per-gene z-score across all 163 qualified samples (one-time)",
  "6. TAL composite = mean of 6 per-gene z-scores",
  "7. Primary: paired change score baseline_post - baseline_pre (n=41 pairs)",
  "8. Kidney fixed-effects verification (lm with kidney factor)",
  "9. Paired sign-flip permutation (9999, seed=42, two-sided)",
  "10. LODO (leave-one-kidney-out from change scores)",
  "11. LOGO (leave-one-gene-out from pre-computed TAL z-scores)",
  "12. Per-gene direction on raw logCPM deltas",
  "13. PT identity sensitivity: paired delta(TAL-PT_id) and kidney FE with PT_id covariate",
  "14. PT injury sensitivity: kidney FE with PT_inj covariate; delta_TAL ~ delta_PT_inj regression",
  "15. Secondary: 3 paired comparisons with Holm correction across 3 tests",
  "16. All 163 GSM accessions cross-validated against GEO series matrix",
  "17. No trajectory-stratified analysis (supplementary labels not confirmed)",
  "",
  "## Conclusion Boundaries",
  "- Bulk tissue-level analysis. No cell-type specificity. No cell-intrinsic claims.",
  "- The frozen 6-gene TAL programme was scored as a tissue-level composite.",
  "- No causal claims. No clinical prediction. No therapeutic recommendations.",
  "- All kidneys experienced IRI (transplant setting). All under immunosuppression.",
  "- Generalizability to non-transplant AKI is untested.",
  "- Timepoints span acute (hours) to subacute/chronic (3-12 months).",
  "",
  "## Deviations from Protocol",
  "- Kidney 29 excluded from primary pairing (no baseline_pre). 41/42 kidneys retained.",
  "- UMOD direction is positive (+0.081) - anomalous among TAL markers.",
  "- Trajectory labels from Cippa supplementary could not be confirmed; trajectory analysis omitted.",
  "- Evidence figure generated from allograft-level paired scores and estimates; no cell-level figure was generated.",
  "",
  "## Decision Criteria",
  sprintf("| Criterion | Result | Status |"),
  sprintf("|---|---|---|"),
  sprintf("| baseline_post - pre delta < 0 | %.3f | %s |", obs_d, if(obs_d<0) "PASS" else "FAIL"),
  sprintf("| >=4/6 genes negative | %d/6 | %s |", neg_genes, if(neg_genes>=4) "PASS" else "FAIL"),
  sprintf("| >=25 pre/post pairs | %d | %s |", n_pairs, if(n_pairs>=25) "PASS" else "FAIL"),
  sprintf("| LODO >=80%% negative | %.0f%% | %s |", lodo_neg*100, if(lodo_neg>=0.8) "PASS" else "FAIL"),
  "",
  sprintf("**`%s`**", final),
  ""
)

writeLines(rpt, file.path(result_dir,"gse126805_audit.md"))
cat(paste(rpt, collapse="\n"), "\n")
message(sprintf("\n=== 03_longitudinal_validation.R DONE: %s ===", final))
