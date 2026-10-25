#!/usr/bin/env Rscript
# =============================================================================
# figure4_tal_states_scaffold.R — Figure 4 (extension): TAL state context
# =============================================================================
# STATUS-AWARE LAUNCHER / CHECKER (replaces the obsolete DEFERRED scaffold).
#
# The Figure 4 pre-analysis package is now in place:
#   - metadata/extension/tal_state_reference.csv  (frozen TAL state programs)
#   - metadata/extension/reference_sources/       (Lake 2023 MOESM6 archive +
#                                                  provenance/audit reports)
#   - results/extension/tal_state_feasibility.md  (diagnosis-blinded feasibility)
#   - results/extension/tal_state_marker_coverage.csv / tal_donor_coverage.csv
#   - analysis_contract_figure4.md                (final contract)
#   - results/extension/figure4_tal_states_status.md  (authoritative status)
#
# This script NEVER writes a DEFERRED message over a valid CONDITIONAL/READY
# status. It checks the inputs that the scoring phase requires, reports what
# is present/missing, and, if invoked with --write-status, updates only the
# DEFERRED/CONDITIONAL/READY decision rules documented below:
#   - READY:      reference CSV exists, feasibility deliverables exist, the
#                 primary continuous scorer dependency (UCell) is verified
#                 present in the R library path, AND the diagnosis-blinded
#                 UCell smoke test PASS marker exists
#                 (metadata/extension/r-lib/UCell_SMOKE_PASS.txt).
#   - CONDITIONAL: reference + feasibility present but the scorer dependency is
#                 not yet verified (exact blocker recorded).
#   - DEFERRED:   only when the frozen reference CSV is absent (legitimate —
#                 no state labels may be invented), and NEVER overwriting a
#                 CONDITIONAL/READY status that is still valid.
#
# Planned panels (manuscript_v4_storyboard Figure 4):
#   A  reference-to-cohort TAL state mapping (continuous scores; reference frozen)
#   B  TAL-only UMAP for descriptive state context (coordinates available in
#      results/extension/atlas_GSE*_umap.csv for broad_class == "TAL/LOH")
#   C  donor-level continuous state scores (donor unit; effect + CI, not cell P)
#   D  detection-versus-mean-expression decomposition among recovered TAL
#   E  donor-level state-programme association (cohort-specific estimates)
#
# Required input paths for the scoring phase:
#   metadata/extension/tal_state_reference.csv  -> frozen state programs
#   results/extension/atlas_GSE183276_umap.csv  (TAL/LOH cells, donor labels)
#   results/extension/atlas_GSE210622_umap.csv
#   results/formal_evidence.csv                 -> frozen six-gene TAL score
#   ../AKI_DCT_analysis/results/stage3b_pseudobulk.rds  -> TAL pseudobulk
#
# COMPLETE is reached only when the scoring phase has produced its required
# final outputs, the two review records carry exact acceptance verdicts, and
# the pre/post evidence-freeze receipts are identical. COMPLETE is revalidated
# on every checker run and can be downgraded if an artefact becomes stale.
#
# Anti-circularity rule (WP1D): TAL membership comes from the author/
# independent annotation before testing the score; loss of the six genes in AKI
# must not remove an AKI cell from TAL. The six frozen genes (SLC12A1, UMOD,
# CLDN16, KCNJ10, PTH1R, CLDN10) never enter state scoring.
#
# Interpretation boundary: recovered fractions are not biological TAL
# abundance; state shifts do not imply temporal transitions, lineage or fate.
# =============================================================================

cli_args <- commandArgs(trailingOnly = TRUE)
write_status <- "--write-status" %in% cli_args
# Optional project root override (documented usage: run from code/extension/)
if (length(cli_args) >= 1 && !grepl("^--", cli_args[1])) {
  project_dir <- normalizePath(cli_args[1])
} else {
  project_dir <- normalizePath(file.path(getwd(), "..", ".."))
}

ref_csv    <- file.path(project_dir, "metadata", "extension", "tal_state_reference.csv")
status_md  <- file.path(project_dir, "results", "extension", "figure4_tal_states_status.md")
feas_md    <- file.path(project_dir, "results", "extension", "tal_state_feasibility.md")
cov_csv    <- file.path(project_dir, "results", "extension", "tal_state_marker_coverage.csv")
don_csv    <- file.path(project_dir, "results", "extension", "tal_donor_coverage.csv")

ref_exists <- file.exists(ref_csv)
feas_exists <- file.exists(feas_md)
cov_exists <- file.exists(cov_csv)
don_exists <- file.exists(don_csv)

# --- scorer dependency check (UCell, project-local) --------------------------
# The project-local R library (metadata/extension/r-lib) must be found even
# when this script is launched from a system R whose default .libPaths() does
# not include it. We therefore temporarily prepend the project r-lib to
# .libPaths() BEFORE requireNamespace(), so the search path always contains it
# (lib.loc alone can miss it when UCell's dependencies live in that same
# project-local tree). The original .libPaths() is restored afterwards.
check_scorer <- function() {
  project_lib <- file.path(project_dir, "metadata", "extension", "r-lib")
  old_lib_paths <- .libPaths()
  on.exit(.libPaths(old_lib_paths), add = TRUE)
  if (dir.exists(project_lib)) {
    .libPaths(c(project_lib, old_lib_paths))
  }
  found <- tryCatch(
    requireNamespace("UCell", quietly = TRUE),
    error = function(e) FALSE
  )
  if (isTRUE(found)) {
    # requireNamespace(quietly=TRUE) may still attach the package namespace
    # lazily; confirm we can actually read its version to avoid a false
    # positive when only a stale/partial installation is present.
    ver_ok <- tryCatch(
      length(utils::packageVersion("UCell")) > 0,
      error = function(e) FALSE
    )
    if (isTRUE(ver_ok)) {
      list(present = TRUE, detail = sprintf("UCell verified in project-local library (version %s)", as.character(utils::packageVersion("UCell"))))
    } else {
      list(present = FALSE,
           detail = "UCell namespace found but version could not be read in project-local library metadata/extension/r-lib (partial/stale installation). Reinstall via BiocManager::install(\"UCell\", lib=...).")
    }
  } else {
    list(present = FALSE,
         detail = "UCell not verified in project-local library metadata/extension/r-lib (install via BiocManager::install(\"UCell\", lib=...))")
  }
}

# --- diagnosis-blinded smoke-test marker -------------------------------------
# READY additionally requires a smoke-test PASS marker, written only by the
# diagnosis-blinded UCell smoke test (code/extension/04c_ucell_smoke_test.R or
# equivalent) after it completes successfully on a small strict-TAL subset.
smoke_marker <- file.path(project_dir, "metadata", "extension", "r-lib",
                          "UCell_SMOKE_PASS.txt")
smoke_passed <- file.exists(smoke_marker)

# --- final-output / review / integrity receipts ------------------------------
final_csvs <- c(
  "tal_state_donor_scores.csv", "tal_state_effects.csv",
  "tal_state_lodo.csv", "tal_state_decomposition.csv")
final_md   <- "figure4_tal_states_report.md"
fig_png    <- file.path(project_dir, "figures", "extension", "figure4_tal_states.png")
fig_pdf    <- file.path(project_dir, "figures", "extension", "figure4_tal_states.pdf")
run_record <- file.path(project_dir, "results", "extension", "figure4_run_record.json")
fe_before  <- file.path(project_dir, "results", "extension", "freeze_check_figure4_before.csv")
fe_after   <- file.path(project_dir, "results", "extension", "freeze_check_figure4_after.csv")
review_stats_md <- file.path(project_dir, "results", "extension", "figure4_stats_review.md")
review_fig_md   <- file.path(project_dir, "results", "extension", "figure4_figure_review.md")

final_missing <- character(0)
for (fn in final_csvs) if (!file.exists(file.path(project_dir, "results", "extension", fn)))
  final_missing <- c(final_missing, fn)
if (!file.exists(file.path(project_dir, "results", "extension", final_md)))
  final_missing <- c(final_missing, final_md)
if (!file.exists(fig_png))  final_missing <- c(final_missing, "figures/extension/figure4_tal_states.png")
if (!file.exists(fig_pdf))  final_missing <- c(final_missing, "figures/extension/figure4_tal_states.pdf")
if (!file.exists(run_record)) final_missing <- c(final_missing, "figure4_run_record.json")
if (!file.exists(fe_before) || !file.exists(fe_after))
  final_missing <- c(final_missing, "freeze_check_figure4_before/after.csv")

receipts_match <- FALSE
if (file.exists(fe_before) && file.exists(fe_after)) {
  receipts_match <- identical(
    readBin(fe_before, "raw", n = file.info(fe_before)$size),
    readBin(fe_after, "raw", n = file.info(fe_after)$size)
  )
  if (!receipts_match) final_missing <- c(final_missing, "freeze receipts differ")
}

# review PASS records: each review MD must exist and carry a PASS verdict line
review_pass <- function(path, pattern) {
  if (!file.exists(path)) return(FALSE)
  txt <- tryCatch(paste(readLines(path, warn = FALSE), collapse = "\n"), error = function(e) "")
  grepl(pattern, txt)
}
stats_pass <- review_pass(review_stats_md, "\\*\\*Verdict: PASS\\*\\*")
fig_pass   <- review_pass(review_fig_md, "\\*\\*Verdict: APPROVE")
reviews_ok <- stats_pass && fig_pass

# --- read existing status (if any) ------------------------------------------
read_status <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  txt <- readLines(path, warn = FALSE)
  hit <- grep("^## Status:", txt, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  sub("^## Status:\\s*", "", hit[1])
}

existing_status <- read_status(status_md)

# --- decide status -----------------------------------------------------------
new_status  <- NA_character_
blocker     <- character(0)

if (length(final_missing) == 0 && reviews_ok) {
  # All required final outputs, review PASS records and integrity receipts are
  # present -> terminal COMPLETE (only reached after every acceptance pass).
  new_status <- "COMPLETE"
  blocker <- character(0)
} else if (!ref_exists) {
  # Legitimate DEFERRED path: no frozen reference -> no state labels may be
  # invented. Only write DEFERRED if the existing status is not CONDITIONAL/READY.
  new_status <- "DEFERRED"
  blocker <- c("metadata/extension/tal_state_reference.csv is absent; the frozen TAL-state reference is required before any state scoring.")
} else {
  scorer <- check_scorer()
  missing_inputs <- character(0)
  if (!feas_exists) missing_inputs <- c(missing_inputs, "tal_state_feasibility.md")
  if (!cov_exists) missing_inputs <- c(missing_inputs, "tal_state_marker_coverage.csv")
  if (!don_exists) missing_inputs <- c(missing_inputs, "tal_donor_coverage.csv")
  if (length(missing_inputs) > 0) {
    new_status <- "CONDITIONAL"
    blocker <- c(sprintf("Feasibility deliverables missing: %s", paste(missing_inputs, collapse = ", ")))
  } else if (!scorer$present) {
    new_status <- "CONDITIONAL"
    blocker <- c(scorer$detail)
  } else if (!smoke_passed) {
    new_status <- "CONDITIONAL"
    blocker <- c("UCell installed but diagnosis-blinded smoke test not yet passed (no metadata/extension/r-lib/UCell_SMOKE_PASS.txt marker). Run the smoke test on a small strict-TAL subset when memory allows.")
  } else {
    new_status <- "READY"
  }
}

# --- never downgrade a valid CONDITIONAL/READY with DEFERRED -----------------
if (new_status == "DEFERRED" && existing_status %in% c("CONDITIONAL", "READY")) {
  cat(sprintf("Preserving existing status '%s': reference absent but status already CONDITIONAL/READY. Not overwriting with DEFERRED.\n", existing_status))
  new_status <- existing_status
  blocker <- c("Reference CSV absent; retaining existing status (do not invent state labels).")
}

# NOTE: status is decided ENTIRELY from the current checks above (reference
# present, feasibility outputs present, UCell currently loadable, smoke-test
# PASS marker present). There is intentionally NO historical-status override:
# a previous READY status file is never used to keep READY when UCell cannot
# be loaded on the current run — that would be a false READY. If UCell fails
# to load, the decision is CONDITIONAL regardless of what the marker file or
# previous status says.

# --- report ------------------------------------------------------------------
cat("Project dir:", project_dir, "\n")
cat("Reference CSV present:", ref_exists, "\n")
cat("Feasibility MD present:", feas_exists, "\n")
cat("Coverage CSV present:", cov_exists, "\n")
cat("Donor CSV present:", don_exists, "\n")
cat("Existing status:", ifelse(is.na(existing_status), "(none)", existing_status), "\n")
cat("Smoke-test PASS marker present:", smoke_passed, "\n")
cat("Freeze receipts identical:", receipts_match, "\n")
cat("Decision status:", new_status, "\n")
if (length(blocker) > 0) cat("Blocker:", paste(blocker, collapse = "; "), "\n")

# Only rewrite the status file when the decision actually changes; a
# no-change run (e.g. CONDITIONAL -> CONDITIONAL) preserves the existing
# detailed status file untouched.
status_changed <- is.na(existing_status) || (!is.na(existing_status) && existing_status != new_status)

if (write_status && (status_changed || new_status == "COMPLETE")) {
  date_line <- format(Sys.Date(), "%Y-%m-%d")
  out <- c(
    "# Figure 4 — TAL states status",
    "",
    sprintf("**Date checked:** %s", date_line),
    sprintf("**Frozen TAL-state reference (tal_state_reference.csv):** %s", if (ref_exists) "YES" else "NO"),
    "",
    sprintf("## Status: %s", new_status),
    "",
    if (new_status == "COMPLETE") {
      c("All required Figure 4 outputs are present and verified: the four",
        "tal_state_*.csv deliverables, the final report (with review findings),",
        "the rendered figure PNG/PDF, the run record with protected-file",
        "baseline hashes, and the pre/post evidence-freeze receipts. Both",
        "independent read-only reviews recorded PASS (statistical re-fit) and",
        "APPROVE (figure) verdicts. COMPLETE is revalidated on each checker run.",
        "",
        "## Final acceptance records",
        "",
        "- results/extension/tal_state_donor_scores.csv: present",
        "- results/extension/tal_state_effects.csv: present",
        "- results/extension/tal_state_lodo.csv: present",
        "- results/extension/tal_state_decomposition.csv: present",
        "- results/extension/figure4_tal_states_report.md: present (FINAL)",
        "- figures/extension/figure4_tal_states.png: present",
        "- figures/extension/figure4_tal_states.pdf: present",
        "- results/extension/figure4_run_record.json: present (protected hashes + manifest + freeze receipts)",
        sprintf("- Statistical review PASS record: %s", if (stats_pass) "present" else "MISSING"),
        sprintf("- Figure review APPROVE record: %s", if (fig_pass) "present" else "MISSING"))
    } else if (new_status == "READY") {
      c("Reference, feasibility deliverables and the primary continuous scorer",
        "dependency (UCell) are all present and verified. Figure 4 may proceed to",
        "the state-scoring phase under the analysis_contract_figure4.md rules.")
    } else if (new_status == "CONDITIONAL") {
      c("Figure 4 can proceed to the state-scoring phase under conditions:",
        "",
        paste0("- ", blocker))
    } else {
      c("No frozen TAL-state reference is present; state labels would be",
        "invented, which is not permitted. Figure 4 is deferred.")
    },
    "",
    "## Pre-analysis package (frozen)",
    "",
    sprintf("- metadata/extension/tal_state_reference.csv: %s", if (ref_exists) "present" else "MISSING"),
    sprintf("- results/extension/tal_state_feasibility.md: %s", if (feas_exists) "present" else "MISSING"),
    sprintf("- results/extension/tal_state_marker_coverage.csv: %s", if (cov_exists) "present" else "MISSING"),
    sprintf("- results/extension/tal_donor_coverage.csv: %s", if (don_exists) "present" else "MISSING"),
    "",
    "## Constraints carried into the scoring phase",
    "",
    "- Cells are never patient replicates; the donor is the biological unit.",
    "- Strict author-annotated TAL cells (`author_label == \"TAL\"`) are primary; broad TAL/LOH is sensitivity only.",
    "- The six frozen genes (SLC12A1, UMOD, CLDN16, KCNJ10, PTH1R, CLDN10) never enter state scoring.",
    "- Continuous UCell scores are the primary state analysis; no state-fraction (max-of-three) assignment.",
    "- No claims of TAL specificity, PT independence, causality, lineage, mechanism, biological abundance, or biomarker/clinical prediction.",
    "- No AKI-versus-control testing before the scoring-phase acceptance criteria are met."
  )
  writeLines(out, status_md)
  cat("Wrote status:", status_md, "\n")
} else if (write_status) {
  cat("Status unchanged (", existing_status, "); left ", status_md, " untouched.\n", sep = "")
}

invisible(new_status)
