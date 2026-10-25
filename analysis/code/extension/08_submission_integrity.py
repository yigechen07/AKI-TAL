#!/usr/bin/env python3
"""Validate and hash the AKI-TAL submission-facing package."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
RESULTS = PROJECT / "results" / "submission"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


summary = {row["metric"]: float(row["value"]) for row in read_rows(
    RESULTS / "figure3_count_assumption_summary.csv"
)}
paired = read_rows(RESULTS / "GSE126805_prepost_only_effects.csv")

checks: list[tuple[str, bool, str]] = []


def check(name: str, passed: bool, detail: str) -> None:
    checks.append((name, bool(passed), detail))


check("frozen_1027_rounded_direction_retention",
      summary["frozen_1027_rounded_logFC_sign_retention"] == 1,
      str(summary["frozen_1027_rounded_logFC_sign_retention"]))
check("frozen_1027_limma_direction_retention",
      summary["frozen_1027_limma_logFC_sign_retention"] == 1,
      str(summary["frozen_1027_limma_logFC_sign_retention"]))
check("rounded_strict_genes_near_frozen",
      summary["rounded_edgeR_strict_cross_cohort_genes"] >= 1000,
      str(int(summary["rounded_edgeR_strict_cross_cohort_genes"])))
check("limma_strict_genes_near_frozen",
      summary["fractional_limma_strict_cross_cohort_genes"] >= 1000,
      str(int(summary["fractional_limma_strict_cross_cohort_genes"])))

paired_map = {(r["endpoint"], r["standardisation"]): r for r in paired}
six = paired_map[("frozen_six_gene", "paired_pre_post_82_samples")]
consensus = paired_map[("signed_1027_gene", "paired_pre_post_82_samples")]
check("paired_only_six_gene_negative",
      float(six["estimate"]) < 0 and float(six["ci_high"]) < 0,
      f"effect={six['estimate']}; CI={six['ci_low']} to {six['ci_high']}")
check("paired_only_consensus_positive",
      float(consensus["estimate"]) > 0 and float(consensus["ci_low"]) > 0,
      f"effect={consensus['estimate']}; CI={consensus['ci_low']} to {consensus['ci_high']}")
check("paired_only_loko_41_of_41",
      float(six["loko_same_direction_fraction"]) == 1 and
      float(consensus["loko_same_direction_fraction"]) == 1,
      "six-gene=1; consensus=1")

required = [
    PROJECT / "manuscript" / "manuscript_v4.md",
    PROJECT / "manuscript" / "evidence_map.md",
    PROJECT / "manuscript" / "supplementary_information.md",
    PROJECT / "manuscript" / "materials_manifest.md",
    PROJECT / "manuscript" / "closeout_report.md",
    PROJECT / "manuscript" / "report_for_hu.md",
    PROJECT / "submission" / "README.md",
    PROJECT / "code" / "extension" / "06_submission_sensitivity.R",
    PROJECT / "code" / "extension" / "07_submission_figure_assembly.py",
]
required += [PROJECT / "figures" / "submission" / f"figure{i}{suffix}"
             for i in (1, 2, 4, 5, 6) for suffix in (".png", ".pdf")]
required += [PROJECT / "figures" / "submission" / f"figure3_complete{suffix}"
             for suffix in (".png", ".pdf")]
required += [PROJECT / "figures" / "supplement" / f"figureS1_submission_sensitivity{suffix}"
             for suffix in (".png", ".pdf")]
required += [PROJECT / "figures" / "supplement" / f"figureS2_state_programme_bridge{suffix}"
             for suffix in (".png", ".pdf")]

for path in required:
    check(f"file:{path.relative_to(PROJECT)}", path.is_file() and path.stat().st_size > 0,
          f"bytes={path.stat().st_size if path.exists() else 0}")

manifest = RESULTS / "submission_manifest.tsv"
with manifest.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["sha256", "bytes", "path"])
    for path in sorted(required):
        writer.writerow([sha256(path), path.stat().st_size, path.relative_to(PROJECT)])

passed = all(item[1] for item in checks)
report = [
    "# Submission integrity report",
    "",
    f"**Overall:** {'PASS' if passed else 'FAIL'}",
    "",
    "| Check | Status | Detail |",
    "|---|---|---|",
]
report += [f"| {name} | {'PASS' if ok else 'FAIL'} | {detail} |"
           for name, ok, detail in checks]
report += [
    "",
    "The historical extension input/output manifest was checked separately and",
    "all frozen files passed. This report covers the new submission-facing",
    "sensitivity, manuscript and figure files.",
]
(RESULTS / "submission_integrity_report.md").write_text("\n".join(report) + "\n")
print("PASS" if passed else "FAIL")
raise SystemExit(0 if passed else 1)
