#!/usr/bin/env python3
"""Write or verify SHA-256 manifests for the AKI–TAL extension package."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
DCT = PROJECT.parent / "AKI_DCT_analysis"
RESULTS = PROJECT / "results" / "extension"

INPUTS = [
    DCT / "data" / "GSE183276_Counts_decomp.RDS",
    DCT / "data" / "GSE183276_Metadata.txt.gz",
    DCT / "data" / "GSE210622_RAW.tar",
    DCT / "data" / "GSE210622_metadata.txt.gz",
    DCT / "results" / "stage2_gse210622_donor_audit.csv",
    DCT / "results" / "stage3b_pseudobulk.rds",
    DCT / "results" / "stage3b_qc.csv",
    DCT / "results" / "stage3c_model_results.csv",
    DCT / "results" / "stage3c_program_scores.csv",
    PROJECT / "results" / "formal_evidence.csv",
    PROJECT / "results" / "formal_effects.csv",
    PROJECT / "results" / "formal_genes.csv",
    PROJECT / "results" / "bulk_scores.csv",
    PROJECT / "results" / "bulk_effects.csv",
    PROJECT / "results" / "bulk_genes.csv",
    PROJECT / "results" / "tal_pt_coupling_results.csv",
    PROJECT / "metadata" / "gse126805_samples.csv",
    PROJECT / "metadata" / "extension" / "tal_state_reference.csv",
    PROJECT / "metadata" / "extension" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt",
    PROJECT / "metadata" / "extension" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt",
    PROJECT / "metadata" / "extension" / "gene_sets" / "c5.go.bp.v2026.1.Hs.symbols.gmt",
    PROJECT / "metadata" / "extension" / "gene_sets" / "provenance.md",
    PROJECT / "results" / "gse126805_results.csv",
]


def output_files() -> list[Path]:
    files: list[Path] = []
    files.extend(sorted((PROJECT / "code" / "extension").glob("*.R")))
    files.extend(sorted((PROJECT / "code" / "extension").glob("*.py")))
    for pattern in ("*.csv", "*.md", "*.rds"):
        files.extend(sorted(RESULTS.glob(pattern)))
    for pattern in ("*.png", "*.pdf"):
        files.extend(sorted((PROJECT / "figures" / "extension").glob(pattern)))
    excluded = {RESULTS / "extension_input_hashes.tsv", RESULTS / "extension_output_hashes.tsv"}
    return [p for p in files if p not in excluded]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT))
    except ValueError:
        return str(Path("..") / path.relative_to(PROJECT.parent))


def write_manifest(path: Path, files: list[Path]) -> None:
    missing = [p for p in files if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(map(str, missing)))
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sha256", "bytes", "path"])
        for item in files:
            writer.writerow([sha256(item), item.stat().st_size, relative(item)])


def check_manifest(path: Path) -> bool:
    ok = True
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            item = (PROJECT / row["path"]).resolve()
            current = sha256(item) if item.is_file() else "MISSING"
            passed = current == row["sha256"] and item.stat().st_size == int(row["bytes"])
            print(f"{'PASS' if passed else 'FAIL'}\t{row['path']}")
            ok &= passed
    return ok


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("write", "check"))
    args = parser.parse_args()
    input_manifest = RESULTS / "extension_input_hashes.tsv"
    output_manifest = RESULTS / "extension_output_hashes.tsv"
    if args.mode == "write":
        write_manifest(input_manifest, INPUTS)
        write_manifest(output_manifest, output_files())
        print(f"wrote {input_manifest}\nwrote {output_manifest}")
        return 0
    return 0 if check_manifest(input_manifest) and check_manifest(output_manifest) else 1


if __name__ == "__main__":
    raise SystemExit(main())
