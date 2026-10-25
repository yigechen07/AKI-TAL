#!/usr/bin/env python3
"""Recreate Figure 1 and Figure 2 from a local original analysis directory."""
from pathlib import Path
import subprocess
import sys
import argparse

root = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--project', type=Path, required=True,
                    help='Local original analysis directory containing the required results/ CSV files')
args = parser.parse_args()
project = args.project.resolve()
if not (project / 'results/formal_effects.csv').is_file():
    parser.error('The local project must contain results/formal_effects.csv; plotting tables are not distributed in this repository.')
subprocess.run([
    sys.executable, str(root / 'scripts/redraw_figures_1_2.py'),
    '--project', str(project), '--out', str(root / 'rebuilt/figures'),
    '--biorender-assets', str(root / 'assets/biorender'),
], check=True)
