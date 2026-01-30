#!/usr/bin/env python3
"""
Generate a sample sheet template for the tRF 5'/3' pipeline.

This script scans BAMs (via --bam-glob) and writes a CSV with columns:
  bam,genotype,rep

It will TRY to infer genotype/rep from filenames like:
  tRNA@ce_male3_sorted.bam  -> genotype=ce_male, rep=3
  tRNA@ce_sperm2_sorted.bam -> genotype=ce_sperm, rep=2

If it can't infer, genotype/rep are left blank for you to fill in.

Usage:
  python make_samplesheet_template.py --bam-glob "*.bam" --out samples_template.csv
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any

import pandas as pd


def infer_genotype_rep(bam_name: str) -> Tuple[Optional[str], Optional[int]]:
    """
    Best-effort inference for filenames like:
      tRNA@ce_male3_sorted.bam
      tRNA@ce_sperm2_sorted.bam
      tRNA@ce_sperm3_sense_sorted.bam

    Returns (genotype, rep) or (None, None) if it can't infer.
    """
    stem = Path(bam_name).stem  # no .bam
    genotype = None
    rep = None

    if "@" not in stem:
        return None, None

    core = stem.split("@", 1)[1]

    for suffix in [
        "_sense_sorted",
        "_antisense_sorted",
        "_sorted",
        "_sense",
        "_antisense",
    ]:
        if core.endswith(suffix):
            core = core[: -len(suffix)]
            break

    # Expect trailing digits for replicate
    m = re.match(r"^(.*?)(\d+)$", core)
    if not m:
        return None, None

    genotype = m.group(1)
    rep = int(m.group(2))
    return genotype, rep


def build_template_rows(bam_paths: List[Path]) -> List[Dict[str, Any]]:
    rows = []
    for p in sorted(bam_paths, key=lambda x: x.name):
        # Skip already split files if user points glob at a mixed directory
        if p.name.startswith(("5p_", "3p_")):
            continue

        genotype, rep = infer_genotype_rep(p.name)

        rows.append(
            {
                "bam": p.name,
                "genotype": genotype if genotype is not None else "",
                "rep": rep if rep is not None else "",
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Generate a sample sheet template CSV from BAM filenames. "
            "Genotype/replicate values are inferred from filenames when possible; "
            "otherwise they are left blank for manual editing."
        )
    )
    ap.add_argument(
        "--bam-glob",
        default="*.bam",
        help=(
            'Glob pattern for input BAM files (default: "*.bam"). ' 
            "Example: \"aligned/*.bam\"" 
        ),
    )
    ap.add_argument(
        "--out",
        default="samples_template.csv",
        help="Output CSV filename (default: samples_template.csv)",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    bam_paths = [Path(p) for p in glob.glob(args.bam_glob)]

    if not bam_paths:
        raise SystemExit(f"No BAMs matched pattern: {args.bam_glob}")

    rows = build_template_rows(bam_paths)
    if not rows:
        raise SystemExit("No input BAMs found (did your glob only match 5p_/3p_ BAMs?)")

    df = pd.DataFrame(rows)
    if df['genotype'].isna().any() or (df['genotype'].astype(str).str.strip() == '').any() or df['rep'].isna().any() or (df['rep'].astype(str).str.strip() == '').any():
        print(
            "Note: Some genotype/rep values could not be inferred from filenames and were left blank. "
            "Please edit the CSV before running downstream analysis."
        )
    df.to_csv(args.out, index=False)
    print(f"[write] {args.out}")
    print("Next: open the CSV, fill/verify genotype and rep, then run the main pipeline with --samples-csv.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

