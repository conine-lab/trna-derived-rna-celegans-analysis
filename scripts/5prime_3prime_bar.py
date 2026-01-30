#!/usr/bin/env python3
"""
tRF 5' / 3' fork + summarize + plot (sample-sheet driven)

Pipeline:
  1) Split BAMs into 5p / 3p by reference_start < cutoff
  2) Count alignments in split BAMs
  3) Assign genotype/rep using a sample sheet (CSV)
  4) Write summary CSV
  5) Plot stacked % (5' bottom, 3' top)

Required sample sheet columns:
  bam, genotype, rep
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path
from typing import Dict, Any, List, Tuple

import pysam
import pandas as pd
import matplotlib.pyplot as plt


# ------------------------------------------------------------
# BAM splitting
# ------------------------------------------------------------

def split_reads_by_position(
    bam_path: Path,
    cutoff: int,
    outdir: Path,
    overwrite: bool
) -> Tuple[Path, Path]:

    outdir.mkdir(parents=True, exist_ok=True)

    five_p = outdir / f"5p_{bam_path.name}"
    three_p = outdir / f"3p_{bam_path.name}"

    for p in (five_p, three_p):
        if p.exists() and not overwrite:
            raise FileExistsError(f"{p} exists (use --overwrite)")

    with pysam.AlignmentFile(str(bam_path), "rb") as bam_in, \
         pysam.AlignmentFile(str(five_p), "wb", header=bam_in.header) as bam_5p, \
         pysam.AlignmentFile(str(three_p), "wb", header=bam_in.header) as bam_3p:

        for read in bam_in:
            if read.is_unmapped or read.reference_start is None:
                continue

            if read.reference_start < cutoff:
                bam_5p.write(read)
            else:
                bam_3p.write(read)

    for p in (five_p, three_p):
        try:
            pysam.index(str(p))
        except Exception:
            pass

    return five_p, three_p


# ------------------------------------------------------------
# Counting
# ------------------------------------------------------------

def count_alignments(bam_path: Path) -> int:
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return sum(1 for _ in bam)


# ------------------------------------------------------------
# Sample sheet loader
# ------------------------------------------------------------

def load_sample_sheet(path: Path) -> Dict[str, Dict[str, Any]]:
    df = pd.read_csv(path)

    required = {"bam", "genotype", "rep"}
    if not required.issubset(df.columns):
        raise ValueError(f"Sample sheet must contain columns: {required}")

    meta = {}
    for _, row in df.iterrows():
        meta[row["bam"]] = {
            "genotype": row["genotype"],
            "rep": int(row["rep"])
        }

    return meta


# ------------------------------------------------------------
# Summary table
# ------------------------------------------------------------

def make_summary_dataframe(
    split_bams: List[Path],
    sample_meta: Dict[str, Dict[str, Any]]
) -> pd.DataFrame:

    rows = []

    for bam in split_bams:
        stem = bam.stem
        original = bam.name.replace("5p_", "").replace("3p_", "")

        if stem.startswith("5p_"):
            half = 5
        elif stem.startswith("3p_"):
            half = 3
        else:
            raise ValueError(f"Could not determine half for {bam.name}")

        if original not in sample_meta:
            raise ValueError(f"{original} not found in sample sheet")

        rows.append({
            "Filename": stem,
            "count": count_alignments(bam),
            "genotype": sample_meta[original]["genotype"],
            "rep": sample_meta[original]["rep"],
            "half": half
        })

    return pd.DataFrame(rows).sort_values(
        by=["genotype", "rep", "half"],
        ascending=[True, True, False]
    ).reset_index(drop=True)


# ------------------------------------------------------------
# Plot (5' bottom, 3' top)
# ------------------------------------------------------------

def plot_stacked_percent(df: pd.DataFrame, out_svg: Path, title: str):
    grp = df.groupby(["genotype", "half"])["count"].sum().unstack(fill_value=0)

    if 3 not in grp.columns:
        grp[3] = 0
    if 5 not in grp.columns:
        grp[5] = 0

    total = grp[3] + grp[5]
    pct5 = (grp[5] / total * 100).fillna(0)
    pct3 = 100 - pct5

    x = range(len(grp.index))

    plt.figure(figsize=(8, 6))
    plt.bar(x, pct5, label="5' tRF")
    plt.bar(x, pct3, bottom=pct5, label="3' tRF")

    plt.xticks(x, grp.index, rotation=45, ha="right")
    plt.ylabel("Percent of alignments (%)")
    plt.title(title)
    plt.ylim(0, 100)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_svg, format="svg")
    plt.close()


# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Fork tRNA BAMs into 5′ and 3′ categories, summarize counts using a sample sheet, "
            "and optionally generate stacked percentage bar plots."
        )
    )
    ap.add_argument(
        "--bam-glob",
        default="*.bam",
        help=(
            "Glob pattern used to locate input BAM files (default: *.bam). "
            "The BAM basenames should match entries in --samples-csv."
        ),
    )
    ap.add_argument(
        "--cutoff",
        type=int,
        default=35,
        help=(
            "Nucleotide position used to split alignments into 5′ and 3′ categories. "
            "Positions ≤ cutoff are classified as 5′; positions > cutoff as 3′. "
            "Coordinates are relative to the reference tRNA. (default: 35)"
        ),
    )
    ap.add_argument(
        "--outdir",
        default="split_bams",
        help=(
            "Output directory for split BAMs and intermediate files (default: split_bams). "
            "Created if it does not exist."
        ),
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing split BAMs if present.",
    )
    ap.add_argument(
        "--samples-csv",
        required=True,
        help="Sample sheet CSV with columns: bam, genotype, rep.",
    )
    ap.add_argument(
        "--summary-csv",
        default="trf_summary.csv",
        help="Output CSV summarizing 5′ and 3′ counts per genotype and replicate (default: trf_summary.csv).",
    )
    ap.add_argument(
        "--plot-svg",
        default="5prime_3prime_comp.svg",
        help="Output SVG filename for the stacked percentage bar plot (default: 5prime_3prime_comp.svg).",
    )
    ap.add_argument(
        "--plot-title",
        default="Comparison of 5′ and 3′ tRFs by Genotype",
        help="Title to display on the generated plot.",
    )
    ap.add_argument(
        "--no-plot",
        action="store_true",
        help="Run summarization only; do not generate plots.",
    )
    return ap.parse_args()


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main() -> int:
    args = parse_args()
    if args.cutoff <= 0:
        print("[error] --cutoff must be a positive integer", file=sys.stderr)
        return 2
    outdir = Path(args.outdir)

    sample_meta = load_sample_sheet(Path(args.samples_csv))

    bams = [Path(p) for p in glob.glob(args.bam_glob)]
    bams = [b for b in bams if not b.name.startswith(("5p_", "3p_"))]

    if not bams:
        print("[error] No input BAMs found", file=sys.stderr)
        return 2

    split_bams: List[Path] = []

    for bam in bams:
        print(f"[split] {bam.name}")
        five_p, three_p = split_reads_by_position(
            bam, args.cutoff, outdir, args.overwrite
        )
        split_bams.extend([five_p, three_p])

    df = make_summary_dataframe(split_bams, sample_meta)
    df.to_csv(args.summary_csv, index=False)
    print(f"[write] {args.summary_csv}")

    if not args.no_plot:
        plot_stacked_percent(df, Path(args.plot_svg), args.plot_title)
        print(f"[plot] {args.plot_svg}")

    print("[done]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


