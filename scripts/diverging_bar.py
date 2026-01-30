#!/usr/bin/env python3
"""dynamic_diverging_bar.py

Create a stacked *diverging* horizontal bar chart comparing two samples across
length bins, stratified by starting nucleotide (A/C/G/T).

This script was refactored from the original notebook:
`dynamic_diverging_bar_V1.ipynb`.

Expected input format (TSV/CSV):
- a row per observation with at least:
    * Length
    * Starting_Nucleotide
    * Average Reads (TPM)   (or another numeric column you specify)

Typical use (example):
    python3 dynamic_diverging_bar.py \
        --left tRNA@ce_maleAVG.tsv --right tRNA@ce_spermAVG.tsv \
        --left-label "Male" --right-label "Sperm" \
        --title "TPM: Male(-) vs Sperm (+) Sense tRF" \
        --length-min 18 --length-max 40 \
        --ytick-every 1 \
        --out male_spermAVG1040.pdf

Notes:
- The *left* sample is plotted as negative values, the *right* sample as positive.
- Bars are stacked by nucleotide within each length bin.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

DEFAULT_NUC_ORDER: Tuple[str, ...] = ("A", "C", "G", "T")

def read_table(path: Path) -> pd.DataFrame:
    """Read TSV/CSV based on file extension."""
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_table(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    # fall back to pandas' delimiter inference for unknown text formats
    return pd.read_table(path, sep=None, engine="python")

def build_nucleotide_matrix(
    df: pd.DataFrame,
    *,
    nuc_col: str,
    length_col: str,
    value_col: str,
    nuc_order: Sequence[str] = DEFAULT_NUC_ORDER,
    length_min: int | None = None,
    length_max: int | None = None,
) -> Tuple[np.ndarray, List[int]]:
    """Return matrix shaped (len(nuc_order), n_lengths) and the length bins.

    The matrix is ordered by nuc_order, and within each nucleotide by ascending
    Length. Missing (nucleotide, length) pairs are filled with 0.
    """
    missing = [c for c in (nuc_col, length_col, value_col) if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns: {missing}. Available columns: {list(df.columns)}"
        )

    # Determine length domain
    lengths = sorted(df[length_col].dropna().astype(int).unique().tolist())
    if not lengths:
        raise ValueError("No length values found.")

    if length_min is None:
        length_min = int(min(lengths))
    if length_max is None:
        length_max = int(max(lengths))

    length_bins = list(range(int(length_min), int(length_max) + 1))

    # Build a pivot table: rows=Length, cols=Nucleotide, values=sum(value_col)
    # sum() is robust if there are multiple rows per (Length, Nuc)
    pivot = (
        df.assign(**{length_col: df[length_col].astype(int)})
        .pivot_table(
            index=length_col,
            columns=nuc_col,
            values=value_col,
            aggfunc="sum",
            fill_value=0.0,
        )
        .reindex(index=length_bins, fill_value=0.0)
    )

    mat = np.vstack([pivot.get(n, pd.Series(0.0, index=length_bins)).to_numpy() for n in nuc_order])
    return mat, length_bins

def plot_diverging_stacked_bars(
    left_mat: np.ndarray,
    right_mat: np.ndarray,
    *,
    length_bins: Sequence[int],
    nuc_order: Sequence[str] = DEFAULT_NUC_ORDER,
    left_label: str = "Left",
    right_label: str = "Right",
    title: str = "",
    outpath: Path | None = None,
    show: bool = False,
    ytick_every: int = 5,
) -> None:
    """Create a diverging stacked horizontal bar plot and save it."""
    if left_mat.shape != right_mat.shape:
        raise ValueError(f"Matrix shapes differ: {left_mat.shape} vs {right_mat.shape}")
    if left_mat.shape[0] != len(nuc_order):
        raise ValueError(
            f"Expected {len(nuc_order)} nucleotide rows, got {left_mat.shape[0]}"
        )
    if left_mat.shape[1] != len(length_bins):
        raise ValueError(
            f"Expected {len(length_bins)} length bins, got {left_mat.shape[1]}"
        )

    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    index = np.arange(len(length_bins))
    bar_height = 1.0

    # Choose colors (stable + readable). We keep it simple & deterministic.
    cmap = plt.cm.viridis
    colors = cmap(np.linspace(0.15, 0.85, len(nuc_order)))

    # Left side: negative
    left_bottom = np.zeros(len(index))
    for i, nuc in enumerate(nuc_order):
        vals = -left_mat[i]
        ax.barh(
            index,
            vals,
            left=left_bottom,
            color=colors[i],
            height=bar_height,
            label=f"{left_label} {nuc}",
        )
        left_bottom += vals

    # Right side: positive
    right_bottom = np.zeros(len(index))
    for i, nuc in enumerate(nuc_order):
        vals = right_mat[i]
        ax.barh(
            index,
            vals,
            left=right_bottom,
            color=colors[i],
            height=bar_height,
            label=f"{right_label} {nuc}",
        )
        right_bottom += vals

    # Y axis: label as lengths
    if ytick_every < 1:
        ytick_every = 1
    major_idx = index[::ytick_every]
    ax.set_yticks(major_idx)
    ax.set_yticklabels([str(length_bins[i]) for i in major_idx])

    # Minor ticks each 1 bin, without labels (helps reading)
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis="y", which="minor", length=2)
    ax.tick_params(axis="y", which="major", length=4)

    ax.axvline(0, linewidth=1)

    ax.set_xlabel("Normalized abundance", fontweight="bold")
    if title:
        ax.set_title(title, fontweight="bold", fontsize=14)

    # Legend: keep but you can disable by setting --no-legend
    ax.legend(loc="best", fontsize=8, frameon=False)

    fig.tight_layout()

    if outpath is not None:
        outpath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outpath, format=outpath.suffix.lstrip(".") or "pdf")

    if show:
        plt.show()
    else:
        plt.close(fig)

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
    description=(
        "Generate a stacked diverging horizontal bar plot comparing two samples across length bins, "
        "with bars stratified by starting nucleotide."
    )
)
    p.add_argument("--left", type=Path, required=True, help="Left sample table (plotted negative).") 
    p.add_argument("--right", type=Path, required=True, help="Right sample table (plotted positive).") 
    p.add_argument("--out", type=Path, required=True, help="Output figure path (e.g., .pdf, .png).")

    p.add_argument("--left-label", default="Left", help="Label prefix for left sample legend.")
    p.add_argument("--right-label", default="Right", help="Label prefix for right sample legend.")
    p.add_argument("--title", default="", help="Plot title.")

    p.add_argument("--nuc-col", default="Starting_Nucleotide", help="Column name for nucleotide.")
    p.add_argument("--length-col", default="Length", help="Column name for length bin.")
    p.add_argument("--value-col", default="Average Reads (TPM)", help="Column name for numeric values.")

    p.add_argument("--length-min", type=int, default=None, help="Minimum length bin (inclusive).") 
    p.add_argument("--length-max", type=int, default=None, help="Maximum length bin (inclusive).") 
    p.add_argument(
        "--nuc-order",
        default=",".join(DEFAULT_NUC_ORDER),
        help="Comma-separated nucleotide order (default: A,C,G,T). ",
    )
    p.add_argument("--ytick-every", type=int, default=5, help="Show every Nth y tick label.") 
    p.add_argument("--show", action="store_true", help="Show the plot window (instead of just saving).") 
    p.add_argument(
        "--no-legend",
        action="store_true",
        help="Disable legend (useful for manuscript assembly). ",
    )
    return p.parse_args()

def main() -> None:
    args = parse_args()

    nuc_order = tuple([x.strip() for x in args.nuc_order.split(",") if x.strip()])
    if not nuc_order:
        raise SystemExit("--nuc-order is empty")

    left_df = read_table(args.left)
    right_df = read_table(args.right)

    left_mat, length_bins = build_nucleotide_matrix(
        left_df,
        nuc_col=args.nuc_col,
        length_col=args.length_col,
        value_col=args.value_col,
        nuc_order=nuc_order,
        length_min=args.length_min,
        length_max=args.length_max,
    )
    right_mat, length_bins2 = build_nucleotide_matrix(
        right_df,
        nuc_col=args.nuc_col,
        length_col=args.length_col,
        value_col=args.value_col,
        nuc_order=nuc_order,
        length_min=args.length_min,
        length_max=args.length_max,
    )

    # Ensure consistent length bins between the two inputs
    if list(length_bins2) != list(length_bins):
        raise SystemExit(
            "Left/right length bins differ. Use --length-min/--length-max to force a common range."
        )

    # Make plot
    if args.no_legend:
        # temporarily disable legend by monkeypatching after plotting
        # simplest: create plot, then remove legend object
        # We'll implement by plotting then removing legend in-place.
        import matplotlib.pyplot as _plt  # local import to keep global tidy

        # Use a custom plotting call to remove legend
        fig, ax = plt.subplots(figsize=(6.5, 6.5))
        index = np.arange(len(length_bins))
        bar_height = 1.0
        cmap = plt.cm.viridis
        colors = cmap(np.linspace(0.15, 0.85, len(nuc_order)))

        left_bottom = np.zeros(len(index))
        for i, nuc in enumerate(nuc_order):
            vals = -left_mat[i]
            ax.barh(index, vals, left=left_bottom, color=colors[i], height=bar_height)
            left_bottom += vals

        right_bottom = np.zeros(len(index))
        for i, nuc in enumerate(nuc_order):
            vals = right_mat[i]
            ax.barh(index, vals, left=right_bottom, color=colors[i], height=bar_height)
            right_bottom += vals

        major_idx = index[::max(1, args.ytick_every)]
        ax.set_yticks(major_idx)
        ax.set_yticklabels([str(length_bins[i]) for i in major_idx])
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
        ax.tick_params(axis="y", which="minor", length=2)
        ax.tick_params(axis="y", which="major", length=4)
        ax.axvline(0, linewidth=1)
        ax.set_xlabel("Normalized abundance", fontweight="bold")
        if args.title:
            ax.set_title(args.title, fontweight="bold", fontsize=14)
        fig.tight_layout()
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out, format=args.out.suffix.lstrip(".") or "pdf")
        if args.show:
            _plt.show()
        else:
            _plt.close(fig)
    else:
        plot_diverging_stacked_bars(
            left_mat,
            right_mat,
            length_bins=length_bins,
            nuc_order=nuc_order,
            left_label=args.left_label,
            right_label=args.right_label,
            title=args.title,
            outpath=args.out,
            show=args.show,
            ytick_every=args.ytick_every,
        )

if __name__ == "__main__":
    main()
