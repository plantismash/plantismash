#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Single stacked histogram for all BGCs:
Each BGC is one bar, stacked by confidence × validation:
[Strong (val), Strong (no), Medium (val), Medium (no), Weak (val), Weak (no)]

Input CSV must have (at least):
Record,Cluster,ClusterStart,ClusterEnd,Product,Confidence,Validated

Usage:
  python tfbs_bgc_stacked_hits.py tfbs_hits_all_bgcs.csv --out all_bgcs.png [--pdf all_bgcs.pdf]
"""

import argparse
import os
from typing import List, Tuple

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

CONF_ORDER = ["Strong", "Medium", "Weak"]
VAL_ORDER  = ["yes", "no"]  # plot validated first (bottom), then not validated
STACK_ORDER: List[Tuple[str, str]] = [
    ("Strong", "yes"),
    ("Strong", "no"),
    ("Medium", "yes"),
    ("Medium", "no"),
    ("Weak", "yes"),
    ("Weak", "no"),
]

def _norm_conf(x: str) -> str:
    if not isinstance(x, str):
        return "Weak"
    s = x.strip().lower()
    if s.startswith("stron"): return "Strong"
    if s.startswith("mediu"): return "Medium"
    if s.startswith("weak"):  return "Weak"
    return "Weak"

def _norm_val(x: str) -> str:
    if isinstance(x, str) and x.strip().lower() in {"yes","y","true","validated"}:
        return "yes"
    return "no"

def load_hits(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    required = {"Record","Cluster","ClusterStart","ClusterEnd","Product","Confidence","Validated"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {', '.join(sorted(missing))}")

    df["Confidence"] = df["Confidence"].map(_norm_conf)
    df["Validated"]  = df["Validated"].map(_norm_val)

    # Clean/grouping keys
    df["Cluster"]       = pd.to_numeric(df["Cluster"], errors="coerce").astype("Int64")
    df["ClusterStart"]  = pd.to_numeric(df["ClusterStart"], errors="coerce")
    df["ClusterEnd"]    = pd.to_numeric(df["ClusterEnd"], errors="coerce")
    df["Product"]       = df["Product"].fillna("-").astype(str)

    # Drop rows lacking key identifiers
    df = df.dropna(subset=["Record","Cluster"])
    return df

def build_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Return a table indexed by BGC (Record,Cluster,Start,End,Product),
       columns=STACK_ORDER, values=counts, missing filled with 0.
    """
    gb = (df
          .groupby(["Record","Cluster","ClusterStart","ClusterEnd","Product","Confidence","Validated"])
          .size()
          .rename("count")
          .reset_index())

    # Pivot to columns for (Confidence, Validated)
    pivot = gb.pivot_table(index=["Record","Cluster","ClusterStart","ClusterEnd","Product"],
                           columns=["Confidence","Validated"],
                           values="count",
                           fill_value=0,
                           aggfunc="sum")

    # Ensure all stacks exist in the requested order
    # Build a flat column index of tuples matching STACK_ORDER
    for c in CONF_ORDER:
        for v in VAL_ORDER:
            if (c,v) not in pivot.columns:
                pivot[(c,v)] = 0

    # Reorder columns to STACK_ORDER
    pivot = pivot.reindex(columns=STACK_ORDER)
    # Sort BGCs by Record, Cluster (then start)
    pivot = pivot.sort_values(by=[("Strong","yes")], ascending=False)  # quick stable sort by strong val
    pivot = pivot.sort_index(level=[0,1,2])  # then overall index order
    return pivot

def make_labels(idx: pd.MultiIndex) -> List[str]:
    labels = []
    for rec, clu, cstart, cend, prod in idx:
        span = f"{int(cstart)}-{int(cend)}"
        prod_txt = prod if prod and prod != "nan" else "-"
        labels.append(f"{rec}\n#{int(clu)} [{span}]")
        # If you prefer shorter: labels.append(f"{rec}-#{int(clu)}")
    return labels

def main():
    ap = argparse.ArgumentParser(description="Single stacked histogram of TFBS hits per BGC (confidence × validation).")
    ap.add_argument("csv", help="Input CSV (e.g., tfbs_hits_all_bgcs.csv)")
    ap.add_argument("--out", default="all_bgcs.png", help="Output image file (PNG/SVG/etc.)")
    ap.add_argument("--pdf", default=None, help="Optional PDF to save the same figure")
    ap.add_argument("--min-hits", type=int, default=1, help="Hide BGCs with < N total hits (default: 1)")
    ap.add_argument("--max-bgcs", type=int, default=0, help="Limit to first N BGCs after sorting (0 = no limit)")
    ap.add_argument("--rotate", type=float, default=90, help="X tick label rotation (default: 90)")
    args = ap.parse_args()

    df = load_hits(args.csv)
    table = build_counts(df)

    # Filter by minimum hits
    totals = table.sum(axis=1)
    table = table.loc[totals >= args.min_hits]
    if args.max_bgcs and args.max_bgcs > 0:
        table = table.iloc[:args.max_bgcs]

    if table.empty:
        raise SystemExit("No BGCs remain after filtering. Nothing to plot.")

    n_bgcs = table.shape[0]
    # Figure width proportional to number of BGCs (readable ticks)
    fig_w = max(10.0, min(0.35 * n_bgcs, 60.0))
    fig_h = 6.5

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)

    bottoms = None
    labels = make_labels(table.index)
    x = range(n_bgcs)

    # Plot in STACK_ORDER (legend entries in the same order)
    legend_labels = [
        "Strong (validated)",
        "Strong (not validated)",
        "Medium (validated)",
        "Medium (not validated)",
        "Weak (validated)",
        "Weak (not validated)",
    ]
    for (conf, val), leg_label in zip(STACK_ORDER, legend_labels):
        vals = table[(conf, val)].values
        ax.bar(x, vals, bottom=bottoms, label=leg_label)
        bottoms = vals if bottoms is None else (bottoms + vals)

    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=args.rotate, ha="right", fontsize=8)

    ax.set_ylabel("TFBS hit count")
    ax.set_title("TFBS hits per BGC (stacked by confidence × validation)")

    ax.legend(ncol=3, frameon=False, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.6)

    # Save
    out_path = os.path.abspath(args.out)
    fig.savefig(out_path, dpi=200)
    if args.pdf:
        with PdfPages(args.pdf) as pdf:
            pdf.savefig(fig)
    plt.close(fig)

    print(f"Wrote: {out_path}")
    if args.pdf:
        print(f"Wrote PDF: {os.path.abspath(args.pdf)}")

if __name__ == "__main__":
    main()