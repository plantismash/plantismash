#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Per-BGC stacked histogram of TFBS hits, showing ONLY Strong + Validated hits,
split by validation Method parsed from the 'Evidence' column.

Input CSV needs at least:
Record,Cluster,ClusterStart,ClusterEnd,Product,Confidence,Validated,Evidence

Usage:
  python tfbs_bgc_validated_methods.py tfbs_hits_all_bgcs.csv \
      --out validated_strong_methods.png \
      [--pdf validated_strong_methods.pdf] \
      [--top-methods 8] \
      [--min-hits 1] \
      [--max-bgcs 0] \
      [--rotate 90]
"""

import argparse
import os
import re
from typing import List

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ------------------------- helpers -------------------------

def _norm_conf(x: str) -> str:
    if not isinstance(x, str): return "Weak"
    s = x.strip().lower()
    if s.startswith("stron"): return "Strong"
    if s.startswith("mediu"): return "Medium"
    if s.startswith("weak"):  return "Weak"
    return "Weak"

def _norm_val(x: str) -> str:
    if isinstance(x, str) and x.strip().lower() in {"yes", "y", "true", "validated"}:
        return "yes"
    return "no"

_METHOD_MAP = {
    "chip-seq": "ChIP-seq", "chip seq": "ChIP-seq", "chipseq": "ChIP-seq",
    "dap-seq": "DAP-seq", "dap seq": "DAP-seq", "dapseq": "DAP-seq", "ampdap": "DAP-seq",
    "y1h": "Y1H", "yeast one hybrid": "Y1H",
    "emsa": "EMSA",
    "pbm": "PBM",
    "grn": "GRN", "gr": "GRN",
}

def _parse_method(evidence: str) -> str:
    """
    Extract 'Method: <...>' from free-text Evidence. Return tidy method or 'Unspecified'.
    """
    if not isinstance(evidence, str) or not evidence.strip():
        return "Unspecified"
    m = re.search(r"(?i)method\s*:\s*([^;]+)", evidence)
    if not m:
        return "Unspecified"
    raw = m.group(1).strip()
    # take first token if multiple listed
    first = re.split(r"[,/|&]+", raw)[0].strip()
    key = first.lower().replace("–", "-").replace("—", "-")
    key = re.sub(r"\s+", " ", key)
    key = key.replace("dap seq", "dap-seq").replace("chip seq", "chip-seq")
    # map common variants
    for k, v in _METHOD_MAP.items():
        if key == k:
            return v
    # fallback to cleaned original casing
    return first

def _labels_from_index(idx: pd.MultiIndex) -> List[str]:
    labels = []
    for rec, clu, cstart, cend, prod in idx:
        labels.append(f"{rec}\n#{int(clu)} [{int(cstart)}-{int(cend)}]")
    return labels


# ------------------------- core -------------------------

def load_and_filter(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    need = {"Record","Cluster","ClusterStart","ClusterEnd","Product","Confidence","Validated","Evidence"}
    miss = need - set(df.columns)
    if miss:
        raise ValueError(f"CSV missing columns: {', '.join(sorted(miss))}")

    df["Confidence"] = df["Confidence"].map(_norm_conf)
    df["Validated"]  = df["Validated"].map(_norm_val)

    # Keep ONLY Strong + Validated
    df = df[(df["Confidence"] == "Strong") & (df["Validated"] == "yes")].copy()
    if df.empty:
        raise SystemExit("No Strong + Validated hits found.")

    # tidy types
    df["Cluster"]      = pd.to_numeric(df["Cluster"], errors="coerce").astype("Int64")
    df["ClusterStart"] = pd.to_numeric(df["ClusterStart"], errors="coerce")
    df["ClusterEnd"]   = pd.to_numeric(df["ClusterEnd"], errors="coerce")
    df["Product"]      = df["Product"].fillna("-").astype(str)

    df = df.dropna(subset=["Record","Cluster"])

    # parse method
    df["Method"] = df["Evidence"].map(_parse_method)
    return df

def pivot_counts(df: pd.DataFrame, top_methods: int) -> pd.DataFrame:
    # count per BGC × Method
    gb = (df.groupby(["Record","Cluster","ClusterStart","ClusterEnd","Product","Method"])
            .size().rename("count").reset_index())

    # limit methods to top N globally, rest -> "Other"
    method_totals = gb.groupby("Method")["count"].sum().sort_values(ascending=False)
    if top_methods and top_methods > 0 and len(method_totals) > top_methods:
        keep = set(method_totals.head(top_methods).index)
        gb["Method"] = gb["Method"].where(gb["Method"].isin(keep), "Other")

        # re-aggregate after bucketing
        gb = (gb.groupby(["Record","Cluster","ClusterStart","ClusterEnd","Product","Method"])
                .agg(count=("count","sum")).reset_index())

    # pivot to columns by Method
    pivot = gb.pivot_table(index=["Record","Cluster","ClusterStart","ClusterEnd","Product"],
                           columns="Method", values="count", fill_value=0, aggfunc="sum")

    # sort BGCs by total validated-strong hits desc
    pivot = pivot.assign(_total=pivot.sum(axis=1)).sort_values("_total", ascending=False).drop(columns="_total")
    return pivot


def main():
    ap = argparse.ArgumentParser(description="Stacked histogram per BGC of Strong+Validated TFBS hits, split by Method.")
    ap.add_argument("csv", help="Input CSV (e.g., tfbs_hits_all_bgcs.csv)")
    ap.add_argument("--out", default="validated_strong_methods.png", help="Output image path")
    ap.add_argument("--pdf", default=None, help="Optional PDF path")
    ap.add_argument("--top-methods", type=int, default=8, help="Keep top N methods, bucket rest as 'Other' (default: 8)")
    ap.add_argument("--min-hits", type=int, default=1, help="Hide BGCs with < N Strong+Validated hits (default: 1)")
    ap.add_argument("--max-bgcs", type=int, default=0, help="Limit to first N BGCs after sorting (0 = no limit)")
    ap.add_argument("--rotate", type=float, default=90, help="X tick label rotation (default: 90)")
    args = ap.parse_args()

    df = load_and_filter(args.csv)
    table = pivot_counts(df, top_methods=args.top_methods)

    # filter by minimum hits
    totals = table.sum(axis=1)
    table = table.loc[totals >= args.min_hits]
    if args.max_bgcs and args.max_bgcs > 0:
        table = table.iloc[:args.max_bgcs]

    if table.empty:
        raise SystemExit("No BGCs remain after filtering. Nothing to plot.")

    methods = list(table.columns)  # stacked series order = as in columns
    n_bgcs = table.shape[0]

    # figure size
    fig_w = max(10.0, min(0.35 * n_bgcs, 60.0))
    fig_h = 6.5

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)

    bottoms = None
    x = range(n_bgcs)
    labels = _labels_from_index(table.index)

    for m in methods:
        vals = table[m].values
        ax.bar(x, vals, bottom=bottoms, label=m)
        bottoms = vals if bottoms is None else (bottoms + vals)

    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=args.rotate, ha="right", fontsize=8)

    ax.set_ylabel("Strong + Validated TFBS hits")
    ax.set_title("Validated methods per BGC (Strong hits only)")
    ax.legend(ncol=3, frameon=False, fontsize=9)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.6)

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