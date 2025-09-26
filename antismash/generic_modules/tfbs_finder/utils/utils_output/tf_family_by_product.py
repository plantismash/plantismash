#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Summarize TFBS hits by gene-cluster product type after adding TF family info.

Inputs:
  1) tfbs_hits_all_bgcs.csv  (columns include: Record,Cluster,Product,CDS,TF_locus,Confidence,Validated,Evidence,Method)
  2) Ath_TF_binding_motifs_information.txt  (columns: Gene_id,Family,Matrix_id,Species,Method,Datasource,Datasource_ID)

Outputs:
  - tf_family_by_product_counts.csv         (# unique TF×CDS pairs per Product×Family, or site counts if --count-sites)
  - tf_family_by_product_props_row.csv      (row-normalized proportions by Product)
  - tf_family_by_product_topN_stacked.png   (stacked bar: top N families per product; + optional PDF)
  - tf_family_by_product_heatmap.png        (heatmap of row-normalized proportions; + optional PDF)
  - tf_family_by_product_examples.csv       (helper table listing top TFs per Product×Family cell)

Defaults:
  - validated-only (toggle with --include-unvalidated)
  - unique TF×CDS pairs (toggle with --count-sites)
  - optionally restrict to strong only (--strong-only)

Usage:
  python tf_family_by_product.py tfbs_hits_all_bgcs.csv Ath_TF_binding_motifs_information.txt \
      --top-families 12 --out-prefix tf_family_by_product


      # validated hits only, unique TF×CDS pairs (recommended), top 12 families
python tf_family_by_product.py tfbs_hits_all_bgcs.csv Ath_TF_binding_motifs_information.txt \
  --top-families 12 --out-prefix tf_family_by_product

# include non-validated hits too:
python tf_family_by_product.py tfbs_hits_all_bgcs.csv Ath_TF_binding_motifs_information.txt \
  --include-unvalidated

# only strong evidence (validated strong), and also save PDFs:
python tf_family_by_product.py tfbs_hits_all_bgcs.csv Ath_TF_binding_motifs_information.txt \
  --strong-only --pdf

"""

import argparse, csv, os, re
from typing import Optional, Tuple, List
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ----------------------------
# Helpers / normalization
# ----------------------------

def _strip_ver(agi: str) -> str:
    """AT1G12345.1 -> AT1G12345"""
    if not isinstance(agi, str):
        return ""
    return agi.split(".", 1)[0].strip()

def _norm_conf(x: str) -> str:
    if not isinstance(x, str): return "Weak"
    s = x.strip().lower()
    if s.startswith("stron"): return "Strong"
    if s.startswith("mediu"): return "Medium"
    if s.startswith("weak"):  return "Weak"
    return "Weak"

def _norm_val(x: str) -> str:
    return "yes" if isinstance(x, str) and x.strip().lower() in {"yes","y","true","validated"} else "no"

def _read_table_any_delim(path: str) -> pd.DataFrame:
    """Try TSV, then CSV."""
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path)

def _safe(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_", ".", " ") else "_" for c in str(s))

# ----------------------------
# Load / merge
# ----------------------------

def load_hits(path: str,
              validated_only: bool = True,
              strong_only: bool = False) -> pd.DataFrame:
    df = pd.read_csv(path)
    needed = {"Record","Cluster","Product","CDS","TF_locus","Confidence","Validated"}
    miss = needed - set(df.columns)
    if miss:
        raise SystemExit(f"Hits file missing columns: {', '.join(sorted(miss))}")

    # Normalize
    df["TF_AGI"]      = df["TF_locus"].astype(str).map(_strip_ver)
    df["Confidence"]  = df["Confidence"].map(_norm_conf)
    df["Validated"]   = df["Validated"].map(_norm_val)
    df["Product"]     = df["Product"].fillna("-").astype(str)

    if validated_only:
        df = df[df["Validated"] == "yes"]
    if strong_only:
        df = df[df["Confidence"] == "Strong"]

    # Drop rows lacking key identifiers
    df = df.dropna(subset=["Record","Cluster","Product","CDS","TF_AGI"])
    return df

def load_tf_families(path: str) -> pd.DataFrame:
    fam = _read_table_any_delim(path)
    # Expected columns per your example:
    # Gene_id, Family, Matrix_id, Species, Method, Datasource, Datasource_ID
    if "Gene_id" not in fam.columns or "Family" not in fam.columns:
        raise SystemExit("Family file must contain columns: Gene_id and Family")
    fam["TF_AGI"] = fam["Gene_id"].astype(str).map(_strip_ver)
    fam["Family"] = fam["Family"].fillna("Unknown").astype(str)
    fam = fam[["TF_AGI","Family"]].drop_duplicates()
    return fam

def attach_family(hits: pd.DataFrame, fam: pd.DataFrame) -> pd.DataFrame:
    merged = hits.merge(fam, on="TF_AGI", how="left")
    merged["Family"] = merged["Family"].fillna("Unknown")
    return merged

# ----------------------------
# Aggregation
# ----------------------------

def build_counts(merged: pd.DataFrame, count_sites: bool) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      counts (Product x Family) of either:
        - unique TF×CDS pairs (default), or
        - raw site counts (if count_sites=True)
      and a long-form helper table with representative TFs per cell.
    """
    if count_sites:
        work = merged.copy()
    else:
        # unique TF×CDS pairs (per cluster, too)
        work = merged.drop_duplicates(subset=["Record","Cluster","CDS","TF_AGI","Family","Product"])

    # counts
    gb = (work
          .groupby(["Product","Family"], dropna=False)
          .size()
          .rename("count")
          .reset_index())

    # helper table: for each Product×Family, list top TF genes contributing
    # (by unique TF×CDS pairs to avoid overcounting multiple sites)
    helper = (merged.drop_duplicates(subset=["Record","Cluster","CDS","TF_AGI","Family","Product"])
                    .groupby(["Product","Family"])["TF_AGI"]
                    .agg(lambda x: ", ".join(pd.Series(x).value_counts().index[:5]))
                    .rename("example_TFs")
                    .reset_index())

    pivot = gb.pivot_table(index="Product", columns="Family", values="count", fill_value=0, aggfunc="sum")
    return pivot, helper

# ----------------------------
# Plots
# ----------------------------

def stacked_bars(pivot: pd.DataFrame, top_families: int, title: str,
                 out_png: str, out_pdf: Optional[str] = None):
    # choose top families by overall sum
    totals = pivot.sum(axis=0).sort_values(ascending=False)
    keep = list(totals.index[:top_families])
    other = pivot.drop(columns=keep, errors="ignore").sum(axis=1)

    plot_df = pivot[keep].copy()
    plot_df["Other"] = other

    n_prod = plot_df.shape[0]
    fig_w = max(8.0, min(0.55 * n_prod, 48.0))
    fig_h = 6.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)

    x = np.arange(n_prod)
    bottoms = np.zeros(n_prod)
    labels = [str(p) for p in plot_df.index]

    for col in plot_df.columns:
        vals = plot_df[col].values
        ax.bar(x, vals, bottom=bottoms, label=col)
        bottoms = bottoms + vals

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Count (unique TF×CDS pairs)")
    ax.set_title(title)
    ax.legend(ncol=3, frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.6)

    fig.savefig(out_png, dpi=200)
    if out_pdf:
        with PdfPages(out_pdf) as pdf:
            pdf.savefig(fig)
    plt.close(fig)

def heatmap_rowprops(pivot: pd.DataFrame, title: str,
                     out_png: str, out_pdf: Optional[str] = None):
    # row-normalize to proportions
    row_sums = pivot.sum(axis=1).replace(0, np.nan)
    props = pivot.div(row_sums, axis=0).fillna(0.0)

    fig_w = max(8.0, min(0.5 * props.shape[1], 48.0))
    fig_h = max(6.0, min(0.5 * props.shape[0], 48.0))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)

    im = ax.imshow(props.values, aspect="auto")
    ax.set_yticks(np.arange(props.shape[0]))
    ax.set_yticklabels([str(p) for p in props.index], fontsize=8)
    ax.set_xticks(np.arange(props.shape[1]))
    ax.set_xticklabels([str(c) for c in props.columns], rotation=90, ha="center", fontsize=8)

    ax.set_title(title)
    ax.set_xlabel("TF family")
    ax.set_ylabel("Cluster product type")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Proportion of family within product")

    fig.savefig(out_png, dpi=200)
    if out_pdf:
        with PdfPages(out_pdf) as pdf:
            pdf.savefig(fig)
    plt.close(fig)

# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser(description="Summarize TF families per cluster product type.")
    ap.add_argument("hits_csv", help="tfbs_hits_all_bgcs.csv")
    ap.add_argument("families_table", help="Ath_TF_binding_motifs_information.txt (Gene_id, Family, ...)")
    ap.add_argument("--out-prefix", default="tf_family_by_product", help="Prefix for outputs")
    ap.add_argument("--include-unvalidated", action="store_true", help="Include non-validated hits")
    ap.add_argument("--strong-only", action="store_true", help="Keep only Confidence=Strong")
    ap.add_argument("--count-sites", action="store_true", help="Count raw site hits instead of unique TF×CDS pairs")
    ap.add_argument("--top-families", type=int, default=12, help="Families to show explicitly in the stacked bars")
    ap.add_argument("--pdf", action="store_true", help="Also write PDF versions of plots")
    args = ap.parse_args()

    hits = load_hits(
        args.hits_csv,
        validated_only=(not args.include_unvalidated),
        strong_only=args.strong_only
    )
    fam = load_tf_families(args.families_table)
    merged = attach_family(hits, fam)

    pivot, helper = build_counts(merged, count_sites=args.count_sites)

    # Save tables
    counts_csv = f"{args.out_prefix}_counts.csv"
    props_csv  = f"{args.out_prefix}_props_row.csv"
    helper_csv = f"{args.out_prefix}_examples.csv"

    pivot.to_csv(counts_csv)
    # row-normalized proportions
    row_props = pivot.div(pivot.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    row_props.to_csv(props_csv)
    helper.to_csv(helper_csv, index=False)

    # Plots
    bars_png = f"{args.out_prefix}_top{args.top_families}_stacked.png"
    heat_png = f"{args.out_prefix}_heatmap.png"
    bars_pdf = f"{args.out_prefix}_top{args.top_families}_stacked.pdf" if args.pdf else None
    heat_pdf = f"{args.out_prefix}_heatmap.pdf" if args.pdf else None

    title_suffix = []
    title_suffix.append("validated only" if not args.include_unvalidated else "all hits")
    title_suffix.append("unique TF×CDS" if not args.count_sites else "site counts")
    if args.strong_only: title_suffix.append("Strong only")
    title = "TF families per cluster product (" + ", ".join(title_suffix) + ")"

    stacked_bars(pivot, args.top_families, title, bars_png, bars_pdf)
    heatmap_rowprops(pivot, "Row-normalized TF family proportions per product", heat_png, heat_pdf)

    # Console summary
    print(f"Wrote: {counts_csv}")
    print(f"Wrote: {props_csv}")
    print(f"Wrote: {helper_csv}")
    print(f"Wrote: {bars_png}" + (f" and {bars_pdf}" if bars_pdf else ""))
    print(f"Wrote: {heat_png}" + (f" and {heat_pdf}" if heat_pdf else ""))

if __name__ == "__main__":
    main()