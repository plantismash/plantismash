#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Usage
# python make_tfbs_hist.py --input path/to/tfbs_hits_all_bgcs.csv --outdir results/result_tfbs/tfbs26


import argparse, os
import pandas as pd
import matplotlib.pyplot as plt

CONF_ORDER = ["Strong", "Medium", "Weak"]

def coerce_conf(x: str) -> str:
    if not isinstance(x, str):
        return "Other"
    x = x.strip().title()
    return x if x in CONF_ORDER else "Other"

def coerce_yesno(x) -> bool:
    if isinstance(x, str):
        return x.strip().lower() in {"yes", "y", "true", "1"}
    if isinstance(x, (int, float)):
        return bool(x)
    return False

def main():
    ap = argparse.ArgumentParser(description="Make TFBS histograms from consolidated CSV.")
    ap.add_argument("--input", required=True, help="Path to tfbs_hits_all_bgcs.csv (or similar).")
    ap.add_argument("--outdir", default=None, help="Where to save plots (default: alongside input).")
    ap.add_argument("--dpi", type=int, default=180, help="Figure DPI (default: 180).")
    args = ap.parse_args()

    outdir = args.outdir or os.path.dirname(os.path.abspath(args.input)) or "."
    os.makedirs(outdir, exist_ok=True)

    # Load
    df = pd.read_csv(args.input)

    # Normalize columns we rely on
    for col in ["Record", "Cluster", "Product", "Confidence", "Validated"]:
        if col not in df.columns:
            raise SystemExit(f"Missing required column: {col}")

    df["Confidence"] = df["Confidence"].map(coerce_conf)
    df["Validated_bool"] = df["Validated"].map(coerce_yesno)

    # Cluster label (Record:cN Product)
    df["BGC"] = df["Record"].astype(str) + ":c" + df["Cluster"].astype(str)
    df["BGC_label"] = df["BGC"] + " (" + df["Product"].astype(str) + ")"

    # ========= Overall: confidence distribution =========
    overall_conf = (
        df.groupby("Confidence")
          .size()
          .reindex(CONF_ORDER + ["Other"])
          .fillna(0)
          .astype(int)
    )
    ax = overall_conf.plot(kind="bar")
    ax.set_title("TFBS hits by confidence (overall)")
    ax.set_xlabel("Confidence")
    ax.set_ylabel("Count")
    for p in ax.patches:
        ax.annotate(int(p.get_height()), (p.get_x() + p.get_width()/2, p.get_height()),
                    ha="center", va="bottom", fontsize=8, rotation=0, xytext=(0,2), textcoords="offset points")
    plt.tight_layout()
    f1 = os.path.join(outdir, "tfbs_hits_by_confidence_overall.png")
    plt.savefig(f1, dpi=args.dpi)
    plt.close()

    # ========= Overall: validated vs not =========
    overall_val = df["Validated_bool"].value_counts().reindex([True, False]).fillna(0).astype(int)
    overall_val.index = ["Validated", "Not validated"]
    ax = overall_val.plot(kind="bar")
    ax.set_title("Validated TFBS hits (overall)")
    ax.set_xlabel("")
    ax.set_ylabel("Count")
    for p in ax.patches:
        ax.annotate(int(p.get_height()), (p.get_x() + p.get_width()/2, p.get_height()),
                    ha="center", va="bottom", fontsize=8, xytext=(0,2), textcoords="offset points")
    plt.tight_layout()
    f2 = os.path.join(outdir, "tfbs_validated_overall.png")
    plt.savefig(f2, dpi=args.dpi)
    plt.close()

    # ========= Per-BGC: confidence distribution (stacked) =========
    pivot_conf = (
        df.pivot_table(index="BGC_label", columns="Confidence", values="Record",
                       aggfunc="count", fill_value=0)
          .reindex(columns=[c for c in CONF_ORDER if c in df["Confidence"].unique()] + (["Other"] if "Other" in df["Confidence"].unique() else []))
          .sort_index()
    )
    ax = pivot_conf.plot(kind="bar", stacked=True, figsize=(max(6, len(pivot_conf) * 0.5), 4))
    ax.set_title("TFBS hits by confidence per BGC")
    ax.set_xlabel("BGC")
    ax.set_ylabel("Count")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    f3 = os.path.join(outdir, "tfbs_hits_by_confidence_per_bgc.png")
    plt.savefig(f3, dpi=args.dpi)
    plt.close()

    # ========= Per-BGC: validated counts =========
    per_bgc_valid = df.groupby("BGC_label")["Validated_bool"].sum().astype(int).sort_index()
    ax = per_bgc_valid.plot(kind="bar", figsize=(max(6, len(per_bgc_valid) * 0.5), 4))
    ax.set_title("Validated TFBS hits per BGC")
    ax.set_xlabel("BGC")
    ax.set_ylabel("Validated count")
    plt.xticks(rotation=45, ha="right")
    for p in ax.patches:
        ax.annotate(int(p.get_height()), (p.get_x() + p.get_width()/2, p.get_height()),
                    ha="center", va="bottom", fontsize=8, xytext=(0,2), textcoords="offset points")
    plt.tight_layout()
    f4 = os.path.join(outdir, "tfbs_validated_per_bgc.png")
    plt.savefig(f4, dpi=args.dpi)
    plt.close()

    # ========= Export summaries =========
    overall_conf.to_csv(os.path.join(outdir, "summary_overall_confidence_counts.csv"), header=["count"])
    per_bgc_valid.to_csv(os.path.join(outdir, "summary_validated_per_bgc.csv"), header=["validated_count"])

    print("Saved:")
    for f in (f1, f2, f3, f4):
        print(" -", f)

if __name__ == "__main__":
    main()
