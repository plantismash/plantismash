#!/usr/bin/env python3

"""
Author: Hannah Augustijn 
University: Wageningen University and Research & Leiden University
Department: Department of Bioinformatics & Institute of Biology Leiden
Date: 02/05/2025

Generate sequence logos from filtered motif JSON.

- Works with BOTH formats:
  * New format: PWM is 4×N (rows = A,C,G,T)
  * Old format: PWM is N×4 (rows = positions, cols = A,C,G,T)

Usage:
  python json_pwm_to_logo.py --json Athaliana_motifs.filtered.json --outdir logos_filtered
"""

import argparse
import json
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt
import logomaker


BASES = ["A", "C", "G", "T"]


def parse_args():
    p = argparse.ArgumentParser(description="Make sequence logos from motif JSON (4×N or N×4).")
    p.add_argument("--json", required=True, help="Input motifs JSON (filtered).")
    p.add_argument("--outdir", default="logos_filtered", help="Output directory for PNG logos.")
    p.add_argument("--dpi", type=int, default=150, help="PNG DPI.")
    p.add_argument("--fmt", default="png", choices=["png", "svg", "pdf"], help="Figure format.")
    return p.parse_args()


def pwm_positions_df(pwm):
    """
    Return a DataFrame with columns ['A','C','G','T'] and one row per position.
    Accepts either:
      - 4×N (rows = A,C,G,T)  -> transpose to N×4
      - N×4 (rows = positions)
    """
    if not pwm or not isinstance(pwm, list) or not isinstance(pwm[0], list):
        raise ValueError("Malformed PWM (expected list of lists).")

    n_rows = len(pwm)
    n_cols = len(pwm[0])

    # Heuristic:
    # - If there are exactly 4 rows AND not exactly 4 columns, assume 4×N (new format)
    # - Else if there are 4 columns, assume N×4 (old format)
    # - If 4×4, prefer 'new format' by default (rows=bases), transpose to N×4
    if n_rows == 4 and (n_cols != 4 or True):
        # 4×N -> transpose to N×4
        pwm_cols = list(zip(*pwm))  # N×4
        df = pd.DataFrame(pwm_cols, columns=BASES)
    elif n_cols == 4:
        # N×4 already
        df = pd.DataFrame(pwm, columns=BASES)
    else:
        # Fallback: try to coerce; raise if inconsistent
        raise ValueError(f"Cannot infer PWM orientation (shape {n_rows}×{n_cols}).")

    # Ensure numeric and non-negative
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    df[df < 0] = 0.0
    return df


def ppm_with_information(df):
    """
    Convert to position probability matrix (rows sum to 1), then weight each row by
    information content: IC = 2 - H, where H = -Σ p log2 p.
    """
    # Row-normalize to probabilities
    row_sums = df.sum(axis=1).replace(0, np.nan)
    ppm = df.div(row_sums, axis=0).fillna(0.0)

    # Information content per position
    with np.errstate(divide='ignore', invalid='ignore'):
        logp = np.where(ppm.values > 0, np.log2(ppm.values), 0.0)
    entropy = -np.sum(ppm.values * logp, axis=1)  # H
    ic = 2.0 - entropy  # log2(4)=2 for DNA

    # Weight each row by IC
    ic_matrix = ppm.mul(ic, axis=0)
    return ic_matrix


def create_sequence_logo(pwm, tf_name, outdir, dpi=150, fmt="png"):
    df = pwm_positions_df(pwm)              # N×4, cols A,C,G,T
    ic_matrix = ppm_with_information(df)    # N×4, weighted by IC

    # Plot width scales with motif length
    motif_len = ic_matrix.shape[0]
    plt.figure(figsize=(max(3, motif_len / 2), 3))
    logo = logomaker.Logo(ic_matrix, shade_below=.5, fade_below=.5, font_name='Arial')
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.ax.set_ylabel('Information (bits)')
    logo.ax.set_title(tf_name)

    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"{tf_name}_logo.{fmt}")
    plt.savefig(outpath, bbox_inches='tight', dpi=dpi)
    plt.close()
    return outpath


def main():
    args = parse_args()

    # Load motifs
    with open(args.json, "r") as f:
        motifs = json.load(f)

    made = 0
    errors = 0
    for name, motif in motifs.items():
        pwm = motif.get("pwm")
        if pwm is None:
            errors += 1
            print(f"[WARN] {name}: missing PWM, skipping", file=sys.stderr)
            continue
        try:
            outpath = create_sequence_logo(pwm, name, args.outdir, dpi=args.dpi, fmt=args.fmt)
            made += 1
        except Exception as e:
            errors += 1
            print(f"[WARN] {name}: failed to render logo: {e}", file=sys.stderr)

    print(f"✅ Generated {made} logo(s) in '{args.outdir}/' ({errors} skipped).")


if __name__ == "__main__":
    main()
