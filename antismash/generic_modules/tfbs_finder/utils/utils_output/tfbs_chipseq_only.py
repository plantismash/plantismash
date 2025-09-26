#!/usr/bin/env python3
# tfbs_chipseq_only.py


"""
Filter TFBS hits table to only those validated by ChIP-seq (optionally only Strong confidence).


# default: CSV output with all confidences
python tfbs_chipseq_only.py tfbs_hits_all_bgcs.csv

# only Strong validated hits
python tfbs_chipseq_only.py tfbs_hits_all_bgcs.csv --strong-only

# write TSV instead
python tfbs_chipseq_only.py tfbs_hits_all_bgcs.csv --tsv --out tfbs_hits_chipseq_only.tsv

# write csv 
"""

# ran like

import argparse, re
import pandas as pd

def _norm_val(x: str) -> str:
    return "yes" if isinstance(x, str) and x.strip().lower() in {"yes","y","true","validated"} else "no"

def _extract_method(evidence: str) -> str:
    """Pull 'Method: ...' value from Evidence, else best-effort guess."""
    if not isinstance(evidence, str):
        return ""
    m = re.search(r"Method:\s*([^;]+)", evidence, flags=re.I)
    if m:
        return m.group(1).strip()
    # fallback: detect common tokens
    s = evidence.lower()
    if "chip-seq" in s or "chip seq" in s or "chipseq" in s: return "ChIP-seq"
    if "dap-seq" in s or "dap seq" in s or "dapseq" in s:   return "DAP-seq"
    return ""

def main():
    ap = argparse.ArgumentParser(description="Filter table to ChIP-seq validated hits.")
    ap.add_argument("csv", help="Input CSV (e.g., tfbs_hits_all_bgcs.csv)")
    ap.add_argument("--out", default="tfbs_hits_chipseq_only.csv",
                    help="Output filename (default: tfbs_hits_chipseq_only.csv)")
    ap.add_argument("--strong-only", action="store_true",
                    help="Keep only Confidence=Strong (default: keep all confidences)")
    ap.add_argument("--tsv", action="store_true",
                    help="Write TSV instead of CSV")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)

    # normalize needed fields
    if "Validated" not in df.columns or "Evidence" not in df.columns:
        raise SystemExit("Input is missing required columns: Validated, Evidence")
    df["Validated"] = df["Validated"].map(_norm_val)
    df["Method"]    = df["Evidence"].map(_extract_method)

    # filter: validated + ChIP-seq
    chip = df[(df["Validated"] == "yes") & (df["Method"].str.contains(r"chip[-\s]?seq", case=False, na=False))]

    # optional filter: strong only
    if args.strong_only:
        chip = chip[chip["Confidence"].astype(str).str.lower().str.startswith("stron")]

    # write
    if args.tsv:
        chip.to_csv(args.out, sep="\t", index=False)
    else:
        chip.to_csv(args.out, index=False)

    # summary
    total = len(df)
    val_yes = (df["Validated"] == "yes").sum()
    print(f"Total rows: {total}")
    print(f"Validated rows: {val_yes}")
    print(f"ChIP-seq validated rows{' (Strong only)' if args.strong_only else ''}: {len(chip)}")
    print(f"Wrote: {args.out}")

if __name__ == "__main__":
    main()