#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Motif filter & converter

- Input: one or more MEME files (--in-meme ...) or a single JSON (--in-json ...)
- Output: filtered JSON with 4×N PWMs (rows A,C,G,T), metadata, and an optional TSV report.

Quality checks (per motif):
  1) Shape/range: PWM is 4×N (after transpose), 0 ≤ p ≤ 1
  2) Column sums ~ 1 (tolerance --col-sum-tol, default 0.05)
  3) Minimum length (--min-length, default 6)
  4) Over-dominance: total probability mass of any single base across all columns
     must be ≤ --max-base-share (default 0.80)
  5) (Optional) p-value threshold viability using lightmotif:
     - If lightmotif is installed, compute score threshold at --pvalue (default 1e-5)
     - By default, negative thresholds do NOT cause failure (set --fail-on-threshold to enforce)

Also:
  - Normalizes old PlantTFDB links to the new domain (https://planttfdb.gao-lab.org) and strips '#...' fragments.
  - Computes consensus (by column-wise argmax) and min/max PWM "scores" (sum of per-column min/max probabilities).

Usage examples:
  From MEME:
    python motif_filter.py --in-meme Ath_TF_binding_motifs.meme \
      --out-json Athaliana_motifs.filtered.json --report filter_report.tsv

  From JSON (already transposed):
    python motif_filter.py --in-json Athaliana_transposed.json \
      --out-json Athaliana_motifs.filtered.json --report filter_report.tsv

"""

import argparse
import json
import math
import os
import sys
from typing import Dict, List, Tuple, Optional

# Optional: lightmotif for threshold computation
try:
    import lightmotif
    _HAS_LIGHTMOTIF = True
except Exception:
    _HAS_LIGHTMOTIF = False


BASES = ["A", "C", "G", "T"]


# --- new: MOODS as preferred threshold engine ---
try:
    from MOODS import tools as moods_tools, parsers
    _HAS_MOODS = True
except Exception:
    _HAS_MOODS = False


# ----------------------------- Parsing ---------------------------------- #

def normalize_tfdb_link(link: str, fallback_id: str = "") -> str:
    """Map old PlantTFDB URLs to the new domain and strip old anchors."""
    if not link:
        return f"https://planttfdb.gao-lab.org/tf.php?sp=Ath&did={fallback_id}.1" if fallback_id else ""
    link = link.replace("http://planttfdb.cbi.pku.edu.cn", "https://planttfdb.gao-lab.org")
    link = link.split("#", 1)[0]
    return link


def transpose_pwm_cols_to_rows(pwm_cols: List[List[float]]) -> List[List[float]]:
    """Convert Nx4 (per-position rows) to 4xN (per-base rows A,C,G,T)."""
    if not pwm_cols:
        return []
    # pwm_cols is a list of [A,C,G,T] per position
    transposed = list(map(list, zip(*pwm_cols)))  # 4 x N
    return transposed


def compute_min_max_scores_from_cols(pwm_cols: List[List[float]]) -> Tuple[float, float]:
    """Min/max 'score' as sum of per-column min/max probabilities (PFM semantics)."""
    if not pwm_cols:
        return (0.0, 0.0)
    min_score = sum(min(col) for col in pwm_cols)
    max_score = sum(max(col) for col in pwm_cols)
    return (round(min_score, 4), round(max_score, 4))


def compute_consensus_from_cols(pwm_cols: List[List[float]]) -> str:
    """Consensus as the base with max probability per column (PFM semantics)."""
    if not pwm_cols:
        return ""
    cons = []
    for col in pwm_cols:
        if len(col) != 4:
            cons.append("N")
            continue
        idx = max(range(4), key=lambda i: col[i])
        cons.append(BASES[idx])
    return "".join(cons)


def parse_meme_file(path_meme: str) -> Dict[str, Dict]:
    """
    Parse a MEME file and return dict of motifs:
    {
      motif_name: {
        "name": ...,
        "description": "Unknown transcription factor",
        "species": "Arabidopsis thaliana",
        "link": "https://planttfdb.gao-lab.org/...",
        "consensus": "ACGT...",
        "max_score": float,
        "min_score": float,
        "pwm": 4xN list (A,C,G,T rows)
      }, ...
    }
    """
    motifs: Dict[str, Dict] = {}
    cur_name: Optional[str] = None
    cur_pwm_cols: List[List[float]] = []
    cur_link: str = ""

    def flush_current():
        nonlocal cur_name, cur_pwm_cols, cur_link
        if cur_name and cur_pwm_cols:
            pwm_rows = transpose_pwm_cols_to_rows(cur_pwm_cols)
            min_s, max_s = compute_min_max_scores_from_cols(cur_pwm_cols)
            consensus = compute_consensus_from_cols(cur_pwm_cols)
            motifs[cur_name] = {
                "name": cur_name,
                "description": "Unknown transcription factor",
                "species": "Arabidopsis thaliana",
                "link": normalize_tfdb_link(cur_link, cur_name),
                "consensus": consensus,
                "max_score": max_s,
                "min_score": min_s,
                "pwm": pwm_rows,  # 4xN
            }
        # reset
        cur_name = None
        cur_pwm_cols = []
        cur_link = ""

    with open(path_meme, "r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("MOTIF"):
                # starting a new motif: flush previous one
                flush_current()
                parts = line.split()
                # Typical: "MOTIF AT1G01060 MP00119" -> use the first token as name
                cur_name = parts[1] if len(parts) >= 2 else None
                cur_pwm_cols = []
                cur_link = ""
                continue

            if line.startswith("URL"):
                # e.g. "URL http://...."
                url = line.split(" ", 1)[1].strip() if " " in line else ""
                cur_link = url
                continue

            # matrix lines follow 'letter-probability matrix: ...' header,
            # but we can simply collect numeric lines with 4 floats
            # (robust to minor format changes)
            if line[0].isdigit() or line[0] == ".":
                parts = line.split()
                try:
                    vals = list(map(float, parts))
                except Exception:
                    vals = []
                if len(vals) == 4:
                    cur_pwm_cols.append(vals)
                continue

    # flush trailing motif
    flush_current()
    return motifs


def load_json_motifs(path_json: str) -> Dict[str, Dict]:
    with open(path_json, "r", encoding="utf-8") as fh:
        motifs = json.load(fh)
    # Ensure link normalization and 4xN
    fixed = {}
    for name, m in motifs.items():
        pwm = m.get("pwm", [])
        # If motif appears Nx4, transpose to 4xN
        if pwm and len(pwm[0]) == 4 and (len(pwm) != 4):
            pwm_rows = transpose_pwm_cols_to_rows(pwm)
        else:
            pwm_rows = pwm
        link = normalize_tfdb_link(m.get("link", ""), name)
        fixed[name] = {
            "name": m.get("name", name),
            "description": m.get("description", "Unknown transcription factor"),
            "species": m.get("species", "Arabidopsis thaliana"),
            "link": link,
            "consensus": m.get("consensus", ""),
            "max_score": float(m.get("max_score", 0.0)),
            "min_score": float(m.get("min_score", 0.0)),
            "pwm": pwm_rows,
        }
    return fixed


# ----------------------------- Quality checks --------------------------- #

def check_shape_and_range(pwm_rows: List[List[float]]) -> bool:
    """4xN, all in [0,1]."""
    if not pwm_rows or len(pwm_rows) != 4:
        return False
    N = len(pwm_rows[0])
    if N == 0:
        return False
    for r in pwm_rows:
        if len(r) != N:
            return False
        for v in r:
            if not (0.0 <= v <= 1.0) or math.isnan(v):
                return False
    return True


def check_col_sums(pwm_rows: List[List[float]], tol: float) -> bool:
    """Columns should sum to ~1 within tolerance."""
    N = len(pwm_rows[0])
    for j in range(N):
        s = pwm_rows[0][j] + pwm_rows[1][j] + pwm_rows[2][j] + pwm_rows[3][j]
        if abs(s - 1.0) > tol:
            return False
    return True


def check_min_length(pwm_rows: List[List[float]], min_len: int) -> bool:
    return len(pwm_rows[0]) >= min_len


def check_over_dominance(pwm_rows: List[List[float]], max_share: float) -> bool:
    """
    Sum each base across all positions and divide by motif length.
    Fail if any base's global share exceeds max_share.
    """
    N = len(pwm_rows[0])
    totals = [sum(row) for row in pwm_rows]  # sums across positions
    # since each column ~1, total mass ~ N
    shares = [t / max(N, 1) for t in totals]
    return max(shares) <= max_share


def try_threshold_lightmotif(
    pwm_rows: List[List[float]],
    pvalue: float
) -> Tuple[Optional[float], Optional[float]]:
    """
    If lightmotif is available, compute:
      - threshold = ScoringMatrix.score(pvalue)
      - max_pssm  = sum of per-column maxima in raw probabilities (PFM-ish max)
    Returns (threshold, max_pssm) or (None, None) if unavailable.
    """
    if not _HAS_LIGHTMOTIF:
        return (None, None)

    # Build dict expected by lightmotif
    mat_dict = {BASES[i]: [float(x) for x in pwm_rows[i]] for i in range(4)}

    scoring = None
    # Try the common signatures across lightmotif versions
    try:
        # Some builds expect (values, alphabet)
        scoring = lightmotif.ScoringMatrix(mat_dict, "ACGT")
    except TypeError:
        try:
            # Others expect (alphabet, values)
            scoring = lightmotif.ScoringMatrix("ACGT", mat_dict)
        except TypeError:
            try:
                # Named arguments for maximum compatibility
                scoring = lightmotif.ScoringMatrix(values=mat_dict, alphabet="ACGT")
            except Exception as e:
                print(f"[WARN] lightmotif.ScoringMatrix signature not recognized: {e}")
                return (None, None)
    except Exception as e:
        print(f"[WARN] Could not construct lightmotif ScoringMatrix: {e}")
        return (None, None)

    try:
        thr = float(scoring.score(pvalue))
    except Exception as e:
        print(f"[WARN] lightmotif score(pvalue={pvalue}) failed: {e}")
        return (None, None)

    # max_pssm comparable baseline in prob space (not identical scales—just a sanity metric)
    N = len(pwm_rows[0]) if pwm_rows and pwm_rows[0] else 0
    max_per_col = [max(pwm_rows[0][j], pwm_rows[1][j], pwm_rows[2][j], pwm_rows[3][j]) for j in range(N)]
    max_pssm = float(sum(max_per_col))
    return (thr, max_pssm)


# ----------------------------- Driver ---------------------------------- #

def build_argparser():
    p = argparse.ArgumentParser(
        description="Filter / convert TFBS motifs (MEME or JSON) and write filtered JSON + report."
    )
    g_in = p.add_mutually_exclusive_group(required=True)
    g_in.add_argument("--in-json", help="Input JSON with motifs (4xN rows A,C,G,T or Nx4 will be auto-transposed).")
    g_in.add_argument("--in-meme", nargs="+", help="One or more MEME files.")

    p.add_argument("--out-json", required=True, help="Output JSON with filtered motifs (4xN).")
    p.add_argument("--report", help="TSV report path (optional).")

    # Quality parameters
    p.add_argument("--pvalue", type=float, default=1e-5, help="P-value for optional threshold computation (default: 1e-5).")
    p.add_argument("--max-base-share", type=float, default=0.80, help="Max global share for any single base across the motif (default: 0.80).")
    p.add_argument("--min-length", type=int, default=6, help="Minimum motif length (default: 6).")
    p.add_argument("--col-sum-tol", type=float, default=0.05, help="Tolerance for per-column sum (~1.0) (default: 0.05).")
    p.add_argument("--fail-on-threshold", action="store_true",
                   help="If set, negative/invalid lightmotif thresholds cause failure (default: do NOT fail).")

    return p

def _pfm_to_lod_moods(pwm_rows, background=(0.25, 0.25, 0.25, 0.25), pseudocount=1e-3):
    """Convert a 4xN probability matrix (PFM) to MOODS log-odds via a temp .pfm file."""
    import tempfile, os
    # PFM string: 4 lines (A,C,G,T), space-separated values, newline at end
    pfm_str = "\n".join(" ".join(f"{v:.6f}" for v in row) for row in pwm_rows) + "\n"
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".pfm", delete=False)
    try:
        tmp.write(pfm_str); tmp.flush(); tmp.close()
        lod = parsers.pfm_to_log_odds(tmp.name, list(background), float(pseudocount))
        return lod.tolist() if hasattr(lod, "tolist") else lod
    finally:
        try: os.unlink(tmp.name)
        except Exception: pass


def try_threshold(pwm_rows, pvalue):
    """
    Prefer MOODS for computing an absolute threshold from p-value.
    Falls back to lightmotif (your existing function) if MOODS is unavailable.
    Returns (threshold, max_pssm, engine_name) where engine_name is 'MOODS', 'lightmotif', or 'none'.
    """
    # 1) MOODS path
    if _HAS_MOODS:
        try:
            lod = _pfm_to_lod_moods(pwm_rows, background=(0.25, 0.25, 0.25, 0.25), pseudocount=1e-3)
            thr = float(moods_tools.threshold_from_p(lod, [0.25, 0.25, 0.25, 0.25], float(pvalue)))
            # Simple “max PSSM” baseline in prob space: sum of per-column maxima
            N = len(pwm_rows[0]) if pwm_rows and pwm_rows[0] else 0
            max_pssm = float(sum(max(pwm_rows[0][j], pwm_rows[1][j], pwm_rows[2][j], pwm_rows[3][j]) for j in range(N)))
            return thr, max_pssm, "MOODS"
        except Exception as e:
            print(f"[WARN] MOODS threshold failed: {e}")

    # 2) Fallback to your lightmotif attempt (if present in your script)
    try:
        thr_lm, max_pssm_lm = try_threshold_lightmotif(pwm_rows, pvalue)  # your existing function
        if thr_lm is not None:
            return thr_lm, max_pssm_lm, "lightmotif"
    except NameError:
        pass

    # 3) No engine available
    return None, None, "none"



def main():
    args = build_argparser().parse_args()

    # ---- Load motifs ----
    motifs: Dict[str, Dict] = {}
    if args.in_json:
        motifs = load_json_motifs(args.in_json)
        print(f"[INFO] Loaded {len(motifs)} motifs from {args.in_json}")
    else:
        total = 0
        for path_meme in args.in_meme:
            part = parse_meme_file(path_meme)
            motifs.update(part)
            total += len(part)
        print(f"[INFO] Parsed {total} motifs from {len(args.in_meme)} MEME file(s)")

    # ---- Filter ----
    failed_shape = 0
    failed_dominance = 0
    failed_colsum = 0
    failed_minlen = 0
    failed_threshold = 0
    passed = 0

    out_motifs: Dict[str, Dict] = {}
    report_rows: List[Tuple[str, str, str]] = []  # (name, status, notes)

    for name, m in motifs.items():
        pwm_rows = m.get("pwm", [])
        notes: List[str] = []
        status = "PASS"

        # 1) shape/range
        if not check_shape_and_range(pwm_rows):
            failed_shape += 1
            status = "FAIL"
            notes.append("bad-shape-or-range")
        else:
            # 2) column sums
            if not check_col_sums(pwm_rows, tol=args.col_sum_tol):
                failed_colsum += 1
                status = "FAIL"
                notes.append("colsum-out-of-range")

            # 3) min length
            if not check_min_length(pwm_rows, args.min_length):
                failed_minlen += 1
                status = "FAIL"
                notes.append("too-short")

            # 4) over-dominance
            if not check_over_dominance(pwm_rows, args.max_base_share):
                failed_dominance += 1
                status = "FAIL"
                notes.append("over-dominant-base")

            # 5) optional threshold sanity (lightmotif)
            thr, max_pssm, thr_engine = try_threshold(pwm_rows, args.pvalue)

            if thr is not None:
                # We only "fail" if user requested strictness.
                if (math.isnan(thr) or math.isinf(thr)) or (args.fail_on_threshold and thr <= 0):
                    failed_threshold += 1
                    status = "FAIL"
                    notes.append(f"threshold-invalid:{thr}")
                # Note: max_pssm is in prob space; thr is a scoring value—don’t hard-fail on comparison.
                # You can still record it:
                notes.append(f"thr={thr:.4f};max_pssm={max_pssm:.4f}")

        if status == "PASS":
            passed += 1
            # Ensure consistent metadata fields, recompute consensus/min/max if missing
            pwm_cols = list(map(list, zip(*pwm_rows)))  # back to Nx4 just to compute
            min_s, max_s = compute_min_max_scores_from_cols(pwm_cols)
            consensus = compute_consensus_from_cols(pwm_cols)
            out_motifs[name] = {
                "name": m.get("name", name),
                "description": m.get("description", "Unknown transcription factor"),
                "species": m.get("species", "Arabidopsis thaliana"),
                "link": normalize_tfdb_link(m.get("link", ""), name),
                "consensus": m.get("consensus", consensus),
                "max_score": float(m.get("max_score", max_s)),
                "min_score": float(m.get("min_score", min_s)),
                "pwm": pwm_rows,  # 4xN
            }

        report_rows.append((name, status, ";".join(notes) if notes else ""))

    # ---- Write outputs ----
    with open(args.out_json, "w", encoding="utf-8") as out_f:
        json.dump(out_motifs, out_f, indent=2)
    print(f"[INFO] Wrote {len(out_motifs)} passing motifs to {args.out_json}")

    if args.report:
        with open(args.report, "w", encoding="utf-8") as rep:
            rep.write("motif\tstatus\tnotes\n")
            for name, status, notes in sorted(report_rows):
                rep.write(f"{name}\t#{status}\t{notes}\n")
        print(f"[INFO] Wrote report to {args.report}")

    # ---- Summary ----
    total_in = len(motifs)
    print("=== Filtering summary ===")
    print(f"Total motifs:                       {total_in}")
    print(f"Failed shape/range/colsums/minlen:  {failed_shape + failed_colsum + failed_minlen}")
    print(f"Failed over-dominance (> {args.max_base_share:.2f}): {failed_dominance}")
    if _HAS_LIGHTMOTIF:
        print(f"Failed threshold (strict):          {failed_threshold}  "
              f"{'(lightmotif present)' if _HAS_LIGHTMOTIF else ''}")
    else:
        print("Threshold check skipped (lightmotif not installed).")
    print(f"Passed all checks:                  {passed}")


if __name__ == "__main__":
    main()
