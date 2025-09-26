# validation.py
# Helpers for TAIR symbols and connectTF validation/evidence

#
# Copyright (C) 2025 Elena Del Pupo
# Wageningen University & Research
# Bioinformatics Group
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from __future__ import annotations
from typing import Dict, Tuple
import os, csv, re, glob
from antismash import utils
import logging

# -------- Paths (data/ next to this module) --------
DATA_DIR = utils.get_full_path(__file__, "data")
if not os.path.isdir(DATA_DIR):  # fallback if utils fails
    DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

CONNECTF_VALIDATED_CSV = os.path.join(DATA_DIR, "connectTF_validated.csv")
CONNECTF_META_DIR = os.path.join(DATA_DIR, "connectTF_metadata")

# -------- Caches --------
_TAIR_SYMBOLS_CACHE: Dict[str, str] = {}
_VALIDATED_MAP: Dict[Tuple[str, str], str] = {}
_META_CACHE: Dict[str, Dict[str, str]] = {}

# -------- Utilities --------
def norm_agi(s: str) -> str:
    """Uppercase AGI locus and drop model suffix like '.1'."""
    if not s:
        return ""
    s = s.strip().upper()
    return re.sub(r"\.\d+$", "", s)

def _extract_symbol_from_primary(text: str) -> str:
    """Return the last (...) chunk from 'Primary Gene Symbol'."""
    if not text:
        return ""
    parts = re.findall(r"\(([^)]+)\)", text)
    return parts[-1].strip() if parts else ""

def _find_tair_tsv_in_datadir() -> str:
    """Pick newest TAIR_GeneDescriptions_*.tsv from data/."""
    candidates = sorted(glob.glob(os.path.join(DATA_DIR, "TAIR_GeneDescriptions_*.tsv")), reverse=True)
    return candidates[0] if candidates else ""

# -------- Public API --------
def load_tair_symbols() -> Dict[str, str]:
    """Return {AGI -> primary symbol (paren chunk)}. Cached."""
    global _TAIR_SYMBOLS_CACHE
    if _TAIR_SYMBOLS_CACHE:
        return _TAIR_SYMBOLS_CACHE

    path = _find_tair_tsv_in_datadir()
    mapping: Dict[str, str] = {}
    if not path:
        logging.warning("validation: no TAIR_GeneDescriptions_*.tsv found in %s", DATA_DIR)
        _TAIR_SYMBOLS_CACHE = mapping
        return mapping

    try:
        with open(path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                agi = norm_agi(row.get("Locus Identifier", ""))
                primary = (row.get("Primary Gene Symbol") or "").strip()
                sym = _extract_symbol_from_primary(primary)
                if agi and sym:
                    mapping[agi] = sym
        logging.info("validation: loaded %d TAIR symbols from %s", len(mapping), os.path.basename(path))
    except Exception as e:
        logging.warning("validation: failed reading TAIR TSV '%s': %s", path, e)

    _TAIR_SYMBOLS_CACHE = mapping
    return mapping

def load_validated_map() -> Dict[Tuple[str, str], str]:
    """Return {(TF, target) -> meta_id} from connectTF_validated.csv. Cached."""
    global _VALIDATED_MAP
    if _VALIDATED_MAP:
        return _VALIDATED_MAP

    mapping: Dict[Tuple[str, str], str] = {}
    if not os.path.isfile(CONNECTF_VALIDATED_CSV):
        logging.info("validation: no connectTF_validated.csv at %s", CONNECTF_VALIDATED_CSV)
        _VALIDATED_MAP = mapping
        return mapping

    try:
        with open(CONNECTF_VALIDATED_CSV, newline="", encoding="utf-8") as fh:
            r = csv.reader(fh)
            header = next(r, None)  # TF,target,metadata
            for row in r:
                if len(row) < 3:
                    continue
                tf = norm_agi(row[0])
                target = norm_agi(row[1])
                meta = os.path.splitext(os.path.basename(row[2]))[0]  # strip .txt
                if tf and target and meta:
                    mapping[(tf, target)] = meta
        logging.info("validation: loaded %d validated TF-target pairs", len(mapping))
    except Exception as e:
        logging.warning("validation: failed to read connectTF_validated.csv: %s", e)

    _VALIDATED_MAP = mapping
    return mapping

def evidence_summary(meta_id: str) -> str:
    """Return concise 'DOI; Source; Type; Method; Edge' string from metadata id."""
    md = _parse_metadata_file(meta_id)
    if not md:
        return ""
    parts = []
    if md.get("DOI"): parts.append(f"DOI: {md['DOI']}")
    if md.get("Data_Source/Publication"): parts.append(f"Source: {md['Data_Source/Publication']}")
    if md.get("Experiment_Type"): parts.append(f"Type: {md['Experiment_Type']}")
    if md.get("Technology/Method"): parts.append(f"Method: {md['Technology/Method']}")
    if md.get("Edge_Type"): parts.append(f"Edge: {md['Edge_Type']}")
    return "; ".join(parts)

# -------- Internal: metadata parsing --------
def _parse_metadata_file(meta_id: str) -> Dict[str, str]:
    global _META_CACHE
    if meta_id in _META_CACHE:
        return _META_CACHE[meta_id]
    d: Dict[str, str] = {}
    path = os.path.join(CONNECTF_META_DIR, f"{meta_id}.txt")
    if not os.path.isfile(path):
        _META_CACHE[meta_id] = d
        return d
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or ":" not in line:
                    continue
                m = re.match(r"^\*?([^:]+):\s*(.*)$", line)
                if not m:
                    continue
                key = m.group(1).strip()
                val = m.group(2).strip()
                d[key] = val
    except Exception as e:
        logging.warning("validation: failed to read metadata '%s': %s", path, e)
    _META_CACHE[meta_id] = d
    return d

# -------- Optional tiny CLI for quick checks --------
if __name__ == "__main__":
    import sys
    if len(sys.argv) >= 3:
        tf, target = norm_agi(sys.argv[1]), norm_agi(sys.argv[2])
        meta = load_validated_map().get((tf, target))
        print("validated:", bool(meta))
        if meta:
            print("evidence:", evidence_summary(meta))
    else:
        print("Usage: python -m tfbs_finder.validation AT1G01060 AT1G01010")
