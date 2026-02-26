#!/usr/bin/env python3

'''
Example usage: 
python antismash/generic_modules/coexpress/utils/count_sra_only_geos.py \
  antismash/generic_modules/coexpress/data/geo_soft \
  --list --show-files
  '''

from __future__ import annotations

import argparse
import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

GSE_RE = re.compile(r"(GSE\d+)", re.IGNORECASE)

# --- Heuristics for supplementary matrix-like files ---
MATRIX_POSITIVE_PATTERNS = [
    r"series_matrix", r"matrix", r"expression", r"expr",
    r"count", r"counts", r"raw_count", r"feature_count",
    r"tpm", r"fpkm", r"rpkm", r"cpm",
    r"normalized", r"normalised", r"norm", r"log2",
    r"abundance",
]
MATRIX_NEGATIVE_PATTERNS = [
    r"readme", r"barcode", r"design", r"protocol", r"annotation",
    r"metadata", r"sample", r"phenodata", r"pheno",
    r"fastq", r"\.sra$", r"\.bam$", r"\.sam$",
]
MATRIX_EXT_OK = (".txt", ".tsv", ".csv", ".mtx", ".xls", ".xlsx", ".rds", ".h5", ".h5ad", ".loom")
COMPRESSED_EXT = (".gz", ".bz2", ".zip")

def is_matrix_like_filename(url: str) -> bool:
    u = url.lower()
    # strip compression suffix for extension checks
    base = u
    for cext in COMPRESSED_EXT:
        if base.endswith(cext):
            base = base[: -len(cext)]
            break

    if not base.endswith(MATRIX_EXT_OK):
        # still allow series_matrix even if extension is odd
        if "series_matrix" not in u:
            return False

    # negative filters
    for pat in MATRIX_NEGATIVE_PATTERNS:
        if re.search(pat, u):
            return False

    # positive signals
    for pat in MATRIX_POSITIVE_PATTERNS:
        if re.search(pat, u):
            return True

    # if it has a good extension but no keyword, treat as unknown -> False
    return False


@dataclass
class GSEInfo:
    gse: str
    soft_path: Path
    n_samples: int = 0
    sample_types: Dict[str, int] = field(default_factory=dict)
    any_sra_samples: bool = False

    # Evidence of "works now"
    any_inline_tables: bool = False  # sample_table blocks OR data_row_count>0
    data_row_count_sum: int = 0

    # Evidence of "works if download"
    series_matrix_urls: List[str] = field(default_factory=list)       # Type 1
    supp_matrix_like_urls: List[str] = field(default_factory=list)    # Type 2
    supp_other_urls: List[str] = field(default_factory=list)

    def classify(self) -> str:
        """
        Returns one of:
          - works_now_inline
          - works_if_type1_series_matrix
          - works_if_type2_supp_matrix
          - unsupported_needs_sra_reprocessing
          - unsupported_unknown
        """
        if self.any_inline_tables:
            return "works_now_inline"
        if self.series_matrix_urls:
            return "works_if_type1_series_matrix"
        if self.supp_matrix_like_urls:
            return "works_if_type2_supp_matrix"
        if self.any_sra_samples:
            return "unsupported_needs_sra_reprocessing"
        return "unsupported_unknown"


def open_text_maybe_gz(path: Path) -> Iterable[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                yield line
    else:
        with path.open("rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                yield line


def parse_family_soft(path: Path) -> Optional[GSEInfo]:
    # Try to infer GSE id from filename first
    m = GSE_RE.search(path.name)
    if not m:
        return None
    gse = m.group(1).upper()

    info = GSEInfo(gse=gse, soft_path=path)

    in_sample = False
    in_table = False

    for raw in open_text_maybe_gz(path):
        line = raw.strip()

        # section boundaries
        if line.startswith("^SAMPLE"):
            info.n_samples += 1
            in_sample = True
            in_table = False
            continue
        if line.startswith("^") and not line.startswith("^SAMPLE"):
            in_sample = False
            in_table = False

        # sample markers
        if in_sample:
            if line.startswith("!Sample_type ="):
                t = line.split("=", 1)[1].strip()
                info.sample_types[t] = info.sample_types.get(t, 0) + 1
                if t.upper() == "SRA":
                    info.any_sra_samples = True

            if line.startswith("!Sample_data_row_count ="):
                try:
                    n = int(line.split("=", 1)[1].strip())
                except ValueError:
                    n = 0
                info.data_row_count_sum += n
                if n > 0:
                    info.any_inline_tables = True

            # table blocks exist even if data_row_count is sometimes 0
            if line.lower().startswith("!sample_table_begin"):
                in_table = True
                info.any_inline_tables = True
            elif line.lower().startswith("!sample_table_end"):
                in_table = False

        # series supplementary / matrix pointers
        if line.startswith("!Series_supplementary_file"):
            url = line.split("=", 1)[1].strip()
            low = url.lower()
            if "series_matrix" in low:
                info.series_matrix_urls.append(url)
            elif is_matrix_like_filename(url):
                info.supp_matrix_like_urls.append(url)
            else:
                info.supp_other_urls.append(url)

    # de-dup URLs
    info.series_matrix_urls = sorted(set(info.series_matrix_urls))
    info.supp_matrix_like_urls = sorted(set(info.supp_matrix_like_urls))
    info.supp_other_urls = sorted(set(info.supp_other_urls))
    return info


def iter_soft_files(dirpath: Path) -> List[Path]:
    paths = []
    for p in dirpath.iterdir():
        if p.is_file() and ("_family.soft" in p.name):
            paths.append(p)
    return sorted(paths)


def choose_one_file_per_gse(paths: List[Path]) -> Dict[str, Path]:
    """
    Prefer .soft.gz if both exist, otherwise whatever exists.
    """
    best: Dict[str, Path] = {}
    for p in paths:
        m = GSE_RE.search(p.name)
        if not m:
            continue
        gse = m.group(1).upper()
        if gse not in best:
            best[gse] = p
        else:
            # Prefer gz (smaller, canonical download)
            if best[gse].suffix != ".gz" and p.suffix == ".gz":
                best[gse] = p
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("soft_dir", type=Path, help="Directory containing *_family.soft or *_family.soft.gz")
    ap.add_argument("--list", action="store_true", help="Print per-GSE classification")
    ap.add_argument("--show-files", action="store_true", help="Also print file pointers for Type1/Type2")
    args = ap.parse_args()

    all_paths = iter_soft_files(args.soft_dir)
    best = choose_one_file_per_gse(all_paths)

    infos: List[GSEInfo] = []
    for gse, p in sorted(best.items()):
        info = parse_family_soft(p)
        if info:
            infos.append(info)

    # counts
    counts = {
        "total_gse": len(infos),
        "works_now_inline": 0,
        "works_if_type1_series_matrix": 0,
        "works_if_type2_supp_matrix": 0,
        "unsupported_needs_sra_reprocessing": 0,
        "unsupported_unknown": 0,
        "is_sra_submission": 0,
    }

    for i in infos:
        counts[i.classify()] += 1
        if i.any_sra_samples:
            counts["is_sra_submission"] += 1

    print("Counts (UNIQUE GSEs):")
    for k in [
        "total_gse",
        "works_now_inline",
        "works_if_type1_series_matrix",
        "works_if_type2_supp_matrix",
        "unsupported_needs_sra_reprocessing",
        "unsupported_unknown",
        "is_sra_submission",
    ]:
        print(f"  {k}: {counts[k]}")

    if args.list:
        print("\nPer-GSE summary:")
        for i in infos:
            cls = i.classify()
            sample_types = ", ".join(f"{k}:{v}" for k, v in sorted(i.sample_types.items()))
            print(f"  {i.gse} ({i.n_samples} samples) -> {cls} | sample_types=[{sample_types}]")

            if args.show_files and cls in ("works_if_type1_series_matrix", "works_if_type2_supp_matrix"):
                if i.series_matrix_urls:
                    print("    Type1 series_matrix pointer(s):")
                    for u in i.series_matrix_urls:
                        print(f"      - {u}")
                if i.supp_matrix_like_urls:
                    print("    Type2 supp matrix-like pointer(s):")
                    for u in i.supp_matrix_like_urls:
                        print(f"      - {u}")

            if args.show_files and cls.startswith("unsupported"):
                # still show non-matrix supplementary files (useful debugging)
                if i.supp_other_urls:
                    print("    Non-matrix supplementary pointer(s) (ignored):")
                    for u in i.supp_other_urls[:6]:
                        print(f"      - {u}")
                    if len(i.supp_other_urls) > 6:
                        print(f"      ... ({len(i.supp_other_urls)-6} more)")

if __name__ == "__main__":
    main()