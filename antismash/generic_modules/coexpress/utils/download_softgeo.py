#!/usr/bin/env python3
"""
Download GEO Series SOFT files in "family" format.

Given one or more GSE accessions (e.g., GSE41935), downloads:
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41935/soft/GSE41935_family.soft.gz

Optionally decompresses .gz to .soft.

Usage examples:
  python download_softgeo.py --gse GSE41935 GSE80744 -o antismash/generic_modules/coexpress/data/geo_soft
  python download_softgeo.py --list antismash/generic_modules/coexpress/data/geo_soft_series.txt -o ... --decompress
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import sys
import urllib.request
from pathlib import Path
from typing import Iterable, List


GSE_RE = re.compile(r"^(GSE)(\d+)$", re.IGNORECASE)


def normalize_gse(gse: str) -> str:
    gse = gse.strip()
    m = GSE_RE.match(gse)
    if not m:
        raise ValueError(f"Invalid GEO Series accession: {gse!r} (expected like 'GSE41935')")
    return f"GSE{m.group(2)}"


def series_bucket(gse: str) -> str:
    """
    GEO series are grouped like:
      GSE41935 -> GSE41nnn
      GSE80744 -> GSE80nnn
      GSE123456 -> GSE123nnn
    i.e., replace last 3 digits with 'nnn'.
    """
    digits = gse[3:]
    if len(digits) < 4:
        # very old/small GSEs are rare; this keeps behavior predictable
        raise ValueError(f"Unexpectedly short GSE number: {gse}")
    return f"GSE{digits[:-3]}nnn"


def soft_url(gse: str) -> str:
    bucket = series_bucket(gse)
    return (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{bucket}/{gse}/soft/{gse}_family.soft.gz"
    )


def read_gse_list(path: Path) -> List[str]:
    ids: List[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        ids.append(normalize_gse(line))
    return ids


def download(url: str, dest: Path, force: bool = False) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        return
    try:
        urllib.request.urlretrieve(url, dest)  # nosec (expected URL)
    except Exception as e:
        if dest.exists():
            # remove partial file
            dest.unlink()
        raise RuntimeError(f"Failed to download {url}: {e}") from e


def decompress_gz(gz_path: Path, out_path: Path, force: bool = False) -> None:
    if out_path.exists() and not force:
        return
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def iter_unique(items: Iterable[str]) -> List[str]:
    seen = set()
    out = []
    for x in items:
        if x not in seen:
            out.append(x)
            seen.add(x)
    return out


def main(argv: List[str]) -> int:
    p = argparse.ArgumentParser(description="Download GEO SOFT family files for GSE accessions.")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--gse", nargs="+", help="One or more GSE IDs (e.g., GSE41935).")
    group.add_argument("--list", type=Path, help="Path to text file with one GSE per line.")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="Output directory for downloads.")
    p.add_argument("--decompress", action="store_true", help="Also write decompressed .soft files.")
    p.add_argument("--force", action="store_true", help="Re-download / overwrite existing files.")
    args = p.parse_args(argv)

    if args.gse:
        gses = [normalize_gse(x) for x in args.gse]
    else:
        gses = read_gse_list(args.list)

    gses = iter_unique(gses)

    ok = 0
    failed = 0

    for gse in gses:
        url = soft_url(gse)
        gz_path = args.outdir / f"{gse}_family.soft.gz"
        soft_path = args.outdir / f"{gse}_family.soft"

        try:
            print(f"[download] {gse} -> {gz_path}")
            download(url, gz_path, force=args.force)

            if args.decompress:
                print(f"[decompress] {gz_path} -> {soft_path}")
                decompress_gz(gz_path, soft_path, force=args.force)

            ok += 1
        except Exception as e:
            failed += 1
            print(f"[ERROR] {gse}: {e}", file=sys.stderr)

    print(f"\nDone. OK={ok}, failed={failed}")
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))