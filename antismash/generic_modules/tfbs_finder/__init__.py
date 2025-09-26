
# Copyright (C) 2025 Hannah E. Augustijn 
# Wageningen University & Research & Leiden University
# Department: Department of Bioinformatics & Institute of Biology Leiden
#
# Copyright (C) 2025 Elena Del Pupo
# Wageningen University & Research
# Bioinformatics Group
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Transcription Factor Binding Site (TFBS) module for plantiSMASH """

from typing import List
from argparse import Namespace
import os
import sys
import time
import hashlib
import logging
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

from Bio.SeqRecord import SeqRecord
from antismash import utils

from .tfbs_detection import (
    PWM_PATH,
    run_tfbs_finder,
    TFBSFinderResults,
)
from . import output  # For HTML rendering

NAME = "tfbs_finder"
SHORT_DESCRIPTION = "Detects transcription factor binding sites (TFBSs)"
REQUIRES = [(PWM_PATH, False)]

# Safe defaults used only when TFBS is enabled
_DEFAULT_TFBS_P = 1e-4
_DEFAULT_TFBS_RANGE = 1000

# ------------------------------------------------------------
# ConnecTF validation file (large, hosted on Zenodo)
# ------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
CONNECT_TF_FILENAME = "connectTF_validated.csv"
CONNECT_TF_PATH = os.path.join(DATA_DIR, CONNECT_TF_FILENAME)

# Direct file endpoint (Zenodo serves with ?download=1)
CONNECT_TF_URL = os.environ.get(
    "PLANTISMASH_TFBS_CONNECT_TF_URL",
    "https://zenodo.org/records/16926250/files/connectTF_validated.csv?download=1",
)

# Expected checksum (MD5). You can override with:
#   1) env var PLANTISMASH_TFBS_CONNECT_TF_MD5
#   2) a sidecar file data/connectTF_validated.csv.md5 whose first token is the hex digest
_DEFAULT_CONNECT_TF_MD5 = "08fc6daa1b4088a6ca341e3b47dca893"
CONNECT_TF_MD5 = os.environ.get("PLANTISMASH_TFBS_CONNECT_TF_MD5", _DEFAULT_CONNECT_TF_MD5).strip().lower()

# Optional: control progress bar
# values: "auto" (default), "on", "off"
PROGRESS_MODE = os.environ.get("PLANTISMASH_TFBS_PROGRESS", "auto").strip().lower()

# Other expected assets in data/
TAIR_TSV = os.path.join(DATA_DIR, "TAIR_GeneDescriptions_2025-08-19.tsv")
CONNECT_TF_META_DIR = os.path.join(DATA_DIR, "connectTF_metadata")


def _enabled(options) -> bool:
    # Only run when the CLI flag --tfbs-detection was set (stored in options.tfbs)
    return bool(getattr(options, "tfbs", False))


def _resolve_p(options) -> float:
    p = getattr(options, "tfbs_pvalue", None)
    if p is None or p <= 0:
        if _enabled(options):
            logging.warning("TFBS: invalid/empty p-value %r; using default %.1e", p, _DEFAULT_TFBS_P)
        return _DEFAULT_TFBS_P
    return p


def _resolve_range(options) -> int:
    r = getattr(options, "tfbs_range", None)
    if r is None or r <= 0:
        if _enabled(options):
            logging.warning("TFBS: invalid/empty scanning range %r; using default %d bp", r, _DEFAULT_TFBS_RANGE)
        return _DEFAULT_TFBS_RANGE
    return r


def _human_size(nbytes: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if nbytes < 1024 or unit == "GB":
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024.0
    return f"{nbytes:.1f} GB"


def _is_tty() -> bool:
    if PROGRESS_MODE == "on":
        return True
    if PROGRESS_MODE == "off":
        return False
    try:
        return sys.stderr.isatty()
    except Exception:
        return False


def _draw_progress_bar(downloaded: int, total: int, width: int = 40) -> None:
    if not total:
        return
    frac = min(max(downloaded / total, 0.0), 1.0)
    filled = int(frac * width)
    bar = "#" * filled + "-" * (width - filled)
    msg = f"\r[ {bar} ] {frac*100:5.1f}%  ({_human_size(downloaded)} / {_human_size(total)})"
    sys.stderr.write(msg)
    sys.stderr.flush()


def _read_expected_md5(dest_path: str) -> str:
    """Return expected MD5 from env, sidecar, or module default."""
    env_val = os.environ.get("PLANTISMASH_TFBS_CONNECT_TF_MD5", "").strip().lower()
    if env_val:
        return env_val
    sidecar = dest_path + ".md5"
    if os.path.exists(sidecar):
        try:
            with open(sidecar, "r", encoding="utf-8") as fh:
                line = fh.readline().strip()
            if not line:
                return CONNECT_TF_MD5
            # allow "<hex>  filename" or just "<hex>"
            first = line.split()[0]
            return first.lower()
        except Exception:
            return CONNECT_TF_MD5
    return CONNECT_TF_MD5


def _md5_file(path: str, chunk: int = 4 * 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            buf = fh.read(chunk)
            if not buf:
                break
            h.update(buf)
    return h.hexdigest()


def _verify_md5_or_raise(path: str) -> None:
    expected = _read_expected_md5(path)
    if not expected:
        logging.warning("TFBS: no MD5 provided for %s; skipping verification", os.path.basename(path))
        return
    actual = _md5_file(path)
    if actual != expected:
        raise RuntimeError(
            f"MD5 mismatch for {os.path.basename(path)}: expected {expected}, got {actual}"
        )
    logging.info("TFBS: checksum OK for %s", os.path.basename(path))


def _download_large_file(url: str, dest_path: str) -> None:
    """Download to dest_path (atomic via .part), with TTY progress bar + INFO logs."""
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    tmp_path = dest_path + ".part"

    logging.info(
        "TFBS: downloading validation file of TFBS–target hits from ConnecTF\n"
        "      %s\n      → %s",
        url, dest_path
    )

    req = Request(url, headers={"User-Agent": "plantiSMASH-tfbs/1.0"})
    start = time.time()
    try:
        with urlopen(req) as resp:
            total = int(resp.headers.get("Content-Length", "0")) or None
            chunk = 1024 * 1024  # 1 MB
            downloaded = 0
            last_log = start
            use_bar = _is_tty() and bool(total)

            with open(tmp_path, "wb") as out:
                while True:
                    buf = resp.read(chunk)
                    if not buf:
                        break
                    out.write(buf)
                    downloaded += len(buf)

                    now = time.time()
                    if use_bar:
                        _draw_progress_bar(downloaded, total)
                    elif now - last_log > 2:
                        if total:
                            pct = 100.0 * downloaded / total
                            logging.info(
                                "TFBS: download progress %s / %s (%.1f%%)",
                                _human_size(downloaded), _human_size(total), pct
                            )
                        else:
                            logging.info("TFBS: download progress %s", _human_size(downloaded))
                        last_log = now

        if _is_tty() and total:
            sys.stderr.write("\n")
            sys.stderr.flush()

        # basic sanity
        if not os.path.exists(tmp_path) or os.path.getsize(tmp_path) < 1024:
            raise IOError("downloaded file is unexpectedly small (<1KB)")

        # atomic replace
        if os.path.exists(dest_path):
            os.remove(dest_path)
        os.rename(tmp_path, dest_path)

        # verify MD5 (will raise on mismatch)
        _verify_md5_or_raise(dest_path)

        dt = time.time() - start
        logging.info(
            "TFBS: ConnecTF validation file ready at %s (%s, %.1fs)",
            dest_path, _human_size(os.path.getsize(dest_path)), dt
        )
    except (HTTPError, URLError, IOError, RuntimeError) as e:
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except Exception:
            pass
        raise RuntimeError(f"Failed to download ConnecTF validation file: {e}")


def _ensure_connect_tf_file() -> None:
    """Ensure the ConnecTF CSV exists in data/; download if missing or bad."""
    if os.path.exists(CONNECT_TF_PATH):
        try:
            _verify_md5_or_raise(CONNECT_TF_PATH)
            return
        except Exception as e:
            logging.warning("TFBS: existing file failed checksum (%s); re-downloading", e)
            try:
                os.remove(CONNECT_TF_PATH)
            except Exception:
                pass
    _download_large_file(CONNECT_TF_URL, CONNECT_TF_PATH)


def _ensure_assets_when_enabled(options) -> List[str]:
    """Ensure required assets exist if TFBS is enabled; return problems for prereq reporting."""
    issues: List[str] = []
    if not _enabled(options):
        return issues

    # PWM JSON must exist locally (do not fetch here)
    if not os.path.exists(PWM_PATH):
        issues.append(f"Missing PWM JSON file at: {PWM_PATH}")

    # TAIR TSV and metadata dir are expected to be present locally
    if not os.path.exists(TAIR_TSV):
        issues.append(f"Missing TAIR gene descriptions TSV: {TAIR_TSV}")
    if not os.path.isdir(CONNECT_TF_META_DIR):
        issues.append(f"Missing ConnecTF metadata directory: {CONNECT_TF_META_DIR}")

    # ConnecTF CSV: try to fetch (and verify) if missing/bad
    try:
        _ensure_connect_tf_file()
    except Exception as e:
        issues.append(str(e))

    return issues


def check_prereqs(options) -> List[str]:
    """Only enforce files when TFBS is enabled. Will attempt to download ConnecTF CSV if absent."""
    return _ensure_assets_when_enabled(options)


def check_options(options) -> List[str]:
    """Validate TFBS-related options (only when enabled)."""
    if not _enabled(options):
        return []
    issues = []
    if getattr(options, "tfbs_pvalue", None) is not None and options.tfbs_pvalue <= 0:
        issues.append("TFBS p-value must be positive.")
    if getattr(options, "tfbs_range", None) is not None and options.tfbs_range <= 0:
        issues.append("TFBS scanning range must be positive.")
    return issues


def run_tfbs_finder_for_record(record: SeqRecord, options) -> None:
    """Record-level interface for TFBS analysis (synchronous)."""
    if not _enabled(options):
        logging.debug("TFBS disabled; skipping record %s", getattr(record, "id", "?"))
        return

    # Ensure assets (download/verify ConnecTF if needed)
    problems = _ensure_assets_when_enabled(options)
    if problems:
        for p in problems:
            logging.error("TFBS: %s", p)
        # proceed with empty results so the pipeline doesn't crash
        results = TFBSFinderResults(
            record_id=record.id,
            pvalue=_resolve_p(options),
            start_overlap=_resolve_range(options),
            hits_by_record={record.id: []},
        )
        options.extrarecord.setdefault(record.id, Namespace())
        ns = options.extrarecord[record.id]
        ns.extradata = getattr(ns, "extradata", {})
        ns.extradata["TFBSFinderResults"] = results
        return

    p = _resolve_p(options)
    r = _resolve_range(options)

    utils.log_status("➡️ Running TFBS analysis for contig #%d [p=%.1e, range=%dbp]"
                     % (getattr(options, "record_idx", 0), p, r))

    try:
        results: TFBSFinderResults = run_tfbs_finder(record, p, r)
        nhits = len(getattr(results, "hits_by_record", {}).get(record.id, []))
        logging.info("TFBS: %s finished with %d hit record(s)", record.id, nhits)
    except Exception:
        logging.exception("💥 TFBS crashed on %s", record.id)
        results = TFBSFinderResults(
            record_id=record.id,
            pvalue=p,
            start_overlap=r,
            hits_by_record={record.id: []},
        )

    # stash results for the output module
    options.extrarecord.setdefault(record.id, Namespace())
    ns = options.extrarecord[record.id]
    ns.extradata = getattr(ns, "extradata", {})
    ns.extradata["TFBSFinderResults"] = results
    logging.info("✅ TFBS finished for record: %s", record.id)


# Output integration
def will_handle(product: str) -> bool:
    return True


def generate_details_div(*args, **kwargs):
    return output.generate_details_div(*args, **kwargs)