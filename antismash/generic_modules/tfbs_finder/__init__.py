# License: GNU Affero General Public License v3 or later
""" Transcription Factor Binding Site (TFBS) module for plantiSMASH """

from typing import Any, Dict, List, Optional
from argparse import Namespace
import os
import time
import logging
import multiprocessing as mp

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

# Optional safety net: set PLANTISMASH_TFBS_TIMEOUT=0 to disable
_DEFAULT_TIMEOUT_SEC = 600


def check_prereqs(options) -> List[str]:
    if not os.path.exists(PWM_PATH):
        return [f"Missing PWM JSON file at: {PWM_PATH}"]
    return []


def check_options(options) -> List[str]:
    issues = []
    if getattr(options, "tfbs_pvalue", 0) <= 0:
        issues.append("TFBS p-value must be positive.")
    if getattr(options, "tfbs_range", 0) <= 0:
        issues.append("TFBS scanning range must be positive.")
    return issues


def _ensure_extradata(options, rec_id: str) -> Namespace:
    ns = options.extrarecord.setdefault(rec_id, Namespace())
    if not hasattr(ns, "extradata") or ns.extradata is None:
        ns.extradata = {}
    return ns


def _run_with_timeout(func, args, timeout_sec: int):
    """Run func(*args) in a subprocess with an optional timeout."""
    if not timeout_sec or timeout_sec <= 0:
        return func(*args)

    q: mp.Queue = mp.Queue(1)

    def _worker(_q, _func, _args):
        try:
            _q.put(("ok", _func(*_args)))
        except Exception as e:  # noqa
            logging.exception("TFBS worker crashed")
            _q.put(("err", e))

    p = mp.Process(target=_worker, args=(q, func, args))
    p.start()
    p.join(timeout_sec)

    if p.is_alive():
        p.terminate()
        p.join(1)
        raise TimeoutError(f"TFBS run exceeded {timeout_sec} s")

    status, payload = q.get()
    if status == "ok":
        return payload
    raise payload  # re-raise original exception


def run_tfbs_finder_for_record(record: SeqRecord, options) -> None:
    """Record-level driver for TFBS analysis with logging, timeout, and safe storage."""
    rec_id = record.id
    utils.log_status(
        "➡️ Running TFBS analysis for contig #%d [p=%.1e, range=%dbp]"
        % (options.record_idx, options.tfbs_pvalue, options.tfbs_range)
    )

    ns = _ensure_extradata(options, rec_id)

    # Quick sanity: any clusters on this record?
    clusters = utils.get_sorted_cluster_features(record)
    logging.info("TFBS: %s has %d cluster(s)", rec_id, len(clusters))
    if not clusters:
        ns.extradata["TFBSFinderResults"] = TFBSFinderResults()  # empty container
        logging.info("TFBS: skipping %s (no clusters)", rec_id)
        return

    # Optional timeout (env or attribute); set PLANTISMASH_TFBS_TIMEOUT=0 to disable
    timeout = getattr(options, "tfbs_timeout", None)
    if timeout is None:
        try:
            timeout = int(os.environ.get("PLANTISMASH_TFBS_TIMEOUT", _DEFAULT_TIMEOUT_SEC))
        except Exception:
            timeout = _DEFAULT_TIMEOUT_SEC

    logging.warning("📥 Starting run_tfbs_finder for record: %s", rec_id)
    t0 = time.time()

    # Run underlying finder robustly
    try:
        results = _run_with_timeout(
            run_tfbs_finder,
            (record, options.tfbs_pvalue, options.tfbs_range),
            timeout,
        )
    except TimeoutError as te:
        logging.error("⏱️ TFBS timed out on %s: %s", rec_id, te)
        results = TFBSFinderResults()
    except Exception:
        logging.exception("💥 TFBS crashed on %s", rec_id)
        results = TFBSFinderResults()

    # Normalize result type just in case
    if not isinstance(results, TFBSFinderResults):
        try:
            results = TFBSFinderResults.from_json(results, record)
        except Exception:
            results = TFBSFinderResults()

    # Persist for HTML + later modules
    ns.extradata["TFBSFinderResults"] = results

    logging.warning(
        "🏁 run_tfbs_finder finished for record: %s in %.2fs",
        rec_id,
        time.time() - t0,
    )
    logging.warning("✅ TFBS finished successfully for record: %s", rec_id)


def run_on_record(record: SeqRecord, results: Optional[TFBSFinderResults], options) -> TFBSFinderResults:
    if results and getattr(results, "record_id", None) == record.id:
        return results
    # Reuse the same robust path
    try:
        return _run_with_timeout(run_tfbs_finder, (record, options.tfbs_pvalue, options.tfbs_range), 0)
    except Exception:
        logging.exception("TFBS run_on_record failed on %s; returning empty results", record.id)
        return TFBSFinderResults()


def regenerate_previous_results(previous: Dict[str, Any], record: SeqRecord, options) -> Optional[TFBSFinderResults]:
    results = TFBSFinderResults.from_json(previous, record)
    if not results:
        return None
    # Only reuse if current settings are at least as permissive/compatible
    if getattr(options, "tfbs_pvalue", 1.0) < getattr(results, "pvalue", 1.0):
        return None
    # If your TFBSFinderResults stores the scan half-range/window under another field,
    # change 'scan_range' accordingly.
    if getattr(options, "tfbs_range", None) != getattr(results, "scan_range", getattr(results, "window", None)):
        return None
    return results


# Output integration
def will_handle(product: str) -> bool:
    return True


def generate_details_div(*args, **kwargs):
    return output.generate_details_div(*args, **kwargs)