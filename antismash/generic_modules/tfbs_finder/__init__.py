# License: GNU Affero General Public License v3 or later
""" Transcription Factor Binding Site (TFBS) module for plantiSMASH """

from typing import List
from argparse import Namespace
import os
import logging

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


def _enabled(options) -> bool:
    # Only run when the CLI flag --tfbs-detection was set (stored in options.tfbs)
    return bool(getattr(options, "tfbs", False))


def _resolve_p(options) -> float:
    p = getattr(options, "tfbs_pvalue", None)
    if p is None or p <= 0:
        if _enabled(options):  # only warn if the feature is actually on
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


def check_prereqs(options) -> List[str]:
    """Only enforce files when TFBS is enabled."""
    if not _enabled(options):
        return []
    if not os.path.exists(PWM_PATH):
        return [f"Missing PWM JSON file at: {PWM_PATH}"]
    return []


def check_options(options) -> List[str]:
    """Validate TFBS-related options (only when enabled)."""
    if not _enabled(options):
        return []
    issues = []
    # We don't fail on None; we default later. Only flag explicitly bad values.
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