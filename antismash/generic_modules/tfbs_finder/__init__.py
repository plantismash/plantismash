# License: GNU Affero General Public License v3 or later

""" Transcription Factor Binding Site (TFBS) module for plantiSMASH """

from typing import Any, Dict, List, Optional
from argparse import Namespace
import os
import logging

from Bio.SeqRecord import SeqRecord
from antismash import utils

from .tfbs_detection import (
    PWM_PATH,
    run_tfbs_finder,
    TFBSFinderResults
)
from . import output  # For HTML rendering


NAME = "tfbs_finder"
SHORT_DESCRIPTION = "Detects transcription factor binding sites (TFBSs)"
REQUIRES = [(PWM_PATH, False)]


def check_prereqs(options) -> List[str]:
    """Check required input files exist (e.g., PWM JSON)."""
    if not os.path.exists(PWM_PATH):
        return [f"Missing PWM JSON file at: {PWM_PATH}"]
    return []


def check_options(options) -> List[str]:
    """Validate TFBS-related options."""
    issues = []
    if options.tfbs_pvalue <= 0:
        issues.append("TFBS p-value must be positive.")
    if options.tfbs_range <= 0:
        issues.append("TFBS scanning range must be positive.")
    return issues


def run_tfbs_finder_for_record(record: SeqRecord, options) -> None:
    """Cluster-level interface for TFBS analysis."""
    utils.log_status("➡️ Running TFBS analysis for contig #%d [p=%.1e, range=%dbp]" %
                     (options.record_idx, options.tfbs_pvalue, options.tfbs_range))
    results = run_tfbs_finder(record, options.tfbs_pvalue, options.tfbs_range)

    logging.warning("✅ TFBS finished successfully for record: %s", record.id)

    options.extrarecord.setdefault(record.id, Namespace())
    options.extrarecord[record.id].extradata = getattr(options.extrarecord[record.id], "extradata", {})
    options.extrarecord[record.id].extradata["TFBSFinderResults"] = results



def run_on_record(record: SeqRecord, results: Optional[TFBSFinderResults], options) -> TFBSFinderResults:
    """Run TFBS prediction unless reusing previous results."""
    if results and results.record_id == record.id:
        return results
    return run_tfbs_finder(record, options.tfbs_pvalue, options.tfbs_range)


def regenerate_previous_results(previous: Dict[str, Any], record: SeqRecord, options) -> Optional[TFBSFinderResults]:
    """Attempt to reuse cached results."""
    results = TFBSFinderResults.from_json(previous, record)
    if not results:
        return None
    if options.tfbs_pvalue < results.pvalue:
        return None
    if options.tfbs_range != results.start_overlap:
        return None
    return results


# Output integration (optional)
def will_handle(product: str) -> bool:
    return True


def generate_details_div(*args, **kwargs):
    return output.generate_details_div(*args, **kwargs)
