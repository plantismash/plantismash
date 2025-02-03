# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Transcription Factor Binding Site (TFBS) Predictor for plantiSMASH
"""

import logging
from os import path
from typing import Any, Dict, List, Optional

from plantismash.common.secmet import Record
from plantismash.config import ConfigType
from plantismash.utils import locate_file

from .tfbs_finder import PWM_PATH, run_tfbs_finder, TFBSFinderResults

NAME = "tfbs_finder"
SHORT_DESCRIPTION = "Detects transcription factor binding sites (TFBSs)"

_required_files = [(PWM_PATH, False)]

def check_prereqs(_options: ConfigType) -> List[str]:
    """ Check if all required files exist """
    failure_messages = []
    if locate_file(PWM_PATH) is None:
        failure_messages.append(f"Failed to locate required file: {PWM_PATH}")
    return failure_messages

def is_enabled(options: ConfigType) -> bool:
    """ Returns True if the module is enabled """
    return options.tfbs

def check_options(options: ConfigType) -> List[str]:
    """ Validates TFBS module options """
    issues = []
    if options.tfbs_pvalue <= 0:
        issues.append(f"TFBS finder p-value must be positive: {options.tfbs_pvalue}")
    if options.tfbs_range <= 0:
        issues.append(f"TFBS finder range must be positive: {options.tfbs_range}")
    return issues

def run_on_record(record: Record, results: Optional[TFBSFinderResults], options: ConfigType
                  ) -> TFBSFinderResults:
    """ Runs TFBS prediction on the given record unless results are already available """
    if results and results.record_id == record.id:
        return results  # Reuse previous results
    logging.info(f"Running TFBS finder on record {record.id}")
    return run_tfbs_finder(record, options.tfbs_pvalue, options.tfbs_range)

def regenerate_previous_results(previous: Dict[str, Any], record: Record, options: ConfigType
                                ) -> Optional[TFBSFinderResults]:
    """ Attempts to reload previous results if valid """
    if not previous:
        return None
    results = TFBSFinderResults.from_json(previous, record)
    if not results:
        return None
    if options.tfbs_pvalue < results.pvalue:
        return None
    if options.tfbs_range != results.start_overlap:
        return None
    return results
