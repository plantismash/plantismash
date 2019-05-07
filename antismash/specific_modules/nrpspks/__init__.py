# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011,2012 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""NRPS/PKS analysis module

"""
import os
import logging
from os import path
from antismash import utils
from .specific_analysis import specific_analysis
from .html_output import will_handle, generate_sidepanel, generate_details_div

name = "nrpspks"

short_description = name.capitalize()

# The tuple is the name of the binary and whether it is an optional requirement
_required_binaries = [
    ('hmmsearch', False),
    ('hmmpress', False),
    ('muscle', False),
    ('java', False),
]

_markov_models = [
    'abmotifs.hmm',
    'dockingdomains.hmm',
    'ksdomains.hmm',
    'nrpspksdomains.hmm'
]

_binary_extensions = [
    '.h3f',
    '.h3i',
    '.h3m',
    '.h3p'
]

def check_prereqs():
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)

    for hmm in _markov_models:
        hmm = utils.get_full_path(__file__, hmm)
        if utils.locate_file(hmm) is None:
            failure_messages.append("Failed to locate file %r" % hmm)
            continue
        for ext in _binary_extensions:
            binary = "%s%s" % (hmm, ext)
            if utils.locate_file(binary) is None:
                _, err, retcode = utils.run_hmmpress(hmm)
                if retcode != 0:
                    failure_messages.append("Failed to hmmpress %r: %r" % (hmm, err))
                break
            else:
                binary_mtime = path.getmtime(binary)
                hmm_mtime = path.getmtime(hmm)
                if hmm_mtime > binary_mtime:
                    try:
                        from glob import glob
                        for f in glob("%s.h3?" % hmm):
                            logging.debug("removing outdated file %s", f)
                            os.remove(f)
                    except OSError as e:
                        failure_messages.append("Failed to remove outdated binary file for %s: %s" % \
                            (hmm, e))
                        break
                    _, err, retcode = utils.run_hmmpress(hmm)
                    if retcode != 0:
                        failure_messages.append("Failed to hmmpress %r: %r" % (hmm, err))
                        import datetime
                        failure_messages.append("HMM binary files outdated. %s (changed: %s) vs %s (changed: %s)" % \
                            (hmm, datetime.datetime.fromtimestamp(hmm_mtime),
                             binary, datetime.datetime.fromtimestamp(binary_mtime)))
                    break



    return failure_messages

def _update_sec_met_entry(feature, results, clustertype):
    result = "; ".join([res.query_id + " (" + str(res.bitscore) + ")" for res in results])

    if not 'sec_met' in feature.qualifiers:
        feature.qualifiers['sec_met'] = [
            "Type: %s" % clustertype,
            "Domains detected: %s" % (result),
            "Kind: biosynthetic"
        ]
    else:
        for ann in feature.qualifiers['sec_met']:
            if not ann.startswith("Domains detected"):
                continue
            ann += "Domains detected: %s" % (result)

__all__ = [ check_prereqs, specific_analysis, will_handle, generate_sidepanel,
            generate_details_div ]
