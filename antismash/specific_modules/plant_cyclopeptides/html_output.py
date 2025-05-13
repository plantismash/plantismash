#!/usr/bin/env python
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
# Copyright (C) 2024 Elena Del Pup 
# Wageningen University & Research, NL
# Bioinformatics Group, Department of Plant Sciences 
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from pyquery import PyQuery as pq
from antismash import utils
from .specific_analysis import Result
from . import output_repeatfinder as output
import logging

def will_handle(product):
    return 'cyclopeptide' in product

def generate_details_div(cluster, seq_record, options, js_domains, details=None):
    logging.info("generating details div")
    """Generate the details div for cyclopeptide clusters."""
    cluster_feature = utils.get_cluster_by_nr(seq_record, cluster['idx'])
    cluster_record = seq_record[cluster_feature.location.start:cluster_feature.location.end]

    details = pq('<div>')
    details.addClass('details')

    header = pq('<h3>')
    header.text('Repeatfinder output')
    details.append(header)

    result_list = gather_results(cluster_record)

    if not result_list:
        details.append(pq("<p>No repeats detected in this cluster.</p>"))
        return details

    output_html = ""
    for r in result_list:
        output_html += output.write_result_summary(r)

    details.append(pq(output_html))

    return details

def generate_sidepanel(cluster, seq_record, options, sidepanel=None):
    """No sidepanel output for cyclopeptide repeatfinder."""
    return None


def gather_results(cluster):
    """Collect cyclopeptide analysis results from cluster features."""
    result_list = []
    for feat in cluster.features:
        if 'cyclopeptide_analysis' in feat.qualifiers:
            result_list.append(Result(feat.qualifiers['cyclopeptide_analysis'][0]))
    return result_list
