# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011-2013 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from pyquery import PyQuery as pq
from antismash import utils
from .specific_analysis import Result

def will_handle(product):
    if product.find('cyclopeptide') > -1:
        return True
    return False

def generate_details_div(cluster, seq_record, options, js_domains, details=None):
    """Generate details div"""
    result = None
    cluster = utils.get_cluster_by_nr(seq_record, cluster['idx']) # use seqrecord.feature
    if "cyclopeptide_analysis" in cluster.qualifiers:
        result = Result(cluster.qualifiers["cyclopeptide_analysis"][0])
    if result != None:
        # write visualization script for sidepanel here
        details = pq('<div>')
        details.text(result.dummy_attribute)

    return details

def generate_sidepanel(cluster, seq_record, options, sidepanel=None):
    """Generate sidepanel div"""
    result = None
    cluster = utils.get_cluster_by_nr(seq_record, cluster['idx']) # use seqrecord.feature
    if "cyclopeptide_analysis" in cluster.qualifiers:
        result = Result(cluster.qualifiers["cyclopeptide_analysis"][0])
    if result != None:
        # write visualization script for sidepanel here
        sidepanel = pq('<div>')
        sidepanel.text(result.dummy_attribute)
        
    return sidepanel