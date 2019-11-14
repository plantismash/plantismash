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
import output_repeatfinder as output
import logging

def will_handle(product):
    if product.find('cyclopeptide') > -1:
        logging.debug("Cyclopeptide found, returning TRUE")
        return True

    logging.debug("Cyclopeptide not found, returning FALSE")
    return False

def generate_details_div(cluster, seq_record, options, js_domains, details=None):
    logging.info("generating details div")
    """Generate details div"""
    result = None
    cluster = utils.get_cluster_by_nr(seq_record, cluster['idx']) # use seqrecord.feature
    text_for_dt = '<div>'
    
    details = pq('<div>')
    if "cyclopeptide_analysis" in cluster.qualifiers:
        logging.debug("cyclopeptide_analysis key in cluster")
        cluster_record = seq_record[cluster.start:cluster.end]
        result_list = gather_results(cluster_record) 
        sidepanel = pq('<div>')
        if len(result_list) > 0:
            logging.debug(" %s results in result list"%(len(result_list)))
            # write visualization script for sidepanel here
            output_html = ""
            for r in result_list:
                output_html.append(output.create_result_output(r))
            
            details.text(output_html)
            logging.debug("writing output")
        else:
            details.text("No repetition found")

    return details

def generate_sidepanel(cluster, seq_record, options, sidepanel=None):
    logging.debug("generating sidepanel")
    """Generate sidepanel div"""
    result_list = None
    cluster = utils.get_cluster_by_nr(seq_record, cluster['idx']) # use seqrecord.feature
    cluster_record = seq_record[cluster.location.start:cluster.location.end]
    result_list = gather_results(cluster_record) 
    sidepanel = pq('<div>')
    if len(result_list) > 0:
    
        # write visualization script for sidepanel here
        output_html = ""
        for r in result_list:
            output_html += output.create_result_output(r)
        
        sidepanel.text(output_html)
    else:
        sidepanel.text("No repetition found")
    return sidepanel

def gather_results(cluster):
    logging.info("type of \"cluster\": %s"%(type(cluster)))
    logging.info("%s"%(cluster))
    result_list = [] 
    for feat in cluster.features:
        logging.info("feat qualifiers:")
        logging.info(feat.qualifiers)
        if 'cyclopeptide_analysis' in feat.qualifiers:
            logging.info("Cyclopeptide analysis found in feature")
            result_list.append(Result(feat.qualifiers['cyclopeptide_analysis'][0]))
    logging.info(result_list)
    return result_list
