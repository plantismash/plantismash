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
        return True

    return False

def generate_details_div(cluster, seq_record, options, js_domains, details=None):
    logging.info("generating details div")
    """Generate details div"""
    cluster = utils.get_cluster_by_nr(seq_record, cluster['idx']) # use seqrecord.feature
    details = pq('<div>')
    details.addClass('details')
    header = pq('<h3>')
    header.text('Repeatfinder output')
    cluster_record = seq_record[cluster.location.start:cluster.location.end]
    result_list = gather_results(cluster_record) 
    sidepanel = pq('<div>')
    if len(result_list) > 0:
        # write visualization script for sidepanel here
        output_html = ""
        for r in result_list:
            output_html += output.create_result_output(r)
        
        details.html(output_html)

    return details

def generate_sidepanel(cluster, seq_record, options, sidepanel=None):
    logging.debug("generating sidepanel")
    """Generate sidepanel div"""
    result_list = None
    cluster = utils.get_cluster_by_nr(seq_record, cluster['idx']) # use seqrecord.feature
    cluster_record = seq_record[cluster.location.start:cluster.location.end]
    result_list = gather_results(cluster_record) 
    sidepanel = pq('<div>')#TODO add class and put it in the details div class
    sidepanel.addClass('sidepanel')
    if len(result_list) > 0:
    
        # write visualization script for sidepanel here
        #output_html = ""
        #for r in result_list:
        #    output_html += output.create_result_output(r)
        
        #sidepanel.html(output_html)
        id_list = []
        for result in result_list:

            if result.cds_id:
                id_list.append(result.cds_id)

            else:
                id_list.append("Region with unknown ID from %s to %s"%(result.position[0],result.position[1]))
                
        sidepanel.html("%s Coding sequences with repeats found:<br> %s"%(len(result_list),"<br>".join(id_list)))
                
        
    else:
        sidepanel.text("No repetition found")
#SIDEPANEL DISABLED
def gather_results(cluster):
    result_list = [] 
    for feat in cluster.features:
        if 'cyclopeptide_analysis' in feat.qualifiers:
            result_list.append(Result(feat.qualifiers['cyclopeptide_analysis'][0]))
        
    return set(result_list)
