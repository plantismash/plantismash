# html_output.py

from pyquery import PyQuery as pq
from antismash import utils
from .specific_analysis import Result
from . import output_repeatfinder as output
import logging

def will_handle(product):
    if product.find('cyclopeptide') > -1:
        return True
    return False

def generate_details_div(cluster, seq_record, options, js_domains, details=None):
    logging.info("generating details div")
    """Generate details div"""
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

    # ➡️ ADD SIDEPANEL INSIDE DETAILS DIV
    sidepanel = pq('<div>')
    sidepanel.addClass('sidepanel')
    id_list = []
    for result in result_list:
        if result.cds_id:
            id_list.append(result.cds_id)
        else:
            id_list.append("Region from {} to {}".format(result.position[0], result.position[1]))
    sidepanel.html("{} coding sequences with repeats found:<br>{}".format(
        len(result_list), "<br>".join(id_list)))
    details.append(sidepanel)

    # ➡️ ADD THE DETAILED OUTPUT
    output_html = ""
    for r in result_list:
        output_html += output.write_result_summary(r)

    details.append(pq(output_html))

    return details


def generate_sidepanel(cluster, seq_record, options, sidepanel=None):
    logging.debug("generating sidepanel")
    cluster_feature = utils.get_cluster_by_nr(seq_record, cluster['idx'])
    sidepanel = pq('<div>')
    sidepanel.addClass('sidepanel')

    cluster_record = seq_record[cluster_feature.location.start:cluster_feature.location.end]
    result_list = gather_results(cluster_record)

    if len(result_list) > 0:
        id_list = []
        for result in result_list:
            if result.cds_id:
                id_list.append(result.cds_id)
            else:
                id_list.append("Region from {} to {}".format(result.position[0], result.position[1]))
        sidepanel.html("{} Coding sequences with repeats found:<br> {}".format(
            len(result_list), "<br>".join(id_list)))

    return sidepanel

def gather_results(cluster):
    result_list = [] 
    for feat in cluster.features:
        if 'cyclopeptide_analysis' in feat.qualifiers:
            result_list.append(Result(feat.qualifiers['cyclopeptide_analysis'][0]))
    return result_list  # 

