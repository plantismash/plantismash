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

import logging
from os import path
from pyquery import PyQuery as pq
from antismash import utils
import re

def will_handle(product):
    if product.find('nrps') > -1:
        return True
    if product.find('t1pks') > -1:
        return True
    if product.find('transatpks') > -1:
        return True

    return False

def generate_details_div(cluster, seq_record, options, js_domains, details=None):
    """Generate details div"""

    cluster_rec = utils.get_cluster_by_nr(seq_record, cluster['idx'])
    if cluster_rec is None:
        return details

    if details is None:
        details = pq('<div>')
        details.addClass('details')

        header = pq('<h3>')
        header.text('Detailed annotation')
        details.append(header)

    js_cluster_domains = { 'id': "cluster-%s-details" % cluster['idx'], 'orfs': [] }
    features = utils.get_cluster_cds_features(cluster_rec, seq_record)
    for feature in features:
        if not 'sec_met' in feature.qualifiers:
            continue

        if 'translation' in feature.qualifiers:
            sequence = feature.qualifiers['translation'][0]
        else:
            sequence = str(utils.get_aa_sequence(feature))

        js_orf = {
            'id': utils.get_gene_id(feature),
            'sequence': sequence,
            'domains': [],
        }

        for qual in feature.qualifiers['sec_met']:
            if not qual.startswith('NRPS/PKS Domain:'):
                continue

            js_domain = _parse_domain(qual, feature, seq_record)
            if len(js_domain) > 0:
                js_orf['domains'].append(js_domain)

        if len(js_orf['domains']) > 0:
            js_cluster_domains['orfs'].append(js_orf)

    if len(js_cluster_domains['orfs']) > 0:
        details_svg = pq('<div>')
        details_svg.addClass('details-svg')
        details_svg.attr('id', '%s-svg' % js_cluster_domains['id'])
        details.append(details_svg)

        js_domains.append(js_cluster_domains)

    return details


def generate_sidepanel(cluster, seq_record, options, sidepanel=None):
    """Generate sidepanel div"""
    cluster_rec = utils.get_cluster_by_nr(seq_record, cluster['idx'])
    if cluster_rec is None:
        return sidepanel

    if sidepanel is None:
        sidepanel = pq('<div>')
        sidepanel.addClass('sidepanel')

    structure = pq('<div>')
    structure.addClass('structure')
    structure_header = pq('<h3>')
    structure_header.text('Predicted core structure')
    structure.append(structure_header)
    a = pq('<a>')
    a.attr('href', _get_structure_image_url(cluster_rec, options.outputfoldername))
    a.attr('target', '_new')
    structure.append(a)
    structure_img = pq('<img>')
    structure_img.attr('src', _get_structure_image_url(cluster_rec, options.outputfoldername))
    a.append(structure_img)
    warning = pq('<div>')
    warning.addClass('as-structure-warning')
    if not 'docking' in options:
        options.docking = {}
    if cluster['idx'] in options.docking and options.docking[cluster['idx']]:
        warning.text('Rough prediction of core scaffold based on assumed '
                    'PKS linker matching; tailoring reactions not taken '
                    'into account')
    else:
        warning.text('Rough prediction of core scaffold based on assumed '
                    'PKS/NRPS colinearity; tailoring reactions not taken '
                    'into account')
    structure.append(warning)
    sidepanel.append(structure)

    details = pq('<div>')
    details.addClass('more-details')
    details_header = pq('<h3>')
    details_header.text('Prediction details')
    details.append(details_header)
    details_list = pq('<dl>')
    details_list.addClass('prediction-text')

    details.append(details_list)
    sidepanel.append(details)
    dt = pq('<dt>')
    dt.text('Monomers prediction:')
    details_list.append(dt)
    dd = pq('<dd>')
    dd.text(_get_monomer_prediction(cluster_rec))
    details_list.append(dd)

    features = utils.get_cluster_cds_features(cluster_rec, seq_record)
    for feature in features:
        if not 'sec_met' in feature.qualifiers:
            continue

        header_printed = False
        per_CDS_predictions = []
        for qual in feature.qualifiers['sec_met']:
            if not qual.startswith('NRPS/PKS Domain:'):
                continue
            # logging.debug("qual: %s" % qual)
            preds = _parse_substrate_predictions(qual)
            
            per_Adomain_predictions = []
            for key, val in preds:
                
                
                if not header_printed:
                    dt = pq('<dt>')
                    dt.text(utils.get_gene_id(feature))
                    details_list.append(dt)
                    header_printed = True
                dd = pq('<dd>')
                dd.html('%s: %s<br>' % (key, val))
                details_list.append(dd)
                if qual.startswith("NRPS/PKS Domain: AMP-binding"):
                    values = _filter_norine_as(val.split(","))
                    if len(values) > 0:
                        per_Adomain_predictions.extend(val.split(","))
            
            if len(preds) > 0:
                if qual.startswith("NRPS/PKS Domain: AMP-binding"):
                    per_Adomains_predictions_unique = list(set(per_Adomain_predictions))
                    per_CDS_predictions.append(per_Adomains_predictions_unique)
                # logging.debug("substrate prediction list: %s" % ",".join(per_Adomains_predictions_unique) )
                dd = pq('<dd>')
                dd.append(pq('<br>'))
                details_list.append(dd)
                
        if len(per_CDS_predictions) > 0:
            url = _get_norine_url_for_specArray(per_CDS_predictions)
            if url:
                dd = pq('<dd>')
                dd.append("Search NORINE for peptide in ")
                a = pq('<a>')
                a.attr('href', url)
                a.attr('target','_new')
                a.text("strict mode")
                dd.append(a)
                dd.append(" // ")
                url = _get_norine_url_for_specArray(per_CDS_predictions, be_strict=False)
                a = pq('<a>')
                a.attr('href', url)
                a.attr('target','_new')
                a.text("relaxed mode")
                dd.append(a)
                dd.append(pq('<br>'))
                dd.append(pq('<br>'))
                details_list.append(dd)
                

    if cluster['type'].find('nrps') > -1:
        cross_refs = pq("<div>")
        refs_header = pq('<h3>')
        refs_header.text('Database cross-links')
        cross_refs.append(refs_header)
        links = pq("<div>")
        links.addClass('prediction-text')

        a = pq("<a>")
        a.attr('href', 'http://bioinfo.lifl.fr/norine/form2.jsp')
        a.attr('target', '_new')
        a.text("Link to NORINE database query form")
        links.append(a)
        links.append("<br>")
        
        a = pq("<a>")
        url = _get_norine_url_for_cluster(cluster_rec)
        logging.debug("NORINE URL string: %s" % url)
        a.attr('href', url)
        a.attr('target', '_new')
        a.text("strict mode")
        links.append("Direct lookup in NORINE database in ")
        links.append(a)
        links.append(" // ")
        url = _get_norine_url_for_cluster(cluster_rec, be_strict=False)
        a = pq("<a>")
        a.attr('href', url)
        a.attr('target', '_new')
        a.text("relaxed mode")
        links.append(a)
        cross_refs.append(links)
        sidepanel.append(cross_refs)

    return sidepanel


def _parse_substrate_predictions(domain):
    "Parse the substrate predictions from the NRPS/PKS domain string"
    predictions = []


    idx = domain.find('Substrate specificity predictions:')
    if idx == -1:
        return []

    specifities = domain[idx + 35:].strip(';').split(', ')

    for spec in specifities:
        try:
            spec = spec.strip()
            substrate, name = spec.split(' ', 1)
            name = name.strip('()')
            predictions.append((name, substrate))
        except ValueError:
            logging.debug(domain)
            logging.debug(spec)
            raise

    return predictions

def _get_structure_image_url(feature, outputfolder):
    "Get the relative url to the structure image"
    for note in feature.qualifiers.get('note', []):
        if not note.startswith("Structure image:"):
            continue

        url = note.split()[-1]
        if not path.exists(path.join(outputfolder, url)):
            logging.debug('No file at %r' % url)
            continue
        return url
    return 'images/nostructure_icon.png'

def _get_monomer_prediction(feature):
    "Get the monomer prediction of the cluster"
    for note in feature.qualifiers.get('note', []):
        if not note.startswith("Monomers prediction:"):
            continue

        monomers = note.split(':')[-1].strip()
        return monomers
    return 'N/A'

def _map_as_names_to_norine(as_name):
    """ Just a dictionary helper function to map antiSMASH amino acid nomenclature to NORINE"""

    as_replacement_dict = {'bht' : 'bOH-Tyr',
                           'dhb' : 'diOH-Bz',
                           'iva' : 'Ival',
                           'pip' : 'Hpr',
                           'sal' : 'diOH-Bz',
                           'nrp' : 'X'}
    if as_replacement_dict.has_key(as_name):
        return as_replacement_dict[as_name].lower()
    else:
        return as_name

def _filter_norine_as(as_list, be_strict=False):
    """ Remove PKS and unknow substrate predictions
        use be_strict = False to also filter nrp/X
    """
    
    filtered_list=[]
    for as_monomer in as_list:
        if not as_monomer in ['pk', 'N/A', 'hydrophilic', 'hydrophobic', 'mal', 'mmal']:
            if be_strict and as_monomer in ['nrp','X']:
                continue
            else:
                filtered_list.append(as_monomer)
    return filtered_list

def _get_norine_url_for_cluster(cluster_rec, be_strict=True):
    """ Get a NORINE URL string for direct querying
        use be_strict=False to add * after each monomer"""
    
    
    monomer_string = _get_monomer_prediction(cluster_rec)
    monomers_per_protein_list = re.findall("\(.*?\)", monomer_string)
    i = 1
    nrpslist = []
    for monomers_per_protein in monomers_per_protein_list:
        monomers = monomers_per_protein[1:-1].split("-")

        if be_strict:
            monomers = [_map_as_names_to_norine(element.lower()) for element in _filter_norine_as(monomers, be_strict=True)]
        else:
            monomers = [_map_as_names_to_norine(element.lower())+"*" for element in _filter_norine_as(monomers)]
        # Norine doesn't allow "x*" as a "relaxed" monomer, so we have to replace this with "x"
        monomers = map(lambda x: "x" if x=="x*" else x, monomers)
        
        if len(monomers) > 0:
            nrpslist.append("nrps" + str(i) + "=" + ",".join(monomers))
        i += 1
    urlstring = "http://bioinfo.lifl.fr/norine/fingerPrintSearch.jsp?"+"&".join(nrpslist)
    logging.debug("NORINE URL from array: %s" % urlstring)
    return urlstring

def _get_norine_url_for_specArray(specArray, be_strict=True):
    """ generate NORINE URL string for direct querying from array of specificity predictions
        use be_strict=False to add * after each monomer"""
    modulelist = []
    if len(specArray)> 0:
        for domain_specificity_list in specArray:
            if len(domain_specificity_list) == 1:
                if be_strict:
                    modulelist.append(_map_as_names_to_norine(domain_specificity_list[0]))
                else:
                    modulelist.append(_map_as_names_to_norine(domain_specificity_list[0])+"*")
            elif len(domain_specificity_list) > 1:
                if be_strict:
                    modulelist.append("["+"|".join([_map_as_names_to_norine(element) for element in _filter_norine_as(domain_specificity_list, be_strict=True)])+"]")
                else:
                    # we have to use be_strict to remove X from the list of predictions, as otherwise consenus: nrp will always match
                    monomers = [_map_as_names_to_norine(element)+"*" for element in _filter_norine_as(domain_specificity_list, be_strict=True)]
                    # Norine doesn't allow "x*" as a "relaxed" monomer, so we have to replace this with "x"
                    monomers = map(lambda x: "x" if x=="x*" else x, monomers)
                    modulelist.append("["+"|".join(monomers)+"]")
            else:
                logging.warn("retrieved emtpy domain list for assembling per protein NORINE link string- this should not happen!!")
        urlstring = "http://bioinfo.lifl.fr/norine/fingerPrintSearch.jsp?nrps1="+",".join(modulelist)
        logging.debug("NORINE URL from array: %s" % urlstring)
        return urlstring
                

def _parse_domain(domain, feature, seq_record):
    "Convert a NRPS/PKS domain string to a dict useable by json.dumps"
    text = domain[17:]
    type_, location, prediction_string = text.split(' ', 2)
    predictions = _parse_substrate_predictions(prediction_string)

    location = location.strip('().')
    coordinates = location.split('-')

    #Create url_link to NaPDoS for C and KS domains
    napdoslink = ""
    domainseq = str(utils.get_aa_sequence(feature))[int(coordinates[0]):int(coordinates[-1])]
    if "PKS_KS" in text:
        napdoslink = "http://napdos.ucsd.edu/cgi-bin/process_request.cgi?query_type=aa&amp;ref_seq_file=all_KS_public_12062011.faa&amp;Sequence=%3EKS_domain_from_antiSMASH%0D" + domainseq
    elif "Condensation" in text:
        napdoslink = "http://napdos.ucsd.edu/cgi-bin/process_request.cgi?query_type=aa&amp;ref_seq_file=all_C_public_12062011.faa&amp;Sequence=%3EC_domain_from_antiSMASH%0D" + domainseq
    blastlink = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + domainseq + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"

    try:
        js_domain = { 'type': type_, 'start': int(coordinates[0]), 'end': int(coordinates[1]),
                      'predictions': predictions, 'napdoslink': napdoslink, 'blastlink': blastlink,
                      'sequence': domainseq }
        return js_domain
    except ValueError:
        logging.debug('%r' % text)
        logging.debug('%r  %r' % (type_, location))
        logging.debug(coordinates)
        raise
