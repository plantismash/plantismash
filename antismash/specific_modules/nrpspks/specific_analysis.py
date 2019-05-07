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
'''
In-depth analysis and annotation of NRPS/PKS gene clusters.
'''

import os
#import logging
from os import path
from antismash import utils
from nrpspksdomains import annotate_pksnrps
from substrates import run_nrps_substr_spec_predictions, run_pks_substr_spec_predictions, parse_substr_spec_predictions
from parsers import calculate_consensus_prediction, generate_domainnamesdict, insert_modified_monomers
from orderfinder import analyse_biosynthetic_order
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re
import logging
import sys

def _map_domaintype(domain_type):
    "This function maps specific domain types inferred from the best matching hmm NAME field to the generic type name."

    #TODO: We should think about putting the mapping in a (auto generated?) config file...

    domain_type_dict = {'Condensation_DCL' : 'Condensation',
                       'Condensation_LCL' : 'Condensation',
                       'Condensation_Dual' : 'Condensation',
                       'Condensation_Starter' : 'Condensation',
                       'CXglyc' : 'Condensation',
                       'Cglyc' : 'Condensation',
                       'cMT' : 'MT',
                       'oMT' : 'MT',
                       'nMT' : 'MT',
                       'Polyketide_cyc' : 'Polyketide_cyc',
                       'Polyketide_cyc2' : 'Polyketide_cyc'
                       }
    if domain_type_dict.has_key(domain_type):
        return domain_type_dict[domain_type]

def create_rawoutput_storage_folder(options):
    #Make output folder for storing raw predictions
    options.raw_predictions_outputfolder = path.abspath(path.join(options.outputfoldername, "nrpspks_predictions_txt"))
    if not os.path.exists(options.raw_predictions_outputfolder):
        os.mkdir(options.raw_predictions_outputfolder)

def create_pksnrpsvars_object():
    #Storage object for all NRPS/PKS data
    pksnrpsvars = utils.Storage()
    #Dictionary, key: gene cluster nr, value: gene cluster type
    pksnrpsvars.nrpspkstypedict = {}
    #Dictionary, key: gene cluster nr, value: monomers string
    pksnrpsvars.compound_pred_dict = {}
    #Dictionary, key: gene ID, value: lists of result lists which each contain [result.hit_id, result.query_start, result.query_end, result.evalue, result.bitscore]
    pksnrpsvars.consensuspred_gene_dict = {}
    #Dictionary, key: gene ID, value: lists of result lists which each contain [result.hit_id, result.query_start, result.query_end, result.evalue, result.bitscore]
    pksnrpsvars.domaindict = {}
    #List of gene cluster nrs with failed structure generation
    pksnrpsvars.failedstructures = []
    #List of gene cluster nrs for which to create docking domain analysis details HTML files
    pksnrpsvars.dockingdomainanalysis = []
    #List of gene IDs of PKS/NRPS core genes
    pksnrpsvars.pksnrpscoregenes = []
    return pksnrpsvars

def specific_analysis(seq_record, options):
    pksnrpsvars = create_pksnrpsvars_object()
    create_rawoutput_storage_folder(options)
    pksnrpsvars = annotate_pksnrps(pksnrpsvars, seq_record, options)

    if len(pksnrpsvars.pksnrpscoregenes) > 0:
        run_nrps_substr_spec_predictions(pksnrpsvars, seq_record, options)
        run_pks_substr_spec_predictions(pksnrpsvars, seq_record, options)
        parse_substr_spec_predictions(pksnrpsvars, options)

        calculate_consensus_prediction(pksnrpsvars, seq_record)
        # print "pksnrpsvars.consensuspreds:"
        # print pksnrpsvars.consensuspreds

        generate_domainnamesdict(pksnrpsvars)
        # print "pksnrpsvars.domainnamesdict:"
        # print pksnrpsvars.domainnamesdict

        insert_modified_monomers(pksnrpsvars, seq_record)
        # print "pksnrpsvars.consensuspreds <- modified:"
        # print pksnrpsvars.consensuspreds

        analyse_biosynthetic_order(pksnrpsvars, seq_record, options)

#    generate_chemical_structure_preds(pksnrpsvars, seq_record, options)

    #write_substr_spec_preds_to_html(options, pksnrpsvars)
    write_data_to_seq_record(pksnrpsvars, seq_record, options)



def write_data_to_seq_record(pksnrpsvars, seq_record, options):
    #Save substrate specificity predictions in NRPS/PKS domain sec_met info of seq_record
    #
    # Workaround to extract positional information for CDS_motifs from the sec_met qualifiers

    for f in utils.get_cluster_features(seq_record):
	cluster_info = f.qualifiers

    for feature in pksnrpsvars.pksnrpscoregenes:
        nrat = 0
        nra = 0
        nrcal = 0
        nrkr = 0
        nrXdom = 0
        secmetqualifiers = feature.qualifiers['sec_met']
        updated_secmetqualifiers = []
        # BiosynML:creating object to add detailed substrate predictions
        updated_secmetqualifiers_predictions = []
        domainFeatures = []
        gene_id = utils.get_gene_id(feature)
        for qualifier in secmetqualifiers:
            if "NRPS/PKS Domain:" not in qualifier:
                updated_secmetqualifiers.append(qualifier)
                updated_secmetqualifiers_predictions.append(qualifier)
            else:
                # extract domain type, start and end position from qualifier string
                match_pos_obj = re.search("NRPS/PKS Domain: ([\w-]+) \((\d+)\-(\d+)\)\. E-value: ([\de\.-]+)\. Score: ([\de\.a-]+);", qualifier)
                if not match_pos_obj:
                    logging.exception("Exception: could not extract domain string from qualifier %s:" % qualifier)
                    sys.exit(1)
                domain_type = match_pos_obj.group(1)
                start_aa = int(match_pos_obj.group(2))
                end_aa = int(match_pos_obj.group(3))
                evalue = float(match_pos_obj.group(4))
                score = float (match_pos_obj.group(5))

                #calculate respective positions based on aa coordinates
                if feature.location.strand==1:
                    start = feature.location.start + ( 3 * start_aa )
                    end = feature.location.start + ( 3* end_aa )
                else:
                    end = feature.location.end - ( 3 * start_aa )
                    start = feature.location.end - ( 3 * end_aa)
                loc = FeatureLocation(start, end, strand=feature.strand)

                # set up new CDS_motif feature
                domainFeature = SeqFeature(loc, type=options.FeatureTags.pksnrpsdomains_tag)
                domainFeature.qualifiers['domain'] = [domain_type]
                if feature.qualifiers.has_key('locus_tag'):
                    domainFeature.qualifiers['locus_tag'] = feature.qualifiers['locus_tag']
                else:
                    domainFeature.qualifiers['locus_tag'] = [gene_id]
                domainFeature.qualifiers['detection'] = ["hmmscan"]
                domainFeature.qualifiers['database'] = ["nrpspksdomains.hmm"]
                domainFeature.qualifiers['evalue'] = [str("{:.2E}".format(float(evalue)))]
                domainFeature.qualifiers['score'] = [score]
                if feature.qualifiers.has_key('transl_table'):
                    [transl_table] = feature.qualifiers['transl_table']
                else:
                    transl_table = 1
                domainFeature.qualifiers['translation'] = [str(domainFeature.extract(seq_record).seq.translate(table=transl_table))]

                domainFeature_specificity = []

                if domain_type == "AMP-binding":
                    nra += 1
                    domainname = gene_id + "_A" + str(nra)
                    domainFeature.qualifiers['label'] = [domainname]
                    domainFeature.qualifiers['asDomain_id'] = ["nrpspksdomains_"+domainname]
                    domainFeature_specificity.append("NRPSpredictor2 SVM: %s" % pksnrpsvars.nrps_svm_preds[domainname])
                    domainFeature_specificity.append("Stachelhaus code: %s" % pksnrpsvars.nrps_code_preds[domainname])
                    domainFeature_specificity.append("Minowa: %s" % pksnrpsvars.minowa_nrps_preds[domainname])
                    domainFeature_specificity.append("consensus: %s" % pksnrpsvars.consensuspreds[domainname])


                    newqualifier = qualifier + " NRPS/PKS Domain: %s; Substrate specificity predictions: %s (NRPSPredictor2 SVM), %s (Stachelhaus code), %s (Minowa), %s (consensus);" % (domainname, pksnrpsvars.nrps_svm_preds[domainname], pksnrpsvars.nrps_code_preds[domainname], pksnrpsvars.minowa_nrps_preds[domainname], pksnrpsvars.consensuspreds[domainname])
                    # BiosynML: appending substrate prediction data into 'newqualifier_detailed'
                    newqualifier_detailed = qualifier + " NRPS/PKS Domain: %s; Substrate specificity predictions: %s (NRPSPredictor2 SVM), %s (Stachelhaus code), %s (Minowa), %s (consensus);" % (domainname,pksnrpsvars.nrps_code_preds_details[domainname], pksnrpsvars.nrps_svm_preds_details[domainname],  pksnrpsvars.minowa_nrps_preds_details[domainname], pksnrpsvars.consensuspreds[domainname])
                    updated_secmetqualifiers.append(newqualifier)
                    updated_secmetqualifiers_predictions.append(newqualifier_detailed)
                elif domain_type == "PKS_AT":
                    nrat += 1
                    domainname = gene_id + "_AT" + str(nrat)
                    domainFeature.qualifiers['label'] = [domainname]
                    domainFeature.qualifiers['asDomain_id'] = ["nrpspksdomains_"+domainname]
                    domainFeature_specificity.append("PKS signature: %s" % pksnrpsvars.pks_code_preds[domainname])
                    domainFeature_specificity.append("Minowa: %s" % pksnrpsvars.minowa_pks_preds[domainname])
                    #For t1pks, t2pks and t3pks
                    if 'transatpks' not in cluster_info['product'][0]:
                        domainFeature_specificity.append("consensus: %s" % pksnrpsvars.consensuspreds[domainname])
                        newqualifier = qualifier + " Substrate specificity predictions: %s (PKS signature), %s (Minowa), %s (consensus);" %(pksnrpsvars.pks_code_preds[domainname], pksnrpsvars.minowa_pks_preds[domainname], pksnrpsvars.consensuspreds[domainname])
                        # BiosynML: appending substrate prediction data into 'newqualifier_detailed'
                        newqualifier_detailed = qualifier + " Substrate specificity predictions: %s (PKS signature), %s (Minowa), %s (consensus);" %(pksnrpsvars.pks_code_preds_details[domainname], pksnrpsvars.minowa_pks_preds_details[domainname], pksnrpsvars.consensuspreds[domainname])
                        updated_secmetqualifiers.append(newqualifier)
                        updated_secmetqualifiers_predictions.append(newqualifier_detailed)
                    #For transatpks
                    elif 'transatpks' in cluster_info['product'][0]:
                        domainFeature_specificity.append("consensus: %s" % pksnrpsvars.consensuspreds_transat[domainname])
                        newqualifier = qualifier + " Substrate specificity predictions: %s (PKS signature), %s (Minowa), %s (consensus);" %(pksnrpsvars.pks_code_preds[domainname], pksnrpsvars.minowa_pks_preds[domainname], pksnrpsvars.consensuspreds_transat[domainname])
                        # BiosynML: appending substrate prediction data into 'newqualifier_detailed'
                        newqualifier_detailed = qualifier + " Substrate specificity predictions: %s (PKS signature), %s (Minowa), %s (consensus);" %(pksnrpsvars.pks_code_preds_details[domainname], pksnrpsvars.minowa_pks_preds_details[domainname], pksnrpsvars.consensuspreds_transat[domainname])

                        updated_secmetqualifiers.append(newqualifier)
                        updated_secmetqualifiers_predictions.append(newqualifier_detailed)
                elif domain_type == "CAL_domain":
                    nrcal += 1
                    domainname = gene_id + "_CAL" + str(nrcal)
                    domainFeature.qualifiers['label'] = [domainname]
                    domainFeature.qualifiers['asDomain_id'] = ["nrpspksdomains_"+domainname]
                    domainFeature_specificity.append("Minowa: %s" % pksnrpsvars.minowa_cal_preds[domainname])
                    newqualifier = qualifier + " Substrate specificity predictions: %s (Minowa);" %(pksnrpsvars.minowa_cal_preds[domainname])
                    # BiosynML: appending substrate prediction data into 'newqualifier_detailed'
                    newqualifier_detailed = qualifier + " Substrate specificity predictions: %s (Minowa);" %(pksnrpsvars.minowa_cal_preds_details[domainname])
                    updated_secmetqualifiers.append(newqualifier)
                    updated_secmetqualifiers_predictions.append(newqualifier_detailed)
                elif domain_type == "PKS_KR":
                    nrkr += 1
                    domainname = gene_id + "_KR" + str(nrkr)
                    domainFeature.qualifiers['label'] = [domainname]
                    domainFeature.qualifiers['asDomain_id'] = ["nrpspksdomains_"+domainname]
                    domainFeature_specificity.append("KR activity: %s" % pksnrpsvars.kr_activity_preds[domainname])
                    domainFeature_specificity.append("KR stereochemistry: %s" % pksnrpsvars.kr_stereo_preds[domainname])
                    newqualifier = qualifier + " Predicted KR activity: %s; Predicted KR stereochemistry: %s;" %(pksnrpsvars.kr_activity_preds[domainname], pksnrpsvars.kr_stereo_preds[domainname])
                    # BiosynML: appending substrate prediction data into 'newqualifier_detailed'
                    newqualifier_detailed = qualifier + " Predicted KR activity: %s; Predicted KR stereochemistry: %s;" %(pksnrpsvars.kr_activity_preds[domainname], pksnrpsvars.kr_stereo_preds[domainname])
                    updated_secmetqualifiers.append(newqualifier)
                    updated_secmetqualifiers_predictions.append(newqualifier_detailed)
                else:
                    nrXdom += 1
                    domainFeature.qualifiers['asDomain_id'] = ["nrpspksdomains_" + gene_id.partition(".")[0] + "_Xdom"+'{:02d}'.format(nrXdom)]
                    updated_secmetqualifiers.append(qualifier)
                domainFeature.qualifiers['specificity'] = domainFeature_specificity
                if _map_domaintype(domain_type):
                    domainFeature.qualifiers['domain_subtype'] = [domain_type]
                    domainFeature.qualifiers['domain'] = [_map_domaintype(domain_type)]
                domainFeatures.append(domainFeature)

        feature.qualifiers['sec_met'] = updated_secmetqualifiers
        # BiosynML: creating new 'sec_met_predictions' qualifier
        #feature.qualifiers['sec_met_predictions'] = updated_secmetqualifiers_predictions
        seq_record.features.extend(domainFeatures)

        if pksnrpsvars.consensuspred_gene_dict.has_key(gene_id):
            feature.qualifiers[options.QualifierTags.product_prediction] = "-".join(pksnrpsvars.consensuspred_gene_dict[gene_id])

    #Save consensus structure + link to structure image to seq_record
    clusters = utils.get_cluster_features(seq_record)
    for cluster in clusters:
        clusternr = utils.get_cluster_number(cluster)
        if pksnrpsvars.compound_pred_dict.has_key(clusternr):
            structpred = pksnrpsvars.compound_pred_dict[clusternr]
            cluster.qualifiers['note'].append("Monomers prediction: " + structpred)
            cluster.qualifiers['note'].append("Structure image: structures/genecluster%s.png" % clusternr)

            # as a demo: save pksnrpsvars storage object in cluster record
            #cluster.qualifiers['note'].append(pksnrpsvars.Storage_to_NoteQualifier('pksnrpsvars'))
