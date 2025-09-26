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

import sys
from antismash import utils
from string import ascii_letters#
import logging
import re
from Bio.SeqFeature import SeqFeature, FeatureLocation

def prepare_visualization(options, seq_record):
    
    # Check, whether (Sub)ClusterBlast data is encoded in source feature
    sourceFeatures = utils.get_all_features_of_type(seq_record, 'source')
    if len(sourceFeatures) == 0:
        loc = FeatureLocation(0, len(seq_record.seq))
        source_feature = SeqFeature(loc, type="source")
        seq_record.features.append(source_feature)
        sourceFeatures = utils.get_all_features_of_type(seq_record, 'source')
    
    if 'extrarecord' in options:
        if seq_record.id in options.extrarecord:
            # As there is only one source feature per record we just can take the first one without cycling through all features
            for key in list(options.extrarecord[seq_record.id].extradata.keys()):
                if key == 'ClusterBlastData':
                    logging.debug("prepare_visualization: Found ClusterBlastData storage object")
                    options.clusterblast = True
                    
                    clusterBlastResults = options.extrarecord[seq_record.id].extradata[key]
                                                                                     
                    seq_record.internalhomologygroupsdict = clusterBlastResults.internalhomologygroupsdict
                    seq_record.known_compound_dict = clusterBlastResults.known_compound_dict
                    seq_record.nrhitgeneclusters = clusterBlastResults.nrhitgeneclusters
                    seq_record.qgeneclusterdata = clusterBlastResults.qgeneclusterdata
                    seq_record.queryclusterdata = clusterBlastResults.queryclusterdata
                    seq_record.pubchem_dict = clusterBlastResults.pubchem_dict
                    seq_record.pubmed_dict = clusterBlastResults.pubmed_dict
        
                elif key == 'SubClusterBlastData':
                    logging.debug("prepare_visualization: Found SubClusterBlastData storage object")
                    options.subclusterblast = True
                   
                    subclusterBlastResults = options.extrarecord[seq_record.id].extradata[key]
                    seq_record.internalhomologygroupsdict = subclusterBlastResults.internalhomologygroupsdict
                    seq_record.sc_nrhitgeneclusters = subclusterBlastResults.sc_nrhitgeneclusters
            #        seq_record.sc_qgeneclusterdata = subclusterBlastResults.sc_qgeneclusterdata
                    seq_record.sc_queryclusterdata = subclusterBlastResults.sc_queryclusterdata
                    seq_record.pubchem_dict = subclusterBlastResults.pubchem_dict
                    seq_record.pubmed_dict = subclusterBlastResults.pubmed_dict
                    
                elif key == 'KnownClusterBlastData':
                    logging.debug("prepare_visualization: Found KnownClusterBlastData storage object")
                    options.knownclusterblast = True
                    
                    knownclusterBlastResults = options.extrarecord[seq_record.id].extradata[key]
                    seq_record.internalhomologygroupsdict = knownclusterBlastResults.internalhomologygroupsdict
                    seq_record.kc_nrhitgeneclusters = knownclusterBlastResults.kc_nrhitgeneclusters
            #        seq_record.kc_qgeneclusterdata = knownclusterBlastResults.sc_qgeneclusterdata
                    seq_record.kc_queryclusterdata = knownclusterBlastResults.kc_queryclusterdata
                    seq_record.pubchem_dict = knownclusterBlastResults.pubchem_dict
                    seq_record.pubmed_dict = knownclusterBlastResults.pubmed_dict
#                 elif key == 'MetabolicModelDataObj':
#                     pass
#                 else:
#                     logging.warn('Found key %s in options.clusterblastdata which does not match the hard coded choices!' % key)
#     # 
    #     load_pubmed_pubchem_links(seq_record)
    #     if options.clusterblast:
    #         load_clusterblast_outputdata(seq_record, options)
    #     if options.subclusterblast:
    #         load_subclusterblast_outputdata(seq_record, options)
    #     if options.knownclusterblast:
    #         load_knownclusterblast_outputdata(seq_record, options)
    #     load_genecluster_info(seq_record, options)
    
        load_genecluster_info(seq_record, options)




def construct_colorgroups(colorgroupsdict, clusternr, blasthitdict, blastdetailsdict, internalhomologygroupsdict, hitclusternr):
    #Make groups of genes for coloring
    colorgroups = []
    internalgroups = internalhomologygroupsdict[clusternr]
    for i in internalgroups:
        querygenes_and_hits = []
        for j in i:
            #Make list of query gene and its hits
            additionalhits = []
            #For each hit, check if it was also hit by another gene; if so, only add it to the group if this hit had the lowest blast score
            queryscore = 0
            if j in blasthitdict:
                for k in blasthitdict[j]:
                    otherscores = []
                    for l in list(blastdetailsdict.keys()):
                        if j == l.partition("_|_|_")[0] and k == l.rpartition("_|_|_")[2]:
                            queryscore = blastdetailsdict[l][1]
                        if k in l and j not in l:
                            otherscores.append(blastdetailsdict[l][1])
                    allscores = otherscores + [queryscore]
                    if int(queryscore) == max([int(m) for m in allscores]):
                        additionalhits.append(k)
                #Add additional hits to the querygenes_and_hits list that will form a colorgroup
                querygenes_and_hits = querygenes_and_hits + additionalhits
                if j not in querygenes_and_hits:
                    querygenes_and_hits.append(j)
        if len(querygenes_and_hits) > 0:
            colorgroups.append(querygenes_and_hits)
    colorgroupsdict[hitclusternr] = colorgroups
    return colorgroupsdict



def test_accession(accession):
    #Test if accession number is probably real GenBank/RefSeq acc nr
    numbers = list(range(0,10))
    letters = [i for i in ascii_letters]
    nrletters = len([i for i in accession if i in ascii_letters])
    nrnumbers = len([i for i in accession if i.isdigit()])
    if "." in accession:
        period_index = accession.index(".")
    else:
        period_index = 10
    if nrnumbers < 3 or nrletters < 1 or len(accession) > 16 or period_index < 4:
        return False
    else:
        return True



def load_genecluster_info(seq_record, options):
    #Gather and store data on each gene cluster
    smcogdict, smcogdescriptions = utils.get_smcog_annotations(seq_record)
    gtrcoglist = ['SMCOG1045','SMCOG1062','SMCOG1102']
    transportercoglist = ['SMCOG1000','SMCOG1005','SMCOG1011','SMCOG1020','SMCOG1029','SMCOG1033','SMCOG1035','SMCOG1044','SMCOG1065','SMCOG1067','SMCOG1069','SMCOG1074','SMCOG1085','SMCOG1096','SMCOG1106','SMCOG1118','SMCOG1131','SMCOG1166','SMCOG1169','SMCOG1184','SMCOG1202','SMCOG1205','SMCOG1214','SMCOG1234','SMCOG1243','SMCOG1245','SMCOG1252','SMCOG1254','SMCOG1288']
    seq_record.qgeneclusterdata = {}
    geneclusters = utils.get_cluster_features(seq_record)
    for genecluster in geneclusters:
        geneclusternr = utils.get_cluster_number(genecluster)
        clustergenes, clustertype, annotations, colors, starts, ends, strands, pksnrpsprots, gtrs, transporters, clustersize = retrieve_gene_cluster_annotations(seq_record, smcogdict, gtrcoglist, transportercoglist, geneclusternr)
        if options.clusterblast:
            hitgeneclusterdata = retrieve_clusterblast_info(seq_record, geneclusternr)
        else:
            hitgeneclusterdata = {}
        pksnrpsprotsnames, pksnrpsdomains, domlist, domsdetails, substrspecnrpspredictordict, substrspecminowadict, substrspecpkssigdict, substrspecconsensusdict, krpredictionsdict, structpred = retrieve_pksnrps_info(seq_record, geneclusternr, pksnrpsprots)
        seq_record.qgeneclusterdata[geneclusternr] = [clustertype, clustersize, clustergenes, annotations, starts, ends, strands, pksnrpsprots, pksnrpsprotsnames, pksnrpsdomains, substrspecnrpspredictordict, substrspecminowadict, substrspecpkssigdict, substrspecconsensusdict, gtrs, transporters, colors, hitgeneclusterdata, structpred, krpredictionsdict]

def retrieve_gene_cluster_annotations(seq_record, smcogdict, gtrcoglist, transportercoglist, geneclusternr):
    allcoregenes = [utils.get_gene_id(cds) for cds in utils.get_secmet_cds_features(seq_record)]
    pksnrpscoregenes = [utils.get_gene_id(cds) for cds in utils.get_pksnrps_cds_features(seq_record)]
    feature_by_id = utils.get_feature_dict(seq_record)
    clustergenes = [utils.get_gene_id(cds) for cds in utils.get_cluster_cds_features(utils.get_cluster_by_nr(seq_record, geneclusternr), seq_record)]
    clustertype = utils.get_cluster_type(utils.get_cluster_by_nr(seq_record, geneclusternr))
    annotations = {}
    colors = []
    starts = []
    ends = []
    strands = []
    pksnrpsprots = []
    gtrs = []
    transporters = []
    for j in clustergenes:
        cdsfeature = feature_by_id[j]
        if 'product' in cdsfeature.qualifiers:
            annotations[j] = cdsfeature.qualifiers['product'][0]
        else:
            annotations[j] = 'Unannotated gene'
        starts.append(cdsfeature.location.start)
        ends.append(cdsfeature.location.end)
        if cdsfeature.strand == -1:
            strands.append("-")
        else:
            strands.append("+")
        if j in allcoregenes:
            colors.append("#810E15")
        else:
            colors.append("grey")
        if j in pksnrpscoregenes:
            pksnrpsprots.append(j)
        if j in smcogdict:
            if len(smcogdict[j]) > 0 and smcogdict[j][0] in gtrcoglist:
                gtrs.append(j)
            if len(smcogdict[j]) > 0 and smcogdict[j][0] in transportercoglist:
                transporters.append(j)
    clustersize = abs(max(starts+ends) - min(starts+ends))
    return clustergenes, clustertype, annotations, colors, starts, ends, strands, pksnrpsprots, gtrs, transporters, clustersize

def retrieve_clusterblast_info(seq_record, geneclusternr):
    hitgeneclusters = list(range(1,(seq_record.nrhitgeneclusters[geneclusternr] + 1)))
    hitgeneclusterdata = {}
    hitgeneclusterdata[geneclusternr] = [hitgeneclusters]
    return hitgeneclusterdata

def retrieve_pksnrps_info(seq_record, geneclusternr, pksnrpsprots):
    pksnrpsprotsnames = [utils.get_gene_id(cds) for cds in utils.get_pksnrps_cds_features(seq_record)]
    domaindict = utils.get_nrpspks_domain_dict(seq_record)
    substr_spec_preds = utils.get_nrpspks_substr_spec_preds(seq_record)
    pksnrpsdomains = {}
    domlist = []
    domsdetails = {}
    substrspecnrpspredictordict = {}
    substrspecminowadict = {}
    substrspecpkssigdict = {}
    substrspecconsensusdict = {}
    krpredictionsdict = {}
    for i in pksnrpsprots:
        domlist = []
        domsdetails = {}
        doms = domaindict[i]
        for j in doms:
            nr = 1
            while j[0] + str(nr) in domlist:
                nr += 1
            domname = j[0] + str(nr)
            domlist.append(domname)
            domsdetails[domname] = [j[1],j[2]]
            if "AMP-binding" in domname or "A-OX" in domname:
                domname2 = i + "_" + "A" + str(nr)
                substrspecminowadict[domname2] = substr_spec_preds.minowa_nrps_preds[i + "_A" + str(nr)]
                substrspecnrpspredictordict[domname2] = [substr_spec_preds.nrps_code_preds[i + "_A" + str(nr)], substr_spec_preds.nrps_svm_preds[i + "_A" + str(nr)]]
                substrspecconsensusdict[domname2] = substr_spec_preds.consensuspreds[i + "_A" + str(nr)]
            if "PKS_AT" in domname:
                domname2 = i + "_" + "AT" + str(nr)
                substrspecminowadict[domname2] = substr_spec_preds.minowa_pks_preds[i + "_AT" + str(nr)]
                substrspecpkssigdict[domname2] = substr_spec_preds.pks_code_preds[i + "_AT" + str(nr)]
                substrspecconsensusdict[domname2] = substr_spec_preds.consensuspreds[i + "_AT" + str(nr)]
            if "CAL_domain" in domname:
                domname2 = i + "_" + "CAL" + str(nr)
                substrspecminowadict[domname2] = substr_spec_preds.minowa_cal_preds[i + "_CAL" + str(nr)]
                substrspecconsensusdict[domname2] = substr_spec_preds.consensuspreds[i + "_CAL" + str(nr)]
            if "CAL_domain" in domname:
                domname2 = i + "_" + "CAL" + str(nr)
                substrspecminowadict[domname2] = substr_spec_preds.minowa_cal_preds[i + "_CAL" + str(nr)]
                substrspecconsensusdict[domname2] = substr_spec_preds.consensuspreds[i + "_CAL" + str(nr)]
            if "PKS_KR" in domname:
                domname2 = i + "_" + "KR" + str(nr)
                krpredictionsdict[domname2] = [substr_spec_preds.kr_activity_preds[i + "_KR" + str(nr)], substr_spec_preds.kr_stereo_preds[i + "_KR" + str(nr)]]
        pksnrpsdomains[i] = [domlist,domsdetails]
    structpred = utils.get_structure_pred(utils.get_cluster_by_nr(seq_record, geneclusternr))
    return pksnrpsprotsnames, pksnrpsdomains, domlist, domsdetails, substrspecnrpspredictordict, substrspecminowadict, substrspecpkssigdict, substrspecconsensusdict, krpredictionsdict, structpred
