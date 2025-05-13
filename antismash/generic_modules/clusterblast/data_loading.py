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

import antismash
from antismash import utils
from argparse import Namespace
import os
from os import path
import logging
from string import ascii_letters
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re
import sys

def prepare_data(seq_record, options, searchtype="general"):
    load_pubmed_pubchem_links(seq_record)
    if searchtype == "general":
        load_clusterblast_outputdata(seq_record, options)
    if searchtype == "subclusters":
        load_subclusterblast_outputdata(seq_record, options)
    if searchtype == "knownclusters":
        load_knownclusterblast_outputdata(seq_record, options)
    load_genecluster_info(seq_record, options, searchtype=searchtype)

def load_pubmed_pubchem_links(seq_record):
    #Read in PubMed / PubChem links of database gene clusters
    seq_record.pubmed_dict = {}
    seq_record.pubchem_dict = {}
    seq_record.known_compound_dict = {}
    seq_record.closestcompounddict = {}
    pubchem_link_file = path.join(
            path.abspath(path.dirname(antismash.__file__)),
            'lib', 'pubmed_pubchem_links.txt')
    pubfile = open(pubchem_link_file, "r")
    pubfile = pubfile.read()
    publines = pubfile.split("\n")
    for i in publines:
        tabs = i.split("\t")
        acc = tabs[0]
        if tabs[1] != "":
            seq_record.pubmed_dict[acc] = tabs[1]
        if tabs[2] != "":
            seq_record.pubchem_dict[acc] = tabs[2]
        if tabs[3] != "":
            seq_record.known_compound_dict[acc] = tabs[3]

def parse_clusterblast_details(options, seq_record, clusternr, details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict, searchtype="general"):
    #For every gene cluster, store hit genes and details
    colorgroupsdict = {}
    hitclusterdata = {}
    hitclusternr = 1
    compound_found = "n"
    allcoregenes = [utils.get_gene_id(cds) for cds in utils.get_secmet_cds_features(seq_record)]
    if searchtype == "general":
        seq_record.nrhitgeneclusters[clusternr] = 0
    elif searchtype == "subclusters":
        seq_record.sc_nrhitgeneclusters[clusternr] = 0
    elif searchtype == "knownclusters":
        seq_record.kc_nrhitgeneclusters[clusternr] = 0
    for i in details:
        hitclustergenes = []
        hitclustergenesdetails = {}
        hits_accessions_dict = {}
        #Only calculate for first ten hit gene clusters
        if hitclusternr <= options.nclusters:
            if searchtype == "general":
                seq_record.nrhitgeneclusters[clusternr] = hitclusternr
            elif searchtype == "subclusters":
                seq_record.sc_nrhitgeneclusters[clusternr] = hitclusternr
            elif searchtype == "knownclusters":
                seq_record.kc_nrhitgeneclusters[clusternr] = hitclusternr
            accession = cb_accessiondict[hitclusternr]
            hitclustergeneslines = ((i.split("Table of genes, locations, strands and annotations of subject cluster:\n")[1]).split("\n\nTable of Blast hits ")[0]).split("\n")
            for j in hitclustergeneslines:
                tabs = j.split("\t")
                hitclustergenes.append(tabs[0])
                hitclustergenesdetails[tabs[0]] = [tabs[2],tabs[3],tabs[4],tabs[5],tabs[1]]
            blasthitslines = ((i.split("%coverage, e-value):\n")[1]).split("\n\n")[0]).split("\n")
            querygeneswithhits = []
            coregeneswithhits = []
            for k in blasthitslines:
                if k.split("\t")[0] not in querygeneswithhits:
                    querygeneswithhits.append(k.split("\t")[0])
                if k.split("\t")[1] in hits_accessions_dict:
                    hits_accessions_dict[k.split("\t")[1]].append(k.split("\t")[0])
                else:
                    hits_accessions_dict[k.split("\t")[1]] = [k.split("\t")[0]]
                if k.split("\t")[0] in allcoregenes and k.split("\t")[0] not in coregeneswithhits:
                    coregeneswithhits.append(k.split("\t")[0])
            if searchtype == "general":
                for k in list(seq_record.known_compound_dict.keys()):
                    if k in i and compound_found == "n" and len(querygeneswithhits) > 2 and len(coregeneswithhits) > 0:
                        seq_record.closestcompounddict[clusternr] = seq_record.known_compound_dict[k]
            compound_found = "y"
            #Create Blast dict
            blasthitdict, blastdetailsdict, querygenes, hitgenes, revblasthitdict = create_blastdicts(blasthitslines)
            #Create colorgroups dict
            colorgroupsdict = construct_colorgroups(colorgroupsdict, clusternr, blasthitdict, blastdetailsdict, seq_record.internalhomologygroupsdict, hitclusternr)
            #Store all data in hitclusterdata dict
            hitclusterdata[hitclusternr] = [colorgroupsdict, hitclustergenes, hitclustergenesdetails, queryclustergenes, queryclustergenesdetails, toptenhitclusters, accession, hits_accessions_dict, blastdetailsdict]
            hitclusternr += 1
        elif hitclusternr > options.nclusters and hitclusternr <= 50:
            blasthitslines = ((i.split("%coverage, e-value):\n")[1]).split("\n\n")[0]).split("\n")
            querygeneswithhits = []
            coregeneswithhits = []
            for k in blasthitslines:
                if k.split("\t")[0] not in querygeneswithhits:
                    querygeneswithhits.append(k.split("\t")[0])
                if k.split("\t")[0] in allcoregenes and k.split("\t")[0] not in coregeneswithhits:
                    coregeneswithhits.append(k.split("\t")[0])
            hitclusternr += 1
    if searchtype == "general":
        seq_record.queryclusterdata[clusternr] = [nrhitclusters,hitclusterdata]
    elif searchtype == "subclusters":
        seq_record.sc_queryclusterdata[clusternr] = [nrhitclusters,hitclusterdata]
    elif searchtype == "knownclusters":
        seq_record.kc_queryclusterdata[clusternr] = [nrhitclusters,hitclusterdata]




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

def create_blastdicts(blasthitslines):
    blasthitdict = {}
    blastdetailsdict = {}
    querygenes = []
    revblasthitdict = {}
    hitgenes = []
    for i in blasthitslines:
        tabs = i.split("\t")
        if len(tabs) <= 1:
            continue
        if tabs[0] in blasthitdict:
            hits = blasthitdict[tabs[0]]
            hits.append(tabs[1])
            blasthitdict[tabs[0]] = hits
            if tabs[1] in revblasthitdict:
                revhits = revblasthitdict[tabs[1]]
                revhits.append(tabs[0])
                revblasthitdict[tabs[1]] = revhits
            else:
                revblasthitdict[tabs[1]] = [tabs[0]]
            blastdetailsdict[tabs[0] + "_|_|_" + tabs[1]] = [tabs[5], tabs[3], tabs[2], tabs[4]]
            if tabs[0] not in querygenes:
                querygenes.append(tabs[0])
            hitgenes.append(tabs[1])
        else:
            blasthitdict[tabs[0]] = [tabs[1]]
            if tabs[1] in revblasthitdict:
                revhits = revblasthitdict[tabs[1]]
                revhits.append(tabs[0])
                revblasthitdict[tabs[1]] = revhits
            else:
                revblasthitdict[tabs[1]] = [tabs[0]]
            blastdetailsdict[tabs[0] + "_|_|_" + tabs[1]] = [tabs[5], tabs[3], tabs[2], tabs[4]]
            if tabs[0] not in querygenes:
                querygenes.append(tabs[0])
            hitgenes.append(tabs[1])
    return blasthitdict, blastdetailsdict, querygenes, hitgenes, revblasthitdict

def generate_Storage_for_cb(options, seq_record, searchtype="ClusterBlastData"):
    """
    This is a very ugly helper function to convert all data stored in "non_standard" lists/dictionaries within the seq_record object
    into a storage object which can be saved in a qualifier...
    
    THIS SHOULD BE REFACTORED so that the information is directly stored in this object instead of this ugly conversion...
    """
    clusterBlastResults = utils.Storage()
    try:
        clusterBlastResults.internalhomologygroupsdict = seq_record.internalhomologygroupsdict
    except AttributeError:
        logging.debug("seq_record.internalhomologygroupsdict does not exist")
    try:
        clusterBlastResults.known_compound_dict = seq_record.known_compound_dict
        del seq_record.known_compound_dict
    except AttributeError:
        logging.debug("seq_record.known_compound_dict does not exist.")
    try:
        clusterBlastResults.pubchem_dict = seq_record.pubchem_dict
        del seq_record.pubchem_dict
    except AttributeError:
        logging.debug("seq_record.pubchem_dict does not exist.")
    try:
        clusterBlastResults.pubmed_dict = seq_record.pubmed_dict
        del seq_record.pubmed_dict
    except AttributeError:
        logging.debug("seq_record.pubmed_dict does not exist.")

    if searchtype == "ClusterBlastData":
        try:
            clusterBlastResults.nrhitgeneclusters = seq_record.nrhitgeneclusters
            del seq_record.nrhitgeneclusters
        except AttributeError:
            logging.debug("seq_record.nrhitgeneclusters does not exist.")
        try:
            clusterBlastResults.qgeneclusterdata = seq_record.qgeneclusterdata
            del seq_record.qgeneclusterdata
        except AttributeError:
            logging.debug("qgeneclusterdata does not exist.")
        try:
            clusterBlastResults.queryclusterdata = seq_record.queryclusterdata
            del seq_record.queryclusterdata
        except AttributeError:
            logging.debug("seq_record.queryclusterdata does not exist.")
            
        if 'pubmed_dict' in seq_record:
            clusterBlastResults.pubmed_dict = seq_record.pubmed_dict
        if 'pubchem_dict' in seq_record:
            clusterBlastResults.pubchem_dict = seq_record.pubchem_dict
        if 'known_compound_dict' in seq_record:
            clusterBlastResults.known_compound_dict = seq_record.known_compound_dict
        if 'closestcompounddict' in seq_record:
            clusterBlastResults.closestcompounddict = seq_record.closestcompounddict
    
    if searchtype == "SubClusterBlastData":
        try:
            clusterBlastResults.sc_nrhitgeneclusters = seq_record.sc_nrhitgeneclusters
            del seq_record.sc_nrhitgeneclusters
        except AttributeError:
            logging.debug("seq_record.sc_nrhitgeneclusters does not exist.")
        try:
            clusterBlastResults.sc_qgeneclusterdata = seq_record.sc_qgeneclusterdata
            del seq_record.sc_qgeneclusterdata
        except AttributeError:
            logging.debug("seq_record.sc_qgeneclusterdata does not exist.")
        try:
            clusterBlastResults.sc_queryclusterdata = seq_record.sc_queryclusterdata
            del seq_record.sc_queryclusterdata
        except AttributeError:
            logging.debug("seq_record.sc_queryclusterdata does not exist.")
            
    if searchtype == "KnownClusterBlastData":
        try:
            clusterBlastResults.kc_nrhitgeneclusters = seq_record.kc_nrhitgeneclusters
            del seq_record.kc_nrhitgeneclusters
        except AttributeError:
            logging.debug("seq_record.kc_nrhitgeneclusters does not exist.")
        try:
            clusterBlastResults.kc_qgeneclusterdata = seq_record.kc_qgeneclusterdata
            del seq_record.sc_qgeneclusterdata
        except AttributeError:
            logging.debug("seq_record.kc_qgeneclusterdata does not exist.")
        try:
            clusterBlastResults.kc_queryclusterdata = seq_record.kc_queryclusterdata
            del seq_record.kc_queryclusterdata
        except AttributeError:
            logging.debug("seq_record.kc_queryclusterdata does not exist.")
            
            
    
        
    if not 'extrarecord' in options:
        options.extrarecord = {}
    if seq_record.id not in options.extrarecord:
        options.extrarecord[seq_record.id] = Namespace()
        
    if not 'extradata' in options.extrarecord[seq_record.id]:
        options.extrarecord[seq_record.id].extradata = {}
    logging.debug("Storing data for %s in storage object" % searchtype)
    options.extrarecord[seq_record.id].extradata[searchtype] = clusterBlastResults


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

def read_clusterblastfile(seq_record, options, clusternr, searchtype="general"):
    
    if searchtype == "general":
        if os.path.exists(options.dbgclusterblast + os.sep + "clusterblast" + os.sep + "cluster" + str(clusternr) + ".txt"):
            clusterblastfile = open(options.dbgclusterblast + os.sep + "clusterblast" + os.sep + "cluster" + str(clusternr) + ".txt", "r")
        else:
            clusterblastfile = open(options.clusterblast_outputfolder + os.sep + "cluster" + str(clusternr) + ".txt","r")
    elif searchtype == "subclusters":
        if os.path.exists(options.dbgclusterblast + os.sep +"subclusterblast" + os.sep + "cluster" + str(clusternr) + ".txt"):
            clusterblastfile = open(options.dbgclusterblast + os.sep +"subclusterblast" + os.sep + "cluster" + str(clusternr) + ".txt", "r")
        else:
            clusterblastfile = open(options.subclusterblast_outputfolder + os.sep + "cluster" + str(clusternr) + ".txt","r")
    elif searchtype == "knownclusters":
        if os.path.exists(options.dbgclusterblast + os.sep +"knownclusterblast" + os.sep + "cluster" + str(clusternr) + ".txt"):
            clusterblastfile = open(options.dbgclusterblast + os.sep +"knownclusterblast" + os.sep + "cluster" + str(clusternr) + ".txt", "r")
        else:
            clusterblastfile = open(options.knownclusterblast_outputfolder + os.sep + "cluster" + str(clusternr) + ".txt","r")
    else:
        logging.exception("Illegal search type %s" % searchtype)
    
    
    clusterblastfile = clusterblastfile.read()
    clusterblastfile = clusterblastfile.replace("\r","\n")
    toptenhitclusters = []
    #Identify top ten hits for visualization
    hitlines = ((clusterblastfile.split("Significant hits: \n")[1]).split("\nDetails:")[0]).split("\n")
    a = 0
    cb_accessiondict = {}
    b = 1
    for i in hitlines:
        if " " in i:
            cb_accessiondict[b] = (i.split("\t")[0]).split(" ")[1]
        if not test_accession(seq_record.id) or (". " + seq_record.id.partition(".")[0] not in i.split("\t")[0] and "_" + seq_record.id.partition(".")[0] not in i.split("\t")[0]):
            b += 1
            if a < options.nclusters:
                if len(i) > 3:
                    percentage_genes_with_hits = get_hitnumbers_from_clusterblastfile(clusterblastfile, i)
                if len(i) < 80 and len(i) > 3:
                    toptenhitclusters.append(i + " (%i%% of genes show similarity)" % percentage_genes_with_hits)
                elif len(i) >= 80:
                    j = i[0:77] + "... (%i%% of genes show similarity)" % percentage_genes_with_hits
                    toptenhitclusters.append(j)
            a += 1
    if test_accession(seq_record.id):
        details = [details for details in (clusterblastfile.split("\nDetails:")[1]).split(">>")[1:] if (seq_record.id.partition(".")[0] not in details.partition("Source:")[0] and "_" + seq_record.id.partition(".")[0] not in details.partition("Source:")[0])]
    else:
        details = (clusterblastfile.split("\nDetails:")[1]).split(">>")[1:]
    nrhitclusters = len(toptenhitclusters)
    #Save query gene cluster data
    querylines = ((clusterblastfile.split("Table of genes, locations, strands and annotations of query cluster:\n")[1]).split("\n\n\nSignificant hits:")[0]).split("\n")
    queryclustergenes = []
    queryclustergenesdetails = {}
    for i in querylines:
        tabs = i.split("\t")
        if len(tabs) > 1:
            queryclustergenes.append(tabs[0])
            queryclustergenesdetails[tabs[0]] = [tabs[1],tabs[2],tabs[3],tabs[4]]
    return details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict

def get_hitnumbers_from_clusterblastfile(clusterblastfile, clusterhit):
    #Get number of genes with homologues
    accession = clusterhit.partition(" ")[2].partition("\t")[0]
    hitdetails = [details for details in (clusterblastfile.split("\nDetails:")[1]).split(">>")[1:] if accession in details][0]
    #Get total number of genes in hitcluster
    nr_hitclustergenes = len(clusterblastfile.partition("of subject cluster:\n")[2].partition("\n\n")[0].split("\n"))
    #Get percentage
    blasthits = hitdetails.partition(" %coverage, e-value):\n")[2].partition("\n\n")[0].split("\n")
    nr_genes_with_hits = len(set([hit.split("\t")[1] for hit in blasthits]))
    percentage_genes_with_hits = int((float(nr_genes_with_hits) / float(nr_hitclustergenes)) * 100)
    return percentage_genes_with_hits

def load_subclusterblast_outputdata(seq_record, options):
    #Read in SubClusterBlast data
    seq_record.sc_queryclusterdata = {}
    seq_record.sc_nrhitgeneclusters = {}
    geneclusters = utils.get_sorted_cluster_features(seq_record)
    for genecluster in geneclusters:
        clusternr = utils.get_cluster_number(genecluster)
        details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict = read_clusterblastfile(seq_record, options, clusternr, searchtype="subclusters")
        parse_clusterblast_details(options, seq_record, clusternr, details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict, searchtype="subclusters")
        genecluster.qualifiers['subclusterblast'] = toptenhitclusters

def load_knownclusterblast_outputdata(seq_record, options):
    #Read in KnownClusterBlast data
    seq_record.kc_queryclusterdata = {}
    seq_record.kc_nrhitgeneclusters = {}
    geneclusters = utils.get_sorted_cluster_features(seq_record)
    for genecluster in geneclusters:
        clusternr = utils.get_cluster_number(genecluster)
        details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict = read_clusterblastfile(seq_record, options, clusternr, searchtype="knownclusters")
        parse_clusterblast_details(options, seq_record, clusternr, details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict, searchtype="knownclusters")
        genecluster.qualifiers['knownclusterblast'] = toptenhitclusters

def load_clusterblast_outputdata(seq_record, options):
    #Read in ClusterBlast data
    seq_record.queryclusterdata = {}
    seq_record.nrhitgeneclusters = {}
    geneclusters = utils.get_sorted_cluster_features(seq_record)
    for genecluster in geneclusters:
        clusternr = utils.get_cluster_number(genecluster)
        details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict = read_clusterblastfile(seq_record, options, clusternr)
        parse_clusterblast_details(options, seq_record, clusternr, details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict)
        genecluster.qualifiers['clusterblast'] = toptenhitclusters
# def read_clusterblastfile(seq_record, options, clusternr, searchtype="general"):
#     logging.debug("Searchtype is %s" % searchtype)
#     if searchtype == "general":
#         if os.path.exists(options.dbgclusterblastres):
#             clusterblastfile = open(options.dbgclusterblastres, "r")
#         else:
#             clusterblastfile = open(options.clusterblast_outputfolder + os.sep + "cluster" + str(clusternr) + ".txt","r")
#     elif searchtype == "subclusters":
#         if os.path.exists(options.dbgsubclusterblastres):
#             clusterblastfile = open(options.dbgsubclusterblastres)
#         else:
#             clusterblastfile = open(options.subclusterblast_outputfolder + os.sep + "cluster" + str(clusternr) + ".txt","r")
#     else:
#         logging.exception("Illegal search type %s" % searchtype)
#     clusterblastfile = clusterblastfile.read()
#     clusterblastfile = clusterblastfile.replace("\r","\n")
#     toptenhitclusters = []
#     #Identify top ten hits for visualization
#     hitlines = ((clusterblastfile.split("Significant hits: \n")[1]).split("\nDetails:")[0]).split("\n")
#     a = 0
#     cb_accessiondict = {}
#     b = 1
#     for i in hitlines:
#         if " " in i:
#             cb_accessiondict[b] = (i.split("\t")[0]).split(" ")[1]
#         if not test_accession(seq_record.id) or (". " + seq_record.id.partition(".")[0] not in i.split("\t")[0] and "_" + seq_record.id.partition(".")[0] not in i.split("\t")[0]):
#             b += 1
#             if a < 10:
#                 if len(i) < 80 and len(i) > 3:
#                     toptenhitclusters.append(i)
#                 elif len(i) >= 80:
#                     j = i[0:77] + "..."
#                     toptenhitclusters.append(j)
#             a += 1
#     if test_accession(seq_record.id):
#         details = [details for details in (clusterblastfile.split("\nDetails:")[1]).split(">>")[1:] if (seq_record.id.partition(".")[0] not in details.partition("Source:")[0] and "_" + seq_record.id.partition(".")[0] not in details.partition("Source:")[0])]
#     else:
#         details = (clusterblastfile.split("\nDetails:")[1]).split(">>")[1:]
#     nrhitclusters = len(toptenhitclusters)
#     #Save query gene cluster data
#     querylines = ((clusterblastfile.split("Table of genes, locations, strands and annotations of query cluster:\n")[1]).split("\n\n\nSignificant hits:")[0]).split("\n")
#     queryclustergenes = []
#     queryclustergenesdetails = {}
#     for i in querylines:
#         tabs = i.split("\t")
#         if len(tabs) > 1:
#             queryclustergenes.append(tabs[0])
#             queryclustergenesdetails[tabs[0]] = [tabs[1],tabs[2],tabs[3],tabs[4]]
#     return details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict
# 
# def load_clusterblast_outputdata(seq_record, options):
#     #Read in ClusterBlast data
#     seq_record.queryclusterdata = {}
#     seq_record.nrhitgeneclusters = {}
#     if "subcblastclusternr" not in options:
#         options.subcblastclusternr = 1
#     geneclusters = utils.get_cluster_features(seq_record)
#     for genecluster in geneclusters:
#         clusternr = utils.get_cluster_number(genecluster)
#         details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict = read_clusterblastfile(seq_record, options, clusternr)
#         logging.debug("load_clusterblast_outputdata; options.subcblastclusternr = %s" % options.subcblastclusternr)
#         options.subcblastclusternr = parse_clusterblast_details(seq_record, clusternr, options.subcblastclusternr, details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict)
#         genecluster.qualifiers['clusterblast'] = toptenhitclusters
# 
# def load_subclusterblast_outputdata(seq_record, options):
#     #Read in SubClusterBlast data
#     seq_record.sc_queryclusterdata = {}
#     seq_record.sc_nrhitgeneclusters = {}
#     if "cblastclusternr" not in options:
#         options.cblastclusternr = 1
#     geneclusters = utils.get_cluster_features(seq_record)
#     for genecluster in geneclusters:
#         clusternr = utils.get_cluster_number(genecluster)
#         details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict = read_clusterblastfile(seq_record, options, clusternr, searchtype="subclusters")
#         logging.debug("load_subclusterblast_outputdata; options.cblastclusternr = %s" % options.cblastclusternr)
#         options.cblastclusternr = parse_clusterblast_details(seq_record, clusternr, options.cblastclusternr, details, toptenhitclusters, nrhitclusters, queryclustergenes, queryclustergenesdetails, cb_accessiondict, searchtype="subclusters")
#         genecluster.qualifiers['subclusterblast'] = toptenhitclusters

def load_genecluster_info(seq_record, options, searchtype="general"):
    #Gather and store data on each gene cluster
    smcogdict, smcogdescriptions = utils.get_smcog_annotations(seq_record)
    gtrcoglist = ['SMCOG1045','SMCOG1062','SMCOG1102']
    transportercoglist = ['SMCOG1000','SMCOG1005','SMCOG1011','SMCOG1020','SMCOG1029','SMCOG1033','SMCOG1035','SMCOG1044','SMCOG1065','SMCOG1067','SMCOG1069','SMCOG1074','SMCOG1085','SMCOG1096','SMCOG1106','SMCOG1118','SMCOG1131','SMCOG1166','SMCOG1169','SMCOG1184','SMCOG1202','SMCOG1205','SMCOG1214','SMCOG1234','SMCOG1243','SMCOG1245','SMCOG1252','SMCOG1254','SMCOG1288']
    seq_record.qgeneclusterdata = {}
    geneclusters = utils.get_sorted_cluster_features(seq_record)
    for genecluster in geneclusters:
        geneclusternr = utils.get_cluster_number(genecluster)
        clustergenes, clustertype, annotations, colors, starts, ends, strands, pksnrpsprots, gtrs, transporters, clustersize = retrieve_gene_cluster_annotations(seq_record, smcogdict, gtrcoglist, transportercoglist, geneclusternr)
        if options.clusterblast:
            hitgeneclusterdata = retrieve_clusterblast_info(seq_record, geneclusternr, searchtype=searchtype)
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
    clustersize = max(ends) - min(starts)
    return clustergenes, clustertype, annotations, colors, starts, ends, strands, pksnrpsprots, gtrs, transporters, clustersize

def retrieve_clusterblast_info(seq_record, geneclusternr, searchtype="general"):
    if searchtype == "general":
        hitgeneclusters = list(range(1,(seq_record.nrhitgeneclusters[geneclusternr] + 1)))
    elif searchtype =="subclusters":
        hitgeneclusters = list(range(1,(seq_record.sc_nrhitgeneclusters[geneclusternr] + 1)))
    elif searchtype == "knownclusters":
        hitgeneclusters = list(range(1,(seq_record.kc_nrhitgeneclusters[geneclusternr] + 1)))
    else:
        logging.exception("unknown searchtype in retrieve_clusterblast_info")
        sys.exit(1)
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
