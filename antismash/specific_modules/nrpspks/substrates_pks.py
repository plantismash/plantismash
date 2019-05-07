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
from antismash import utils
import os
from os import path
from Minowa import minowa_CAL, minowa_AT
from kr_analysis import kr_analysis
from pkssignatures import PKS_analysis
from helperlibs.wrappers.io import TemporaryDirectory


def count_pks_genes(pksnrpscoregenes, domaindict, seq_record):
    pkscount = 0
    for feature in pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        domaindetails = domaindict[locus]
        nr = 0
        for tab in domaindetails:
            if tab[0] == "PKS_AT" or tab[0] == "CAL_domain" or tab[0] == "PKS_KR":
                pkscount += 1
    return pkscount

def extract_pks_genes(pksnrpscoregenes, domaindict, seq_record):
    pksnames = []
    pksseqs = []
    for feature in pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        domaindetails = domaindict[locus]
        nr = 0
        for tab in domaindetails:
            if tab[0] == "PKS_AT":
                nr += 1
                start = int(tab[1])
                end = int(tab[2])
                seq = str(utils.get_aa_sequence(feature))[start:end]
                name = locus + "_AT" + str(nr)
                pksnames.append(name)
                pksseqs.append(seq)
    return pksnames, pksseqs

def run_minowa_predictor_pks_at(pksnames, pksseqs, options):
    #Predict PKS AT domain specificities with Minowa et al. method and PKS code (NP searcher / ClustScan / own?)
    utils.writefasta(pksnames, pksseqs,
               options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_pksseqs.fasta")
    #Run PKS signature analysis
    logging.info("Predicting PKS AT domain substrate specificities by Yadav et al. PKS signature sequences")
    with TemporaryDirectory(change=True):
        PKS_analysis.run_pkssignature_analysis(
                path.join(options.raw_predictions_outputfolder, "ctg" + str(options.record_idx) + "_pksseqs.fasta"),
                path.join(options.raw_predictions_outputfolder, "ctg" + str(options.record_idx) + "_pkssignatures.txt"))

    #Minowa method: run Minowa_AT
    logging.info("Predicting PKS AT domain substrate specificities by Minowa et al. method")
    with TemporaryDirectory(change=True):
        minowa_AT.run_minowa_at(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_pksseqs.fasta", options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_minowa_pkspredoutput.txt")

def run_minowa_predictor_pks_cal(pksnrpscoregenes, domaindict, seq_record, options):
    calnames = []
    calseqs = []
    #Predict PKS CAL domain specificities with Minowa et al. method
    logging.info("Predicting CAL domain substrate specificities by Minowa et al. method")
    for feature in pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        domaindetails = domaindict[locus]
        nr = 0
        for tab in domaindetails:
            if tab[0] == "CAL_domain":
                nr += 1
                start = int(tab[1])
                end = int(tab[2])
                seq = str(utils.get_aa_sequence(feature))[start:end]
                name = locus + "_CAL" + str(nr)
                calnames.append(name)
                calseqs.append(seq)
    if len(calnames) > 0:
        utils.writefasta(calnames, calseqs, options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_calseqs.fasta")
        with TemporaryDirectory(change=True):
            minowa_CAL.run_minowa_cal(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_calseqs.fasta", options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_minowa_calpredoutput.txt")
    return calnames, calseqs

def run_kr_stereochemistry_predictions(pksnrpscoregenes, domaindict, seq_record, options):
    #Predict PKS KR domain stereochemistry using pattern as published in ClustScan
    krnames = []
    krseqs = []
    logging.info("Predicting PKS KR activity and stereochemistry using KR " \
        "fingerprints from Starcevic et al.")
    for feature in pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        domaindetails = domaindict[locus]
        nr = 0
        for tab in domaindetails:
            if tab[0] == "PKS_KR":
                nr += 1
                start = int(tab[1])
                end = int(tab[2])
                seq = str(utils.get_aa_sequence(feature))[start:end]
                name = locus + "_KR" + str(nr)
                krnames.append(name)
                krseqs.append(seq)
    if len(krnames) > 0:
        utils.writefasta(krnames, krseqs, options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_krseqs.fasta")
        with TemporaryDirectory(change=True):
            kr_analysis.run_kr_analysis(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_krseqs.fasta", options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_krpredoutput.txt")
    return krnames, krseqs
