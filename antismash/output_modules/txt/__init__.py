# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2014 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011-2014 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# Copyright (C) 2014 Tilmann Weber
# The Novo Nordisk Foundation Center for Biosustainability
# Technical University of Denmark
# Section: Metabolic Engineering for Natural Compounds / New Bioactive Compounds
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""TXT output format module

"""

import os
from antismash import utils
from os import path
from . import write_tables
import logging

name = "txt"
short_description = "txt_database_output"
priority = 10000

def write(seq_records, options):
    logging.debug("Exporting antiSMASH information as txt tables")
    #Don't store TXT tables for protein input
    if options.input_type == 'prot':
        return
    #Localize output folder, create TXT subdirectory
    txt_outfolder = options.full_outputfolder_path + os.sep + "txt"
    if not os.path.exists(txt_outfolder):
        os.mkdir(txt_outfolder)
    #Define table names
    tables = "genome", "BGC", "signature_gene_info", "gene", "NRPS_PKS", "smCOG", "RiPP", "transltable"
    #For each gene cluster, write out info to TXT files
    for seq_record in seq_records:
        if len(utils.get_cluster_features(seq_record)) > 0:
            #Open up TXT files
            txt_files = {}
            for table in tables:
                txt_files[table] = open(path.join(txt_outfolder, "%s_%s.txt" % (seq_record.id.partition(".")[0], table)),"w")
            #Gather all information
            info = utils.Storage()
            info.clustertypes, info.clustergenes, info.accessions, info.cdsmotifs, info.clusternrs = {}, {}, {}, {}, []
            clusters = utils.get_cluster_features(seq_record)
            for cluster in clusters:
                clusternr = utils.get_cluster_number(cluster)
                info.clusternrs.append(clusternr)
                info.clustertypes[clusternr] = utils.get_cluster_type(cluster)
                info.clustergenes[clusternr] = [utils.get_gene_id(cds) for cds in utils.get_cluster_cds_features(cluster, seq_record)]
                info.accessions[clusternr] = [utils.get_gene_acc(cds) for cds in utils.get_cluster_cds_features(cluster, seq_record)]
                info.cdsmotifs[clusternr] = utils.get_all_features_of_type(seq_record, ["CDS_motif"])
            info.seq_record = seq_record
            #Write information to tables
            for table in tables:
                getattr(write_tables, 'write_' + table)(txt_files[table], info, options)
            for table in tables:
                txt_files[table].close()
