# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#clusterblastStorage = utils.Storage()
# Copyright (C) 2011,2012 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os
from antismash import utils
from antismash.generic_modules.clusterblast.clusterblast import *
from helperlibs.wrappers.io import TemporaryDirectory

def perform_subclusterblast(options, seq_record, clusters, proteinlocations, proteinstrands, proteinannotations, proteintags):
    #Run BLAST on gene cluster proteins of each cluster and parse output
    logging.info("Running NCBI BLAST+ subcluster searches..")
    geneclusters = utils.get_sorted_cluster_features(seq_record)
    with TemporaryDirectory(change=True):
        for genecluster in geneclusters:
            clusternumber = utils.get_cluster_number(genecluster)
            if options.debug and os.path.exists(options.dbgclusterblast + os.sep + "subclusterblast" + os.sep + "cluster" + str(clusternumber) + ".txt"):
                logging.debug ("Skipping SubClusterblast calculations, using results from %s instead" % options.dbgclusterblast + os.sep + "subclusterblast" + os.sep + "cluster" + str(clusternumber) + ".txt")
            else:
                logging.info("   Gene cluster " + str(clusternumber))
                queryclusternames, queryclusterseqs, queryclusterprots = create_blast_inputs(genecluster, seq_record)
                write_clusterblast_inputfiles(options, queryclusternames, queryclusterseqs)
                run_clusterblast_processes(options, searchtype="subclusters")
                blastoutput = read_clusterblast_output(options)
                write_raw_clusterblastoutput(options.full_outputfolder_path, blastoutput, searchtype="subclusters")
                logging.info("   Blast search finished. Parsing results...")
                minseqcoverage = 40
                minpercidentity = 45
                blastdict, querylist, hitclusters = parse_blast(blastoutput, seq_record, minseqcoverage, minpercidentity)
                querylist = remove_queries_without_hits(querylist, blastdict)
                allcoregenes = [utils.get_gene_acc(cds) for cds in utils.get_secmet_cds_features(seq_record)]
                rankedclusters, rankedclustervalues, hitclusterdict, hitclusterdata = score_clusterblast_output(blastdict, querylist, hitclusters, clusters, allcoregenes)
                
                # store all clusterblast related data in a utils.Storage object and serialize it
                subclusterblastStorage = utils.Storage()
                subclusterblastStorage.clusternumber = clusternumber
                subclusterblastStorage.queryclusterprots = queryclusterprots
                subclusterblastStorage.clusters = clusters
                subclusterblastStorage.hitclusterdata = hitclusterdata
                subclusterblastStorage.rankedclusters = rankedclusters
                subclusterblastStorage.rankedclustervalues = rankedclustervalues
                subclusterblastStorage.proteintags = proteintags
                subclusterblastStorage.proteinlocations = proteinlocations
                subclusterblastStorage.proteinannotations = proteinannotations
                subclusterblastStorage.proteinstrands = proteinstrands
                    
                write_clusterblast_output(options, seq_record, subclusterblastStorage, searchtype="subclusters")
