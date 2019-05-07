# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2014 Marnix H. Medema
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
import os
from antismash import utils
from antismash.generic_modules.clusterblast.clusterblast import *
from helperlibs.wrappers.io import TemporaryDirectory

def perform_knownclusterblast(options, seq_record, clusters, proteinlocations, proteinstrands, proteinannotations, proteintags):
    #Run BLAST on gene cluster proteins of each cluster and parse output
    logging.info("Running DIAMOND knowncluster searches..")
    geneclusters = utils.get_sorted_cluster_features(seq_record)
    with TemporaryDirectory(change=True) as tempdir:
        for genecluster in geneclusters:
            clusternumber = utils.get_cluster_number(genecluster)
            if options.debug and os.path.exists(options.dbgclusterblast + os.sep + "knwonclusterblast" + os.sep + "cluster" + str(clusternumber) + ".txt"):
                logging.debug ("Skipping SubClusterblast calculations, using results from %s instead" % options.dbgclusterblast + os.sep + "knownclusterblast" + os.sep + "cluster" + str(clusternumber) + ".txt")
            else:
                
                logging.info("   Gene cluster " + str(clusternumber))
                queryclusternames, queryclusterseqs, queryclusterprots = create_blast_inputs(genecluster, seq_record)
                utils.writefasta([qcname.replace(" ","_") for qcname in queryclusternames], queryclusterseqs, "input.fasta")
                out, err, retcode = run_diamond("input.fasta", path.join(options.knownclusterblastdir, 'knownclusterprots'), tempdir, options)
                if retcode != 0:
                    logging.debug("out: %r, err: %r, retcode: %s", out, err, retcode)
                convert_to_tabular(tempdir)
                with open("input.out", 'r') as fh:
                    blastoutput = fh.read()
                write_raw_clusterblastoutput(options.full_outputfolder_path, blastoutput, searchtype="knownclusters")
                logging.info("   DIAMOND search finished. Parsing results...")
                minseqcoverage = 40
                minpercidentity = 45
                blastdict, querylist, hitclusters = parse_blast(blastoutput, seq_record, minseqcoverage, minpercidentity)
                querylist = remove_queries_without_hits(querylist, blastdict)
                allcoregenes = [utils.get_gene_id(cds) for cds in utils.get_secmet_cds_features(seq_record)]
                rankedclusters, rankedclustervalues, hitclusterdict, hitclusterdata = score_clusterblast_output(blastdict, querylist, hitclusters, clusters, allcoregenes)
                
                # store all clusterblast related data in a utils.Storage object and serialize it
                knownclusterblastStorage = utils.Storage()
                knownclusterblastStorage.clusternumber = clusternumber
                knownclusterblastStorage.queryclusterprots = queryclusterprots
                knownclusterblastStorage.clusters = clusters
                knownclusterblastStorage.hitclusterdata = hitclusterdata
                knownclusterblastStorage.rankedclusters = rankedclusters
                knownclusterblastStorage.rankedclustervalues = rankedclustervalues
                knownclusterblastStorage.proteintags = proteintags
                knownclusterblastStorage.proteinlocations = proteinlocations
                knownclusterblastStorage.proteinannotations = proteinannotations
                knownclusterblastStorage.proteinstrands = proteinstrands
                
                write_clusterblast_output(options, seq_record, knownclusterblastStorage, searchtype="knownclusters")
