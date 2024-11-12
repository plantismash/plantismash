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
# Copyright (C) 2024 Ziqiang Luo
# Wageningen University & Research, NL
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""ClusterBlast comparative gene cluster analysis"""

import logging
import subprocess
from os import path, makedirs
from antismash import utils
from clusterblast import load_clusterblast_database, internal_homology_blast, perform_clusterblast, filter_overlap
from data_loading import prepare_data, generate_Storage_for_cb

name = "clusterblast"
short_description = name.capitalize()
priority = 10000


# Tuple is ( binary_name, optional)
_required_binaries = [
    ('blastp', False),
    ('makeblastdb', False),
    ('diamond', False),
]

_required_files = [
    # ('geneclusterprots.dmnd', False),
    # ('geneclusterprots.fasta', False),
    # ('geneclusters.txt', False),
    ('plantgeneclusterprots.dmnd', False),
    ('plantgeneclusterprots.fasta', False),
    ('plantgeneclusters.txt', False),
]

def check_prereqs(options):
    "Check if all required applications are around"
    if options.clusterblastdir =="":
        options.clusterblastdir = utils.get_full_path(__file__, '')

    failure_messages = []
    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    for file_name, optional in _required_files:
        if utils.locate_file(path.join(options.clusterblastdir, file_name)) is None and not optional:
            if file_name == "plantgeneclusterprots.dmnd":
                # Run diamond makedb command
                diamond_command = ["diamond", "makedb", "--in", path.join(options.clusterblastdir, "plantgeneclusterprots.fasta"), "-d",
                                   path.join(options.clusterblastdir, "plantgeneclusterprots")]
                subprocess.check_call(diamond_command)
            else:
                failure_messages.append("Failed to locate file: %r" % file_name)

    return failure_messages

def run_clusterblast(seq_record, options):
    logging.info('Running ClusterBlast')
    clusters, proteinlocations, proteinstrands, proteinannotations, proteintags = load_clusterblast_database(seq_record)
    seq_record.internalhomologygroupsdict = internal_homology_blast(seq_record)
    perform_clusterblast(options, seq_record, clusters, proteinlocations, proteinstrands, proteinannotations, proteintags)
    prepare_data(seq_record, options, searchtype="general")
    generate_Storage_for_cb(options, seq_record)


def make_geneclusterprots(seq_records, options):
    """make gene cluster proteins fasta file for clusterblast"""
    names = []
    seqs = []
    for seq_record in seq_records:
        geneclusters = utils.get_sorted_cluster_features(seq_record)
        for genecluster in geneclusters:
            queryclusternames = []
            queryclusterseqs = []
            strand_start_ends = []
            queryclusterprots = filter_overlap(utils.get_cluster_cds_features(genecluster, seq_record))
            # only completely overlapping filtered

            for cds in queryclusterprots:
                if cds.strand == 1:
                    strand = "+"
                else:
                    strand = "-"
                start = str(cds.location.start).replace(">", "").replace("<", "")
                end = str(cds.location.end).replace(">", "").replace("<", "")
                strand_start_end = (strand, start, end)

                if strand_start_end not in strand_start_ends:
                    # todo: Incompletely overlapping splicing should be treated as one gene if their sequences are 50% similar？
                    strand_start_ends.append(strand_start_end)
                    annotation = utils.get_gene_annotation(cds)

                    fullname = "|".join([seq_record.id, "c" + str(utils.get_cluster_number(genecluster)),
                                         start + "-" + end,
                                         strand, utils.get_gene_id(cds), annotation.replace(' ', '_') , utils.get_gene_acc(cds)])
                    # scaffold|cluster number|location|strand|locustag|annotation|protein or accession
                    queryclusternames.append(fullname)
                    queryclusterseqs.append(str(utils.get_aa_sequence(cds)))

            for i in range(len(queryclusternames)):
                names.append(queryclusternames[i])
                seqs.append(queryclusterseqs[i])

    # put the fasta file in output folder
    outputname = path.join(path.abspath(options.outputfoldername), "plantgeneclusterprots.fasta")
    utils.writefasta(names, seqs, outputname)

def where_is_clusterblast():
    return utils.get_full_path(__file__, '')

