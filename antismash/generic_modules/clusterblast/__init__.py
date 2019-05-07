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

"""ClusterBlast comparative gene cluster analysis"""

import logging
from os import path
from antismash import utils
from clusterblast import load_clusterblast_database, internal_homology_blast, perform_clusterblast
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
    ('geneclusterprots.dmnd', False),
    ('geneclusterprots.fasta', False),
    ('geneclusters.txt', False),
    ('plantgeneclusterprots.dmnd', False),
    ('plantgeneclusterprots.fasta', False),
    ('plantgeneclusters.txt', False),
]

def check_prereqs(options):
    "Check if all required applications are around"
    if 'clusterblastdir' not in options:
        options.clusterblastdir = utils.get_full_path(__file__, '')

    failure_messages = []
    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    for file_name, optional in _required_files:
        if utils.locate_file(path.join(options.clusterblastdir, file_name)) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % file_name)

    return failure_messages

def run_clusterblast(seq_record, options):
    logging.info('Running ClusterBlast')
    clusters, proteinlocations, proteinstrands, proteinannotations, proteintags = load_clusterblast_database(seq_record)
    seq_record.internalhomologygroupsdict = internal_homology_blast(seq_record)
    perform_clusterblast(options, seq_record, clusters, proteinlocations, proteinstrands, proteinannotations, proteintags)
    prepare_data(seq_record, options, searchtype="general")
    generate_Storage_for_cb(options, seq_record)
