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

"""ClusterBlast comparative gene cluster analysis for known clusters"""

import logging
from os import path
from antismash import utils
from antismash.generic_modules.clusterblast.clusterblast import internal_homology_blast, load_clusterblast_database
from antismash.generic_modules.clusterblast.data_loading import prepare_data, generate_Storage_for_cb
from .knownclusterblast import perform_knownclusterblast

# Tuple is ( binary_name, optional)
_required_binaries = [
    ('blastp', False),
    ('makeblastdb', False),
    ('diamond', False),
]

_required_files = [
    ('knownclusterprots.fasta', False),
    ('knownclusterprots.dmnd', False),
    ('knownclusters.txt', False)
]

def check_prereqs(options):
    "Check if all required applications are around"
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    for file_name, optional in _required_files:
        if utils.locate_file(path.join(utils.get_full_path(__file__, ''), file_name)) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % file_name)

    return failure_messages

def run_knownclusterblast(seq_record, options):
    logging.info('Running known cluster search')
    knownclusterblastvars = utils.Storage()
    knownclusterblastvars.internalhomologygroupsdict = {}
    knownclusterblastvars.clusterblastpositiondata = {}
    knownclusterblastvars.queryclusterdata = {}
    clusters, proteinlocations, proteinstrands, proteinannotations, proteintags = load_clusterblast_database(seq_record, searchtype="knownclusters")
    if not options.clusterblast:
        seq_record.internalhomologygroupsdict = internal_homology_blast(seq_record)
    perform_knownclusterblast(options, seq_record, clusters, proteinlocations, proteinstrands, proteinannotations, proteintags)
    prepare_data(seq_record, options, searchtype="knownclusters")
    generate_Storage_for_cb(options, seq_record, searchtype="KnownClusterBlastData")


