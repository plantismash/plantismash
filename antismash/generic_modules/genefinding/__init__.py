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

"""Gene finding using Glimmer / GlimmerHMM

"""

from os import path
from os import listdir
from antismash import utils
from antismash.utils import locate_executable
from .all_orfs import find_all_orfs
from .run_glimmerhmm import run_glimmerhmm
from .run_glimmer import run_glimmer
from .run_prodigal import run_prodigal

name = "genefinding"
short_description = name.capitalize()
priority = 1000

_required_binaries = [
    ('long-orfs', False),
    ('extract', False),
    ('build-icm', False),
    ('glimmer3', False),
    ('glimmerhmm', False),
    ('prodigal', False)
]

def find_genes(seq_record, options):
    "Find genes in a seq_record"
    if options.genefinding == 'none':
        raise ValueError("Called find_genes, but genefinding disabled")
    if options.eukaryotic:
        if options.all_orfs:
            find_all_orfs(seq_record, options)
        run_glimmerhmm(seq_record, options)
    else:
        if options.all_orfs:
            find_all_orfs(seq_record, options)
        if options.genefinding == "prodigal" or options.genefinding == "prodigal-m":  #PRODIGAL CHANGE IMPLEMENTED HERE
            run_prodigal(seq_record, options)
        else:
            run_glimmer(seq_record, options)

def check_prereqs(options):
    "Check if all required applications are around"
    failure_messages = []
    if options.genefinding == 'none':
        return failure_messages
    basedir = utils.get_genefinding_basedir(options)
    for binary_name, optional in _required_binaries:
        new_binary_name = path.join(basedir, binary_name)
        if locate_executable(new_binary_name) is None and not optional:
            if locate_executable(binary_name) is None and not optional:
                failure_messages.append("Failed to locate executable for %r" % binary_name)

    return failure_messages


def get_supported_training_models():
    "Get a list of all supported training models"
    result = []
    for fname in listdir(path.dirname(path.abspath(__file__))):
        if path.isdir(path.join(path.dirname(path.abspath(__file__)), fname)):
            if fname.startswith("train_"):
                result.append(fname)
    return result
