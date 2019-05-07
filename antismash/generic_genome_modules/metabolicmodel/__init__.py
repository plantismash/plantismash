# vim: set fileencoding=utf-8 :
#

#
# Copyright (C) 2014 Tilmann Weber
# The Novo Nordisk Foundation Center for Biosustainability
# Technical University of Denmark
# Section: Metabolic Engineering for Natural Compounds / New Bioactive Compounds
#
# Copyright (C) 2014 Hyun Uk Kim
# The Novo Nordisk Foundation Center for Biosustainability
# Technical University of Denmark
# Section: Metabolic Engineering for Natural Compounds / New Bioactive Compounds
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""
Generate homology metabolic model of input Sequence
"""


import logging
import os

from antismash import utils

try:
    from antismash.generic_genome_modules.metabolicmodel.automodel import run_automodel
    libImport = True
except (ImportError, ImportWarning):
    libImport = False

name = "MetabolicModel"
short_description = "Metabolic Model generates a homology based metabolic model"
priority = 90000


def check_prereqs():
    "Check if all required files and applications are around"

    # Tuple is ( binary_name, optional)
    _required_binaries = [
        ('blastp', False),
    ]

    failure_messages = []

    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    try:
        import cobra
        logging.debug("Found cobra version %s", cobra.__version__)
    except (ImportError, ImportWarning) as err:
        failure_messages.append(str(err))

    try:
        import libsbml
        logging.debug("Found libsmbl version %s", libsbml.getLibSBMLVersion())
    except (ImportError, ImportWarning) as err:
        failure_messages.append(str(err))


    if not libImport:
        failure_messages.append("Failed to import automodel")

    return failure_messages


def run_analyses(seq_record, options):
    "Wrapper to calculate metabolic model"

    if options.modeling == "none":
        return False


    model_dirname = "metabolicModel"

    # set up output directory
    basename = options.outputfoldername
    options.metabolicmodeldir = os.path.join(basename, model_dirname)
    logging.debug("Writing metabolic models to %r", options.metabolicmodeldir)
    if not os.path.exists(options.metabolicmodeldir):
        os.mkdir(options.metabolicmodeldir)

    result = run_automodel(seq_record, options)

    if not result:
        options.modeling = "none"
        options.modeling_successful = False
    else:
        options.modeling_successful = True

    return result
