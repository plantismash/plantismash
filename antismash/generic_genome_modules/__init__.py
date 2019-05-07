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

# This directory contains modules that are executed in genome-wide context


from antismash.generic_genome_modules import (
       metabolicmodel
    )

def check_prereqs(options):
    failure_msgs = []

    if options.modeling != "none":
        failure_msgs.extend(metabolicmodel.check_prereqs())

    return failure_msgs
