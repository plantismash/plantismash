# vim: set fileencoding=utf-8 :
#
#
# Copyright (C) 2014 Tilmann Weber, Hyun Uk Kim
# The Novo Nordisk Foundation Center for Biosustainability
# Technical University of Denmark
# Section: Metabolic Engineering for Natural Compounds / New Bioactive Compounds
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Predict EC number for CDS features
"""

import logging
import straight.plugin
from antismash import utils

name = "ecpredictor"
short_description = name.capitalize()
priority = 10000

def run(seq_record, options):
    logging.debug('Predicting EC numbers')
    ecpred_plugins = list(straight.plugin.load('antismash.generic_modules.ecpredictor'))
    
    logging.debug('ECpredictor: Found plugins: %s' % ", ".join([plugin.short_description for plugin in ecpred_plugins]))
    
    for ecpred_plugin in ecpred_plugins:
        ecpred_plugin.getECs(seq_record, options)
    
