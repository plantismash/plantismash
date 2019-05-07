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

import logging
from substrates_nrps import *
from substrates_pks import *
from parsers import *
from orderfinder import *

def run_nrps_substr_spec_predictions(pksnrpsvars, seq_record, options):
    #Predict NRPS A domain specificities with NRPSPredictor and Minowa et al. method
    logging.info("Predicting NRPS A domain substrate specificities by " \
        "NRPSPredictor")
    pksnrpsvars.nrpsnames, pksnrpsvars.nrpsseqs = extract_nrps_genes(pksnrpsvars.pksnrpscoregenes, pksnrpsvars.domaindict, seq_record, extra_aa=120)
    if len(pksnrpsvars.nrpsnames) > 0:
        run_nrpspredictor(seq_record, pksnrpsvars.nrpsnames, pksnrpsvars.nrpsseqs, options)
        run_minowa_predictor_nrps(pksnrpsvars.pksnrpscoregenes, pksnrpsvars.domaindict, seq_record, options)

def run_pks_substr_spec_predictions(pksnrpsvars, seq_record, options):
    pksnrpsvars.pksnames, pksnrpsvars.pksseqs = extract_pks_genes(pksnrpsvars.pksnrpscoregenes, pksnrpsvars.domaindict, seq_record)
    if count_pks_genes(pksnrpsvars.pksnrpscoregenes, pksnrpsvars.domaindict, seq_record) > 0:
        run_minowa_predictor_pks_at(pksnrpsvars.pksnames, pksnrpsvars.pksseqs, options)
        pksnrpsvars.calnames, pksnrpsvars.calseqs = run_minowa_predictor_pks_cal(pksnrpsvars.pksnrpscoregenes, pksnrpsvars.domaindict, seq_record, options)
        pksnrpsvars.krnames, pksnrpsvars.krseqs = run_kr_stereochemistry_predictions(pksnrpsvars.pksnrpscoregenes, pksnrpsvars.domaindict, seq_record, options)
    else:
        pksnrpsvars.calnames, pksnrpsvars.calseqs, pksnrpsvars.krnames, pksnrpsvars.krseqs = [], [], [], []

def parse_substr_spec_predictions(pksnrpsvars, options):
    #Read and parse all substrate specificity prediction output files
    parse_nrps_preds(options, pksnrpsvars)
    parse_pks_preds(options, pksnrpsvars)
    parse_kr_activity_preds(options, pksnrpsvars)
