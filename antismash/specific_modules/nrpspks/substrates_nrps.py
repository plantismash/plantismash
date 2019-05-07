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
from antismash import utils
import os
from os import path
import sys
import shutil
from NRPSPredictor2 import nrpscodepred
from Minowa import minowa_A
from helperlibs.wrappers.io import TemporaryDirectory

def extract_nrps_genes(pksnrpscoregenes, domaindict, seq_record, extra_aa=0):
    nrpsnames = [] 
    nrpsseqs = []
    for feature in pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        domaindetails = domaindict[locus]
        nr = 0
        for tab in domaindetails:
            if tab[0] == "AMP-binding" or tab[0] == "A-OX":
                nr += 1
                start = int(tab[1])
                end = int(tab[2]) + extra_aa
                seq = str(utils.get_aa_sequence(feature))[start:end]
                name = locus + "_A" + str(nr)
                nrpsnames.append(name)
                nrpsseqs.append(seq)
    return nrpsnames, nrpsseqs

def run_nrpspredictor(seq_record, nrpsnames, nrpsseqs, options):
    #NRPSPredictor: extract AMP-binding + 120 residues N-terminal of this domain, extract 8 Angstrom residues and insert this into NRPSPredictor
    with TemporaryDirectory(change=True):
        nrpsseqs_file = "nrpsseqs.fasta"
        NRPSPredictor2_dir = utils.get_full_path(__file__, "NRPSPredictor2")
        utils.writefasta(nrpsnames, nrpsseqs, nrpsseqs_file)
        #Get NRPSPredictor2 code predictions, output sig file for input for NRPSPredictor2 SVMs
        nrpscodepred.run_nrpscodepred(options)
        #Run NRPSPredictor2 SVM
        datadir = path.join(NRPSPredictor2_dir, 'data')
        libdir = path.join(NRPSPredictor2_dir, 'lib')
        jarfile = path.join(NRPSPredictor2_dir, 'build', 'NRPSpredictor2.jar')
        classpath = [ jarfile,
                     '%s/java-getopt-1.0.13.jar' % libdir,
                     '%s/Utilities.jar' % libdir,
                     '%s/libsvm.jar' % libdir
                    ]
        if sys.platform == ("linux2") or sys.platform == ("darwin"):
            java_separator = ":"
        elif sys.platform == ("win32"):
            java_separator = ";"
        commands = ['java', '-Ddatadir=%s' % datadir, '-cp', java_separator.join(classpath),
                    'org.roettig.NRPSpredictor2.NRPSpredictor2', '-i', 'input.sig',
                    '-r', path.join(options.raw_predictions_outputfolder, "ctg" + str(options.record_idx) + '_nrpspredictor2_svm.txt'),
                    '-s', '1', '-b', options.eukaryotic and '1' or '0']
        out, err, retcode = utils.execute(commands)
        if err != '':
            logging.debug('running nrpspredictor2 gave error %r' % err)
        #Copy NRPSPredictor results and move back to original directory
        try:
            os.remove(path.join(options.raw_predictions_outputfolder, "ctg" + str(options.record_idx) + "_nrpspredictor2_codes.txt"))
        except:
            pass
        shutil.move("ctg" + str(options.record_idx) + "_nrpspredictor2_codes.txt", options.raw_predictions_outputfolder)

def run_minowa_predictor_nrps(pksnrpscoregenes, domaindict, seq_record, options):
    #Minowa method: extract AMP-binding domain, and run Minowa_A
    logging.info("Predicting NRPS A domain substrate specificities by Minowa " \
        "et al. method")
    nrpsnames2, nrpsseqs2 = extract_nrps_genes(pksnrpscoregenes, domaindict, seq_record, extra_aa=0)
    #Make Minowa output folder
    utils.writefasta(nrpsnames2, nrpsseqs2,
            path.join(options.raw_predictions_outputfolder, "ctg" + str(options.record_idx) + "_nrpsseqs.fasta"))
    with TemporaryDirectory(change=True):
        minowa_A.run_minowa_a(path.join(options.raw_predictions_outputfolder, "ctg" + str(options.record_idx) + "_nrpsseqs.fasta"),
                              path.join(options.raw_predictions_outputfolder, "ctg" + str(options.record_idx) + "_minowa_nrpspredoutput.txt"))
