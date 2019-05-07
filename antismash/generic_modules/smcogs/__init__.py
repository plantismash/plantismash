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

"""Full genome PFAM anotation"""

import logging
from os import path
import os
from antismash import utils
from smcogs import smcog_analysis
from antismash.lib.hmmscanparser import parse_hmmscan_results
from multiprocessing import Process
import time
from helperlibs.wrappers.io import TemporaryDirectory

name = "smcogs"
short_description = name.capitalize()
priority = 10000


# Tuple is ( binary_name, optional)
_required_binaries = [
    ('muscle', False),
    ('hmmscan', False),
    ('hmmpress', False),
    ('fasttree', False),
    ('java', False)
]

_markov_models = [
    'smcogs.hmm',
]

_binary_extensions = [
    '.h3f',
    '.h3i',
    '.h3m',
    '.h3p'
]



def check_prereqs(options):
    "Check if all required applications are around"
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    for hmm in _markov_models:
        hmm = utils.get_full_path(__file__, hmm)
        if utils.locate_file(hmm) is None:
            failure_messages.append("Failed to locate file %r" % hmm)
            continue
        for ext in _binary_extensions:
            binary = "%s%s" % (hmm, ext)
            if utils.locate_file(binary) is None:
                command = ['hmmpress', hmm]
                try:
                    out, err, retcode = utils.execute(command)
                except OSError as e:
                    retcode = 1
                    err = str(e)
                if retcode != 0:
                    failure_messages.append("Failed to hmmpress %r: %r" % (hmm, err))
                break


    return failure_messages

def run_smcog_analysis(seq_record, options):
    #run_smcog_analysis(opts, globalvars, geneclustervars, pksnrpscoregenes)
    logging.info('Running smCOG analysis')
    smcogvars = utils.Storage()
    smcogvars.smcogtreedict = {}
    smcogvars.smcogdict = {}
    geneclustergenes = utils.get_withincluster_cds_features(seq_record)
    pksnrpscoregenes = utils.get_pksnrps_cds_features(seq_record)
    logging.info("Performing smCOG analysis")
    smcogs_fasta = utils.get_specific_multifasta(geneclustergenes)
    smcogs_opts = ["-E", "1E-6"]
    smcogs_results = utils.run_hmmscan(utils.get_full_path(__file__, "smcogs.hmm"), smcogs_fasta, smcogs_opts)
    hmmlengthsdict = utils.hmmlengths(utils.get_full_path(__file__, "smcogs.hmm"))
    smcogvars.smcogdict = parse_hmmscan_results(smcogs_results, hmmlengthsdict)
    #Write output
    options.smcogsfolder = path.abspath(path.join(options.outputfoldername, "smcogs"))
    if not os.path.exists(options.smcogsfolder):
        os.mkdir(options.smcogsfolder)
    originaldir = os.getcwd()
    os.chdir(options.smcogsfolder)
    smcogfile = open("smcogs.txt","w")
    pksnrpscoregenenames = [utils.get_gene_id(feature) for feature in pksnrpscoregenes]
    for feature in geneclustergenes:
        k = utils.get_gene_id(feature)
        if k not in pksnrpscoregenenames:
            if smcogvars.smcogdict.has_key(k):
                l = smcogvars.smcogdict[k]
                smcogfile.write(">> " + k + "\n")
                smcogfile.write("name\tstart\tend\te-value\tscore\n")
                smcogfile.write("** smCOG hits **\n")
                for i in l:
                    smcogfile.write(str(i[0]) + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + str(i[3]) + "\t" + str(i[4]) + "\n")
                smcogfile.write("\n\n")
    smcogfile.close()
    #smCOG phylogenetic tree construction
    logging.info("Calculating and drawing phylogenetic trees of cluster genes "
        "with smCOG members")
    with TemporaryDirectory(change=True):
        smcoganalysisgenes = []
        for feature in geneclustergenes:
            k = utils.get_gene_id(feature)
            if k not in pksnrpscoregenenames:
                smcoganalysisgenes.append(feature)
        smcogsets = []
        equalpartsizes = int(len(smcoganalysisgenes)/options.cpus)
        for i in range(options.cpus):
            if i == 0:
                geneslist = smcoganalysisgenes[:equalpartsizes]
            elif i == (options.cpus - 1):
                geneslist = smcoganalysisgenes[(i*equalpartsizes):]
            else:
                geneslist = smcoganalysisgenes[(i*equalpartsizes):((i+1)*equalpartsizes)]
            smcogsets.append(geneslist)
        processes = []
        z = 0
        for k in smcogsets:
            processes.append(Process(target=smcog_analysis,
                                     args=[k, z, seq_record,
                                        smcogvars.smcogdict, options.smcogsfolder]))
            z += 1
        for k in processes:
            k.start()
        time.sleep(1)
        while True:
            processrunning = "n"
            for k in processes:
                if k.is_alive():
                    processrunning = "y"
            if processrunning == "y":
                time.sleep(5)
            else:
                break
        for k in processes:
            k.join()
    os.chdir(options.smcogsfolder)
    dircontents = os.listdir(os.getcwd())
    for k in dircontents:
        if ".png" in k:
            tag = k.split(".png")[0]
            smcogvars.smcogtreedict[tag] = tag + ".png"
    os.chdir(originaldir)
    _annotate(geneclustergenes, smcogvars, options)

def _annotate(geneclustergenes, smcogvars, options):
    #Annotate smCOGS in CDS features
    for feature in geneclustergenes:
        gene_id = utils.get_gene_id(feature)
        if smcogvars.smcogdict.has_key(gene_id):
            detailslist = smcogvars.smcogdict[gene_id]
            if not feature.qualifiers.has_key('note'):
                feature.qualifiers['note'] = []
            if len(detailslist) > 0:
                feature.qualifiers['note'].append("smCOG: " + detailslist[0][0] + " (Score: " + str(detailslist[0][4]) + "; E-value: " + str(detailslist[0][3]) + ");")
        if smcogvars.smcogtreedict.has_key(gene_id):
            if not feature.qualifiers.has_key('note'):
                feature.qualifiers['note'] = []
            feature.qualifiers['note'].append("smCOG tree PNG image: smcogs/%s"  % smcogvars.smcogtreedict[gene_id])


