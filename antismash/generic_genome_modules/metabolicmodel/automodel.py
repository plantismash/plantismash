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
Automodel core pipeline script
"""

from .prunPhase import (
    get_targetGenomeInfo,
    make_blastDB,
    run_blastp,
    parseBlaspResults,
    makeBestHits_dict,
    getBBH,
    get_nonBBH,
    labelRxnToRemove,
    pruneModel,
    swap_locusTag_tempModel
)

from .augPhase import (
    get_targetGenome_locusTag_ec_nonBBH_dict,
    make_all_rxnInfo_fromRefSeq,
    get_mnxr_list_from_modelPrunedGPR,
    check_existing_rxns,
    get_mnxr_using_kegg,
    extract_rxn_mnxm_coeff,
    add_nonBBH_rxn
)
import cobra
import pickle
import os
import logging

from antismash import utils

def run_automodel(seq_records, options):


    #List of input (static) files as pickles
    ###################################################################
    #For model pruning phase
    #Choose "eco" or "sco"
    #root = os.path.dirname(utils.get_full_path('__file__', 'input1'))
    if not cobra.__version__ == "0.2.1":
        logging.error("The modeling pipeline is only compatible wit COBRApy version 0.2.1; your insstalled version is %s",
                      cobra.__version__)
        return False

    root = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'input1'
    temp_fasta = options.metabolicmodeldir

    #Get template model-specific pickles
    logging.debug("[metabolicmodel] set up model dir as %s and output dir as %s", root, temp_fasta)

    model = pickle.load(open(root+os.sep+options.modeling+os.sep+'model.p','rb'))
    tempModel_biggRxnid_locusTag_dict = pickle.load(open(root+os.sep+options.modeling+os.sep+'tempModel_biggRxnid_locusTag_dict.p','rb'))
    tempModel_exrxnid_flux_dict = pickle.load(open(root+os.sep+options.modeling+os.sep+'tempModel_exrxnid_flux_dict.p','rb'))

    #Template model-independent pickles for model augmentation phase
    logging.debug("loading pickle files of the parsed template model and its relevant genbank data..")
    bigg_mnxr_dict = pickle.load(open(root+os.sep+'bigg_mnxr_dict.p','rb'))
    kegg_mnxr_dict = pickle.load(open(root+os.sep+'kegg_mnxr_dict.p','rb'))
    mnxr_kegg_dict = pickle.load(open(root+os.sep+'mnxr_kegg_dict.p','rb'))
    mnxr_rxn_dict = pickle.load(open(root+os.sep+'mnxr_rxn_dict.p','rb'))
    bigg_mnxm_compound_dict = pickle.load(open(root+os.sep+'bigg_mnxm_compound_dict.p','rb'))
    mnxm_bigg_compound_dict = pickle.load(open(root+os.sep+'mnxm_bigg_compound_dict.p','rb'))
    kegg_mnxm_compound_dict = pickle.load(open(root+os.sep+'kegg_mnxm_compound_dict.p','rb'))
    mnxm_kegg_compound_dict = pickle.load( open(root+os.sep+'mnxm_kegg_compound_dict.p','rb'))
    mnxm_compoundInfo_dict = pickle.load(open(root+os.sep+'mnxm_compoundInfo_dict.p','rb'))
    ###################################################################

    logging.debug("pruning phase starting..")
    ###################################################################

    logging.debug("reading genbank file of the target genome.."    )
    targetGenome_locusTag_ec_dict, targetGenome_locusTag_prod_dict, target_fasta = get_targetGenomeInfo(seq_records, options)

    if len(targetGenome_locusTag_ec_dict) == 0:
        logging.error("Error: no EC_number in sequence record; skipping modeling")
        return False

    logging.debug("generating a DB for the genes from the target genome..")
    make_blastDB(query_fasta=target_fasta, options=options)

    logging.debug("running BLASTP #1: genes in the target genome against genes in the template model..")
    run_blastp(target_fasta=options.metabolicmodeldir+os.sep+'targetGenome_locusTag_aaSeq.fa', \
               blastp_result=options.metabolicmodeldir+os.sep+'blastp_targetGenome_against_tempGenome.txt',\
               db_dir=root+os.sep+options.modeling+os.sep+'tempBlastDB', evalue=1e-30)

    logging.debug("running BLASTP #2: genes in the template model against genes in the target genome..")
    run_blastp(target_fasta=root+os.sep+options.modeling+os.sep+'tempModel_locusTag_aaSeq.fa', \
               blastp_result=options.metabolicmodeldir+os.sep+'blastp_tempGenome_against_targetGenome.txt',\
               db_dir = options.metabolicmodeldir+os.sep+'targetBlastDB', evalue=1e-30)

    logging.debug("parsing the results of BLASTP #1..")
    blastpResults_dict1 = parseBlaspResults(options.metabolicmodeldir+os.sep+'blastp_targetGenome_against_tempGenome.txt', \
                                            options.metabolicmodeldir+os.sep+'blastp_targetGenome_against_tempGenome_parsed.txt')

    logging.debug("parsing the results of BLASTP #2..")
    blastpResults_dict2 = parseBlaspResults(options.metabolicmodeldir+os.sep+'blastp_tempGenome_against_targetGenome.txt', \
                                            options.metabolicmodeldir+os.sep+'blastp_tempGenome_against_targetGenome_parsed.txt')

    logging.debug("selecting the best hits for BLASTP #1..")
    bestHits_dict1 = makeBestHits_dict(options.metabolicmodeldir+os.sep+'blastp_targetGenome_against_tempGenome_parsed.txt')

    logging.debug("selecting the best hits for BLASTP #2..")
    bestHits_dict2 = makeBestHits_dict(options.metabolicmodeldir+os.sep+'blastp_tempGenome_against_targetGenome_parsed.txt')

    logging.debug("selecting the bidirectional best hits..")
    targetBBH_list, temp_target_BBH_dict = getBBH(bestHits_dict1, bestHits_dict2)


    logging.debug("selecting genes that are not bidirectional best hits..")
    nonBBH_list = get_nonBBH(targetGenome_locusTag_ec_dict, targetBBH_list)
    ###################################################################

    ###################################################################
    logging.debug("labeling reactions with nonhomologous genes to remove from the template model..")
    rxnToRemove_dict = labelRxnToRemove(model, temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict)

    logging.debug("removing reactions with nonhomologous genes from the template model..")
    modelPruned, rxnToRemoveEssn_dict, rxnRemoved_dict, rxnRetained_dict = pruneModel(model, rxnToRemove_dict, options.automodel.solver)

    logging.debug("correcting GPR associations in the template model..")
    modelPrunedGPR = swap_locusTag_tempModel(modelPruned, temp_target_BBH_dict)
    ###################################################################


    logging.debug("augmentation phase starting..")
    ###################################################################
    logging.debug("creating various dictionary files for the nonBBH gene-associted reactions...")

    targetGenome_locusTag_ec_nonBBH_dict = get_targetGenome_locusTag_ec_nonBBH_dict(targetGenome_locusTag_ec_dict, nonBBH_list)

    rxnid_info_dict, rxnid_locusTag_dict = make_all_rxnInfo_fromRefSeq(targetGenome_locusTag_ec_nonBBH_dict, options)

    modelPrunedGPR_mnxr_list = get_mnxr_list_from_modelPrunedGPR(modelPrunedGPR, bigg_mnxr_dict)
    ###################################################################

    ###################################################################
    logging.debug("adding the nonBBH gene-associated reactions...")
    rxnid_to_add_list = check_existing_rxns(kegg_mnxr_dict, modelPrunedGPR_mnxr_list, rxnid_info_dict)

    mnxr_to_add_list = get_mnxr_using_kegg(rxnid_to_add_list, kegg_mnxr_dict)

    rxnid_mnxm_coeff_dict = extract_rxn_mnxm_coeff(mnxr_to_add_list, mnxr_rxn_dict, mnxm_bigg_compound_dict, mnxm_kegg_compound_dict, mnxr_kegg_dict)

    target_model = add_nonBBH_rxn(modelPrunedGPR, rxnid_info_dict, rxnid_mnxm_coeff_dict, rxnid_locusTag_dict, bigg_mnxm_compound_dict, kegg_mnxm_compound_dict, mnxm_compoundInfo_dict, targetGenome_locusTag_prod_dict, tempModel_exrxnid_flux_dict, options)
    ###################################################################

    #Output on screen
    model = pickle.load(open(root+os.sep+options.modeling+os.sep+'model.p','rb'))
    logging.debug("Number of genes in template and pruned models: %s / %s", len(model.genes), len(modelPruned.genes))
    logging.debug("Number of reactions in template and pruned models: %s / %s", len(model.reactions), len(modelPruned.reactions))
    logging.debug("Number of metabolites in template and pruned models: %s / %s", len(model.metabolites), len(modelPruned.metabolites))

    # Set up extrarecord data structure within options, if not already set
    if "extrarecord" not in options:
        options.extrarecord = {}

    # store model data in seq_records[0]
    seq_record = seq_records[0]
    if seq_record.id not in options.extrarecord:
        options.extrarecord[seq_record.id] = utils.Storage()
    if "extradata" not in options.extrarecord[seq_record.id]:
        options.extrarecord[seq_record.id].extradata = {}

    # as the cobra model object does not provide an own serialization, let's try with pickle...
    options.extrarecord[seq_record.id].extradata["MetabolicModelDataObj"] = pickle.dumps(target_model)

    if 'MetabolicModelDataObj' in options.extrarecord[seq_record.id].extradata:
        logging.debug("Generate options.extrarecord entry")
    else:
        logging.warning("Could not generate options.extrarecord for %s", seq_record.id)

    return True
