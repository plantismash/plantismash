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
Automodel augmentation phase
"""

from Bio import SeqIO
from cobra import Reaction, Metabolite
from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file
from cobra.manipulation.delete import prune_unused_metabolites

import copy
import os
import sys
try:
    import cPickle as pickle
except ImportError:
    import pickle

import re
import urllib2
import logging


#Retrieve a list of reaction IDs using their EC numbers from KEGG
#Input: E.C number in string form (e.g., 4.1.3.6)
#Output: reactionID in list form (e.g., ['R00362'])
def get_rxnid_from_ECNumber(enzymeEC):
    url = "http://rest.kegg.jp/get/enzyme:%s"%(enzymeEC)
    try:
        ecinfo_text = urllib2.urlopen(url).read()
    except (urllib2.URLError, urllib2.HTTPError):
        logging.exception("can't connect to KEGG with API query %s", url)
        # FIXME: modules should never call sys.exit()
        sys.exit(1)

    #Original line also extracted genes in other organisms: R50912; R50345 (NOT rxnid)
    #The HTTP error was solved by putting "\\b" only at the end (not at the front) in order to also
    #include reaction ID followed by "other" in KEGG
    rxnid_set = re.findall(r'\s+R[0-9]+[0-9]+[0-9]+[0-9]+[0-9]'+'\\b', ecinfo_text)
    rxnid_list = []
    for each_set in rxnid_set:
        rxnid = re.findall('R[0-9]+[0-9]+[0-9]+[0-9]+[0-9]+', each_set)
        rxnid_list+=rxnid

        #Removes redundancy
        #FIXME: Check runtime implications
        rxnid_list = list(set(rxnid_list))
    return rxnid_list


#Get reaction information using its ID from KEGG
#Input: KEGG rxnid in string form (e.g., R00362)
#Output: Reaction information for 'Name', 'Definition', and 'Equation' as dictionary
#form {'NAME': 'citrate oxaloacetate-lyase', 'DEFINITION': Citrate <=> Acetate + Oxaloacetate,
#      'EQUATION': C00158 <=> C00033 + C00036}
def get_rxnInfo_from_rxnid(rxnid):
    url = "http://rest.kegg.jp/get/rn:%s"%(rxnid)
    try:
        reaction_info_text = urllib2.urlopen(url).read()
    except urllib2.URLError:
        logging.exception("can't connect to KEGG with API query %s", url)
        sys.exit(1)
    split_text = reaction_info_text.strip().split('\n')
    NAME = ''
    DEFINITION = ''
    EQUATION = ''
    ENZYME = ''
    PATHWAY = ''

    for line in split_text:
        sptlist = line.split()
        if sptlist[0].strip() == 'NAME':
            NAME = ' '.join(sptlist[1:])
        if sptlist[0].strip() == 'DEFINITION':
            DEFINITION = ' '.join(sptlist[1:])
        if sptlist[0].strip() == 'EQUATION':
            EQUATION = ' '.join(sptlist[1:])
        if sptlist[0].strip() == 'ENZYME':
            ENZYME = ' '.join(sptlist[1:])

        #Considers only reactions mapped in pathways
        #Otherwise reactions have unspecified molecules having R groups
        if sptlist[0].strip() == 'PATHWAY':
            PATHWAY = ' '.join(sptlist[1:])
            return {'NAME':NAME, 'DEFINITION':DEFINITION, 'EQUATION':EQUATION, 'ENZYME':ENZYME, 'PATHWAY':PATHWAY}



def get_targetGenome_locusTag_ec_nonBBH_dict(targetGenome_locusTag_ec_dict, nonBBH_list):
    targetGenome_locusTag_ec_nonBBH_dict = {}

    for locusTag in nonBBH_list:
        if locusTag in targetGenome_locusTag_ec_dict.keys():
            targetGenome_locusTag_ec_nonBBH_dict[locusTag] = targetGenome_locusTag_ec_dict[locusTag]
    return targetGenome_locusTag_ec_nonBBH_dict


#Two nested function calling four functions above
def make_all_rxnInfo_fromRefSeq(targetGenome_locusTag_ec_nonBBH_dict, options):
    rxnid_info_dict ={}
    rxnid_locusTag_dict = {}

    logging.debug("Querying KEGG to get reaction ID's for %s gene products..., this step can take a while",
                  len(targetGenome_locusTag_ec_nonBBH_dict.keys()))

    KEGG_EC_rxnCache = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'input1' + os.sep + 'KEGG_EC_rxnCache.p'
    KEGG_rxnInfoCache = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'input1' + os.sep + 'KEGG_rxnInfoCache.p'

    if "automodel" in options:
        if "ecrxncachefile" in options.automodel:
            KEGG_EC_rxnCache = options.automodel.ecrxncachefile

        if "rxninfocache" in options.automodel:
            KEGG_rxnInfoCache = options.automodel.rxninfocache

    # Prepare cache of EC-rxn hash
    ECrxnMapping = {}
    try:
        with open(KEGG_EC_rxnCache, "rb") as fh:
            ECrxnMapping = pickle.load(fh)

    except pickle.UnpicklingError as e:
        logging.warning("could not read KEGG EC-rxn cache file %s: %s", KEGG_EC_rxnCache, e)
    except IOError as e:
        logging.warning("Can't open %s in rb mode: %s", KEGG_EC_rxnCache, e)

    # Prepare cache for rxn info
    rxnid_info_dictCache = {}
    try:
        with open(KEGG_rxnInfoCache, "rb") as fh:
            rxnid_info_dictCache = pickle.load(fh)
    except pickle.UnpicklingError as e:
        logging.warning("could not read KEGG rxn-Info cache file %s: %s", KEGG_rxnInfoCache, e)
    except IOError as e:
        logging.warning("Can't open %s in rb mode: %s", KEGG_rxnInfoCache, e)

    for locusTag in targetGenome_locusTag_ec_nonBBH_dict.keys():
        for enzymeEC in targetGenome_locusTag_ec_nonBBH_dict[locusTag]:
            logging.debug("EC_number for locusTag: %s: %s", locusTag, enzymeEC)

            #KEGG REST does not accept unspecific EC_number: e.g., 3.2.2.-
            if '-' not in enzymeEC:
                if enzymeEC in ECrxnMapping:
                    rxnid_list = ECrxnMapping[enzymeEC]
                    if len(rxnid_list) == 0:
                        logging.debug("ECrxnmapping for %s resulted in empty list", enzymeEC)

                else:
                    rxnid_list = get_rxnid_from_ECNumber(enzymeEC)
                    ECrxnMapping[enzymeEC] = rxnid_list
                    logging.debug("%s: %s", enzymeEC, ",".join(rxnid_list))

                for rxnid in rxnid_list:

                    if rxnid in rxnid_info_dictCache:
                        rxnid_info_dict[rxnid] = rxnid_info_dictCache[rxnid]
                        #if rxnid_info_dict[rxnid] == None:
                        #    raise TypeError("NoneType retrieved from cache for rxnid %s" % rxnid)
                        logging.debug("getting reaction information for %s from local cache", rxnid)
                    else:
                        logging.debug("querying KEGG to get reaction information for %s", rxnid)
                        rxnid_info_dict[rxnid] = get_rxnInfo_from_rxnid(rxnid)

                        # Don't store NoneTypes in Cache dictionary
                        if rxnid_info_dict[rxnid]:
                            rxnid_info_dictCache[rxnid] = rxnid_info_dict[rxnid]
                        else:
                            logging.warn("Skipping rxnid %s for adding to cache as it equals None", rxnid)

                    if rxnid not in rxnid_locusTag_dict.keys():
                        rxnid_locusTag_dict[rxnid] = [(locusTag)]
                    #Append additional different genes to the same reaction ID
                    else:
                        rxnid_locusTag_dict[rxnid].append((locusTag))

    # There could be a file locking conflict, if two processes simultaneously want to write...
    # So if write fails just ignore it => the data will be retrieved again
    try:
        with open(KEGG_EC_rxnCache,"wb") as fh:
            pickle.dump(ECrxnMapping, fh)
    except pickle.PicklingError as e:
        logging.warning('Error in serializing ECrxnMapping data; error: %s', e)
    except IOError as e:
        logging.warning("Can't open %s in wb mode: %s", KEGG_EC_rxnCache, e)

    try:
        with open(KEGG_rxnInfoCache,"wb") as fh:
            pickle.dump(rxnid_info_dictCache, fh)
    except pickle.PicklingError as e:
        logging.warning('Error in serializing rxnInfoCache data; error: %s', e)
    except IOError as e:
        logging.warning("Can't open %s in wb mode: %s", KEGG_rxnInfoCache, e)

    return rxnid_info_dict, rxnid_locusTag_dict


#Output: a list of MNXRs available in the modelPrunedGPR
def get_mnxr_list_from_modelPrunedGPR(modelPrunedGPR, bigg_mnxr_dict):
    modelPrunedGPR_mnxr_list = []

    index_last = len(modelPrunedGPR.reactions)
    index_last = index_last - 1
    index = 0

    while index <= index_last:
        rxn = modelPrunedGPR.reactions[index].id
        if rxn in bigg_mnxr_dict.keys():
            modelPrunedGPR_mnxr_list.append(bigg_mnxr_dict[rxn])
        index+=1

    return modelPrunedGPR_mnxr_list


def check_existing_rxns(kegg_mnxr_dict, modelPrunedGPR_mnxr_list, rxnid_info_dict):
    rxnid_to_add_list =[]

    for rxnid in rxnid_info_dict.keys():
        #Consider only reactions mapped in pathways
        if rxnid in kegg_mnxr_dict.keys():
            kegg_mnxr = kegg_mnxr_dict[rxnid]

            #Check with reactions in the template model through MNXref
            if kegg_mnxr not in modelPrunedGPR_mnxr_list and rxnid not in rxnid_to_add_list:
                rxnid_to_add_list.append(rxnid)

    rxnid_to_add_list = list(sorted(set(rxnid_to_add_list)))
    return rxnid_to_add_list


#Output: MNXR for the reactions to add, converted from KEGG rxnid
def get_mnxr_using_kegg(rxnid_to_add_list, kegg_mnxr_dict):
    mnxr_to_add_list = []
    for rxnid in rxnid_to_add_list:
        if rxnid in kegg_mnxr_dict.keys():
            mnxr_to_add_list.append(kegg_mnxr_dict[rxnid])
    return mnxr_to_add_list


#Check multiple presence of a metabolite as substrates or products
def get_correct_metab_coeff(converted_metab_id, metab_coeff, metab_type, mnxm_coeff_dict, mnxm_metab_list):

    #If the same metabolite appears multiple times as either substrates or products,
    #their stoichiometric coeff's are all added up
    if converted_metab_id in mnxm_metab_list:
        overlap_metab_coeff = float(mnxm_coeff_dict[converted_metab_id])
        mnxm_coeff_dict[converted_metab_id] = overlap_metab_coeff+float(metab_coeff)*-1
    else:
        if metab_type == 'substrate':
            mnxm_coeff_dict[converted_metab_id] = float(metab_coeff)*-1
        elif metab_type == 'product':
            mnxm_coeff_dict[converted_metab_id] = float(metab_coeff)

    return mnxm_coeff_dict


#Check if the same metabolite appears as a substrate and a product
def check_overlap_subs_prod(mnxm_subs_list, mnxm_prod_list):

    for each_substrate in mnxm_subs_list:
        if each_substrate in mnxm_prod_list:
            overlap_check = True
            break
        else:
            overlap_check = False

    return overlap_check


#Creating: e.g., {'R03232': {'f1p': -1.0, 'C04261': 1.0, 'fru': 1.0, 'C00615': -1.0}}
#Metabolites are presented primarily with bigg, otherwise with KEGG
def extract_rxn_mnxm_coeff(mnxr_to_add_list, mnxr_rxn_dict, mnxm_bigg_compound_dict, mnxm_kegg_compound_dict, mnxr_kegg_dict):
    rxnid_mnxm_coeff_dict = {}
    mnxm_coeff_dict = {}

    for mnxr in mnxr_to_add_list:
        unparsed_equation = mnxr_rxn_dict[mnxr]
        logging.debug(unparsed_equation)

        #"substrates" and "products" contain stoichiometric coeff of each compound
        sptReaction = unparsed_equation.split('=')
        substrates = sptReaction[0].strip()
        products = sptReaction[1].strip()

        #Discard polymerization reactions with undefinite coeff's
        #e.g., 1 MNXM9 + (n+2) MNXM90033 = 1 MNXM5617 + (n) MNXM90033
        if '(' not in substrates and '(' not in products:
            #Creating: e.g., {bigg compoundID:(-1)coeff}, {kegg compoundID:(-1)coeff} or {mnxm:(-1)coeff}
            substrates = substrates.split(' + ')
            mnxm_coeff_dict = {}
            mnxm_subs_list = []
            mnxm_prod_list = []

            for substrate in substrates:
                metab_type = 'substrate'
                substrate = substrate.split()

                if substrate[1] in mnxm_bigg_compound_dict.keys():
                    mnxm_coeff_dict = get_correct_metab_coeff(mnxm_bigg_compound_dict[substrate[1]], substrate[0], metab_type, mnxm_coeff_dict, mnxm_subs_list)
                    mnxm_subs_list.append(mnxm_bigg_compound_dict[substrate[1]])

                elif substrate[1] in mnxm_kegg_compound_dict.keys():
                    mnxm_coeff_dict = get_correct_metab_coeff(mnxm_kegg_compound_dict[substrate[1]], substrate[0], metab_type, mnxm_coeff_dict, mnxm_subs_list)
                    mnxm_subs_list.append(mnxm_kegg_compound_dict[substrate[1]])

                else:
                    mnxm_coeff_dict = get_correct_metab_coeff(substrate[1], substrate[0], metab_type, mnxm_coeff_dict, mnxm_subs_list)
                    mnxm_subs_list.append(substrate[1])

            #Creating: e.g., {bigg compoundID:coeff}, {kegg compoundID:coeff} or {mnxm:coeff}
            products = products.split(' + ')
            for product in products:
                metab_type = 'product'
                product = product.split()

                if product[1] in mnxm_bigg_compound_dict.keys():
                    mnxm_coeff_dict = get_correct_metab_coeff(mnxm_bigg_compound_dict[product[1]], product[0], metab_type, mnxm_coeff_dict, mnxm_prod_list)
                    mnxm_prod_list.append(mnxm_bigg_compound_dict[product[1]])

                elif product[1] in mnxm_kegg_compound_dict.keys():
                    mnxm_coeff_dict = get_correct_metab_coeff(mnxm_kegg_compound_dict[product[1]], product[0], metab_type, mnxm_coeff_dict, mnxm_prod_list)
                    mnxm_prod_list.append(mnxm_kegg_compound_dict[product[1]])
                else:
                    mnxm_coeff_dict = get_correct_metab_coeff(product[1], product[0], metab_type, mnxm_coeff_dict, mnxm_prod_list)
                    mnxm_prod_list.append(product[1])

            #Check overlapping metabolites as a substrate and a product
            #e.g., ATP + ADP <=> ADP + ATP
            overlap_check = check_overlap_subs_prod(mnxm_subs_list, mnxm_prod_list)
            if overlap_check == True:
                continue
            else:
                #Creating:
                #e.g., {'R03232': {'f1p': -1.0, 'C04261': 1.0, 'fru': 1.0, 'C00615': -1.0}}
	        rxnid_mnxm_coeff_dict[mnxr_kegg_dict[mnxr]] = mnxm_coeff_dict

    return rxnid_mnxm_coeff_dict


def add_nonBBH_rxn(modelPrunedGPR, rxnid_info_dict, rxnid_mnxm_coeff_dict, rxnid_locusTag_dict, bigg_mnxm_compound_dict, kegg_mnxm_compound_dict, mnxm_compoundInfo_dict, targetGenome_locusTag_prod_dict, tempModel_exrxnid_flux_dict, options):

    for rxnid in rxnid_mnxm_coeff_dict.keys():
        logging.debug("%s being examined for a new nonBBH reaction..", rxnid)
        if rxnid_info_dict[rxnid] != None:
            #ID
	    rxn = Reaction(rxnid)

            #Name
            #Some reaction IDs do not have NAME despite the presence of PATHWAY
            rxn.name = rxnid_info_dict[rxnid]['NAME']

            #Reversibility / Lower and upper bounds
            rxn.lower_bound = -1000
            rxn.uppwer_bound = 1000

            #Metabolites and their stoichiometric coeff's
            for metab in rxnid_mnxm_coeff_dict[rxnid]:
                metab_compt = '_'.join([metab,'c'])

                #Add metabolites already in the model
                if metab_compt in modelPrunedGPR.metabolites:
                    rxn.add_metabolites({modelPrunedGPR.metabolites.get_by_id(metab_compt):rxnid_mnxm_coeff_dict[rxnid][metab]})

                #Add metabolites with bigg compoundID, but not in the model
                elif metab in bigg_mnxm_compound_dict.keys():
                    mnxm = bigg_mnxm_compound_dict[metab]
                    metab_compt = Metabolite(metab, formula = mnxm_compoundInfo_dict[mnxm][1], name = mnxm_compoundInfo_dict[mnxm][0], compartment='c')
                    rxn.add_metabolites({metab_compt:rxnid_mnxm_coeff_dict[rxnid][metab]})

                #Add metabolites with MNXM and not in the model
                else:
                    mnxm = kegg_mnxm_compound_dict[metab]
                    metab_compt = Metabolite(mnxm, formula = mnxm_compoundInfo_dict[mnxm][1], name = mnxm_compoundInfo_dict[mnxm][0], compartment='c')
                    rxn.add_metabolites({metab_compt:rxnid_mnxm_coeff_dict[rxnid][metab]})

            #GPR association
            if rxnid not in rxnid_locusTag_dict:
                logging.warning("rxnid_locusTag_dict error: has no key %s; content: %s", rxnid, str(rxnid_locusTag_dict))
            if len(rxnid_locusTag_dict[rxnid]) == 1:
                gpr = '( %s )' %(rxnid_locusTag_dict[rxnid][0])
            else:
                count = 1
                for locusTag in rxnid_locusTag_dict[rxnid]:

                    #Check the submitted gbk file contains "/product" for CDS
                    if locusTag in targetGenome_locusTag_prod_dict:

                        #Consider "and" relationship in the GPR association
                        if 'subunit' in targetGenome_locusTag_prod_dict[locusTag]:
                            count += 1
                if count == len(rxnid_locusTag_dict[rxnid]):
                    gpr = ' and '.join(rxnid_locusTag_dict[rxnid])
                else:
                    gpr = ' or '.join(rxnid_locusTag_dict[rxnid])
                gpr = '( %s )' %(gpr)
            rxn.add_gene_reaction_rule(gpr)

            #Subsystem
            rxn.subsystem = rxnid_info_dict[rxnid]['PATHWAY']

            #E.C. number: not available feature in COBRApy
            #Objective coeff: default
            rxn.objective_coefficient = 0

            #Add a reaction to the model if it does not affect Exchange reaction flux direction
            modelPrunedGPR.add_reaction(rxn)

            write_cobra_model_to_sbml_file(modelPrunedGPR, options.metabolicmodeldir+os.sep+"modelPrunedGPR.xml")
            modelPrunedGPR = create_cobra_model_from_sbml_file(options.metabolicmodeldir+os.sep+"modelPrunedGPR.xml")

            target_exrxnid_flux_dict = get_exrxnid_flux(modelPrunedGPR, tempModel_exrxnid_flux_dict)
            exrxn_flux_change_list = check_exrxn_flux_direction(tempModel_exrxnid_flux_dict, target_exrxnid_flux_dict)

            if 'F' in exrxn_flux_change_list:
	        modelPrunedGPR.remove_reactions(rxn)

                write_cobra_model_to_sbml_file(modelPrunedGPR, options.metabolicmodeldir+os.sep+"modelPrunedGPR.xml")
                modelPrunedGPR = create_cobra_model_from_sbml_file(options.metabolicmodeldir+os.sep+"modelPrunedGPR.xml")

    prune_unused_metabolites(modelPrunedGPR)
    target_model = copy.deepcopy(modelPrunedGPR)
    return target_model


#Output: a dictionary file for major Exchange reactions {Exchange reaction ID:flux value}
def get_exrxnid_flux(model, tempModel_exrxnid_flux_dict):

    target_exrxnid_flux_dict = {}
    model.optimize()

    for exrxn_id in tempModel_exrxnid_flux_dict.keys():
        if exrxn_id in model.solution.x_dict:
            target_exrxnid_flux_dict[exrxn_id] = model.solution.x_dict[exrxn_id]
        else:
            continue
    return target_exrxnid_flux_dict


#Output: a list file having either T or F for major Exchange reactions
def check_exrxn_flux_direction(tempModel_exrxnid_flux_dict, target_exrxnid_flux_dict):

    exrxn_flux_change_list = []

    for exrxn_id in tempModel_exrxnid_flux_dict.keys():
        if exrxn_id in target_exrxnid_flux_dict.keys():
            template_exrxn_flux = tempModel_exrxnid_flux_dict[exrxn_id]
            target_exrxn_flux = target_exrxnid_flux_dict[exrxn_id]
            ratio_exrxn_flux = float(target_exrxn_flux)/float(template_exrxn_flux)

            #Similar species are allowed to uptake nutrients within a decent range
            if float(target_exrxn_flux)*float(template_exrxn_flux) > 0.0 and 0.2 < ratio_exrxn_flux and ratio_exrxn_flux < 2.0:
                exrxn_flux_change_list.append('T')

            #Cause drastic changes in Exchange reaction fluxes (direction and/or magnitude)
            else:
                exrxn_flux_change_list.append('F')

    return exrxn_flux_change_list

