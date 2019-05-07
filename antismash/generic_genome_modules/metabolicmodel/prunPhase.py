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
Automodel pruning phase
"""

from cobra.flux_analysis import single_deletion
from cobra.manipulation.delete import prune_unused_metabolites
import copy
import os

import sys
from antismash import utils
import logging


#Looks for .xml and .gb(k) files in the pre-defined folder
def get_temp_fasta(orgName):
    for root, _, files in os.walk('./input1/%s/' %(orgName)):
        for f in files:
            if f.endswith('.fa'):
                tempFasta = os.path.join(root, f)
                return root, tempFasta
            else:
                #FIXME: Don't use sys.exit
                sys.exit(1)


def get_targetGenomeInfo(seq_records, options):
    targetFastaFilename = options.metabolicmodeldir + os.sep + 'targetGenome_locusTag_aaSeq.fa'
    fp = open(targetFastaFilename, 'w')
    targetGenome_locusTag_aaSeq_dict = {}
    targetGenome_locusTag_ec_dict = {}
    targetGenome_locusTag_prod_dict = {}

    counter_for_temp_locusTags = 1
    # Reads GenBank file
    for seq_record in seq_records:
        logging.debug('[MetabolicModeling] processing sequence id "%s" out of %s sequences', seq_record.id, len(seq_records))
        for feature in utils.get_cds_features(seq_record):

            # Retrieving "locus_tag (i.e., ORF name)" for each CDS
            locusTag = feature.qualifiers.get('locus_tag',['-'])[0]
            logging.debug("Found locus_tag %s for feature %s", locusTag, utils.get_gene_id(feature))
            # Assign own locus tag, if not set:
            if locusTag == "-":
                if not utils.get_gene_id(feature) == "no_tag_found":
                    locusTag = utils.get_gene_id(feature)
                else:
                    locusTag = "automodelorf{0:05d}".format(counter_for_temp_locusTags)
                feature.qualifiers['locus_tag'] = [locusTag]
                logging.debug("replaced locus tag to %s for %s.", locusTag, utils.get_gene_id(feature))
                counter_for_temp_locusTags += 1

            # Some locus_tag's have multiple same qualifiers (e.g., EC_number)
            for item in feature.qualifiers:

                # Note that the numbers of CDS and "translation" do not match.
                # There are occasions that CDS does not have "translation".
                if item == 'translation':

                    # Retrieving "translation (i.e., amino acid sequences)" for each CDS
                    translation = feature.qualifiers.get('translation')
                    targetGenome_locusTag_aaSeq_dict[locusTag] = translation[0]
                    print >> fp, '>%s\n%s' % (str(locusTag), str(translation[0]))

                # Used to find "and" relationship in the GPR association
                if item == 'product':
                    product = feature.qualifiers.get('product')[0]
                    targetGenome_locusTag_prod_dict[locusTag] = product

                # Watch multiple EC_number's
                if item == 'EC_number':
                    ecnum = feature.qualifiers.get('EC_number')
                    targetGenome_locusTag_ec_dict[locusTag] = ecnum

    # Check if the gbk file has EC_number
    # Additional conditions should be given upon setup of in-house EC_number assigner
    logging.debug("len(targetGenome_locusTag_ec_dict.keys):")
    logging.debug(len(targetGenome_locusTag_ec_dict))
    logging.debug("len(targetGenome_locusTag_prod_dict.keys):")
    logging.debug(len(targetGenome_locusTag_prod_dict))

    fp.close()
    return targetGenome_locusTag_ec_dict, targetGenome_locusTag_prod_dict, targetFastaFilename


#making database files using fasta files
def make_blastDB(query_fasta, options):
    db_dir = options.metabolicmodeldir + os.sep + 'targetBlastDB'
    DBprogramName = utils.locate_executable('makeblastdb')

    utils.execute([DBprogramName,'-in',query_fasta,'-out',db_dir,'-dbtype','prot'])

    #Checks if DB is properly created; otherwise shutdown
    if os.path.isfile(options.metabolicmodeldir+os.sep+'targetBlastDB.psq') == False:
        logging.exception("error in make_blastDB: blast DB not created")
        #FIXME: don't use sys.exit
        sys.exit(1)


#Output: b0002,ASPK|b0002,0.0,100.00,820
#"1e-30" is set as a threshold for bidirectional best hits
#Four blastp output files need to be removed before running the modeling pipeline through a locally installed antiSMASH
def run_blastp(target_fasta = '', blastp_result = '', db_dir = '', evalue = 1e-30):
    BLASTPprogramName =  utils.locate_executable('blastp')
    # Execute blast if output file is not present
    if not os.path.isfile(blastp_result):
        args = [
            BLASTPprogramName,
            '-query', target_fasta,
            '-out', blastp_result,
            '-db', db_dir,
            '-evalue', str(evalue),
            '-outfmt', "10 qseqid sseqid evalue score length pident"
        ]
        out, err, retcode = utils.execute(args)
        if retcode != 0:
            logging.debug("out: %r, err: %r, retcode: %s", out, err, retcode)
    else:
        logging.warn("Found blast file %s, skipping new caluclation", blastp_result)


#Input: Results file from "run_blastp"
#Output: '\t' inserted between each element of the input
def parseBlaspResults(inputFile, outputFile):
    blastpResults_dict = {}
    fp = open(inputFile,'r')
    fp2 = open(outputFile,'w')
    itemnum=0
    for line in fp:
        key = itemnum
        sptList = line.strip().split(',')
        qseqid = sptList[0].strip()
        sseqid = sptList[1].strip()
        evalue = sptList[2].strip()
        score = float(sptList[3].strip())
        length = int(sptList[4].strip())
        pident = float(sptList[5].strip())

        blastpResults_dict[key] = {"query_locusTag": qseqid, "db_locusTag": sseqid, "evalue": evalue, "score": score, "length": length, "identity": pident}
        print >>fp2, '%s\t%s\t%s\t%f\t%d\t%f' % (qseqid, sseqid, evalue, score, length, pident)
        itemnum += 1
    fp.close()
    fp2.close()
    return blastpResults_dict


#Search the best hit of a particular gene to a target genome
#Input: Results file from "parseBlaspResults"
#Output: query_locusTag '\t' db_locusTag
def makeBestHits_dict(inputFile):
    bestHits_dict = {}
    fp1 = open(inputFile,'r')

    for line in fp1:
        sptList = line.strip().split('\t')
        query_locusTag = sptList[0].strip()
        db_locusTag = sptList[1].strip()

        if query_locusTag not in bestHits_dict.keys():
            #Value is in List to enable "append" below
            bestHits_dict[query_locusTag] = ([db_locusTag])

        #This additional condition is necessary because blastp strangely produces a redundant set of gene pairs:
        #e.g., SCO5892    ['SAV_7184', 'SAV_419', 'SAV_419', 'SAV_419', 'SAV_419', ...]
        elif '%s' %(db_locusTag) not in bestHits_dict[query_locusTag]:
            bestHits_dict[query_locusTag].append((db_locusTag))

    fp1.close()
    return bestHits_dict


#Finding bidirectional best hits
#Input: two dict data from "selectBestHits" (e.g.,bestHits_dict)
def getBBH(dic1,dic2):
    targetBBH_list = []
    temp_target_BBH_dict = {}

    for target_locusTag in dic1.keys():
        for temp_locusTag in dic1[target_locusTag]:
            if temp_locusTag in dic2.keys():
                for target_locusTag2 in dic2[temp_locusTag]:
                    # The BBH case
                    if target_locusTag == target_locusTag2:
                        if target_locusTag not in targetBBH_list:
                            targetBBH_list.append(target_locusTag)

                        # Some genes in template model have more than one homologous gene in a target genome
                        if temp_locusTag not in temp_target_BBH_dict.keys():
                            temp_target_BBH_dict[temp_locusTag] = ([target_locusTag])
                        else:
                            temp_target_BBH_dict[temp_locusTag].append((target_locusTag))
    return targetBBH_list, temp_target_BBH_dict


#A set of locusTag not included in BBH_list were considered nonBBH_list.
#Their respective reactions, if available, are added to the model in augPhase.
def get_nonBBH(targetGenome_locusTag_ec_dict, targetBBH_list):
    nonBBH_list = []

    for locusTag in targetGenome_locusTag_ec_dict.keys():
        if locusTag not in targetBBH_list:
            nonBBH_list.append(locusTag)

    nonBBH_list = sorted(set(nonBBH_list))
    return nonBBH_list


def calcBoolean(booleanList):
    booleanList2 = copy.deepcopy(booleanList)
    finalList = []
    threshold = 1

    if len(booleanList) == 0:
        return False

    for i in range(len(booleanList)):
        if type(booleanList[i])==list:
            for j in range(len(booleanList[i])):

                #Threshold created to differentiate BBH and nonBBH genes.
                if float(booleanList[i][j]) >= threshold:
                    booleanList[i][j] = 0
                if float(booleanList[i][j]) <= threshold:
                    booleanList2[i][j] = 1
            value=1
            for j in range(len(booleanList[i])):
                value = value * booleanList2[i][j]
            finalList.append(value)
        else:
            #Threshold created to differentiate BBH and nonBBH genes.
            if float(booleanList[i]) >= threshold:
                booleanList[i] = 0
            if float(booleanList[i]) <= threshold:
                booleanList2[i] = 1
            finalList.append(booleanList2[i])

    value=0
    for i in range(len(finalList)):
        value = value + finalList[i]

    if value == 0:
        return False
    else:
        return True


#Output: e.g., [[2, '0'], ['0', 2]]
#Now considers nonBBH genes without removing them in the Boolean list.
#For ( A and B), if one of them is nonBBH, its rxn should be False (subject to removal).
def makeBooleanFormat(temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict):
    booleanList = tempModel_biggRxnid_locusTag_dict
    valueList = copy.deepcopy(booleanList)

    for i in range(len(booleanList)):
        if type(booleanList[i])==list:
            for j in range(len(booleanList[i])):
                geneid = booleanList[i][j]
                if geneid in temp_target_BBH_dict.keys():
                    value = 0
                #Now considers nonBBH genes without removing them in the Boolean list
                else:
                    value = 2 #For nonBBH genes
                    booleanList[i][j]='nonBBH'
                valueList[i][j] = value
        else:
            geneid = booleanList[i]
            if geneid in temp_target_BBH_dict.keys():
                value = 0
            else:
                value = 'na'
                booleanList[i]='na'
            valueList[i] = value

    while 'na' in booleanList or [] in booleanList:
        if 'na' in booleanList:
            booleanList.pop( booleanList.index('na') )
        if [] in booleanList:
            booleanList.pop( booleanList.index([]) )

    newbooleanList = []
    for i in range(len(booleanList)):
        if type(booleanList[i]) != list:
            newbooleanList.append(booleanList[i])
    newbooleanList = list(set(newbooleanList))

    for i in range(len(booleanList)):
        if type(booleanList[i]) == list:
            tmpList = booleanList[i]
            tmpList = list(set(tmpList))
            newbooleanList.append(tmpList)

    booleanList = newbooleanList
    valueList2 = copy.deepcopy(booleanList)

    for i in range(len(booleanList)):
        if type(booleanList[i])==list:
            for j in range(len(booleanList[i])):
                geneid = booleanList[i][j]
                if 'nonBBH' not in geneid:
                    value = 0
                else:
                    value = 2
                valueList2[i][j] = value
        else:
            geneid = booleanList[i]
            value = 0
            valueList2[i] = value
    return valueList2


def labelRxnToRemove(model, temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict):
    rxnToRemove_dict = {}

    for biggRxnid in tempModel_biggRxnid_locusTag_dict.keys():
        rxn = model.reactions.get_by_id(biggRxnid)
        #Prevent removal of transport reactions from the template model
        if 'Transport' not in rxn.name and 'transport' not in rxn.name and 'Exchange' not in rxn.name and 'exchange' not in rxn.name:
            booleanList = makeBooleanFormat(temp_target_BBH_dict, tempModel_biggRxnid_locusTag_dict[biggRxnid])
            rxnToRemove_dict[biggRxnid] = calcBoolean(booleanList)

    return rxnToRemove_dict


def pruneModel(model, rxnToRemove_dict, solver_arg):
    rxnToRemoveEssn_dict = {}
    rxnRemoved_dict = {}
    rxnRetained_dict = {}

    for rxnid in rxnToRemove_dict.keys():

        #Single reaction deletion is performed only for reactions labelled as "False"
        if rxnToRemove_dict[rxnid] == False:
            growth_rate_dict, solution_status_dict, problem_dict = single_deletion(model, list([rxnid]), element_type='reaction', solver=solver_arg)

            #Checks optimality first.
            if str(solution_status_dict.values()[0]) == 'optimal':

                #Full list of reactions and predicted growth rates upon their deletions
                rxnToRemoveEssn_dict[rxnid] = float(growth_rate_dict.values()[0])

                #Checks growth rate upon reaction deletion
                if float(growth_rate_dict.values()[0]) >= 0.01:
                    model.remove_reactions(rxnid)
                    #List of reactions removed from the template model
                    rxnRemoved_dict[rxnid] = float(growth_rate_dict.values()[0])
                    logging.debug("Removed reaction: %s, %s, %s, %s", rxnid, growth_rate_dict.values()[0], len(model.reactions), len(model.metabolites))
                else:
                    #List of reactions retained in the template model
                    rxnRetained_dict[rxnid] = float(growth_rate_dict.values()[0])
                    logging.debug("Retained reaction: %s, %s, %s, %s", rxnid, growth_rate_dict.values()[0], len(model.reactions), len(model.metabolites))

    #Removing metabolites that are not used in the reduced model
    prune_unused_metabolites(model)
    modelPruned = copy.deepcopy(model)

    return modelPruned, rxnToRemoveEssn_dict, rxnRemoved_dict, rxnRetained_dict


def get_gpr_fromString_toList(line):
    calcNewList = []
    line = line.strip()
    calcList = line.split('or')
    for c in calcList:
        c = c.replace('(','')
        c = c.replace(')','')
        c = c.replace(' ','')
        c = c.strip()
        if 'and' in c:
            newlist = c.split('and')
            newlist = list(set(newlist))
            newlist.sort()
            calcNewList.append(newlist)
        else:
            geneid=c.strip()
            if geneid not in calcNewList:
                calcNewList.append(geneid)

    return calcNewList

def swap_locusTag_tempModel(modelPruned, temp_target_BBH_dict):

    #Retrieves reactions associated with each homologous gene in template model
    for BBHrxn in modelPruned.reactions:
        booleanList = []
        #Retrieves all the genes associated with a reaction having the homologous gene and transforms String to List
        booleanList = get_gpr_fromString_toList(BBHrxn.gene_reaction_rule)

        modified_booleanList = []
        for tempLocusTag in booleanList:

            #Checks if the element itself is List.
            #If the element is not List, then gene in template model is directly replaced with genes in target genome
            if type(tempLocusTag) != list:
                if tempLocusTag in temp_target_BBH_dict:
                    booleanList.pop(booleanList.index(tempLocusTag))
                    for targetLocusTag in temp_target_BBH_dict[tempLocusTag]:
                        modified_booleanList.append(targetLocusTag)
                else:
                    modified_booleanList.append( tempLocusTag )

            #This is the case the element is List
            else:
                temp_gpr_list = []
                for eachLocusTag in tempLocusTag:
                    if eachLocusTag in temp_target_BBH_dict:
                        for targetLocusTag in temp_target_BBH_dict[eachLocusTag]:
                            temp_gpr_list.append(targetLocusTag)
                    else:
                        temp_gpr_list.append(eachLocusTag)
                #This case was not generated, but just in case
                if len(temp_gpr_list)==1:
                    logging.debug(temp_gpr_list)
                    modified_booleanList.append(temp_gpr_list[0])
                elif len(temp_gpr_list) > 1:
                    modified_booleanList.append( temp_gpr_list )


        #Converts GPR in List to String:
        booleanList = copy.deepcopy( modified_booleanList )
        stringlist = []
        for i in range(len(booleanList)):
            booleanstring = ''
            if type(booleanList[i])==list:
                booleanstring=' and '.join(booleanList[i])
                booleanstring='('+booleanstring+')'
                stringlist.append(booleanstring)
            else:
                booleanstring=booleanList[i]
                stringlist.append(booleanstring)

        booleanstring2 = ''
        if len(stringlist) > 0:
            booleanstring2 = '('+' or '.join(stringlist)+')'
        BBHrxn.add_gene_reaction_rule(booleanstring2)

    modelPrunedGPR = copy.deepcopy(modelPruned)
    return modelPrunedGPR
