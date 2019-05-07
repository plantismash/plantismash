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


import shutil
import os
import sys
from os import path
import logging
from antismash import utils
from indigo import Indigo
from indigo_renderer import IndigoRenderer
from helperlibs.wrappers.io import TemporaryDirectory


def _update_sec_met_entry(clusterfeature, smiles_string):
    clusterfeature.qualifiers['structure'] = [smiles_string]


def generate_chemical_structure_preds(pksnrpsvars, seq_record, options):
    #Create directory to store structures
    options.structuresfolder = path.abspath(path.join(options.outputfoldername, "structures"))
    if not os.path.exists(options.structuresfolder):
        os.mkdir(options.structuresfolder)

    #Combine predictions into a prediction of the final chemical structure and generate images
    geneclusters = utils.get_cluster_features(seq_record)

    for genecluster in geneclusters:
        geneclusternr = utils.get_cluster_number(genecluster)
        smiles_string = ""
        if pksnrpsvars.compound_pred_dict.has_key(geneclusternr):

            #print "output_modules/html/pksnrpsvars.compound_pred_dict:"
            #print pksnrpsvars.compound_pred_dict

            residues = pksnrpsvars.compound_pred_dict[geneclusternr].replace("(","").replace(")","").replace(" + "," ").replace("-"," ")


            #Now generates SMILES of predicted secondary metabolites without NP.searcher
            residuesList = residues.split(" ")

            #Counts the number of malonate and its derivatives in polyketides
            mal_count = 0
            for i in residuesList:
                if "mal" in i:
                    mal_count += 1

            nrresidues = len(residuesList)

            #Reflecting reduction states of ketide groups starting at beta carbon of type 1 polyketide
            if "pk" in residuesList and "mal" in residuesList[-1]:
                residuesList.pop(residuesList.index('pk')+1)
                residuesList.append('pks-end1')
            elif mal_count == len(residuesList):
                if residuesList[0] == "mal":
                    residuesList[0] = "pks-start1"
                if residuesList[-1] == "ccmal":
                    residuesList.append('pks-end2')

            if nrresidues > 1:
                #Conventionally used aaSMILES was used;
                #chirality expressed with "@@" causes indigo error
                smiles_monomer = open(os.path.dirname(os.path.realpath(__file__)) + os.sep + 'aaSMILES.txt','r')
                smiles = smiles_monomer.readline()
                smiles = smiles_monomer.readline()

                aa_smiles_dict = {}
                while smiles:
                    smiles = smiles.split()
                    if len(smiles) > 1:
                        smiles[0] = smiles[0].strip()
                        smiles[1] = smiles[1].strip()
                        aa_smiles_dict[smiles[0]] = smiles[1]
                    smiles = smiles_monomer.readline()
                smiles_monomer.close()

                for monomer in residuesList:
                    if monomer in aa_smiles_dict.keys():
                        smiles_string += aa_smiles_dict[monomer]
                logging.debug("Cluster %s: smiles_string: %s", geneclusternr, smiles_string)
                with TemporaryDirectory(change=True):
                    smilesfile = open("genecluster" + str(geneclusternr) + ".smi", "w")
                    smilesfile.write(smiles_string)
                    smilesfile.close()
                    depictstatus = depict_smile(geneclusternr, options.structuresfolder)
                if depictstatus == "failed":
                    pksnrpsvars.failedstructures.append(geneclusternr)
        elif utils.get_cluster_type(genecluster) == "ectoine":
            smiles_string = "CC1=NCCC(N1)C(=O)O"
            with TemporaryDirectory(change=True):
                smilesfile = open("genecluster" + str(geneclusternr) + ".smi", "w")
                smilesfile.write(smiles_string)
                smilesfile.close()
                depictstatus = depict_smile(geneclusternr, options.structuresfolder)
            if depictstatus == "failed":
                pksnrpsvars.failedstructures.append(geneclusternr)
            elif genecluster in pksnrpsvars.failedstructures:
                del pksnrpsvars.failedstructures[pksnrpsvars.failedstructures.index(geneclusternr)]
            pksnrpsvars.compound_pred_dict[geneclusternr] = "ectoine"
        _update_sec_met_entry(genecluster, smiles_string)

def depict_smile(genecluster, structuresfolder):
    indigo = Indigo()
    renderer = IndigoRenderer(indigo)
    query = indigo.loadMoleculeFromFile("genecluster" + str(genecluster) + ".smi")
    indigo.setOption("render-coloring", True)
    renderer.renderToFile(query, "genecluster" + str(genecluster) + ".png")

    indigo.setOption("render-image-size", 200, 150)
    renderer.renderToFile(query, "genecluster" + str(genecluster) + "_icon.png")
    dircontents = os.listdir(os.getcwd())
    geneclusterstring = "genecluster" + str(genecluster) + ".png"
    if geneclusterstring in dircontents:
        shutil.copy("genecluster" + str(genecluster) + ".png", structuresfolder)
        shutil.copy("genecluster" + str(genecluster) + "_icon.png", structuresfolder)
        shutil.copy("genecluster" + str(genecluster) + ".smi", structuresfolder)
        os.remove("genecluster" + str(genecluster) + ".png")
        os.remove("genecluster" + str(genecluster) + "_icon.png")
        os.remove("genecluster" + str(genecluster) + ".smi")
        smiles_input = path.join('SMILES', 'input')
        if path.exists(smiles_input):
            os.remove(smiles_input)
        return "success"
    else:
        return "failed"
