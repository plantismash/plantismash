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
from antismash import utils

def _update_sec_met_entry(clusterfeature, smiles_string):
    clusterfeature.qualifiers['structure'] = [smiles_string]

def generate_chemical_structure_preds(pksnrpsvars, seq_record, options):
    #Create directory to store structures
    options.structuresfolder = path.abspath(path.join(options.outputfoldername, "structures"))
    if not os.path.exists(options.structuresfolder):
        os.mkdir(options.structuresfolder)
    originaldir = os.getcwd()
    structure_drawing_dir = utils.get_full_path(__file__, '') + os.sep + "NRPeditor"
    os.chdir(structure_drawing_dir)
    #Combine predictions into a prediction of the final chemical structure and generate images
    geneclusters = utils.get_cluster_features(seq_record)
    for genecluster in geneclusters:
        smiles_string = "N/A"
        geneclusternr = utils.get_cluster_number(genecluster)
        if geneclusternr in pksnrpsvars.compound_pred_dict:
            # if product is ectoine generate predefined SMILE string and generate structure
            if pksnrpsvars.compound_pred_dict[geneclusternr] == "ectoine":
                smiles_string = "CC1=NCCC(N1)C(=O)O"
                smilesfile = open("genecluster" + str(geneclusternr) + ".smi","w")
                smilesfile.write(smiles_string)
                smilesfile.close()
                depictstatus = depict_smile(geneclusternr,options.structuresfolder)
                if depictstatus == "failed":
                    pksnrpsvars.failedstructures.append(geneclusternr)
                elif genecluster in pksnrpsvars.failedstructures:
                    del pksnrpsvars.failedstructures[pksnrpsvars.failedstructures.index(geneclusternr)]
            else:
                # use information on peptide / polyketide sequence to gernerate structure image
                residues = pksnrpsvars.compound_pred_dict[geneclusternr].replace("(","").replace(")","").replace(" + "," ").replace("-"," ")
                nrresidues = len(residues.split(" "))
                if nrresidues > 1:
                    if sys.platform == ('win32') or sys.platform == ('darwin'):
                        structcommand = 'main input 100 4000 1000 AA DDV DIM ' + str(nrresidues + 1) + ' "'
                    elif sys.platform == ('linux2'):
                        structcommand = './main input 100 4000 1000 AA DDV DIM ' + str(nrresidues + 1) + ' "'
                    for i in [res for res in residues.split(" ") if len(res) > 1]:
                        structcommand = structcommand + i + " "
                    structcommand = structcommand + 'TE"'
                    smilesinfo = os.popen(structcommand)
                    smilesinfo = smilesinfo.read()
                    smiles_string = (smilesinfo.split("core peptide: ")[1]).split("\ntermintype")[0]
                    if sys.platform == ('linux2') or sys.platform == ('darwin'):
                        smiles_string.replace("[X]","[*:X]")
                        smiles_string2 = ""
                        a = 1
                        for k in smiles_string:
                            if k == "X":
                                smiles_string2 = smiles_string2 + str(a)
                                a += 1
                            else:
                                smiles_string2 = smiles_string2 + k
                        smiles_string = smiles_string2
                    smilesfile = open("genecluster" + str(geneclusternr) + ".smi","w")
                    smilesfile.write(smiles_string)
                    smilesfile.close()
                    depictstatus = depict_smile(geneclusternr, options.structuresfolder)
                    if depictstatus == "failed":
                        pksnrpsvars.failedstructures.append(geneclusternr)
        _update_sec_met_entry(genecluster, smiles_string)
    os.chdir(originaldir)

def depict_smile(genecluster,structuresfolder):
    if sys.platform == ('win32') or sys.platform == ('darwin'):
        indigo_depict_command1 = "indigo-depict genecluster" + str(genecluster) + ".smi " + "genecluster" + str(genecluster) + "_icon.png -query -w 200 -h 150"
        indigo_depict_command2 = "indigo-depict genecluster" + str(genecluster) + ".smi " + "genecluster" + str(genecluster) + ".png -query"
    elif sys.platform == ('linux2'):
        indigo_depict_command1 = "./indigo-depict genecluster" + str(genecluster) + ".smi " + "genecluster" + str(genecluster) + "_icon.png -query -w 200 -h 150"
        indigo_depict_command2 = "./indigo-depict genecluster" + str(genecluster) + ".smi " + "genecluster" + str(genecluster) + ".png -query"
    os.system(indigo_depict_command1)
    os.system(indigo_depict_command2)
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
