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

import os
from antismash import utils

def parse_nrps_preds(options, pksnrpsvars):
    pksnrpsvars.minowa_nrps_preds = {}
    pksnrpsvars.minowa_nrps_preds_details = {}
    pksnrpsvars.nrps_svm_preds = {}
    pksnrpsvars.nrps_svm_preds_details = {}
    pksnrpsvars.nrps_code_preds = {}
    pksnrpsvars.nrps_code_preds_details = {}
    substratetransdict2 = {'pipecolate':'pip','fOHOrn':'orn','beta-Lys':'blys','5NhOrn':'orn','OHOrn':'orn','Aad':'Aaa','bOHTyr':'bht'}
    if len(pksnrpsvars.nrpsnames) > 0:
        minowa_a_file = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_minowa_nrpspredoutput.txt","r")
        minowa_a_file = minowa_a_file.read()
        minowa_a_file = minowa_a_file.replace("\r","\n")
        parts = minowa_a_file.split("\\\\\n")[1:]
        for i in parts:
            partlines = i.split("\n")
            acc = partlines[0]
            tophit = partlines[2].split("\t")[0]
            if tophit in substratetransdict2.keys():
                tophit = substratetransdict2[tophit]
            pksnrpsvars.minowa_nrps_preds[acc] = tophit.lower()
            pksnrpsvars.minowa_nrps_preds_details[acc] = "<b>Minowa HMM method A-domain<br>Substrate specificity prediction top hits:</b><br>\n" + partlines[1] + "<br>\n" + partlines[2] + "<br>\n" + partlines[3] + "<br>\n" + partlines[4] + "<br><br>\n\n"
        nrpspredictorfile1 = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_nrpspredictor2_svm.txt","r")
        nrpspredictorfile2 = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_nrpspredictor2_codes.txt","r")
        nrpspredictorfile1 = nrpspredictorfile1.read()
        nrpspredictorfile1 = nrpspredictorfile1.replace("\r","\n")
        lines = nrpspredictorfile1.split("\n")[1:-1]
        for k in lines:
            tabs = k.split("\t")
            if tabs[6] != "N/A":
                pksnrpsvars.nrps_svm_preds[tabs[0]] = tabs[6]
            elif tabs[5] != "N/A":
                pksnrpsvars.nrps_svm_preds[tabs[0]] = tabs[5]
            elif tabs[4] != "N/A":
                pksnrpsvars.nrps_svm_preds[tabs[0]] = tabs[4]
            else:
                pksnrpsvars.nrps_svm_preds[tabs[0]] = tabs[3]
            pksnrpsvars.nrps_svm_preds_details[tabs[0]] = "<b> NRPSPredictor2 SVM prediction details:</b><br>\n8 Angstrom 34 AA code:<br>\n" + tabs[1] + "<br>\nPredicted physicochemical class:<br>\n" + tabs[3] + "<br>\nLarge clusters prediction:<br>\n" + tabs[4] + "<br>\nSmall clusters prediction:<br>\n" + tabs[5] + "<br>\nSingle AA prediction:<br>\n" + tabs[6] + "<br><br>\n\n"
        nrpspredictorfile2 = nrpspredictorfile2.read()
        nrpspredictorfile2 = nrpspredictorfile2.replace("\r","\n")
        lines = nrpspredictorfile2.split("\n")[:-1]
        for k in lines:
            tabs = k.split("\t")
            pksnrpsvars.nrps_code_preds[tabs[0]] = tabs[1]
            pksnrpsvars.nrps_code_preds_details[tabs[0]] = "<b> NRPSPredictor2 Stachelhaus code prediction:</b><br>\n" + tabs[1] + "<br><br>\n\n"

def parse_pks_preds(options, pksnrpsvars):
    pksnrpsvars.minowa_pks_preds_details = {}
    pksnrpsvars.minowa_pks_preds = {}
    pksnrpsvars.pks_code_preds ={}
    pksnrpsvars.pks_code_preds_details ={}
    pksnrpsvars.minowa_cal_preds = {}
    pksnrpsvars.minowa_cal_preds_details = {}
    substratetransdict = {'Malonyl-CoA':'mal','Methylmalonyl-CoA':'mmal','Methoxymalonyl-CoA':'mxmal','Ethylmalonyl-CoA':'emal','Isobutyryl-CoA':'isobut','2-Methylbutyryl-CoA':'2metbut','trans-1,2-CPDA':'trans-1,2-CPDA','Acetyl-CoA':'Acetyl-CoA','Benzoyl-_CoA':'benz','Propionyl-CoA':'prop','3-Methylbutyryl-CoA':'3metbut','Ethylmalonyl-CoA':'Ethyl_mal','CE-Malonyl-CoA':'cemal','2-Rhyd-Malonyl-CoA':'2Rhydmal','CHC-CoA':'CHC-CoA','inactive':'inactive'}
    if len(pksnrpsvars.pksnames) > 0:
        minowa_at_file = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_minowa_pkspredoutput.txt","r")
        minowa_at_file = minowa_at_file.read()
        minowa_at_file = minowa_at_file.replace("\r","\n")
        parts = minowa_at_file.split("\\\\\n")[1:]
        for i in parts:
            partlines = i.split("\n")
            acc = partlines[0]
            if substratetransdict.has_key(partlines[2].split("\t")[0]):
                tophit = substratetransdict[partlines[2].split("\t")[0]]
            else:
                tophit = "pk"
            pksnrpsvars.minowa_pks_preds[acc] = tophit
            pksnrpsvars.minowa_pks_preds_details[acc] = "<b>Minowa HMM method AT-domain<br>Substrate specificity prediction top hits:</b><br>\n" + partlines[1] + "<br>\n" + partlines[2] + "<br>\n" + partlines[3] + "<br>\n" + partlines[4] + "<br><br>\n\n"
        pkssignaturefile = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_pkssignatures.txt","r")
        pkssignaturefile = pkssignaturefile.read()
        pkssignaturefile = pkssignaturefile.replace("\r","\n")
        parts = pkssignaturefile.split("//\n")[1:]
        for i in parts:
            partlines = i.split("\n")
            partlines2 = []
            for j in partlines:
                if j != "":
                    partlines2.append(j)
            partlines = partlines2
            acc = partlines[0].split("\t")[0]
            if len(partlines) > 2:
                tophit = (partlines[1].split("\t")[0]).split("__")[1]
                pksnrpsvars.pks_code_preds[acc] = tophit
                codes = []
                prots = []
                scores = []
                for i in partlines[1:4]:
                    codes.append(i.split("\t")[0])
                    prot = i.split("\t")[1]
                    prot = prot.replace("_AT"," (AT")
                    prot = prot.replace("__","): ")
                    prots.append(prot)
                    scores.append(i.split("\t")[2])
                if len(prots) >= 3:
                    pksnrpsvars.pks_code_preds_details[acc] = "<b>PKS Active Site Signature method<br>AT-domain substrate specificity prediction top hits:</b><br>\nCode:" + partlines[0].split("\t")[1] + "<br>\n" + codes[0] + " - " + prots[0] + " : (" + scores[0] + "% identity)<br>\n" + codes[1] + " - " + prots[1] + " : (" + scores[1] + "% identity)<br>\n" + codes[2] + " - " + prots[2] + " : (" + scores[2] + "% identity)<br><br>\n\n"
                elif len(prots) == 2:
                    pksnrpsvars.pks_code_preds_details[acc] = "<b>PKS Active Site Signature method<br>AT-domain substrate specificity prediction top hits:</b><br>\nCode:" + partlines[0].split("\t")[1] + "<br>\n" + codes[0] + " - " + prots[0] + " : (" + scores[0] + "% identity)<br>\n" + codes[1] + " - " + prots[1] + " : (" + scores[1] + "% identity)<br><br>\n\n"
                elif len(prots) == 1:
                    pksnrpsvars.pks_code_preds_details[acc] = "<b>PKS Active Site Signature method<br>AT-domain substrate specificity prediction top hits:</b><br>\nCode:" + partlines[0].split("\t")[1] + "<br>\n" + codes[0] + " - " + prots[0] + " : (" + scores[0] + "% identity)<br><br>\n\n"
            else:
                pksnrpsvars.pks_code_preds[acc] = "N/A"
                pksnrpsvars.pks_code_preds_details[acc] = "<b>PKS Active Site Signature method<br>No AT-domain substrate specificity prediction hits above 40% identity.<br>\n\n"
    if len(pksnrpsvars.calnames) > 0:
        minowa_cal_file = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_minowa_calpredoutput.txt","r")
        minowa_cal_file = minowa_cal_file.read()
        minowa_cal_file = minowa_cal_file.replace("\r","\n")
        parts = minowa_cal_file.split("\\\\\n")[1:]
        for i in parts:
            partlines = i.split("\n")
            acc = partlines[0]
            tophit = partlines[2].split("\t")[0]
            pksnrpsvars.minowa_cal_preds[acc] = tophit
            pksnrpsvars.minowa_cal_preds_details[acc] = "<b>Minowa HMM method<br>CAL-domain substrate specificity prediction top hits:</b><br>\n" + partlines[1] + "<br>\n" + partlines[2] + "<br>\n" + partlines[3] + "<br>\n" + partlines[4] + "<br><br>\n\n"

def parse_kr_activity_preds(options, pksnrpsvars):
    pksnrpsvars.kr_activity_preds = {}
    pksnrpsvars.kr_stereo_preds = {}
    if len(pksnrpsvars.krnames) > 0:
        krfile = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_krpredoutput.txt","r")
        krfile = krfile.read()
        krfile = krfile.replace("\r","\n")
        krlines = krfile.split("\n")[:-1]
        for i in krlines:
            tabs = i.split("\t")
            pksnrpsvars.kr_activity_preds[tabs[0]] = tabs[1]
            pksnrpsvars.kr_stereo_preds[tabs[0]] = tabs[2]


def calculate_consensus_prediction(pksnrpsvars, seq_record):
    #Combine substrate specificity predictions into consensus prediction
    pksnrpsvars.consensuspreds = {}
    pksnrpsvars.consensuspreds_transat = {}
    available_smiles_parts = ['GLY','ALA','VAL','LEU','ILE','MET','PRO','PHE','TRP','SER','THR','ASN','GLN','TYR','CYS','LYS','ARG',
    'HIS','ASP','GLU','MPRO','ORN','PGLY','DAB','BALA','AEO','DHA','PIP','BMT','gly','ala','val','leu','ile','met','pro','phe','trp','ser',
    'thr','asn','gln','tyr','cys','lys','arg','his','asp','glu','aaa','mpro','dhb','2hiva','orn','pgly','dab','bala','aeo','4mha','pico','phg',
    'dha','scy','pip','bmt','adds','aad','abu','hiv','dhpg','bht','3-me-glu','4pPro','ala-b','ala-d','dht','Sal','tcl','lys-b','hpg','hyv-d',
    'iva','vol','mal','mmal','ohmal', 'redmal', 'mxmal','emal','nrp','pk','Gly','Ala','Val','Leu','Ile','Met','Pro','Phe','Trp','Ser','Thr','Asn','Gln','Tyr',
    'Cys','Lys','Arg','His','Asp','Glu','Mpro','23Dhb','34Dhb','2Hiva','Orn','Pgly','Dab','Bala','Aeo','4Mha','Pico','Aaa','Dha','Scy','Pip',
    'Bmt','Adds','DHpg','DHB','nrp','pk']

    #Extracting gene cluster type (e.g., "transatpks")
    for f in utils.get_cluster_features(seq_record):
	cluster_info = f.qualifiers

    for feature in pksnrpsvars.pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        nra = 0
        nrat = 0
        nrcal = 0
	nrtransat = 0
        j = pksnrpsvars.domaindict[locus]

        for k in j:
            if 'transatpks' not in cluster_info['product'][0]:
                if k[0] == "PKS_AT":
                    nrat += 1
                    preds = []
                    preds.append(pksnrpsvars.minowa_pks_preds[locus + "_AT" + str(nrat)])
                    preds.append(pksnrpsvars.pks_code_preds[locus + "_AT" + str(nrat)])
                    cpred = "n"
                    for l in preds:
                        if preds.count(l) > 1:
                            if l in available_smiles_parts:
                                pksnrpsvars.consensuspreds[locus + "_AT" + str(nrat)] = l
                            else:
                                pksnrpsvars.consensuspreds[locus + "_AT" + str(nrat)] = "pk"
                            cpred = "y"
                    if cpred == "n":
                        pksnrpsvars.consensuspreds[locus + "_AT" + str(nrat)] = "pk"
	    elif 'transatpks' in cluster_info['product'][0]:
                if k[0] == "PKS_AT":
                    nrat += 1
                    preds = []
                    preds.append(pksnrpsvars.minowa_pks_preds[locus + "_AT" + str(nrat)])
                    preds.append(pksnrpsvars.pks_code_preds[locus + "_AT" + str(nrat)])
                    cpred = "n"

		    #Only for the writing purpose in sec_record (i.e., trans-AT)
                    for l in preds:
                        if preds.count(l) > 1:
                            if l in available_smiles_parts:
                                pksnrpsvars.consensuspreds_transat[locus + "_AT" + str(nrat)] = l
                            else:
                                pksnrpsvars.consensuspreds_transat[locus + "_AT" + str(nrat)] = "pk"
                            cpred = "y"
                    if cpred == "n":
                        pksnrpsvars.consensuspreds_transat[locus + "_AT" + str(nrat)] = "pk"
		#For chemical display purpose for chemicals from trans-AT PKS gene cluster
		#mal is always assumed for trans-AT
		if k[0] == "PKS_KS":
		    nrtransat += 1
		    pksnrpsvars.consensuspreds[locus + "_KS" + str(nrtransat)] = "mal"
		    cpred = "y"
            if k[0] == "AMP-binding" or k[0] == "A-OX":
                nra +=1
                preds = []
                preds.append(pksnrpsvars.minowa_nrps_preds[locus + "_A" + str(nra)])
                preds.append(pksnrpsvars.nrps_svm_preds[locus + "_A" + str(nra)])
                preds.append(pksnrpsvars.nrps_code_preds[locus + "_A" + str(nra)])
                cpred = "n"
                for l in preds:
                    if preds.count(l) > 1:
                        if l in available_smiles_parts:
                            pksnrpsvars.consensuspreds[locus + "_A" + str(nra)] = l
                        else:
                            pksnrpsvars.consensuspreds[locus + "_A" + str(nra)] = "nrp"
                        cpred = "y"
                if cpred == "n":
                    pksnrpsvars.consensuspreds[locus + "_A" + str(nra)] = "nrp"
            if k[0] == "CAL_domain":
                nrcal += 1
                if pksnrpsvars.minowa_cal_preds[locus + "_CAL" + str(nrcal)] in available_smiles_parts:
                    pksnrpsvars.consensuspreds[locus + "_CAL" + str(nrcal)] = pksnrpsvars.minowa_cal_preds[locus + "_CAL" + str(nrcal)]
                else:
                    pksnrpsvars.consensuspreds[locus + "_CAL" + str(nrcal)] = "pk"


def generate_domainnamesdict(pksnrpsvars):
    pksnrpsvars.domainnamesdict = {}
    for feature in pksnrpsvars.pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        j = pksnrpsvars.domaindict[locus]
        domainnames = []
        for k in j:
            domainnames.append(k[0])
        pksnrpsvars.domainnamesdict[locus] = domainnames

def find_duplicate_position(domainList,item):
    start_at = -1
    locs = []
    while True:
        try:
	    loc = domainList.index(item,start_at+1)
	except ValueError:
	    break
	else:
	    locs.append(loc)
	    start_at = loc
    return locs

def insert_modified_monomers(pksnrpsvars, seq_record):
    locusTag_domain = []

    #Extracting gene cluster type (e.g., "transatpks")
    for f in utils.get_cluster_features(seq_record):
	cluster_info = f.qualifiers

    #pksnrpsvars.domainnamesdict = {'CRYAR_RS43165': ['PKS_KS', 'PKS_AT',...]}
    #Get a unique set of genes having ATs
    for key in pksnrpsvars.domainnamesdict.keys():
        if key not in locusTag_domain:
            locusTag_domain.append(key)

    locusTag_domain =  sorted(set(locusTag_domain))

    for locusTag in locusTag_domain:
	at_list = find_duplicate_position(pksnrpsvars.domainnamesdict[locusTag], 'PKS_AT')
        #For transatpks
        ks_list = find_duplicate_position(pksnrpsvars.domainnamesdict[locusTag], 'PKS_KS')
        kr_list = find_duplicate_position(pksnrpsvars.domainnamesdict[locusTag], 'PKS_KR')
        dh_list = find_duplicate_position(pksnrpsvars.domainnamesdict[locusTag], 'PKS_DH')
        er_list = find_duplicate_position(pksnrpsvars.domainnamesdict[locusTag], 'PKS_ER')

        if 'transatpks' not in cluster_info['product'][0]:
	    for at_idx in range(len(at_list)):
                #Monomer change caused by only KR
	        for kr_idx in range(len(kr_list)):
		    if at_idx + 1 <= len(at_list) - 1:
		        if kr_list[kr_idx] > at_list[at_idx] and kr_list[kr_idx] < at_list[at_idx + 1]:
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "mal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ohmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "mmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ohmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "mxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ohmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "emal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ohemal"
		    if at_idx + 1 > len(at_list) - 1:
		        if kr_list[kr_idx] > at_list[at_idx]:
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "mal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ohmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "mmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ohmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "mxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ohmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "emal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ohemal"

                #Monomer change caused by KR and DH
	        for dh_idx in range(len(dh_list)):
		    if at_idx + 1 <= len(at_list) - 1:
		        if dh_list[dh_idx] > at_list[at_idx] and dh_list[dh_idx] < at_list[at_idx + 1]:
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ohmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ccmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ohmmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ccmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ohmxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ccmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ohemal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ccemal"
		    if at_idx + 1 > len(at_list) - 1:
		        if dh_list[dh_idx] > at_list[at_idx]:
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ohmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ccmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ohmmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ccmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ohmxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ccmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ohemal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "ccemal"

                #Monomer change caused by KR, DH and ER
	        for er_idx in range(len(er_list)):
		    if at_idx + 1 <= len(at_list) - 1:
		        if er_list[er_idx] > at_list[at_idx] and er_list[er_idx] < at_list[at_idx + 1]:
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ccmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "redmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ccmmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "redmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ccmxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "redmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ccemal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "redemal"
		    if at_idx + 1 > len(at_list) - 1:
		        if er_list[er_idx] > at_list[at_idx]:
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ccmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "redmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ccmmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "redmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ccmxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "redmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] == "ccemal":
			        pksnrpsvars.consensuspreds[locusTag+"_AT"+str(at_idx+1)] = "redemal"

        if 'transatpks' in cluster_info['product'][0]:
	    for ks_idx in range(len(ks_list)):
                #Monomer change caused by only KR
	        for kr_idx in range(len(kr_list)):
		    if ks_idx + 1 <= len(ks_list) - 1:
		        if kr_list[kr_idx] > ks_list[ks_idx] and kr_list[kr_idx] < ks_list[ks_idx + 1]:
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "mal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ohmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "mmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ohmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "mxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ohmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "emal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ohemal"
		    if ks_idx + 1 > len(ks_list) - 1:
		        if kr_list[kr_idx] > ks_list[ks_idx]:
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "mal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ohmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "mmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ohmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "mxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ohmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "emal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ohemal"

                #Monomer change caused by KR and DH
	        for dh_idx in range(len(dh_list)):
		    if ks_idx + 1 <= len(ks_list) - 1:
		        if dh_list[dh_idx] > ks_list[ks_idx] and dh_list[dh_idx] < ks_list[ks_idx + 1]:
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ohmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ccmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ohmmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ccmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ohmxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ccmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ohemal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ccemal"
		    if ks_idx + 1 > len(ks_list) - 1:
		        if dh_list[dh_idx] > ks_list[ks_idx]:
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ohmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ccmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ohmmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ccmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ohmxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ccmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ohemal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "ccemal"

                #Monomer change caused by KR, DH and ER
	        for er_idx in range(len(er_list)):
		    if ks_idx + 1 <= len(ks_list) - 1:
		        if er_list[er_idx] > ks_list[ks_idx] and er_list[er_idx] < ks_list[ks_idx + 1]:
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ccmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "redmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ccmmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "redmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ccmxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "redmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ccemal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "redemal"
		    if ks_idx + 1 > len(ks_list) - 1:
		        if er_list[er_idx] > ks_list[ks_idx]:
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ccmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "redmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ccmmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "redmmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ccmxmal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "redmxmal"
			    if pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] == "ccemal":
			        pksnrpsvars.consensuspreds[locusTag+"_KS"+str(ks_idx+1)] = "redemal"


def write_substr_spec_preds_to_html(options, pksnrpsvars):
    #Make output folder for storing HTML files with substrate specificities
    originaldir = os.getcwd()
    options.substrspecsfolder = options.full_outputfolder_path + os.sep + "html"
    if not os.path.exists(options.substrspecsfolder):
        os.mkdir(options.substrspecsfolder)
    #Write all prediction details to HTML files for each gene to be used as pop-up window
    for feature in pksnrpsvars.pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        if "PKS_AT" in pksnrpsvars.domainnamesdict[locus] or "AMP-binding" in pksnrpsvars.domainnamesdict[locus] or "A-OX" in pksnrpsvars.domainnamesdict[locus] or "CAL_domain" in pksnrpsvars.domainnamesdict[locus]:
            j = pksnrpsvars.domaindict[locus]
            nrat = 0
            nra = 0
            nrcal = 0
            for k in j:
                if k[0] == "PKS_AT":
                    nrat += 1
                    domainname = locus + "_AT" + str(nrat)
                    htmloutfile = open(options.substrspecsfolder + os.sep + domainname + ".html","w")
                    htmloutfile.write('<html>\n<head>\n<title>Prediction details</title>\n<STYLE type="text/css">\nbody{\n  text-align:left;\n  background-color:white;\n  font-family: Tahoma, sans-serif;\n  font-size: 0.8em;\n  color: #810E15;\n}\n</STYLE>\n</head>\n<body>')
                    htmloutfile.write(pksnrpsvars.minowa_pks_preds_details[domainname])
                    htmloutfile.write(pksnrpsvars.pks_code_preds_details[domainname])
                    htmloutfile.write("<b><i>Consensus Predictions: " + pksnrpsvars.consensuspreds[domainname] + "</b></i>")
                    htmloutfile.write('\n</body>\n</html>')
                    htmloutfile.close()
                if k[0] == "AMP-binding" or k[0] == "A-OX":
                    nra += 1
                    domainname = locus + "_A" + str(nra)
                    htmloutfile = open(options.substrspecsfolder + os.sep + domainname + ".html","w")
                    htmloutfile.write('<html>\n<head>\n<title>Prediction details</title>\n<STYLE type="text/css">\nbody{\n  text-align:left;\n  background-color:white;\n  font-family: Tahoma, sans-serif;\n  font-size: 0.8em;\n  color: #810E15;\n}\n</STYLE>\n</head>\n<body>')
                    htmloutfile.write(pksnrpsvars.nrps_svm_preds_details[domainname])
                    htmloutfile.write(pksnrpsvars.nrps_code_preds_details[domainname])
                    htmloutfile.write(pksnrpsvars.minowa_nrps_preds_details[domainname])
                    htmloutfile.write("<b><i>Consensus Prediction: '" + pksnrpsvars.consensuspreds[domainname] + "'</b></i>")
                    htmloutfile.write('\n</body>\n</html>')
                    htmloutfile.close()
                if k[0] == "CAL_domain":
                    nrcal += 1
                    domainname = locus + "_CAL" + str(nrcal)
                    htmloutfile = open(options.substrspecsfolder + os.sep + domainname + ".html","w")
                    htmloutfile.write('<html>\n<head>\n<title>Prediction details</title>\n<STYLE type="text/css">\nbody{\n  text-align:left;\n  background-color:white;\n  font-family: Tahoma, sans-serif;\n  font-size: 0.8em;\n  color: #810E15;\n}\n</STYLE>\n</head>\n<body>')
                    htmloutfile.write(pksnrpsvars.minowa_cal_preds_details[domainname])
                    htmloutfile.write('\n</body>\n</html>')
                    htmloutfile.close()
