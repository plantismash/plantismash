# Copyright (C) 2010-2012 Srikanth Duddela, Daniel Krug
# Helmholtz Institute for Pharmaceutical Research Saarland
# License: GNU Affero General Public License v3 or later

"""BiosynML output format module

"""
import logging
import warnings
from helperlibs.bio import seqio
from os import path
import os
import sys
import json
import math
import re
import shutil
from antismash import utils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
import copy
import xml.etree.ElementTree as etree
from xml.dom.minidom import Document
from time import gmtime, strftime

name = "BiosynML"
short_description = "BiosynML output"
priority = 1

def generate_template():
    #Generate BiosynML template
    doc = Document()
    list = []
    pseq=""
    name_gene = ""
    doc1 = Document()
    # Create the <wml> base element
    wml = doc.createElement("root")
    doc.appendChild(wml)
    filef=""
    qualifiers4biosynML = ""
    # Create the main <card> element
    maincard = doc.createElement("Header")
    dbxref=""
    genelist_motif=[]
    wml.appendChild(maincard)
    domainlist = doc.createElement("domainlist")
    wml.appendChild(domainlist) 

    maincard1 = doc.createElement("genelist")
    wml.appendChild(maincard1)
    cdsmotif = doc.createElement("motiflist")
    wml.appendChild(cdsmotif)   
    Sequencelist = doc.createElement("Sequencelist")
    wml.appendChild(Sequencelist)   
    kl=0;
    # writing header
    paragraph1 = doc.createElement("system")
    maincard.appendChild(paragraph1)
    ptext = doc.createTextNode("BiosynML")

    paragraph1.appendChild(ptext)

    paragraph1 = doc.createElement("version")
    maincard.appendChild(paragraph1)
    ptext = doc.createTextNode("1.0")
    paragraph1.appendChild(ptext)

    paragraph = doc.createElement("date")
    maincard.appendChild(paragraph)
    paragraph.setAttribute("action", "created")
    ptext = doc.createTextNode(strftime("%d.%m.%Y", gmtime()))
    paragraph.appendChild(ptext)

    paragraph = doc.createElement("date")
    maincard.appendChild(paragraph)
    paragraph.setAttribute("action", "modified")
    ptext = doc.createTextNode(strftime("%d.%m.%Y", gmtime()))
    paragraph.appendChild(ptext)
    paragraph1 = doc.createElement("author")
    maincard.appendChild(paragraph1)
    ptext = doc.createTextNode("antismash")
    paragraph1.appendChild(ptext)

    paragraph1 = doc.createElement("description")
    maincard.appendChild(paragraph1)
    ptext = doc.createTextNode("Output from antismash analysis")
    paragraph1.appendChild(ptext)
    return doc, wml, maincard, domainlist, maincard1, cdsmotif, Sequencelist, paragraph1, ptext

def write(seq_records, options):
    doc, wml, maincard, domainlist, maincard1, cdsmotif, Sequencelist, paragraph1, ptext = generate_template()
    seq_name = ""
    seq_record_copy = copy.deepcopy(seq_records)
    output_dir = options.outputfoldername
    sequenceListCounter= 0
    CDSProteinSequence=""
    motifCounter = 0
    geneCounter= 0
    modelCounter = 0
    nodeCount = 0
    CDSfeatureStart= 0
    domainPresence = 0
    # reading domains functions from xml to have uniform domain functions
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    data_file = (os.path.join(__location__,'domains.xml'))
    docbb = etree.parse((os.path.join(__location__,'buildingblocks.xml')))
    doc1 = etree.parse(data_file)
    root1 = doc1.getroot()
    seq_name = seq_record_copy[0].id
    logging.info("Started Writing BiosynML at "+strftime("%H:%M:%S", gmtime()))
    for rec in seq_record_copy:
        sequenceListCounter += 1
        withinclusterfeatures = []
        cdsmotifdict = {}
        clustera = []
        # reading CDS motifs data
        cdsmotifs = utils.get_all_features_of_type(rec, ["CDS_motif"])
        clusters = utils.get_cluster_features(rec)
        
        # copying prediction details
        substr_spec_preds = get_nrpspks_substr_spec_preds(rec)
        
        for cluster in clusters:
            smcog_descr = ""
            if cluster not in clustera:
                #transferring cluster data to models
                modelNumber = str(cluster.qualifiers['note'][0]).split(':')
                cluster_type = str(cluster.qualifiers['product'][0])
                clusterStartPosition = cluster.location.start
                clusterEndPosition = cluster.location.end
                modelCounter = modelCounter+1
                maincard14 = doc.createElement("model")
                wml.appendChild(maincard14)
                maincard14.setAttribute("id", str(modelCounter))    
                title = doc.createElement("title")
                maincard14.appendChild(title)
                ptext = doc.createTextNode("Biosynthetic pathway"+str(modelCounter))
                title.appendChild(ptext)
                confidence = doc.createElement("confidence")
                maincard14.appendChild(confidence)
                ptext = doc.createTextNode("Sequence-based prediction")
                confidence.appendChild(ptext)
                generator = doc.createElement("generator")
                maincard14.appendChild(generator)
                ptext = doc.createTextNode("antiSMASH")
                generator.appendChild(ptext)    
                label = doc.createElement("label")
                maincard14.appendChild(label)
                ptext = doc.createTextNode("model"+str(modelNumber[1]))
                label.appendChild(ptext)
                organism = doc.createElement("organism")
                maincard14.appendChild(organism)
                na3 = doc.createElement("name")
                organism.appendChild(na3)
                ptext = doc.createTextNode("")
                na3.appendChild(ptext)
                strain1 = doc.createElement("strain")
                organism.appendChild(strain1)
                ptext = doc.createTextNode(" ")
                strain1.appendChild(ptext)
                identifier1 = doc.createElement("identifier")
                organism.appendChild(identifier1)
                identifier1.setAttribute("source", " ")
                ptext = doc.createTextNode(" ")
                identifier1.appendChild(ptext)
                status1 = doc.createElement("status")
                organism.appendChild(status1)
                ptext = doc.createTextNode("")
                status1.appendChild(ptext)
                
                compound = doc.createElement("compound")
                maincard14.appendChild(compound)
                na2 = doc.createElement("name")
                compound.appendChild(na2)
                ptext = doc.createTextNode(" ")
                na2.appendChild(ptext)
                identifier2= doc.createElement("identifier")
                compound.appendChild(identifier2)
                identifier2.setAttribute("source", "")
                ptext = doc.createTextNode(" ")
                identifier2.appendChild(ptext)
                status = doc.createElement("status")
                compound.appendChild(status)
                ptext = doc.createTextNode(" ")
                status.appendChild(ptext)
                citation2 = doc.createElement("citation")
                compound.appendChild(citation2)
                citation2.setAttribute("type", "")
                ptext = doc.createTextNode("")
                citation2.appendChild(ptext)
                main1 = doc.createElement("genecluster")
                maincard14.appendChild(main1)    
                
                name_element = doc.createElement("name")
                main1.appendChild(name_element)
                ptext = doc.createTextNode("Biosynthetic genecluster")
                name_element.appendChild(ptext) 
                
                shortname = doc.createElement("shortname")
                main1.appendChild(shortname)
                ptext = doc.createTextNode("")
                shortname.appendChild(ptext)
                status = doc.createElement("status")
                main1.appendChild(status)
                ptext = doc.createTextNode("")
                status.appendChild(ptext)
                type1 = doc.createElement("type")
                main1.appendChild(type1)
                if(len(cluster_type)>2):
                    ptext = doc.createTextNode(cluster_type)
                    type1.appendChild(ptext)
                else:
                    ptext = doc.createTextNode(cluster_type)
                    type1.appendChild(ptext)
                identifier = doc.createElement("identifier")
                main1.appendChild(identifier)
                identifier.setAttribute("source", "")
                ptext = doc.createTextNode(" ")
                identifier.appendChild(ptext)
                
                citation = doc.createElement("citation")
                main1.appendChild(citation)
                citation.setAttribute("type", "")
                ptext = doc.createTextNode("")
                citation.appendChild(ptext)
                
                sequencelist = doc.createElement("sequence")
                main1.appendChild(sequencelist)
                sequencelist.setAttribute("source", "sequencelist")
                ptext = doc.createTextNode(str(sequenceListCounter))
                sequencelist.appendChild(ptext)     
                
                region = doc.createElement("region")
                main1.appendChild(region)
                start4xml1 = doc.createElement("begin") 
                region.appendChild(start4xml1)
                if(clusterStartPosition == 0):
                    clusterStartPosition= 1
                ptext = doc.createTextNode(str(clusterStartPosition))
                start4xml1.appendChild(ptext)  
                end4xml1 = doc.createElement("end") 
                region.appendChild(end4xml1)
                ptext = doc.createTextNode(str(clusterEndPosition))
                end4xml1.appendChild(ptext)
                #creating  Nodelist under models
                main = doc.createElement("nodelist")
                maincard14.appendChild(main)
                clustera.append(cluster)   
                # reading geneslist which are related to cluster 
                cdsfeatures = utils.get_cluster_cds_features(cluster,rec)
                for cdsfeature in cdsfeatures:
                    domianListForMotif= {}
                    countingA= 0
                    countingAT= 0
                    countingKR = 0
                    try:
                        geneName = utils.get_gene_id(cdsfeature)
                    except:
                        s= 1   
                    
                    cds_domains = []
                    o1 = 1
                    #adding genes to BiosynML file
                    geneCounter = geneCounter+1
                    paragrap11 = doc.createElement("gene")
                    maincard1.appendChild(paragrap11)
                 
                    paragrap11.setAttribute("id", str(geneCounter))
                    clas1 = doc.createElement("gene_name")
                    paragrap11.appendChild(clas1)
                    ptext = doc.createTextNode(geneName)
                    clas1.appendChild(ptext)
                    clas1 = doc.createElement("sequence")
                    paragrap11.appendChild(clas1)
                    clas1.setAttribute("source", "sequencelist")
                    ptext = doc.createTextNode(str(sequenceListCounter))
                    clas1.appendChild(ptext)
                    CDSfeatureStart= cdsfeature.location.start+1
                    clas11 = doc.createElement("gene_location")
                    paragrap11.appendChild(clas11)
                    CdsFeaturestart = CDSfeatureStart
                    CdsFeatureend = cdsfeature.location.end
                    temp = 0
                    if(CDSfeatureStart == 0):
                        CdsFeaturestart = 1
                    else:
                        CdsFeaturestart = CDSfeatureStart
                    if(cdsfeature.location.strand == -1):
                        temp = CdsFeaturestart
                        CdsFeaturestart = CdsFeatureend
                        CdsFeatureend = temp
                    start4xml1 = doc.createElement("begin") 
                    clas11.appendChild(start4xml1)
                    ptext = doc.createTextNode(str(CdsFeaturestart))
                    start4xml1.appendChild(ptext)  
                    end4xml1 = doc.createElement("end") 
                    clas11.appendChild(end4xml1)
                    ptext = doc.createTextNode(str(CdsFeatureend))
                    end4xml1.appendChild(ptext)  
                    gene_qualifiers = doc.createElement("gene_qualifiers")
                    paragrap11.appendChild(gene_qualifiers)
                    smcog_dict=[]
                    smc =""
                    domainPresence =0
                    smcog_score = ""
                    smcog_eval = ""
                    
                    for key in cdsfeature.qualifiers:
                        for value in cdsfeature.qualifiers[key]:
                            #search for genes with smCOG id and if the gene doesn't have domains, use this ID to generate domains from domainslist.xml
                            if('smCOG' in value):
                                if value not in smcog_dict:
                                    if('SMCOG' in value.replace("smCOG:","")):
                                        smcog_descr = value.partition("smCOG: ")[2].partition(":")[2].partition("(Score:")[0]
                                        smcog_score = value.partition("(Score: ")[2].partition(";")[0]
                                        smcog_eval = value.partition("E-value: ")[2].partition(");")[0]
                                        qualifier = doc.createElement("qualifier")
                                        gene_qualifiers.appendChild(qualifier)
                                        smc = value.replace("smCOG:","")
                                        ptext = doc.createTextNode(value.replace("smCOG:","").strip())
                                        qualifier.appendChild(ptext)
                                        qualifier.setAttribute("ori", "auto-annotation")
                                        qualifier.setAttribute("style", "genbank") 
                                        qualifier.setAttribute("name", "smCOG")
                                        smcog_dict.append(value)
                                        domainPresence = 1
                            else:
                                qualifier = doc.createElement("qualifier")
                                gene_qualifiers.appendChild(qualifier)
                                ptext = doc.createTextNode(str(value))
                                qualifier.appendChild(ptext)
                                qualifier.setAttribute("ori", "auto-annotation")
                                qualifier.setAttribute("style", "genbank") 
                                qualifier.setAttribute("name", str(key))
                                if 'translation' in str(key):
                                    CDSProteinSequence = value
                    
                    
                    clas1 = doc.createElement("operon")
                    paragrap11.appendChild(clas1)
                    ptext = doc.createTextNode("1")
                    clas1.appendChild(ptext)
                    sec_metQualifier = ""
                    
                    if(domainPresence == 1):
                        if 'sec_met' not in cdsfeature.qualifiers:
                            try:     
                                for child1 in root1:
                                    if(len(splitted)>1): ## splitted is undefined here!
                                        delimiters = "_"," ","-"
                                        regexPattern = '|'.join(map(re.escape, delimiters))
                                        keywords = re.split(regexPattern, smcog_descr)
                                        for func in keywords:
                                            if(len(func.strip())>1):
                                                if(func.lower().strip() in child1.find('chemistry').text.lower().strip()):
                                                    sec_metQualifier = child1.find('function').text.strip()
                                                    break  
                                                elif(func.lower().strip() in child1.find('substrate').text.lower().strip()):
                                                    sec_metQualifier = child1.find('function').text.strip()
                                                    break
                                                elif(func.lower().strip() in child1.find('name')):
                                                    sec_metQualifier = child1.find('function').text.strip()
                                                    break
                                                elif(func.lower().strip() in  child1.find('keywords').text.lower().strip()):
                                                    sec_metQualifier = child1.find('function').text.strip()
                                                    break
                                        if(len(sec_metQualifier)>0):
                                            break
                            except:
                                s=1
                            #modify local copy of the seq object to include domains generated from smCOG
                            if(len(sec_metQualifier)>0):
                                cdsfeature.qualifiers['sec_met'] =["NRPS/PKS Domain: "+sec_metQualifier+" (1-"+str(abs(CdsFeaturestart-CdsFeatureend)/3)+"). E-value: "+smcog_eval+". Score: "+smcog_score+";"]
                    #extract domains and their meta information   
                    if 'sec_met' in cdsfeature.qualifiers:
                        cds_domains = [qual for qual in cdsfeature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qual]
                        for domain in cds_domains:
                            
                            hit_id =  domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[0]
                            domstart = domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[2].partition("-")[0]
                            domend = domain.partition("NRPS/PKS Domain: ")[2].partition("). ")[0].rpartition("-")[2]
                           
                            d = hit_id.split('_')
                            func=""
                            if(len(d) == 1):
                                func = hit_id
                            else:
                                if(d[1].strip().lower() == "docking"):
                                    func = "Dck-pks"
                                else:
                                    func = d[1]
                            
                            if(func.strip() == "domain"):
                                func = d[0]
                            if(func.strip() == "TD"):
                                func = "Red" 
                            geb = ""
                            domainPresence =1
                            domainStart = int(domstart.strip())
                            domainEnd = int(domend.strip())
                            pdomainStart = int(domstart.strip())
                            pdomainEnd = int(domend.strip())
                            domainProteinSeq = CDSProteinSequence[pdomainStart:pdomainEnd]
                            cla = ""
                            che = ""
                            subs = ""
                            subtype = " "
                            dom =""
                            
                           
                            for child1 in root1:
                                geb = ""
                                cla = ""
                                che = ""
                                subs = ""
                                k2 = ""
                                subfunction=" "
                                get_it = ""
                                try:
                                    get_it = child1.find('keywords').text
                                except:
                                    get_it = ""
                                if(get_it == ""):
                                    c =  ""
                                else:
                                    splitted = []
                                    try:
                                        splitted1 = child1.find('keywords').text.split(',')
                                        splitted = [ i for i in splitted1 ]
                                    except:
                                        et_it = ""
                                    stt=""
                                    for i in splitted:
                                        if(stt.strip() == ""):
                                            stt = i.lower().strip()+"=="+func.lower().strip()
                                        else:
                                            stt += " or "+i.lower().strip()+"=="+func.lower().strip()
                                    if(len(splitted)>1):
                                        
                                        if (func.strip() in splitted):
                                            
                                            geb = child1.find('default_context').text.strip()
                                            cla = child1.find('class').text.strip()
                                            che = child1.find('chemistry').text.strip()
                                            subs = child1.find('substrate').text.strip()
                                            k2 = child1.find('function').text.strip()
                                            try:
                                            
                                                subfunction = child1.find('subtypes').text.strip()
                                            except:
                                                c =  "" ""
                                            break
                                        
                                        elif(child1.find('chemistry').text.lower().strip() == func.lower().strip()):    
                                            geb = child1.find('default_context').text.strip()
                                            cla = child1.find('class').text.strip()
                                            che = child1.find('chemistry').text.strip()
                                            subs = child1.find('substrate').text.strip()
                                            k2 = child1.find('function').text.strip()
                                            try:
                                            
                                                subfunction = child1.find('subtypes').text.strip()
                                            except:
                                                c =  ""
                                            break                             
            
                                        elif(child1.find('substrate').text.lower().strip() == func.lower().strip()):
                                            geb = child1.find('default_context').text.strip()
                                            cla = child1.find('class').text.strip()
                                            che = child1.find('chemistry').text.strip()
                                            subs = child1.find('substrate').text.strip()
                                            k2 = child1.find('function').text.strip()
                                            try:
                                            
                                                subfunction = child1.find('subtypes').text.strip()
                                            except:
                                                c =  ""
                                            break
                                        elif(child1.find('name').text.lower().strip() == func.lower().strip()):    
                                            geb = child1.find('default_context').text.strip()
                                            cla = child1.find('class').text.strip()
                                            che = child1.find('chemistry').text.strip()
                                            subs = child1.find('substrate').text.strip()
                                            k2 = child1.find('function').text.strip()
                                            try:
                                            
                                                subfunction = child1.find('subtypes').text.strip()
                                            except:
                                                c =  ""
                                            break
                                        elif(child1.find('function').text.lower().strip() == func.lower().strip()):
                                            geb = child1.find('default_context').text.strip()
                                            cla = child1.find('class').text.strip()
                                            che = child1.find('chemistry').text.strip()
                                            subs = child1.find('substrate').text.strip()
                                            k2 = child1.find('function').text.strip()
                                            try:
                                            
                                                subfunction = child1.find('subtypes').text.strip()
                                            except:
                                                c =  ""
                                            break                                
                                                 
                                    else:
                                        if(child1.find('keywords').text == func.strip()):
                                            
                                            geb = child1.find('default_context').text
                                            cla = child1.find('class').text.strip()
                                            che = child1.find('chemistry').text.strip()
                                            subs = child1.find('substrate').text.strip()
                                            k2 = child1.find('function').text.strip()
                                            try:
                                                subfunction = child1.find('subtypes').text.strip()
                                            except:
                                                c =  ""
                                            break
                                        
                                        elif(child1.find('chemistry').text.lower().strip() == func.lower().strip()):
                                            geb = child1.find('default_context').text.strip()
                                            cla = child1.find('class').text.strip()
                                            che = child1.find('chemistry').text.strip()
                                            subs = child1.find('substrate').text.strip()
                                            k2 = child1.find('function').text.strip()
                                            try:
                                            
                                                subfunction = child1.find('subtypes').text.strip()
                                            except:
                                                c =  ""
                                            break                             
            
                                        elif(child1.find('substrate').text.lower().strip() == func.lower().strip()):
                                            geb = child1.find('default_context').text.strip()
                                            cla = child1.find('class').text.strip()
                                            che = child1.find('chemistry').text.strip()
                                            subs = child1.find('substrate').text.strip()
                                            k2 = child1.find('function').text.strip()
                                            try:
                                            
                                                subfunction = child1.find('subtypes').text.strip()
                                            except:
                                                c =  ""
                                            break
                                        elif(child1.find('name').text.lower().strip() == func.lower().strip()):    
                                            geb = child1.find('default_context').text.strip()
                                            cla = child1.find('class').text.strip()
                                            che = child1.find('chemistry').text.strip()
                                            subs = child1.find('substrate').text.strip()
                                            k2 = child1.find('function').text.strip()
                                            try:
                                            
                                                subfunction = child1.find('subtypes').text.strip()
                                            except:
                                                c =  ""
                                            break
                                        
                                        elif(child1.find('function').text.strip() == func.strip()):
                                            geb = child1.find('default_context').text.strip()
                                            cla = child1.find('class').text.strip()
                                            che = child1.find('chemistry').text.strip()
                                            subs = child1.find('substrate').text.strip()
                                            k2 = child1.find('function').text.strip()
                                            try:
                                            
                                                subfunction = child1.find('subtypes').text.strip()
                                            except:
                                                c =  ""
                                            break
                          
                            if("Nterm" in func.strip()):
                                pass
                            else:
                                if not subfunction:
                                    subfunction = " "
                                nodeCount = nodeCount+1
                                paragrap1 = doc.createElement("node")
                                main.appendChild(paragrap1)
                                paragrap1.setAttribute("id", str(nodeCount))
                                clas = doc.createElement("class")
                                paragrap1.appendChild(clas)
                                ptext = doc.createTextNode(cla)
                                clas.appendChild(ptext)
                                clas = doc.createElement("context")
                                paragrap1.appendChild(clas)
                                ptext = doc.createTextNode(str(geb))
                                clas.appendChild(ptext)
                                #Adding   domainlist
                                domain = doc.createElement("domain")
                                domainlist.appendChild(domain)
                                domain.setAttribute("id", str(nodeCount))
                                nodeid = doc.createElement("nodeid") 
                                domain.appendChild(nodeid)
                                ptext = doc.createTextNode(str(nodeCount))
                                nodeid.appendChild(ptext)
                                function = doc.createElement("function") 
                                domain.appendChild(function)
                                ptext = doc.createTextNode(str(k2))
                                function.appendChild(ptext)
                                if str(k2.strip()) != "KR":
                                    subtypes = doc.createElement("subtype") 
                                    domain.appendChild(subtypes)
                                    lis = subfunction.split(',')
                                    if(len(lis) <= 1):
                                        ptext = doc.createTextNode(str(subfunction))
                                        subtypes.appendChild(ptext)
                                activity5 = doc.createElement("dstatus") 
                                domain.appendChild(activity5)
                                ptext = doc.createTextNode("active")
                                activity5.appendChild(ptext) 
                                if (func.strip() == "AT"):
                                    countingAT = countingAT +1
                                    dom = "AT"+str(countingAT)
                                elif (func.strip() == "AMP-binding"):
                                    countingA = countingA +1
                                    dom = "A"+str(countingA)
                                elif (func.strip() == "KR"):
                                    countingKR = countingKR +1
                                    dom = "KR"+str(countingKR)                                    
                                else:
                                    dom = str(k2)
                                label = doc.createElement("label") 
                                domain.appendChild(label)
                                ptext = doc.createTextNode(str(k2))
                                label.appendChild(ptext)
                                
                                chemistry = doc.createElement("chemistry") 
                                domain.appendChild(chemistry)
                                ptext = doc.createTextNode(che)
                                chemistry.appendChild(ptext)
                                
                                substrate = doc.createElement("substrate") 
                                domain.appendChild(substrate)
                                ptext = doc.createTextNode(subs)
                                substrate.appendChild(ptext)
                                location = doc.createElement("location")
                                domain.appendChild(location) 
                                gene1 = doc.createElement("gene") 
                                location.appendChild(gene1)   
                                name_element = doc.createElement("geneid") 
                                gene1.appendChild(name_element)
                                name_element.setAttribute("source", "genelist")
                                ptext = doc.createTextNode(str(geneCounter))
                                name_element.appendChild(ptext)
                                #mapping domain positions related to gene positions
                                if(cdsfeature.location.strand == -1):
                                    if(int(domainStart) == 0):
                                        domainStart = 1
                                    else:
                                        domainStart = domainStart*3
                                    domainEnd = domainEnd*3    
                                    tempDomain = 0
                                    tempPDomain = 0
                                    tempDomain =  domainStart
                                    domainStart = domainEnd
                                    domainEnd = tempDomain
                                    tempPDomain =  pdomainStart
                                    pdomainStart = pdomainEnd
                                    pdomainEnd = tempPDomain
                                    positionOfDomain = str(domainEnd)+"-"+str(domainStart)
                                    domianListForMotif[str(nodeCount)] = positionOfDomain
                                else:
                                    domainStart = domainStart*3
                                    domainEnd = domainEnd*3    
                                    positionOfDomain = str(domainStart)+"-"+str(domainEnd)
                                    domianListForMotif[str(nodeCount)] = positionOfDomain  
                                try:
                                    position = doc.createElement("position") 
                                    gene1.appendChild(position)
                                    
                                    start4xml1 = doc.createElement("begin") 
                                    position.appendChild(start4xml1)
                                    ptext = doc.createTextNode(str(domainStart))
                                    start4xml1.appendChild(ptext)  
                                    end4xml1 = doc.createElement("end") 
                                    position.appendChild(end4xml1)
                                    ptext = doc.createTextNode(str(domainEnd))
                                    end4xml1.appendChild(ptext)
                                     
                                    protein1 = doc.createElement("protein") 
                                    location.appendChild(protein1)
                                    name11 = doc.createElement("name") 
                                    protein1.appendChild(name11)
                                    ptext = doc.createTextNode(geneName.strip().upper()+"_"+str(k2))
                                    name11.appendChild(ptext)
                                        
                                    sequence = doc.createElement("sequence") 
                                    protein1.appendChild(sequence)
                                    ptext = doc.createTextNode(domainProteinSeq)
                                    sequence.appendChild(ptext)
                                    
                                    position1 = doc.createElement("position") 
                                    protein1.appendChild(position1)
                                    
                                    
                                    start4xml = doc.createElement("begin") 
                                    position1.appendChild(start4xml)
                                    ptext = doc.createTextNode(str(pdomainStart)  )
                                    start4xml.appendChild(ptext)  
                                    end4xml = doc.createElement("end") 
                                    position1.appendChild(end4xml)
                                    ptext = doc.createTextNode(str(pdomainEnd))
                                    end4xml.appendChild(ptext)
                                except:
                                    c =  ""
                                predictionFileData = ""
                                                            
                                #reading predictions
                                domainname = geneName.strip()+"_"+dom
                                domainIdentifier = " ".join(re.findall("[a-zA-Z]+", dom))
                                kr_activity_preds = ""
                                kr_stereo_preds = ""
                                # reading predictions based on the "domainnames" as keys. When domains are inserted by biosynML based on smCOG analysis, these domains doesn't have predictions, hence the error is caught
                                try:
                                    if "KR" == domainIdentifier.strip():
                                        kr_activity_preds= substr_spec_preds.kr_activity_preds[domainname]
                                        kr_stereo_preds= substr_spec_preds.kr_stereo_preds[domainname]
                                        subtypes = doc.createElement("subtype") 
                                        domain.appendChild(subtypes)
                                        ptext = doc.createTextNode(str(kr_stereo_preds.strip()))
                                        subtypes.appendChild(ptext)                                    
                                        predictionFileData = kr_activity_preds+"\n"+ kr_stereo_preds
                                    if "AT" == domainIdentifier.strip():
                                        substrspecminowa = substr_spec_preds.minowa_pks_preds[domainname]
                                        substrspecpkssig = substr_spec_preds.pks_code_preds[domainname]
                                        substrspecconsensus = substr_spec_preds.consensuspreds[domainname]
                                        
                                        predictionFileData = substrspecminowa+"\n"+substrspecpkssig+"\n"+"Consensus Prediction:"+substrspecconsensus
                                        #print  predictionFileData
                                    if "A" == domainIdentifier.strip():
                                        substrspecminowa = substr_spec_preds.minowa_nrps_preds[domainname]
                                        substrspecnrpspredictor = [substr_spec_preds.nrps_code_preds[domainname], substr_spec_preds.nrps_svm_preds[domainname]]
                                        substrspecconsensus = substr_spec_preds.consensuspreds[domainname]
                                        predictionFileData = substrspecminowa+"\n"+",".join(substrspecnrpspredictor)+"\n"+"Consensus Prediction:"+substrspecconsensus
                                    if "CAL_domain" in domainIdentifier.strip():
                                        substrspecminowa = substr_spec_preds.minowa_cal_preds[domainname]
                                        substrspecconsensus = substr_spec_preds.consensuspreds[domainname]
                                        predictionFileData = substrspecminowa+"\n"+"Consensus Prediction:"+substrspecconsensus
                                except KeyError:
                                    logging.info("No predictions found for "+domainname)
                                if(len(predictionFileData)>3):
                                    predictionData = re.sub('<[^<]+?>', '', predictionFileData)
                                    predictions = predictionData.rstrip().split('\n')
                                    skip = 0
                                    NRPSPredictor2_SVM=0
                                    Minowa_HMM_method_A = 0
                                    NRPSPredictor2_Stachelhaus = 0
                                    PKS_Active_Signature = 0
                                    Minowa_HMM_method_AT=0
                                    Consensus_Prediction =0
                                    predictionCounter = 0
                                    for i in predictions:
                                        if(i.strip() == "NRPSPredictor2 SVM prediction details:"):
                                            NRPSPredictor2_SVM = 1   
                                            Minowa_HMM_method_A = 0
                                            NRPSPredictor2_Stachelhaus =0
                                            PKS_Active_Signature = 0
                                            Minowa_HMM_method_AT=0
                                            Consensus_Prediction =0
                                            predictionCounter = predictionCounter+1
                                            prediction = doc.createElement("prediction")
                                            domain.appendChild(prediction)
                                            prediction.setAttribute("id", str(predictionCounter))
                                            prediction.setAttribute("type", "substrate")
                                            clas1 = doc.createElement("Method")
                                            prediction.appendChild(clas1)
                                            ptext = doc.createTextNode(i.strip().replace(':',''))
                                            clas1.appendChild(ptext)
                                        elif ('Minowa HMM method A-domain' in i.strip()):
                                            Minowa_HMM_method_A = 1
                                            NRPSPredictor2_SVM = 0
                                            NRPSPredictor2_Stachelhaus =0
                                            PKS_Active_Signature = 0
                                            Minowa_HMM_method_AT=0
                                            Consensus_Prediction =0
                                            predictionCounter = predictionCounter+1
                                            prediction = doc.createElement("prediction")
                                            domain.appendChild(prediction)
                                            prediction.setAttribute("id", str(predictionCounter))
                                            prediction.setAttribute("type", "substrate")
                                            clas1 = doc.createElement("Method")
                                            prediction.appendChild(clas1)
                                            ptext = doc.createTextNode("Minowa HMM method A-domain")
                                            clas1.appendChild(ptext)
                                        elif ('Consensus Prediction' in i.strip()):
                                            Minowa_HMM_method_A = 0
                                            NRPSPredictor2_SVM = 0
                                            NRPSPredictor2_Stachelhaus =0
                                            PKS_Active_Signature = 0
                                            Minowa_HMM_method_AT=0
                                            Consensus_Prediction =1
                                            predictionCounter = predictionCounter+1
                                            prediction = doc.createElement("prediction")
                                            domain.appendChild(prediction)
                                            prediction.setAttribute("id", str(predictionCounter))
                                            prediction.setAttribute("type", "substrate")
                                            clas1 = doc.createElement("Method")
                                            prediction.appendChild(clas1)
                                            ptext = doc.createTextNode("Consensus Predictions")
                                            clas1.appendChild(ptext)
                                            Consensus_Predictions = i.split(':')
                                            function = doc.createElement("readout") 
                                            prediction.appendChild(function)
                                            ptext = doc.createTextNode("prediction")
                                            function.appendChild(ptext)
                                            function1 = doc.createElement("result") 
                                            prediction.appendChild(function1)
                                            ptext = doc.createTextNode(str(Consensus_Predictions[1].replace("'", "")))
                                            function1.appendChild(ptext)
                                            # Consensus Predictions for nodelist
                                            buildingblock = doc.createElement("buildingblock")
                                            paragrap1.appendChild(buildingblock)   
                                            moiety = doc.createElement("moiety") 
                                            buildingblock.appendChild(moiety)
                                            moiety.setAttribute("ratio", "1")
                                            name_element = doc.createElement("name") 
                                            moiety.appendChild(name_element)
                                            ptext = doc.createTextNode(str(Consensus_Predictions[1].strip().replace("'", "")))
                                            name_element.appendChild(ptext)
                                            category = doc.createElement("category") 
                                            moiety.appendChild(category)
                                            ptext = doc.createTextNode(geb.replace(',','!'))
                                            category.appendChild(ptext)
                                            evidence = doc.createElement("evidence") 
                                            moiety.appendChild(evidence)
                                            ptext = doc.createTextNode("prediction")
                                            evidence.appendChild(ptext)
                                            confidence = doc.createElement("confidence") 
                                            moiety.appendChild(confidence)
                                            ptext = doc.createTextNode("Sequence-based prediction")
                                            confidence.appendChild(ptext)
                                            parent_name=""
                                            
                                            rootbb = docbb.getroot()
                                            try:
                                                for childbb in rootbb:
                                                    if(childbb.find('code').text.lower().strip() == str(Consensus_Predictions[1].strip())):
                                                        try:
                                                            parent_name = childbb.find('parent').text.strip()
                                                            
                                                        except:
                                                            pass
                                                            
                                                            
                                            except:
                                                pass
                                                
                                            
                                            parent = doc.createElement("parent") 
                                            moiety.appendChild(parent)
                                            par_name = doc.createElement("name") 
                                            parent.appendChild(par_name)
                                            ptext = doc.createTextNode(parent_name) 
                                            par_name.appendChild(ptext)
                                            par_trans = doc.createElement("transform") 
                                            parent.appendChild(par_trans)
                                            ptext = doc.createTextNode("") 
                                            par_trans.appendChild(ptext)   
                                        elif ('Minowa HMM method AT-domain' in i.strip()):
                                            Minowa_HMM_method_A = 0
                                            NRPSPredictor2_SVM = 0
                                            NRPSPredictor2_Stachelhaus =0
                                            PKS_Active_Signature = 0
                                            Consensus_Prediction =0
                                            Minowa_HMM_method_AT=1
                                            predictionCounter = predictionCounter+1
                                            prediction = doc.createElement("prediction")
                                            domain.appendChild(prediction)
                                            prediction.setAttribute("id", str(predictionCounter))
                                            prediction.setAttribute("type", "substrate")
                                            clas1 = doc.createElement("Method")
                                            prediction.appendChild(clas1)
                                            ptext = doc.createTextNode("Minowa HMM method AT-domain")
                                            clas1.appendChild(ptext)
                                        elif ('PKS Active Site Signature method' in i.strip()):
                                            Minowa_HMM_method_A = 0
                                            NRPSPredictor2_SVM = 0
                                            NRPSPredictor2_Stachelhaus =0
                                            PKS_Active_Signature = 1
                                            Consensus_Prediction =0
                                            Minowa_HMM_method_AT=0
                                            predictionCounter = predictionCounter+1
                                            prediction = doc.createElement("prediction")
                                            domain.appendChild(prediction)
                                            prediction.setAttribute("id", str(predictionCounter))
                                            prediction.setAttribute("type", "substrate")
                                            clas1 = doc.createElement("Method")
                                            prediction.appendChild(clas1)
                                            ptext = doc.createTextNode("PKS Active Site Signature method")
                                            clas1.appendChild(ptext)
                                        elif(i.strip() == "NRPSPredictor2 Stachelhaus code prediction:"):
                                            NRPSPredictor2_Stachelhaus =1
                                            Minowa_HMM_method_A = 0
                                            NRPSPredictor2_SVM = 0
                                            PKS_Active_Signature = 0
                                            Consensus_Prediction =0
                                            Minowa_HMM_method_AT=0
                                            predictionCounter = predictionCounter+1
                                            prediction = doc.createElement("prediction")
                                            domain.appendChild(prediction)
                                            prediction.setAttribute("id", str(predictionCounter))
                                            prediction.setAttribute("type", "substrate")
                                            clas1 = doc.createElement("Method")
                                            prediction.appendChild(clas1)
                                            ptext = doc.createTextNode(i.strip().replace(':',''))
                                            clas1.appendChild(ptext)
                                        elif(NRPSPredictor2_SVM == 1):
                                            if(len(i.strip())>1):
                                                if i.find(":") == -1:
                                                    function1 = doc.createElement("result") 
                                                    prediction.appendChild(function1)
                                                    ptext = doc.createTextNode(str(i.strip().replace(',',';')))
                                                    function1.appendChild(ptext) 
                                                else:
                                                    function = doc.createElement("readout") 
                                                    prediction.appendChild(function)
                                                    ptext = doc.createTextNode(i.strip().replace(':',''))
                                                    function.appendChild(ptext)
                                        elif(Minowa_HMM_method_A == 1):
                                            if(len(i.strip())>1):
                                                if i.find(":") == -1:
                                                    function = doc.createElement("readout") 
                                                    prediction.appendChild(function)
                                                    ptext = doc.createTextNode("Substrate: "+"".join(re.findall("[a-zA-Z]+", i.strip())))
                                                    function.appendChild(ptext)
                                                    function1 = doc.createElement("result") 
                                                    prediction.appendChild(function1)
                                                    ptext = doc.createTextNode("Score: "+str('.'.join(re.findall('\d+', i.strip()))))
                                                    function1.appendChild(ptext)
                                        elif(NRPSPredictor2_Stachelhaus == 1):
                                            if(len(i.strip())>1):
                                                function = doc.createElement("readout") 
                                                prediction.appendChild(function)
                                                ptext = doc.createTextNode("prediction")
                                                function.appendChild(ptext)
                                                function1 = doc.createElement("result") 
                                                prediction.appendChild(function1)
                                                ptext = doc.createTextNode(str(i).strip())
                                                function1.appendChild(ptext)
                                        elif(PKS_Active_Signature == 1):
                                            if(len(i.strip())>1):
                                                if(i.find('AT-domain substrate specificity prediction top hits')==-1):
                                                    
                                                    kk = i.split(':')
                                                    if(len(kk))>1:
                                                        function = doc.createElement("readout") 
                                                        prediction.appendChild(function)
                                                        ptext = doc.createTextNode(kk[0].strip())
                                                        function.appendChild(ptext)
                                                        function1 = doc.createElement("result") 
                                                        prediction.appendChild(function1)
                                                        ptext = doc.createTextNode(str(kk[1].strip()))
                                                        function1.appendChild(ptext)
                                        elif(Minowa_HMM_method_AT == 1):
                                            if(len(i.strip())>1):
                                                if i.find(":") == -1:
                                                    function = doc.createElement("readout") 
                                                    prediction.appendChild(function)
                                                    ptext = doc.createTextNode("Substrate: "+"".join(re.findall("[a-zA-Z]+", i.strip())))
                                                    function.appendChild(ptext)
                                                    function1 = doc.createElement("result") 
                                                    prediction.appendChild(function1)
                                                    ptext = doc.createTextNode("Score: "+str('.'.join(re.findall('\d+', i.strip()))))
                                                    function1.appendChild(ptext)                                
                                
                        #Adding   Motiflist                  
                    for cds_motif in cdsmotifs:
                        domainIDForMotif = 0
                        PositionStart = 0
                        PositionEnd = 0
                        MotifPositionStart = 0
                        MotifPositionEnd = 0                        
                        if cdsfeature.location.start <= cds_motif.location.start <= cdsfeature.location.end:
                            motifCounter = motifCounter+1
                            try:
                                if(cdsfeature.location.strand == -1):
                                    MotifPositionEnd = abs(cds_motif.location.start- cdsfeature.location.end)
                                    MotifPositionStart = abs(cds_motif.location.end- cdsfeature.location.end)
                                    for key, value in domianListForMotif.iteritems():
                                        positionsOfDomain = value.split('-')
                                        PositionStart = CDSfeatureStart+int(positionsOfDomain[0].strip())+10
                                        PositionEnd = CDSfeatureStart+int(positionsOfDomain[1].strip())
                                        if PositionStart <= cds_motif.location.start and PositionEnd >= cds_motif.location.end:
                                            domainIDForMotif = key
                                else:
                                    MotifPositionStart = abs(cds_motif.location.start- CDSfeatureStart)
                                    MotifPositionEnd = abs(cds_motif.location.end- CDSfeatureStart)
                                    for key, value in domianListForMotif.iteritems():
                                        positionsOfDomain = value.split('-')
                                        PositionStart = CDSfeatureStart+int(positionsOfDomain[0].strip())+10
                                        PositionEnd = CDSfeatureStart+int(positionsOfDomain[1].strip())
                                        if PositionStart <= cds_motif.location.start and PositionEnd >= cds_motif.location.end:
                                            domainIDForMotif = key
                            except:
                                no_domains= "-1"
                            CDS_motif = doc.createElement("motif")
                            cdsmotif.appendChild(CDS_motif)
                            CDS_motif.setAttribute("id", str(motifCounter))
                            CDS_motif.setAttribute("geneID", str(geneCounter))
                            CDS_motif.setAttribute("domainID", str(domainIDForMotif))
                            motif_name = doc.createElement("motif_name")
                            CDS_motif.appendChild(motif_name)
                            motifname = cds_motif.qualifiers['note'][0].split('(')
                            
                            getMotifName =""
                            if(len(motifname)>1):
                                motifName = motifname[0].split(':')
                                if(len(motifName)>1):
                                    getMotifName  = motifName[1].strip()
                                else:
                                    getMotifName  = motifname[0].strip()
                            else:
                                getMotifName= cds_motif.qualifiers['note'][0]
                            
                            ptext = doc.createTextNode(str(getMotifName.replace(',',';')))
                            motif_name.appendChild(ptext)
                            if CDSfeatureStart<= cds_motif.location.start and cdsfeature.location.end>= cds_motif.location.end:
                                motif_type = doc.createElement("motif_type")
                                CDS_motif.appendChild(motif_type)
                                ptext = doc.createTextNode("BiosynML motif")
                                motif_type.appendChild(ptext)
                            else:
                                motif_type = doc.createElement("motif_type")
                                CDS_motif.appendChild(motif_type)
                                ptext = doc.createTextNode("motif")
                                motif_type.appendChild(ptext)
                                
                            motif_seq = doc.createElement("sequence")
                            CDS_motif.appendChild(motif_seq)    
                            ptext = doc.createTextNode(str(sequenceListCounter))#TODO
                            motif_seq.appendChild(ptext)
                            motif_seq.setAttribute("source", "sequencelist")
                            motif_location = doc.createElement("motif_location")
                            CDS_motif.appendChild(motif_location)
                            if(cds_motif.location.strand == -1):
                                tempMotif = MotifPositionStart
                                MotifPositionStart = MotifPositionEnd
                                MotifPositionEnd = tempMotif
                            begin = doc.createElement("begin")
                            motif_location.appendChild(begin)
                            ptext = doc.createTextNode(str(MotifPositionStart+1))
                            begin.appendChild(ptext)
                            end = doc.createElement("end")
                            motif_location.appendChild(end)
                            ptext = doc.createTextNode(str(MotifPositionEnd))
                            end.appendChild(ptext)
                            motif_qualifiers = doc.createElement("motif_qualifiers")
                            CDS_motif.appendChild(motif_qualifiers)   
                            for key in cds_motif.qualifiers:
                                for value in cds_motif.qualifiers[key]:
                                    qualifier = doc.createElement("qualifier")
                                    motif_qualifiers.appendChild(qualifier)
                                    ptext = doc.createTextNode(value.replace(',',';'))
                                    qualifier.appendChild(ptext)
                                    qualifier.setAttribute("ori", "auto-annotation")
                                    qualifier.setAttribute("style", "genbank")
                                    qualifier.setAttribute("name", "note")
                                    qualifier = doc.createElement("qualifier")
                            qualifier = doc.createElement("qualifier")
                            motif_qualifiers.appendChild(qualifier)
                            ptext = doc.createTextNode("CDS_motif")
                            qualifier.appendChild(ptext)
                            qualifier.setAttribute("ori", "auto-annotation")
                            qualifier.setAttribute("style", "genbank")
                            qualifier.setAttribute("name", "NCBI Feature Key")
        #Adding  seq information  
        dbxref = get_db_xref_feature(rec)
        sourceXml = doc.createElement("source")
        Sequencelist.appendChild(sourceXml)
        sourceXml.setAttribute("id", str(sequenceListCounter))
        paragrap1 = doc.createElement("label")
        sourceXml.appendChild(paragrap1)
        ptext = doc.createTextNode(str(rec.id))
        paragrap1.appendChild(ptext)
        try:
            if(len(dbxref)>1):
                paragrap1 = doc.createElement("db_xref")
                sourceXml.appendChild(paragrap1)
                ptext = doc.createTextNode(str(dbxref))
                paragrap1.appendChild(ptext)
        except:
            ds = ""
                
        clas = doc.createElement("Sequence")
        sourceXml.appendChild(clas)
        ptext = doc.createTextNode(str(rec.seq))
        clas.appendChild(ptext)
    logging.info("Finished writing BiosynML "+strftime("%H:%M:%S", gmtime()))
    try:
        f = open(output_dir +"/"+"biosynML.xml", "w")
        f.write(doc.toprettyxml(indent="  "))
    except:
        got_exception =  ""
    

def get_source_features(seq_record):
    '''Return all source features for a seq_record'''
    return utils.get_all_features_of_type(seq_record, 'source')


def get_db_xref_feature(seq_record):
    '''Return all db_xref features for a seq_record'''
    cds_list1 = get_source_features(seq_record)
    dbxrefInfomation = ''
    for feature in cds_list1:
        if feature.qualifiers.has_key('db_xref'):
            for db_xref in feature.qualifiers.get('db_xref'):
                if 'myxobase' in db_xref:
                    dbxrefInfomation = db_xref
    return dbxrefInfomation

def get_nrpspks_substr_spec_preds(seq_record):
    substr_spec_preds = utils.Storage()
    substr_spec_preds.consensuspreds = {}
    substr_spec_preds.nrps_svm_preds = {}
    substr_spec_preds.nrps_code_preds = {}
    substr_spec_preds.minowa_nrps_preds = {}
    substr_spec_preds.pks_code_preds = {}
    substr_spec_preds.minowa_pks_preds = {}
    substr_spec_preds.minowa_cal_preds = {}
    substr_spec_preds.kr_activity_preds = {}
    substr_spec_preds.kr_stereo_preds = {}
    features = utils.get_cds_features(seq_record)
    for feature in features:
        nrat, nra, nrcal, nrkr = 0, 0, 0, 0
        if 'sec_met_predictions' in feature.qualifiers:
            domains = [qualifier for qualifier in feature.qualifiers['sec_met_predictions'] if "NRPS/PKS Domain: " in qualifier]
            for domain in domains:
                if "AMP-binding" in domain or "A-OX" in domain:
                    nra += 1
                    domainname = utils.get_gene_id(feature) + "_A" + str(nra)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    nrps_svm_pred = predictionstext.partition(" (NRPSPredictor2 SVM)")[0]
                    nrps_code_pred = predictionstext.partition(" (NRPSPredictor2 SVM), ")[2].partition(" (Stachelhaus code)")[0]
                    minowa_nrps_pred = predictionstext.partition("(Stachelhaus code), ")[2].partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.nrps_svm_preds[domainname] = nrps_svm_pred
                    substr_spec_preds.nrps_code_preds[domainname] = nrps_code_pred
                    substr_spec_preds.minowa_nrps_preds[domainname] = minowa_nrps_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "PKS_AT" in domain:
                    nrat += 1
                    domainname = utils.get_gene_id(feature) + "_AT" + str(nrat)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    pks_code_pred = predictionstext.partition(" (PKS signature)")[0]
                    minowa_pks_pred = predictionstext.partition("(PKS signature), ")[2].partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.pks_code_preds[domainname] = pks_code_pred
                    substr_spec_preds.minowa_pks_preds[domainname] = minowa_pks_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "CAL_domain" in domain:
                    nrcal += 1
                    domainname = utils.get_gene_id(feature) + "_CAL" + str(nrcal)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    minowa_cal_pred = predictionstext.partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.minowa_cal_preds[domainname] = minowa_cal_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "PKS_KR" in domain:
                    nrkr += 1
                    domainname = utils.get_gene_id(feature) + "_KR" + str(nrkr)
                    activityprediction = domain.partition("Predicted KR activity: ")[2].partition(";")[0]
                    stereoprediction = domain.partition("Predicted KR stereochemistry: ")[2].partition(";")[0]
                    substr_spec_preds.kr_activity_preds[domainname] = activityprediction
                    substr_spec_preds.kr_stereo_preds[domainname] = stereoprediction
   
    return substr_spec_preds