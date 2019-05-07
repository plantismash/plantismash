# vim: set fileencoding=utf-8 :
#
#
# Copyright (C) 2014 Tilmann Weber, Hyun Uk Kim
# The Novo Nordisk Foundation Center for Biosustainability
# Technical University of Denmark
# Section: Metabolic Engineering for Natural Compounds / New Bioactive Compounds
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.


"""
Predict EC numbers by querying KEGG genes
"""

import logging
import sys
import urllib2
import re
import pickle
from urllib2 import URLError
from antismash import utils


name = "kegg"
short_description = name.capitalize()
priority = 10000

# URLs for KEGG Rest API
QueryGI = "http://rest.kegg.jp/conv/genes/ncbi-gi:%s"
QueryGeneID = "http://rest.kegg.jp/conv/genes/ncbi-geneid:%s"
QueryKEGG = "http://rest.kegg.jp/get/%s"
ChunkSizeMapping = 200
# Note: For the QueryKEGG API call KEGG returns A MAXIMUM of 10 hits, so don't use chunks>10!
ChunkSizeDataRetrieval = 10

def _getKEGG_speciesLocusTag(CDSFeatureDict):
    "Query KEGG rest API to get species:LocusTag"
     
    #TODO: Check with Hyun Uk whether we reallly need this function. Couldn't we just make one query to 
    #get the species ID and then use the NCBI Locus tags to assemble the Kegg Species:LocusTag...
     
     
    # speciesLocusTagDict[ID] = Kegg-LocusTag
    speciesLocusTagDict = {}
    keylist = CDSFeatureDict.keys()
     
    # Get mapping data from KEGG in chunks of $ChunkSize
    for chunkIterator in range (0, len(keylist), ChunkSizeMapping):
        queryGIList = []
        queryGeneIDList = []
         
        logging.debug("ECpredictor/kegg: doing mapping for CDS %s - %s of %s" % (chunkIterator, chunkIterator + ChunkSizeMapping, len(keylist)))
        #GItoIDDict[KEGG ncbi-gi or ncbi-geneid] = ID
        GItoIDDict = {}
        for i in range (0,ChunkSizeMapping):
             
            keyindex = i +  chunkIterator
            #print keyindex
             
            # Break loop when reaching last list index
            if keyindex >= len(keylist):
                break
                 
            Feature = CDSFeatureDict[keylist[keyindex]]
             
            # Get gi or GeneID from db_xref tag of CDS feature
            if Feature.qualifiers.has_key('db_xref'):
                db_xrefs = Feature.qualifiers['db_xref']
                 
                for db_xref in db_xrefs:
                    KEGGQuery=""
                    if "GI:" in db_xref:
                        queryGIList.append(db_xref.split(':')[1])
                        GItoIDDict["ncbi-gi:"+db_xref.split(':')[1]] = keylist[keyindex]
                        break
                    elif "GeneID:" in db_xref:
                        queryGeneIDList.append(db_xref.split(':')[1])
                        GItoIDDict["ncbi-geneid:"+db_xref.split(':')[1]] = keylist[keyindex]
                         
        KEGGQuery = ""
        # if list of queryGIList is not empty, do KEGG REST API call for ncbi-gi
        if len(queryGIList) > 0:
            try:
                #print QueryGI % "+ncbi-gi:".join(queryGIList)
                KEGGQuery += urllib2.urlopen(QueryGI % "+ncbi-gi:".join(queryGIList)).read()
                print
            except URLError, e:
                logging.warn("Could not retrieve GI - KEGG results for species:LocusTag query for entries %s - %s \n \
                              Error message %s" % (chunkIterator, chunkIterator + i, e.reason))
         
        # if list of queryGeneIDList is not empty, do KEGG REST API call for ncbi-geneid
        if len(queryGeneIDList)>0:
            try:
                KEGGQuery += urllib2.urlopen(QueryGeneID % "+ncbi-geneid:".join(queryGeneIDList)).read()
            except URLError, e:
                logging.warn("Could not retrieve GeneID - KEGG result for species:LocusTag query for features %s - %s \n \
                              Error message %s" % (chunkIterator, chunkIterator + i, e.reason))
                 
        # Process results from KEGG API
        for line in KEGGQuery.strip().split("\n"):
            [GI, KEGGspeciesLocusTag] = line.strip().split()
            #print "GI / KEGGspeciesLocusTag: %s / %s\n" % (GI, KEGGspeciesLocusTag)
             
            antiSMASHID = GItoIDDict[GI]
            speciesLocusTagDict[antiSMASHID] = KEGGspeciesLocusTag
             
            # logging.debug("Found %s - %s" % (antiSMASHID, KEGGspeciesLocusTag))
    
    #DEBUG: Load array from serialized text file generated with pickle.dump(open("filename","w",speciesLocusDict) during debug session
    #f = open("spLocusTagDict.p","r")
    
    #speciesLocusTagDict = pickle.load(f)
    return speciesLocusTagDict
        
def _get_ECNumberDict(KEGGspeciesLocusTagDict):
    "Query KEGG REST API to retrieve EC numbers"
    
        # speciesLocusTagDict[ID] = Kegg-LocusTag
    KeggID_ECDict = {}
    antiSMASHID_ECDict = {}
    keylist = KEGGspeciesLocusTagDict.keys()
    
    KEGGquery = ""
    # Get mapping data from KEGG in chunks of $ChunkSize
    
    for chunkIterator in range (0, len(keylist), ChunkSizeDataRetrieval):
        logging.debug("ECpredictor/kegg: getting results %s - %s of %s" % (chunkIterator, chunkIterator + ChunkSizeDataRetrieval, len(keylist)))
        KEGGLocusList = []
        for i in range (0,ChunkSizeDataRetrieval):
            keyindex = i +  chunkIterator
            
            # Break loop when reaching last list index
            if keyindex >= len(keylist):
                break
            
            KEGGLocusList.append(KEGGspeciesLocusTagDict[keylist[keyindex]])
        logging.debug("+".join(KEGGLocusList))
        KEGGQueryResult = urllib2.urlopen(QueryKEGG % "+".join(KEGGLocusList)).read()
        logging.debug("found %s hits" % len(re.findall("///", KEGGQueryResult)))
        KEGGquery += KEGGQueryResult
        
    
    split_text = KEGGquery.strip().split('\n')
    
    KEGGOrgID = ""
    KEGGLocusTag =""
    ECNumbers = []
    for line in split_text:
        sptlist = re.split('\s+', line)
        
        if sptlist[0] == 'ENTRY':
            KEGGLocusTag = sptlist[1]
        if sptlist[0] == 'ORGANISM':
            KEGGOrgID = sptlist[1]
        if sptlist[0] == 'ORTHOLOGY':
            ECNumbers = re.findall('\d+\.\d+\.\d+\.\d+', line)
        if sptlist[0] == '///':
            # As this always is the last line of the entry, we can now assemble the data
            if len(ECNumbers)>0:
                if ((KEGGLocusTag == "") or (KEGGOrgID == "")):
                    logging.warn('Error parsing KEGG result; could not parse Locus tag or Organism')
                KeggID_ECDict[KEGGOrgID+":"+KEGGLocusTag] = ECNumbers
            KEGGLocusTag = ""
            KEGGOrgID = ""
            ECNumbers = []
    
    for antiSMASHID, KEGGspeciesLocusTag in KEGGspeciesLocusTagDict.items():
        if KeggID_ECDict.has_key(KEGGspeciesLocusTag):
            antiSMASHID_ECDict[antiSMASHID] = KeggID_ECDict[KEGGspeciesLocusTag]
    return antiSMASHID_ECDict

def getECs(seq_record, options):
    
    if not name in options.ecpred:
        logging.debug("ECprediction %s not selected, returning..." % name)
        return
    
    CDSFeatureDict = utils.get_feature_dict(seq_record)
    logging.debug("Predicting EC numbers using KEGG online queries")
    KEGGspeciesLocusTagDict = _getKEGG_speciesLocusTag(CDSFeatureDict)
    ECDict = _get_ECNumberDict(KEGGspeciesLocusTagDict)
        
    notes = []
    # logging.debug("Found %s EC predictions" % len(ECDict.keys()))
    for key in ECDict.keys():
        Feature = CDSFeatureDict[key]
        if Feature.qualifiers.has_key('note'):
            notes = Feature.qualifiers['note']
        
        if len(ECDict[key]) > 0:
            logging.debug("Found EC numbers: %s" % ", ".join(ECDict[key]))
            notes.append('EC number prediction based on KEGG query: %s' % ECDict[key])
            Feature.qualifiers['note'] = notes
            if Feature.qualifiers.has_key('EC_number'):
                logging.warn('ECpredictor[kegg]: Overwriting existing EC annotation: %s  with %s' % \
                             (", ".join(Feature.qualifiers['EC_number']), ", ".join(ECDict[key])))
            
            Feature.qualifiers['EC_number'] = ECDict[key]
        else:
            logging.warn('ECpredictor[KEGG]: Could not find EC number for %s' % utils.get_gene_id(Feature))
        
    
    