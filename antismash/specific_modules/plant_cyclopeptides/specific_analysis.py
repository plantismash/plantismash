# vim: set fileencoding=utf-8 :
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
'''
'''


import logging
#workaround for the import error
from antismash import utils
#import utils
#TODO FIX THIS BEFORE RELEASE
import repeatfinder

def specific_analysis(seq_record, options):
    clusters = utils.get_cluster_features(seq_record)
    for cluster in clusters:
        if 'product' not in cluster.qualifiers or \
           'cyclopeptide' not in cluster.qualifiers['product'][0]:
            continue
        logging.debug("Finding repeats within BURP region (%s)" % (cluster.qualifiers['note'][0]))
        #cluster_record = seq_record[cluster.start:cluster.end]
        find_repeats(seq_record)
         


def fbk_output_to_result(seq_record):
    nr_of_repeat_seqs = 0
    for feat in seq_record.features:
        result = Result()
        if "seq_met" in feat.qualifiers:
            
            feature_type = feat.qualifiers["seq_met"][0]
        if "pattern" in feat.qualifiers:
            if len(feat.qualifiers["table"][0]) > 5 \
                    and len(feat.qualifiers["top_kmer_hits"]) >= 3:
                has_repeat = True
                result.pattern = feat.qualifiers["pattern"]
                result.instances = feat.qualifiers["table"]#top_kmer_hits for locations
                result.sequence = feat.qualifiers["translation"][0]
                result.feature_type = feat.type
                result.position = (feat.location.start,feat.location.end)
                feat.qualifiers['cyclopeptide_analysis'] = [result.encode()]
   
###########  

class Result:
    
    pattern = ""
    instances = 0
    sequence = ""
    feature_type = ""
    position = None
    def __init__(self, qualifier = None):
        if qualifier != None:
            self.decode(qualifier)

    def encode(self):
        # put here your encoding algorithm
        qualifier = ""
        qualifier += self.pattern
        qualifier += "//" + "|".join(str(s) for s in self.instances)
        qualifier += "//" + self.sequence
        qualifier += "//" + self.feature_type
        qualifier += "//" + str(self.position[0]) + "|" + str(self.position[1])
        #qualifier = [self.pattern,str(self.instances),self.sequence,\
        #        self.feature_type ,self.position]
        return qualifier
        
    def decode(self, qualifier):
        self.__init__() # re-initialize
        # put here your decoding algorithm
        print qualifier
        q = qualifier.split("//")
        self.pattern = q[0]
        self.instances = q[1].split("|")
        self.sequence = q[2]
        self.feature_type = q[3]
        self.position = tuple(q[4].split("|"))

def find_repeats(seq_record):
    result = None 

    # put here your repeat-finding procedure
    repeatfinder.run_fbk(seq_record)
    fbk_output_to_result(seq_record) 
    return result
