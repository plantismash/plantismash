# vim: set fileencoding=utf-8 :
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
'''
'''
from antismash import utils
import repeatfinder

def specific_analysis(seq_record, options):
    clusters = utils.get_cluster_features(seq_record)
    for cluster in clusters:
        if 'product' not in cluster.qualifiers or \
           'cyclopeptide' not in cluster.qualifiers['product'][0]:
            continue
        find_repeats(seq_record)

def fbk_output_to_result(seq_record):
    nr_of_repeat_seqs = 0
    for feat in seq_record.features:
        result = Result()
        if "seq_met" in feat.qualifiers:
            
            feature_type = feat.qualifiers["seq_met"][0]
        has_ripp = False
        if "pattern" in feat.qualifiers:
            if (len(feat.qualifiers["table"][0]) > 5 and len(feat.qualifiers["top_kmer_hits"]) >= 3): 
                has_ripp = True
                result.pattern = feat.qualifiers["pattern"]
                result.instances = feat.qualifiers["table"]#top_kmer_hits for locations
        if "ripp_evidence" in feat.qualifiers and \
            len(feat.qualifiers["ripp_evidence"]) > 0:
            has_ripp = True
                       
            result.evidence = feat.qualifiers["ripp_evidence"]

        if has_ripp:
            result.sequence = feat.qualifiers["translation"][0]
            result.feature_type = feat.type
            result.position = (feat.location.start,feat.location.end)
            if 'gene' in feat.qualifiers:
                result.cds_id = feat.qualifiers['gene'][0]
            elif 'locus_tag' in feat.qualifiers:
                result.cds_id = feat.qualifiers['locus_tag'][0]
            elif 'db_xref' in feat.qualifiers:
                result.CDS_id = feat.qualifiers['db_xref'][0]
            feat.qualifiers['cyclopeptide_analysis'] = [result.encode()]
###########  

class Result:
    
    pattern = ""
    instances = []
    sequence = ""
    feature_type = ""
    position = None
    cds_id = None
    evidence = {}
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
        qualifier += "//" + self.cds_id
        qualifier += "//" + "|".join(self.evidence.keys())
        qualifier += "//" + "|".join([str(x) for x in self.evidence.values()])
        

        return qualifier
        
    def decode(self, qualifier):
        self.__init__() # re-initialize
        # put here your decoding algorithm
        print(qualifier)
        q = qualifier.split("//")
        self.pattern = q[0]
        self.instances = q[1].split("|")
        self.sequence = q[2]
        self.feature_type = q[3]
        self.position = tuple(q[4].split("|"))
        self.cds_id = q[5]#TODO removed length if statement. if it breaks this is probably the reason
        keys = q[6].split("|")
        values_strings = q[7].split("|")
        vallist = []
        for v in values_strings:
            vallist.append(list(v[2:-2].split(",")))
        self.evidence = dict(zip(keys,vallist))
        has_evidence = True#TODO CHANGE BACK OR REMOVE THIS PIECE OF CODE 

        #for vallist in vallist_list:
         #   if len(vallist) > 1:
          #      has_evidence = True
        #if has_evidence:    
        #    self.evidence = dict(zip(evidence_keys,vallist_list)) 
        #else: 
        #    self.evidence = None

def find_repeats(seq_record):
    result = None 

    # put here your repeat-finding procedure
    repeatfinder.run_fbk(seq_record)
    fbk_output_to_result(seq_record) 
    return result
