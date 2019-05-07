# vim: set fileencoding=utf-8 :
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
'''
'''

import logging
from antismash import utils

def specific_analysis(seq_record, options):
    clusters = utils.get_cluster_features(seq_record)
    for cluster in clusters:
        if 'product' not in cluster.qualifiers or \
           'cyclopeptide' not in cluster.qualifiers['product'][0]:
            continue
        logging.debug("Finding repeats within BURP region (%s)" % (cluster.qualifiers['note'][0]))
        result = find_repeats(cluster)
        if result != None:
            cluster.qualifiers['cyclopeptide_analysis'] = [result.encode()]

###########  
class Result:
    dummy_attribute = ""

    def __init__(self, qualifier = None):
        if qualifier != None:
            self.decode(self, qualifier)

    def encode(self):
        qualifier = ""
        # put here your encoding algorithm
        qualifier = self.dummy_attribute
        return qualifier
        
    def decode(self, qualifier):
        self = self.__init__(self) # re-initialize
        # put here your decoding algorithm
        self.dummy_attribute = qualifier

def find_repeats(cluster):
    result = None

    # put here your repeat-finding procedure
    result = Result()
    result.dummy_attribute = "DUMMY ATTRIBUTE"

    return result