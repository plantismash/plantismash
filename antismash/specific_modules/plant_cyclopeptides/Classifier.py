from antismash.utils import get_full_path
import re


def check_known_classes(feat):
    #known classes are checked for by comparing the repeat containing sequences to these regex snippets
    #known classes are stored as regex patterns in a fasta like txt file
    known_classes = get_known_classes_from_file()
    qualdict = {}
    for key in known_classes:
        pattern = known_classes[key]
        p = re.compile(pattern)
        pattern_instances = []
        matches = p.finditer(feat.qualifiers["translation"][0])
        mlen = 0 
        for m in matches:
            mlen += 1 
            pattern_instances.append(m.start())
        
        #print "pattern occurrences %s and instances %s, key is %s" %(str(pattern_occurrences),str(pattern_instances),key) 
        if mlen is not 0: 
            qualdict[key+": "+pattern] = pattern_instances
        
    if len(qualdict) is not 0: 
        feat.qualifiers["ripp_evidence"] = qualdict
    else:
        feat.qualifiers["ripp_evidence"] = {}


def get_known_classes_from_file():
    filepath = get_full_path(__file__, "known_motifs.txt")
    known_classes = {}
    infile = open(filepath, "r")

    current_key = ""

    for li in infile:
        line = re.sub("\\n", "", li)
        line = re.sub("\\r", "", line)
        if line.startswith(">"):
            current_key = line[1:]

        else:
            known_classes[current_key] = line

    infile.close()

    return known_classes
