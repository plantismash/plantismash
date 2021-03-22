#import IO.seqparse as sp
import re
import sys
print(sys.path)
from antismash import utils 
#The filepath for the known motifs file
#filepath = "./plantismash/antismash/specific_modules/plant_cyclopeptides/known_motifs.txt"
filepath = utils.get_full_path(__file__,"known_motifs.txt")

def get_properties_table():

    properties_dict = {'A': {"name": "Alanine","three_letter":"Ala","hydro":"hydrophobic","type": "Aliphatic"},
    'I':{"name": "Isoleucine","three_letter":"Ile","hydro":"hydrophobic","type": "Aliphatic"},
    'L':{"name": "Leucine","three_letter":"Leu","hydro":"hydrophobic","type": "Aliphatic"},
    'M':{"name": "Methionine","three_letter":"Met","hydro":"hydrophobic","type": "Aliphatic"},
    'V':{"name": "Valine","three_letter":"Val","hydro":"hydrophobic","type": "Aliphatic"},
    'F':{"name": "Phenylalanine","three_letter":"Phe","hydro":"hydrophobic","type": "Aromatic"},
    'W':{"name": "Tryptophan","three_letter":"Trp","hydro":"hydrophobic","type": "Aromatic"},
    'Y':{"name": "Tyrosine","three_letter":"Tyr","hydro":"hydrophobic","type": "Aromatic"},
    'N':{"name": "Asparagine","three_letter":"Asp","hydro":"hydrophilic","type": "polar-neutral"},
    'C':{"name": "Cysteine","three_letter":"Cys","hydro":"hydrophilic","type": "polar-neutral"},
    'Q':{"name": "Glutamine","three_letter":"Gln","hydro":"hydrophilic","type": "polar-neutral"},
    'S':{"name": "Serine","three_letter":"Ser","hydro":"hydrophilic","type": "polar-neutral"},
    'T':{"name": "Threonine","three_letter":"Thr","hydro":"hydrophilic","type": "polar-neutral"},
    'D':{"name": "Aspartate","three_letter":"Asp","hydro":"hydrophilic","type": "positively charged"},
    'E':{"name": "Glutamate","three_letter":"Glu","hydro":"hydrophilic","type": "positively charged"},
    'R':{"name": "Arginine","three_letter":"Arg","hydro":"hydrophilic","type": "negatively charged"},
    'H':{"name": "Histidine","three_letter":"His","hydro":"hydrophilic","type": "negatively charged"},
    'K':{"name": "Lysine","three_letter":"Lys","hydro":"hydrophilic","type": "negatively charged"},
    'G':{"name": "Glycine","three_letter":"Gly","hydro":"hydrophilic","type": "unique"},
    'P':{"name": "Proline","three_letter":"Pro","hydro":"hydrophilic","type": "unique"},
    'X':{"name": "Unknown", "three_letter": "Unk","type": "unknown"}}

    return properties_dict

def get_aa_properties(aa):

    aa_dict = get_properties_table()

    return aa_dict[aa.upper()]

def get_aa_type(aa):

    aa_dict = get_properties_table()
    if aa.upper() in aa_dict:
        return aa_dict[aa.upper()]["type"]
    else:
        return "amino acid unknown"

def get_aa_property(aa, property_key):

    aa_dict = get_properties_table()

    if aa.upper() in aa_dict and property_key in aa_dict[aa.upper()]:

        return aa_dict[aa.upper()][property_key]
    else:
        return "property not found"

def get_glutamine_seqs(feat):

    seq = feat.qualifiers["translation"][0]
    glutamine_seqs = []
    for i in range(len(seq)-8):

        if get_aa_property(seq[i],"name") == "Glutamine":

            glutamine_seqs.append(seq[i:i+8])

    return glutamine_seqs

def check_aromatic_aa(seq):

    aromat_counter = 0

    for aa in seq:

        if get_aa_type(aa) == "Aromatic":

            aromat_counter += 1

    return aromat_counter

def check_glutamine_aromat_structure(feat):

    seq = feat.qualifiers["translation"][0]

    glutamine_seqs = get_glutamine_seqs(feat)
    feat.qualifiers["glutamine_aromat_structures"] = [glutamine_seqs]+ check_for_aromatic_structure(glutamine_seqs)


def check_for_aromatic_structure(seqs,cutoff = 0.5):

    aromatic_conc = []
    is_aromat_list = []
    aromats_in_seq_list = []
    for seq in seqs:
        aromats_in_seq_list.append(check_aromatic_aa(seq))
        conc = aromats_in_seq_list[-1]/len(seq)
        aromatic_conc.append(conc)
        if conc > 0.5:
            is_aromat_list.append(True)

        else:
            is_aromat_list.append(False)
    return [aromats_in_seq_list,aromatic_conc,is_aromat_list]

def check_instances_aromat_structure(feat):

    instances = feat.qualifiers["table"]
    feat.qualifiers["instances_aromat_structures"] = check_for_aromatic_structure(instances)

def check_known_classes(feat):
    #known clases are checked for by comparing the repeat containing sequences to these regex snippets
    #known classes are stored as regex patterns in a fasta like txt file
    known_classes = get_known_classes_from_file(filepath)
    #feat.qualifiers["ripp_evidence"] = {}
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
            qualdict[pattern] = pattern_instances
        
    if len(qualdict) is not 0: 
        feat.qualifiers["ripp_evidence"] = qualdict
    else:
        feat.qualifiers["ripp_evidence"] = {}

def get_known_classes_from_file(filepath):
    known_classes = {}
    infile = open(filepath,"r")

    current_key = ""

    for li in infile:
        line = re.sub("\\n", "",li)
        line = re.sub("\\r","",line)
        if line.startswith(">"):
            current_key = line[1:]

        else:
            known_classes[current_key] = line

    infile.close()

    return known_classes

def classify_ripp(feat):

    print 'p'

    #define AA into classes DONE
    #check if it is a ripp SEMI_DONE

    #check if it can be specified into lyciumin or others


def classify_feat(feat):
    """"this function checks the repeats for aromatic AA content, it also searches for a general RiPP core peptide pattern,
    as well as checking the sequence for hits of known patterns, all results are written als feat qualifiers"""
    check_glutamine_aromat_structure(feat)
    check_instances_aromat_structure(feat)
    check_known_classes(feat)
    #check spacing

    #align either interrepeat region or align between kmer

    #maybe add coverage to this module


if __name__ == "__main__":


    known_dict = get_known_classes_from_file("known_motifs.txt")

    for k in known_dict:
        print k
        print known_dict[k]
