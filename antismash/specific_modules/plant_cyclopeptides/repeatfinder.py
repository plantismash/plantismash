import os
import Bio
#import seqparse as sp
#import IO.utils as ut
from . import scorer as sc
#import IO.output_repeatfinder as output
from . import Classifier

def get_qualifier_translation(feat):

    return feat.qualifiers["translation"][0]

def print_gbk_file(gbk_file):

    print(repr(gbk_file))
    print(gbk_file.seq)
    print(len(gbk_file.seq))
    print(gbk_file.features[0])
    print(ut.get_aa_translation(gbk_file,gbk_file.features[0]))
    print(len(ut.get_aa_translation(gbk_file,gbk_file.features[0])))

    print_features(gbk_file)


def check_has_burp_repeat(feat):
    """Returns True if the feature has both a repeat and a BURP domain annotation."""
    if feat.qualifiers.get("has_repeat", False):
        domain_info = feat.qualifiers.get("domain_record", [""])[0]  # domain_record is usually a list
        if "plants/BURP" in domain_info:
            return True
    return False


def run_fbk(seq_record):
    repeat_regions = []
    for feat in seq_record.features:
        if feat.type == "CDS":
            scan_seqfeature_translation(seq_record,feat)
            Classifier.check_known_classes(feat)
            find_top_kmer(feat)
            if len(feat.qualifiers["top_kmer_hits"])>0:
                build_table(feat)

                expand_table(feat)
                seq_record.annotations["repetitive_sequence_number"] = 0
                if len(feat.qualifiers["table"][0])>3:#TODO replace this with assesment module as bool check
                    seq_record.annotations["repetitive_sequence_number"] += 1
                    if "gene" in feat.qualifiers:
                        repeat_regions.append(feat.qualifiers["gene"])
                    else:
                        repeat_regions.append("unknown region")
                    make_pattern(feat)
                    feat.qualifiers["has_repeat"] = True
                    seq_record.annotations["repeat_regions"] = repeat_regions
                    if check_has_burp_repeat(feat):
                        feat.qualifiers["has_burp_repeat"] = True
                    should_delete = False
                    for t in feat.qualifiers["table"]:
                        if len(set(list(t))) <=2:
                            should_delete = True
                    if should_delete:
                        if should_delete:
                            feat.qualifiers["table"] = None
                            if "has_repeat" in feat.qualifiers:
                                feat.qualifiers["has_repeat"] = False
                            if "pattern" in feat.qualifiers:
                                del feat.qualifiers["pattern"]

                    #output.add_html_output(seq_record)



def assess_cluster(seq_record):
    """""this will run the assessment of the cluster, looking at coverage, spacing, classification,
    identity inter repeat regions and complexity/compressibility"""
    
    print("placeholder")
    #TODO: write the actual function

def print_features(seq_record):

    for feat in seq_record.features:

        print(feat)


def scan_seqfeature_translation(seq_record,feature, k = 3):
    #TODO make scanning only CDS's and option
    kmer_dict = {}
    seq = None
    #TODO rewrite this to just give a translation when none exists
    if feature.qualifiers["translation"]:
       seq = feature.qualifiers["translation"][0]#feature.qualifiers[key] gives a list of a single string, this gives the translation
    else:
        seq = ut.get_aa_translation(seq_record,feature)

    for i in range(len(seq)-k):

        word = str(seq[i:i+k])
        if word in kmer_dict:

            kmer_dict[word].append(i)

        else:

            kmer_dict[word] = [i]

    feature.qualifiers["kmer_dict"] = kmer_dict

def find_top_kmer(feature):

    kmers = feature.qualifiers["kmer_dict"]
    top_kmer_word = ""
    top_kmer_hits = []
    n_hits = 0
    for key in kmers:
        if len(kmers[key]) > n_hits:

            top_kmer_word = key
            top_kmer_hits = kmers[key]
            n_hits = len(kmers[key])

    feature.qualifiers["top_kmer_word"] = top_kmer_word
    feature.qualifiers["top_kmer_hits"] = top_kmer_hits

def build_table(feature):
    table = [feature.qualifiers["top_kmer_word"]]*len(feature.qualifiers["top_kmer_hits"])

    feature.qualifiers["table"] = table

def expand_table(feature,cutoff = 2):

    seq = get_qualifier_translation(feature)
    table = feature.qualifiers["table"]
    hits = feature.qualifiers["top_kmer_hits"]

    cont = True#continuator boolean

    while cont:

        cont = expand_forward(seq,hits,table,cutoff)
        #print table
    #remove duplicate repeats
        
    #feature.qualifiers["table"] = table


def expand_forward(seq,hits,pattern_table, cutoff = 2.5):

    next_residue_list = []
    max_len = 15

    #print "seq length: %s, hit index: %s, pattern table length: %s --statement from find_by_kmer.expand_forward"%(len(seq),hits[0],len(pattern_table[0]))

    for i,hit in enumerate(hits):
        #print pattern_table
        if hit+len(pattern_table[i]) < len(seq) and len(pattern_table[i]) < max_len:
            next_residue_list.append(seq[hit+len(pattern_table[i])])
        else:
            return False
    if check_cutoff(next_residue_list, cutoff):

        for i,char in enumerate(next_residue_list):

            pattern_table[i] += char


        return True

    else:
        return False

def check_cutoff(nrl, cutoff= 0.4):
    #old code
    #max_res = most_common(nrl)
    #n_mcr = nrl.count(max_res)
    #return n_mcr/len(nrl) > cutoff

    max_res = most_common(nrl)
    order, matrix = sc.blosum62()
    score_list = []
    for res in nrl:
        current_score=sc.score(max_res, res, order, matrix)
        score_list.append(current_score)
    if len(score_list)> 0:
        avg_score = sum(score_list)/len(score_list)
    else:
        return False
    #print avg_score
    if avg_score > cutoff:

        return True
    else:

        return False

def most_common(lst):
    try:
        return max(set(lst), key=lst.count)
    except:
        "print error occurred in most common function in repeatfinder.py"
def make_pattern(feat):
    #this function can stay the same. no changes needed
    #function creates a regex from an expanded kmer table
    pattern_string = ""
    current_residue_list = []

    pattern_table = feat.qualifiers["table"]

    for i in range(len(pattern_table[0])):
        for pattern in pattern_table:
            current_residue_list.append(pattern[i])

        unique_residues = "".join(set(current_residue_list))

        if len(unique_residues) > 1:

            pattern_string += "[" + unique_residues + "]"


        else:
            pattern_string += unique_residues
        current_residue_list = []

    feat.qualifiers["pattern"] = pattern_string

def mass_run(in_path,out_path,output_type = "txt",teiresias = False):

    output_folder_name = "/results"
    filepath_list = sp.make_filepath_list(in_path)
    outpath_list = sp.generate_outpath_list(filepath_list,out_path,output_folder_name)

    run_fbk_per_cluster(filepath_list,outpath_list,output_type,teiresias)

def run_fbk_per_cluster(filepath_list,outpath_list,output_type = "txt",teiresias= False):


    for i,fp in enumerate(filepath_list):
        cluster_nr = fp.split("0")[-1].split(".")[0]


        if not os.path.isdir("/".join(outpath_list[i].split("/")[:-1])):
            os.mkdir("/".join(outpath_list[i].split("/")[:-1]))
        if not os.path.isdir(outpath_list[i]):
            os.mkdir(outpath_list[i])
        #TODO make the below loop more efficient
        repeat_features = []
        seq_record = sp.open_gbk(fp)
        run_fbk(seq_record)
        #summary_file = open("/".join(outpath_list[i].split("/")[:-1])+ "/output_summary.txt","a")
        if output_type == "html":
            out_file = open(outpath_list[i]+ "/fbk_output.html","w+")
            out_file.write("Cluster found in organism %s \n"%(fp.split("/")[-1]))
            if "html_output" in seq_record.annotations:
                out_file.write("<br>"+seq_record.annotations["html_output"])
            for feature in seq_record.features:
                if "html_output" in feature.qualifiers:
                    out_file.write(feature.qualifiers["html_output"])
        if output_type == "genbank":

            sp.write_seq_record_to_genbank(seq_record, outpath_list[i])
def main():

    mass_run("C:/users/ezmod/data/results","C:/users/ezmod/data/algorithm_output_v6",output_type="genbank")

if __name__ == "__main__":

    main()

