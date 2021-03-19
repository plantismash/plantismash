import scorer as sc
import Classifier
from antismash.utils import get_aa_translation


def get_qualifier_translation(feat):

    return feat.qualifiers["translation"][0]


def run_fbk(seq_record):
    repeat_regions = []
    for feat in seq_record.features:
        if feat.type == "CDS":
            scan_seqfeature_translation(seq_record, feat)
            Classifier.check_known_classes(feat)
            find_top_kmer(feat)
            if len(feat.qualifiers["top_kmer_hits"]) > 0:
                build_table(feat)

                expand_table(feat)
                seq_record.annotations["repetitive_sequence_number"] = 0
                if len(feat.qualifiers["table"][0]) > 3:
                    seq_record.annotations["repetitive_sequence_number"] += 1
                    if "gene" in feat.qualifiers:
                        repeat_regions.append(feat.qualifiers["gene"])
                    else:
                        repeat_regions.append("unknown region")
                    make_pattern(feat)
                    feat.qualifiers["has_repeat"] = True
                    seq_record.annotations["repeat_regions"] = repeat_regions
                    should_delete = False
                    for t in feat.qualifiers["table"]:
                        if len(set(list(t))) <= 2:
                            should_delete = True
                    if should_delete:
                        feat.qualifiers["table"] = None
                        feat.qualifiers["has_repeat"] = False
                        del(feat.qualifiers["pattern"]) 


def scan_seqfeature_translation(seq_record, feature, k=3):
    kmer_dict = {}
    seq = None
    if feature.qualifiers["translation"]:
        seq = feature.qualifiers["translation"][0]  #feature.qualifiers[key] gives a list of a single string, this gives the translation
    else:
        seq = get_aa_translation(seq_record, feature)

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


def expand_table(feature, cutoff=2):

    seq = get_qualifier_translation(feature)
    table = feature.qualifiers["table"]
    hits = feature.qualifiers["top_kmer_hits"]

    cont = True  #continuator boolean

    while cont:

        cont = expand_forward(seq, hits, table, cutoff)


def expand_forward(seq, hits, pattern_table, cutoff=2.5):

    next_residue_list = []
    max_len = 15

    for i, hit in enumerate(hits):
        if hit+len(pattern_table[i]) < len(seq) and len(pattern_table[i]) < max_len:
            next_residue_list.append(seq[hit+len(pattern_table[i])])
        else:
            return False
    if check_cutoff(next_residue_list, cutoff):

        for i, char in enumerate(next_residue_list):

            pattern_table[i] += char

        return True

    else:
        return False


def check_cutoff(nrl, cutoff=0.4):

    max_res = most_common(nrl)
    order, matrix = sc.blosum62()
    score_list = []
    for res in nrl:
        current_score = sc.score(max_res, res, order, matrix)
        score_list.append(current_score)
    if len(score_list) > 0:
        avg_score = sum(score_list)/len(score_list)
    else:
        return False
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
