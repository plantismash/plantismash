import re


class Effects:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


class SeqEl:
    s = None  #start
    e = None  #end
    k = ""  #kind
    t = ""  #text

    def __init__(self, start, end, kind, text):
        self.s = start
        self.e = end
        self.k = kind
        self.t = text


class html:
    bold = "<b>"
    bstyle = "<b style=\""
    bclose = "</b>"
    underscore = "<u>"
    uclose = "</u>"
    style = "style=\""
    stclose = "\">"

    class color:
        red = "color:#ff0000"
        blue = "color:#0000ff"
        green = "color:#3cb371"
        pink = "color:#ee82ee"
        yellow = "color:#ffa500"
        purple = "color:#6a5acd"

def generate_html_seq(seq,pattern_a,pattern_b = None):
    #instantiate matches
    SeqEl_list = matches_to_SeqEl(pattern_a,seq,"a")
    if pattern_b:
        SeqEl_list += matches_to_SeqEl(pattern_b,seq,"b")
    #sort list
    SeqEl_list = se_sort(SeqEl_list)
    #fill list
    SeqEl_list = fill_Seqel_list(SeqEl_list,seq)
    #style code
    apply_styles(SeqEl_list)

    #combine elements into text
    total_text = ""
    for el in SeqEl_list:
        total_text += el.t
    return total_text

def matches_to_SeqEl(pattern,seq,kind):
    SeqEl_list = []
    if isinstance(pattern,list):
        print(pattern)
        for p in pattern:
            SeqEl_list += make_SeqEl_list(p,seq,kind)

    else:
        SeqEl_list = make_SeqEl_list(pattern,seq,kind)

    return SeqEl_list

def make_SeqEl_list(pattern,seq,kind):
    SeqEl_list = []
    mts = re.finditer(pattern,seq)

    for m in mts:
        SeqEl_list.append(SeqEl(m.start(), m.end(), kind, m.group(0)))

    return SeqEl_list

def se_sort(SeqEl_list):
    sorted_list = []
    low_el = None
    l = len(SeqEl_list)
    for i in range(l):
        num = 10000
        for el in SeqEl_list:
            if el.s < num:
                num = el.s
                low_el = el

        sorted_list.append(low_el)
        SeqEl_list.remove(low_el)

    return sorted_list

def fill_Seqel_list(SeqEl_list, seq):
    f = SeqEl_list[0].s
    start_el = SeqEl(0,f,"n",seq[0:f])
    new_el_list = [start_el]
    for i,el in enumerate(SeqEl_list[:-1]):
        el1 = el
        el2 = SeqEl_list[i+1]
        if el2.e > el1.e > el2.s and el1.k != el2.k:
            cs = el2.s
            ce = el1.e
            ol = ce - cs
            ct = el2.t[:ol]
            el3 = SeqEl(cs,ce,"c",ct)
            new_el_list.append(el3)
            el1.e = cs #TODO these might have to be incremented
            orig_el2_end = el2.s
            el2.s = ce
            el1.t = el1.t[:el1.e-el1.s]
            el2.t = el2.t[ce- orig_el2_end:]
        elif el1.e > el2.e and el1.k != el2.k:
            el2.k = "c"
            cs = el2.e
            ce = el1.e
            ck = el1.k
            ct = seq[cs:ce]
            new_el_list.append(SeqEl(cs,ce,ck,ct))
            el1.e = el2.s
            el1.t = seq[el1.s:el1.e]

        elif el1.e < el2.s:
            ns = el1.e
            ne = el2.s
            nk = "n"
            nt = seq[ns:ne]
            new_el_list.append(SeqEl(ns,ne,nk,nt))

    SeqEl_list = se_sort(SeqEl_list+new_el_list)
    l = SeqEl_list[-1].e
    SeqEl_list.append(SeqEl(l,len(seq),"n",seq[l:len(seq)]))
    return SeqEl_list

def apply_styles(Seqel_list):
    for el in Seqel_list:
        if el.k == "a":
            el.t = html.bstyle+html.color.red+ html.stclose+ el.t + html.bclose
        if el.k == "b":
            el.t = html.bstyle+html.color.blue+ html.stclose+ el.t + html.bclose
        if el.k == "c":
            el.t = html.bstyle+html.color.blue+ html.stclose+ el.t + html.bclose

def highlight_seq(sequence, instances, color = "#FF0000"):
    """"This function highlights the repeats in a certain sequence"""
    formatted_string = sequence
    for instance in set(instances):

        b_start = "<b style=\"color:" + color +"\">"
        b_end = "</b>"
        replacement = b_start+instance+b_end
        formatted_string = re.sub(instance,replacement,formatted_string)

    return formatted_string

def underscore_multi(sequence, pattern_dict):

    for k in pattern_dict:
       sequence = underscore_seq(sequence,pattern_dict[k])

    return sequence

def underscore_seq(sequence,instances):
    u_start = "<u>"
    u_end = "</b>"
    formatted_string = ""
    for instance in set(instances):
        replacement = u_start + instance + u_end
        formatted_string = re.sub(instance, replacement, sequence)

    return formatted_string

def format_string(sequence, instances,effects_list, stringlength):

    formatted_string = sequence
    str_effects = ""
    for e in effects_list:
        str_effects += e

    for instance in instances:
        replacement = str_effects + instance + Effects.END
        formatted_string = re.sub(instance,replacement,formatted_string)
    
    return formatted_string

def write_result_summary(result,instance_line_nr = 5):
    br = "<br>"#CHANGED FOR PLAINTEXT VERSION, <br> is html
    hr = "<hr>"#CHANGED FOR PLAINTEXT VERSION, <hr> is html
    summary ="" 
    seqprint = False
    repeat_found = False
    known_motif_found = False
    if len(result.instances) > 1:
        repeat_found = True
        summary += hr
        if result.cds_id:
            summary += "Repeat found in %s <br>"%(result.cds_id)
        summary += "Repeat occurs %s times in a sequence of %s amino acids"%(len(result.instances),len(result.sequence))+br
        summary += "Location between %s and %s %s" % (result.position[0], result.position[1], br)
        coverage = round((float(len(result.instances)*len(result.instances[0]))/len(result.sequence))*100,2)
        summary += "Coverage of %s %s"%(coverage,"%")+br
        summary += "Instances:<br> %s%s"%(format_instances(result.instances,instance_line_nr),br)
        summary += "pattern: %s%s"%(result.pattern,br)
        seqprint = True

    if result.evidence:
        summary_text_printed = False
        for ek in result.evidence:
            
            if len(result.evidence[ek]) > 1:
                known_motif_found = True
                if summary == "":
                    summary += hr
                if not summary_text_printed:
                    if repeat_found:
                        summary += "The following known motifs were found: <br>"
                    elif result.cds_id:
                        summary += "The following known motifs were found in CDS %s <br>"%(result.cds_id)
                        summary += "Location between %s and %s %s" % (result.position[0], result.position[1], br)
                    else:
                        summary += "The following known motifs were found between locations %s and %s <br>" % (result.position[0], result.position[1], br)

                    summary_text_printed = True

                summary += "%s was found %s times in this sequence<br>"%(ek,len(result.evidence[ek]))
                seqprint = True
    if seqprint:
        html_seq = ""
        if repeat_found:
            #repeat
            if known_motif_found:
                #repeat and known motifs
                known_motifs = list(result.evidence.keys())
                print(known_motifs)
                print((type(known_motifs)))
                html_seq = generate_html_seq(result.sequence,result.pattern,known_motifs)
            else:
                #repeat but no known motifs
                html_seq = generate_html_seq(result.sequence,result.pattern)
        elif known_motif_found:
           # no repeat, known motif
            known_motifs = list(result.evidence.keys())
            html_seq = generate_html_seq(result.sequence,known_motifs)

            summary += "Sequence: <br>"  # <br> in html
        summary += format_linebreaks(html_seq, 64)
        return summary

    return ""

def write_feature_summary(feature,instance_line_nr = 5):
    """"This function takes all the ancillary data from a repeat feature and stylizes it with HTML
    repeat_feature is a repeat feature object found in the misc package
    returns a string readable by a HTML parser"""
    br = "<br>"
    hr = "<hr>"
    summary = ""
    summary += hr
    summary += "Repeat found between %s and %s %s"%(feature.location.start,feature.location.end,br)
    summary += "Repeat was found %s times in a sequence of %s amino acids"%(len(feature.qualifiers["table"]),len(feature.qualifiers["translation"]))+br
    coverage = round((float(len(feature.qualifiers["table"])*len(feature.qualifiers["table"][0]))/len(feature.qualifiers["translation"][0]))*100,2)
    summary += "Coverage of %s %s"%(coverage,"%")+br
    summary += "Instances:<br> %s%s"%(format_instances(feature.qualifiers["table"],instance_line_nr),br)
    summary += "pattern: %s%s"%(feature.qualifiers["pattern"],br)

    return summary
def write_evidence_summary(feat):
    br = "<br>"
    hr = "<hr>"
    summary = "the following evidence was found for RiPPs " + br +hr
    known_ripp_found = False
    for key in feat.qualifiers["ripp_evidence"]:
        current_evidence_len = len(feat.qualifiers["ripp_evidence"][key][0])
        if current_evidence_len >0:
            substrings = key.split("|")
            known_ripp_str = substrings[0]
            pattern = substrings[1]

            summary += "the pattern %s, associated with %s was found %s times %s"%(pattern,known_ripp_str,current_evidence_len,br)
            known_ripp_found = True
    if known_ripp_found == False:
        summary += "no known ripp patterns found in this sequence" + br

    summary += "Putative glutamine core peptides found:" + br
    glut_struct_instances = feat.qualifiers["glutamine_aromat_structures"][0]

    summary += format_instances(glut_struct_instances,5)
    summary += br
    #summary += feat.qualifiers["html_seq_glutamine"]
    return summary


def format_linebreaks(line,max_length):
    break_index_list = []
    counter = 0
    should_count = True
    for i,char in enumerate(line):

        if char == "<":
            should_count = False

        if should_count ==True:
            counter += 1

        if char ==">":
            should_count = True

        if counter == max_length:

            break_index_list.append(i)

            counter = 0

    correction = 0
    for index in break_index_list:
        index = index + correction
        line = line[:index] + "<br>" + line[index:]#CHANGE BACK TO HTML MAYBE

        correction += 4
    return line

def create_output(feature):
    """"Wrapper function to create stylized html output per repeat feature"""
    highlight_seq(feature.qualifiers["translation"][0],instances=feature.qualifiers["table"],output_key="html_seq")
    highlight_seq(feature.qualifiers["translation"][0],instances=feature.qualifiers["glutamine_aromat_structures"][0],output_key="html_seq_glutamine")
    output_line = write_feature_summary(feature)
    output_line += "Sequence: <br>"
    output_line += format_linebreaks(feature.qualifiers["html_seq"],64)
    output_line += write_evidence_summary(feature)
    output_line += "<br>glutamine sequences highlighted <br>"
    output_line += format_linebreaks(feature.qualifiers["html_seq_glutamine"],64)
    feature.qualifiers["html_output"] = output_line

def format_instances(instances,nr_on_line):
    """formats the amount of instances you see on one line.
    takes instances as a list and nr_on_ line is the amount of instances you can see on one line
    return a formatted html string of the instances"""

    counter = 0
    output_string = ""
    for instance in instances:
        output_string += instance + " | "
        counter += 1

        if counter == nr_on_line:
            output_string =output_string[:-2] + "<br>"
            counter = 0

    return output_string

def add_html_output(seqrecord):
    seqrecord.annotations["html_output"] = "Cluster %s from %s with %s potential repeat regions"%(seqrecord.id,seqrecord.description, len(seqrecord.annotations["repeat_regions"]))

    for feat in seqrecord.features:
        if "has_repeat" in feat.qualifiers:
            print((feat.qualifiers["table"]))
            print((feat.qualifiers["ripp_evidence"]))
            if len(feat.qualifiers["table"]) > 0 or len(feat.qualifiers["ripp_evidence"]) >0: 
                create_output(feat)
                print("condition reached")
                exit(0)
