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
    #remove duplicates
    #SeqEl_list = remove_duplicates(SeqEl_list)
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
    to_remove = [] 
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
        elif el1.s == el2.s and el1.e == el2.e:
            to_remove.append(el1)
    for el in to_remove:
        SeqEl_list.remove(el)

    SeqEl_list = se_sort(SeqEl_list+new_el_list)
    l = SeqEl_list[-1].e
    SeqEl_list.append(SeqEl(l,len(seq),"n",seq[l:len(seq)]))
    return SeqEl_list

def remove_duplicates(SeqEl_list):
    to_remove = [] 
    for i,el in enumerate(SeqEl_list[:-1]):
        el1 = el
        el2 = SeqEl_list[i+1]
        if el1.s == el2.s and el1.e == el2.e and el1.t == el2.t:
            to_remove.append(el1)
    for el in to_remove:
        SeqEl_list.remove(el)
    return SeqEl_list


def apply_styles(Seqel_list):
    for el in Seqel_list:
        if el.k == "a":
            el.t = html.bstyle+html.color.red+ html.stclose+ el.t + html.bclose
        if el.k == "b":
            el.t = html.bstyle+html.color.blue+ html.stclose+ el.t + html.bclose
        if el.k == "c":
            el.t = html.bstyle+html.color.blue+ html.stclose+ el.t + html.bclose


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
            
            if len(result.evidence[ek]) >= 1 and len(ek) > 0:
                known_motif_found = True
                if summary == "":
                    summary += hr
                if not summary_text_printed:
                    if result.cds_id:
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
                known_motifs = populate_known_motifs(list(result.evidence.keys()))
                html_seq = generate_html_seq(result.sequence,result.pattern,known_motifs)
            else:
                #repeat but no known motifs
                html_seq = generate_html_seq(result.sequence,result.pattern)
        elif known_motif_found:
           # no repeat, known motif
            known_motifs = populate_known_motifs(list(result.evidence.keys()))
            html_seq = generate_html_seq(result.sequence,known_motifs)

            summary += "Sequence: <br>"  # <br> in html
        summary += format_linebreaks(html_seq, 64)
        return summary

    return ""

def populate_known_motifs(keylist):
    outlist = []
    for k in keylist:
        outlist = k.split(": ")[-1]
    return outlist


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
