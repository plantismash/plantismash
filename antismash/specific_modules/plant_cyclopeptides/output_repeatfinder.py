import re
from .specific_analysis import Result
class effects:
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
    
def highlight_seq(sequence, instances, color = "#FF0000"):
    """"This function highlights the repeats in a certain sequence"""

    br = "<br>"
    hr = "<hr>"
    formatted_string = sequence 
    for instance in instances:

        b_start = "<b style=\"color:" + color +"\">"
        b_end = "</b>"
        replacement = b_start+instance+b_end
        formatted_string = re.sub(instance,replacement,formatted_string)

    return formatted_string

def format_string(sequence, instances,effects_list, stringlength):

    formatted_string = sequence
    str_effects = ""
    for e in effects_list:
        str_effects += e

    for instance in instance:
        replacement = str_effects + instance + effects.END 
        formatted_string = re.sub(instance,replacement,formatted_string)
    
    return formatted_string

def write_result_summary(result,instance_line_nr = 5):
    br = "<br>"#CHANGED FOR PLAINTEXT VERSION, <br> is html
    hr = "<hr>"#CHANGED FOR PLAINTEXT VERSION, <hr> is html
    summary ="" 
    seqprint = False
    if len(result.instances) > 1:
        summary += hr
        if result.cds_id:
            summary += "Repeat found in %s <br>"%(result.cds_id)
        summary += "Location between %s and %s %s"%(result.position[0],result.position[1],br)
        summary += "Repeat was found %s times in a sequence of %s amino acids"%(len(result.instances),len(result.sequence))+br
        coverage = round((float(len(result.instances)*len(result.instances[0]))/len(result.sequence))*100,2)
        summary += "Coverage of %s %s"%(coverage,"%")+br
        summary += "Instances:<br> %s%s"%(format_instances(result.instances,instance_line_nr),br)
        summary += "pattern: %s%s"%(result.pattern,br)
        seqprint = True
    if result.evidence:
        summary_text_printed = False
        for ek in result.evidence:
            
            if len(result.evidence[ek]) > 1: 
                if summary == "":
                    summary += hr
                if not summary_text_printed:
                    summary += "the following known Motifs were found: <br>"
                    summary_text_printed = True
                    
                summary += "%s was found %s times in this sequence<br>"%(ek,len(result.evidence[ek]))
                seqprint = True
    if seqprint:
        html_seq = highlight_seq(result.sequence,result.instances)
        summary += "Sequence: <br>"#<br> in html
        summary += format_linebreaks(html_seq,64)
    return summary

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
    print "line at end of formatting is %s"%(line)
    return line

def create_result_output(result):
    #highlight_seq(feature,=feature.qualifiers["glutamine_aromat_structures"][0],output_key="html_seq_glutamine")
    output_line = write_result_summary(result)
    #output_line += write_evidence_summary()
    #output_line += "<br>glutamine sequences highlighted <br>"
    #output_line += format_linebreaks(feature.qualifiers["html_seq_glutamine"],64)
    return output_line

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
    print seqrecord.annotations
    seqrecord.annotations["html_output"] = "Cluster %s from %s with %s potential repeat regions"%(seqrecord.id,seqrecord.description, len(seqrecord.annotations["repeat_regions"]))

    for feat in seqrecord.features:
        if "has_repeat" in feat.qualifiers:
            print(feat.qualifiers["table"])
            print(feat.qualifiers["ripp_evidence"])
            if len(feat.qualifiers["table"]) > 0 or len(feat.qualifiers["ripp_evidence"]) >0: 
                create_output(feat)
                print("condition reached")
                exit(0)
