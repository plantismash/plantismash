# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011,2012 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

def purifydomainlist(domainlist, hmmlengthsdict):
    #Purify domain list to remove overlapping domains, only keeping those with the highest scores 
    domainlist2 = [domainlist[0]]
    for i in domainlist[1:]:
        maxoverlap = 0.20 * max([hmmlengthsdict[i[0]], hmmlengthsdict[domainlist2[-1][0]]])
        if i[1] < (domainlist2[-1][2] - maxoverlap):
            if i[4] < domainlist2[-1][4]:
                pass
            elif i[4] > domainlist2[-1][4]:
                del domainlist2[-1]
                domainlist2.append(i)
        else:
            domainlist2.append(i)
    return domainlist2

def remove_incomplete_domains(domainlist,hmmlengthsdict):
    domainlist2 = []
    for i in domainlist:
        alilength = int(i[2]) - int(i[1])
        domainlength = hmmlengthsdict[i[0]]
        if alilength > (0.5 * domainlength):
            domainlist2.append(i)
    return domainlist2

def mergedomainlist(domainlist,hmmlengthsdict):
    domainlist2 = [domainlist[0]]
    for i in domainlist[1:]:
        alilength1 = int(domainlist2[-1][2]) - int(domainlist2[-1][1])
        alilength2 = int(i[2]) - int(i[1])
        domainlength = hmmlengthsdict[i[0]]
        if i[0] == domainlist2[-1][0] and (alilength1 < (0.75 * domainlength) or alilength2 < (0.75 * domainlength)) and (alilength1 + alilength2) < (1.5 * domainlength):
            name = i[0]
            start = domainlist2[-1][1]
            end = i[2]
            evalue = str(float(domainlist2[-1][3]) * float(i[3]))
            score = str(float(domainlist2[-1][4]) + float(i[4]))
            del domainlist2[-1]
            domainlist2.append([name,start,end,evalue,score])
        else:
            domainlist2.append(i)
    return domainlist2

def parse_hmmscan_results(hmmscan_results, hmmlengthsdict):
    hmmscandict = {}
    results_by_id = {}
    for results in hmmscan_results:
        for hsp in results.hsps:
            if not results_by_id.has_key(hsp.query_id):
                results_by_id[hsp.query_id] = [hsp]
            else:
                if not hsp in results_by_id[hsp.query_id]:
                    results_by_id[hsp.query_id].append(hsp)
    for cds in results_by_id.keys():
        domainlist = []
        results = results_by_id[cds]
        for result in results:
            domainlist.append([result.hit_id, result.query_start, result.query_end, result.evalue, result.bitscore])
        domainlist.sort(lambda a, b: cmp(a[1], b[1]))
        #Only keep best hits for overlapping domains
        if len(domainlist) > 1:
            domainlist = purifydomainlist(domainlist, hmmlengthsdict)
        #Merge domain fragments which are really one domain
        if len(domainlist) > 1:
            domainlist = mergedomainlist(domainlist,hmmlengthsdict)
        #Remove incomplete domains (covering less than 60% of total domain hmm length)
        if len(domainlist) > 1:
            domainlist = remove_incomplete_domains(domainlist,hmmlengthsdict)
        if len(domainlist) > 0:
            hmmscandict[cds] = domainlist
    return hmmscandict