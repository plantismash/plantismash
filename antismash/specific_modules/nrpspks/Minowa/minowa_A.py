#!/usr/bin/python2
## Author: Marnix Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre

import os
import sys
from os import path
from antismash import utils

##Functions used in this program
# Function that reads the fasta file into a dictionary
def fastadict(fasta):
    file = open(fasta,"r")
    text = file.read()
    text = text.replace("\r","\n")
    text = text.strip()
    #Replaces all spaces with "_" to avoid problems
    text = text.replace(' ','_')
    text = text.split()
    dictseq = {}
    for a in text:
        if ">" in a[0]:
            f = str()
            d = a[1:64]
        else:
            e = a
            f += e
            dictseq[d] = f
    return dictseq

# Function that extracts all sequence names from the fasta dictionary
def lnames(fastadict):
    items = fastadict.items()
    items.sort()
    return [names for names, seqs in items]


# Function that extracts all sequences from the fasta dictionary
def lseqs(fastadict):
    items = fastadict.items()
    items.sort()
    return [seqs for names, seqs in items]

def fastanames(fasta):
    names = []
    file = open(fasta,"r")
    text = file.read()
    text = text.replace("\r","\n")
    text = text.strip()
    #Replaces all spaces with "_" to avoid problems
    text = text.replace(' ','_')
    text = text.split()
    for a in text:
        if ">" in a[0]:
            f = str()
            d = a[1:]
            names.append(d)
    return names

def fastaseqs(names,fastadict):
    seqs = []
    for i in names:
        seq = fastadict[i]
        seqs.append(seq)
    return seqs

def writefasta(names,seqs,file):
    e = 0
    f = len(names) - 1
    out_file = open(file,"w")
    while e <= f:
        out_file.write(">")
        out_file.write(names[e])
        out_file.write("\n")
        out_file.write(seqs[e])
        out_file.write("\n")
        e += 1
    out_file.close()

def hmmsearch(fasta,hmm):
    lsname = fastanames(fasta)[0]
    text, err, retcode = utils.execute(["hmmsearch", "--noali", hmm, fasta])
    text = text.replace("\r","\n")
    start = text.find('Domain annotation for each sequence:')
    end = text.find('Internal pipeline statistics summary:')
    lines = []
    ls_names = []
    ls_domain_nrs = []
    ls_starts = []
    ls_ends = []
    ls_scores = []
    ls_evalues = []
    lines = text[start:end].split('\n')
    if "[No targets detected" in text:
        ls_names.append(lsname)
        ls_scores.append(str(0))
    else:
        lines = lines[4:-4]
        for i in lines:
            tabs = i.split(" ")
            tabs2 = []
            for i in tabs:
                if i != "":
                    tabs2.append(i)
            ls_names.append(lsname)
            ls_domain_nrs.append(tabs2[0])
            ls_starts.append(tabs2[6])
            ls_ends.append(tabs2[7])
            ls_scores.append(tabs2[2])
            ls_evalues.append(tabs2[4])
    dicthmm = {}
    for i in ls_names:
        j = ls_names.index(i)
        dicthmm[i] = ls_scores[j]
    return dicthmm

def hmmnames(hmmdict):
    items = hmmdict.items()
    return [names for names, scores in items]

def hmmscores(hmmdict):
    items = hmmdict.items()
    return [scores for names, scores in items]


def sortdictkeysbyvalues(dict):
    items = [(value, key) for key, value in dict.items()]
    items.sort()
    items.reverse()
    return [key for value, key in items]

def run_minowa_a(infile2, outfile):
    ## Core
    infile = utils.get_full_path(__file__, "A_domains_muscle.fasta")
    muscle_file = "muscle.fasta"
    out_file = open(outfile, "w")
    dict2 = fastadict(infile2)
    namesa = fastanames(infile2)
    seqsa = fastaseqs(namesa,dict2)
    startpos = 65
    namesb = namesa
    seqsb = seqsa

    for i in namesb:
        seq = seqsb[namesb.index(i)]
        writefasta([i],[seq],"infile.fasta")
        infile2 = "infile.fasta"
        out_file.write("\\\\" + "\n" + i + "\n")
        refsequence = "P0C062_A1"
        namesa = [i]
        seqsa = [seq]

        #Run muscle and collect sequence positions from file
        utils.execute(["muscle", "-profile", "-quiet", "-in1", infile, "-in2", infile2, "-out", "muscle.fasta"])
        file = open(utils.get_full_path(__file__, "Apositions.txt"), "r")
        text = file.read()
        text = text.replace("\r","\n")
        text = text.strip()
        text = text.replace(' ','_')
        positions = text.split("\t")
        positions2 = []
        for i in positions:
            pos = int(i)
            pos = pos - startpos
            positions2.append(pos)
        positions = positions2

        #Count residues in ref sequence and put positions in list
        muscle_dict = fastadict(muscle_file)
        muscle_seqs = lseqs(muscle_dict)
        muscle_names = lnames(muscle_dict)
        refseqnr = muscle_names.index(refsequence)
        refseq = muscle_seqs[refseqnr]
        poslist = []
        a = 0
        b = 0
        c = 0
        while refseq != "":
            i = refseq[0]
            if c in positions and i != "-":
                poslist.append(b)
            if i != "-":
                c += 1
            b += 1
            refseq = refseq[1:]

        #Extract positions from query sequence and create fasta file to use as input for hmm searches
        query = namesa[0]
        query_seqnr = muscle_names.index(query)
        query_seq = muscle_seqs[query_seqnr]
        seq = ""
        for j in poslist:
            aa = query_seq[j]
            if aa == "-":
                aa = "X"
            seq = seq + aa
        query_names = []
        query_names.append(query)
        query_seqs = []
        query_seqs.append(seq)
        writefasta(query_names,query_seqs,"hmm_infile.fasta")
        #- then use list to extract positions from every sequence -> HMMs (one time, without any query sequence)


        #Compare scores and output prediction
        hmm_names = []
        hmm_scores = []
        a_hmms_dir = utils.get_full_path(__file__, 'A_HMMs')
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "2-3-diaminopropionate.hmm"))
        hmmname = "2-3-diaminoproprionate"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "3mGlu.hmm"))
        hmmname = "3mGlu"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "5NhOrn.hmm"))
        hmmname = "5NhOrn"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Abu.hmm"))
        hmmname = "Abu"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Ahp.hmm"))
        hmmname = "Ahp"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "alaninol.hmm"))
        hmmname = "alaninol"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Asn.hmm"))
        hmmname = "Asn"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Asp.hmm"))
        hmmname = "Asp"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "beta-Lys.hmm"))
        hmmname = "beta-Lys"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "bOHTyr.hmm"))
        hmmname = "bOHTyr"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Cys.hmm"))
        hmmname = "Cys"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "DHB.hmm"))
        hmmname = "DHB"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "fOHOrn.hmm"))
        hmmname = "fOHOrn"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Glu.hmm"))
        hmmname = "Glu"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "guanidinoacetic_acid.hmm"))
        hmmname = "guanidinoacetic_acid"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "His.hmm"))
        hmmname = "His"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Hpg2Cl.hmm"))
        hmmname = "Hpg2Cl"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Ile.hmm"))
        hmmname = "Ile"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Leu.hmm"))
        hmmname = "Leu"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "MeAsp.hmm"))
        hmmname = "MeAsp"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "OHOrn.hmm"))
        hmmname = "OHOrn"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Orn.hmm"))
        hmmname = "Orn"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Phenylacetate.hmm"))
        hmmname = "Phenylacetate"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Pro.hmm"))
        hmmname = "Pro"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Qna.hmm"))
        hmmname = "Qna"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Sar.hmm"))
        hmmname = "Sar"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Thr-4-Cl.hmm"))
        hmmname = "Thr-4-Cl"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Trp.hmm"))
        hmmname = "Trp"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Val.hmm"))
        hmmname = "Val"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "3-HPA.hmm"))
        hmmname = "3-HPA"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "4-MHA.hmm"))
        hmmname = "4-MHA"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Aad.hmm"))
        hmmname = "Aad"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Aeo.hmm"))
        hmmname = "Aeo"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Ala.hmm"))
        hmmname = "Ala"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Arg.hmm"))
        hmmname = "Arg"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "B-Ala.hmm"))
        hmmname = "B-Ala"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Bmt.hmm"))
        hmmname = "Bmt"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "capreomycidine.hmm"))
        hmmname = "capreomycidine"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Dab.hmm"))
        hmmname = "Dab"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "DHpg.hmm"))
        hmmname = "DHpg"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Gln.hmm"))
        hmmname = "Gln"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Gly.hmm"))
        hmmname = "Gly"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "hAsn.hmm"))
        hmmname = "hAsn"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "homoTyr.hmm"))
        hmmname = "homoTyr"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Hpg.hmm"))
        hmmname = "Hpg"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Kyn.hmm"))
        hmmname = "Kyn"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Lys.hmm"))
        hmmname = "Lys"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "mPro.hmm"))
        hmmname = "mPro"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "OmAsp.hmm"))
        hmmname = "OmAsp"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Phe.hmm"))
        hmmname = "Phe"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "pipecolate.hmm"))
        hmmname = "pipecolate"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "QA.hmm"))
        hmmname = "QA"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Sal.hmm"))
        hmmname = "Sal"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Ser.hmm"))
        hmmname = "Ser"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Thr.hmm"))
        hmmname = "Thr"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])
        hmmresults = hmmsearch("hmm_infile.fasta", path.join(a_hmms_dir, "Tyr.hmm"))
        hmmname = "Tyr"
        hmm_names.append(hmmname)
        hmmscore = hmmscores(hmmresults)
        hmm_scores.append(hmmscore[0])

        #Sort names & scores by scores:
        scoredict = {}
        a = 0
        for i in hmm_names:
            score = hmm_scores[a]
            scoredict[i] = float(score)
            a += 1
        hmm_names = sortdictkeysbyvalues(scoredict)
        hmm_scores = []
        for i in hmm_names:
            score = str(scoredict[i])
            hmm_scores.append(score)

        out_file.write("Substrate:")
        out_file.write("\t")
        out_file.write("Score:")
        out_file.write("\n")
        for i in hmm_names:
            out_file.write(i)
            out_file.write("\t")
            j = hmm_names.index(i)
            score = hmm_scores[j]
            out_file.write(score)
            out_file.write("\n")
        out_file.write("\n")
