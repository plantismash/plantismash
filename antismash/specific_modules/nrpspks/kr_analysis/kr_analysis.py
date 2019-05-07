#!/usr/bin/env python
## Author: Marnix Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre

import os
from antismash import utils

##Functions used in this program
# Function that reads the fasta file into a dictionary
def fastadict(fasta):
  file = open(fasta,"r")
  text = file.read()
  text = text.strip()
  #Replaces all spaces with "_" to avoid problems
  text = text.replace(' ','_')
  text = text.split()
  dictseq = {}
  for a in text:
    if ">" in a[0]:
      f = str()
      d = a[1:68]
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

def run_kr_analysis(infile2, out_file):
    ##Core script
    #Extract activity and stereochemistry signatures from KR domains
    infile = utils.get_full_path(__file__, "KRdomains_muscle.fasta")
    muscle_file = "muscle.fasta"
    dict2 = fastadict(infile2)
    namesb = fastanames(infile2)
    seqsb = fastaseqs(namesb,dict2)
    #startpos = 2
    querysignames = []
    querysigseqs_act = []
    querysigseqs_ste = []
    for i in namesb:
      seq = seqsb[namesb.index(i)]
      querysignames.append(i)
      writefasta([i],[seq],"infile.fasta")
      infile2 = "infile.fasta"
      refsequence = "MAPSI|PKS|CAM00062.1|Erythromycin_synthase_modules_1_and_2|Sacc_KR1"
      namesa = [i]
      #Run muscle and collect sequence positions from file
      utils.execute(["muscle", "-profile", "-quiet", "-in1", infile, "-in2", infile2, "-out", "muscle.fasta"])
      positions_act = [110,134,147,151]
      positions_ste = [90,91,92,139,144,147,149,151]
      #Count residues in ref sequence and put positions in list
      muscle_dict = fastadict(muscle_file)
      muscle_seqs = lseqs(muscle_dict)
      muscle_names = lnames(muscle_dict)
      refseqnr = muscle_names.index(refsequence)
      #Extract activity signature
      refseq = muscle_seqs[refseqnr]
      poslist_act = []
      b = 0
      c = 0
      while refseq != "":
        i = refseq[0]
        if c in positions_act and i != "-":
          poslist_act.append(b)
        if i != "-":
          c += 1
        b += 1
        refseq = refseq[1:]
      #Extract stereochemistry signature
      refseq = muscle_seqs[refseqnr]
      poslist_ste = []
      b = 0
      c = 0
      while refseq != "":
        i = refseq[0]
        if c in positions_ste and i != "-":
          poslist_ste.append(b)
        if i != "-":
          c += 1
        b += 1
        refseq = refseq[1:]
      #Extract positions from query sequence
      query = namesa[0]
      query_seqnr = muscle_names.index(query)
      query_seq = muscle_seqs[query_seqnr]
      seq_act = ""
      seq_ste = ""
      for j in poslist_act:
        aa = query_seq[j]
        seq_act = seq_act + aa
      querysigseqs_act.append(seq_act)
      for j in poslist_ste:
        aa = query_seq[j]
        seq_ste = seq_ste + aa
      querysigseqs_ste.append(seq_ste)

    #Check activity
    activitydict = {}
    for i in querysignames:
      querysigseq_act = querysigseqs_act[querysignames.index(i)]
      activity = "inactive"
      if querysigseq_act[0] == "K" and (querysigseq_act[1] == "S" or querysigseq_act[1] == "A" or querysigseq_act[1] == "G") and querysigseq_act[2] == "Y" and querysigseq_act[3] == "N":
        activity = "active"
      if querysigseq_act[0] == "E" and (querysigseq_act[1] == "S" or querysigseq_act[1] == "A" or querysigseq_act[1] == "G") and querysigseq_act[2] == "H" and querysigseq_act[3] == "H":
        activity = "active"
      if querysigseq_act[0] == "K" and (querysigseq_act[1] == "S" or querysigseq_act[1] == "A" or querysigseq_act[1] == "G") and querysigseq_act[2] == "Y" and (querysigseq_act[3] == "N" or querysigseq_act[3] == "G"):
        activity = "active"
      activitydict[i] = activity

    #Predict stereochemistry
    stereodict = {}
    for i in querysignames:
      querysigseq_ste = querysigseqs_ste[querysignames.index(i)]
      if querysigseq_ste[0:3] != "LDD" and querysigseq_ste[3] == "W" and querysigseq_ste[4] != "H" and querysigseq_ste[5:] == "YAN":
        stereochemistry = "A1"
      elif querysigseq_ste[0:3] != "LDD" and querysigseq_ste[3:] == "WHYAN":
        stereochemistry = "A2"
      elif querysigseq_ste[0:3] == "LDD" and querysigseq_ste[5] == "Y" and querysigseq_ste[6] != "P" and querysigseq_ste[7] == "N":
        stereochemistry = "B1"
      elif querysigseq_ste[0:3] == "LDD" and querysigseq_ste[5:] == "YPN":
        stereochemistry = "B2"
      elif querysigseq_ste[5] != "Y":
        stereochemistry = "C1"
      elif querysigseq_ste[5] == "Y" and querysigseq_ste[7] != "N":
        stereochemistry = "C2"
      else:
        stereochemistry = "?"
      stereodict[i] = stereochemistry

    #Output to file
    outfile = open(out_file,"w")
    for i in querysignames:
      outfile.write(i + "\t" + activitydict[i] + "\t" + stereodict[i] + "\n")
