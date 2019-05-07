#!/usr/bin/env python
## Author: Marnix Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre

##Functions used in this program
# Function that reads the fasta file into a dictionary

import sys
import os
from antismash import utils

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

def sortdictkeysbyvalues(dict):
    items = [(value, key) for key, value in dict.items()]
    items.sort()
    items.reverse()
    return [key for value, key in items]

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

def run_nrpscodepred(options):
    #Extract 10 / 34 AA NRPS signatures from A domains
    infile = utils.get_full_path(__file__, "A_domains_muscle.fasta")
    infile2 = "nrpsseqs.fasta"
    out_file1 = "input.sig"
    out_file2 = "ctg" + str(options.record_idx) + "_nrpspredictor2_codes.txt"
    muscle_file = "muscle.fasta"
    dict = fastadict(infile)
    dict2 = fastadict(infile2)
    seqs = lseqs(dict)
    names = lnames(dict)
    namesb = fastanames(infile2)
    seqsb = fastaseqs(namesb,dict2)
    startpos = 66
    querysignames = []
    querysigseqs = []
    querysig34codes = []
    illegalcharacters = """!@#$%^&*(){}:"<>?/.,';][`~1234567890*-+-=_\|"""
    for i in namesb:
      seq = seqsb[namesb.index(i)]
      for char in seq:
        if char in illegalcharacters:
            seq = seq.replace(char, "X")
        if len(seq) == 0:
            continue
      querysignames.append(i)
      writefasta([i],[seq],"infile.fasta")
      infile2 = "infile.fasta"
      refsequence = "P0C062_A1"
      namesa = [i]
      seqsa = [seq]
      #Run muscle and collect sequence positions from file
      os.system('muscle -profile -quiet -in1 "' + infile + '" -in2 ' + infile2 + ' -out muscle.fasta')
      file = open(utils.get_full_path(__file__, "A34positions.txt"), "r")
      text = file.read()
      text = text.strip()
      text = text.replace(' ','_')
      angpositions = text.split("\t")
      angpositions2 = []
      for i in angpositions:
        angpos = int(i)
        angpos = angpos - startpos
        angpositions2.append(angpos)
      angpositions = angpositions2
      file = open(utils.get_full_path(__file__, "Apositions.txt"), "r")
      text = file.read()
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
      #Extract positions from query sequence
      query = namesa[0]
      query_seqnr = muscle_names.index(query)
      query_seq = muscle_seqs[query_seqnr]
      seq = ""
      for j in poslist:
        aa = query_seq[j]
        k = j
        l = j
        if aa == "-":
          k += 1
          l = l - 1
          if l not in poslist:
            aa = query_seq[l]
          elif (j + 1) not in poslist:
            aa = query_seq[k]
        seq = seq + aa
      #Add fixed lysine 517
      seq = seq + "K"
      querysigseqs.append(seq)
      #Do the same to extract 34 AA codes
      refseqnr = muscle_names.index(refsequence)
      refseq = muscle_seqs[refseqnr]
      poslist = []
      a = 0
      b = 0
      c = 0
      while refseq != "":
        i = refseq[0]
        if c in angpositions and i != "-":
          poslist.append(b)
        if i != "-":
          c += 1
        b += 1
        refseq = refseq[1:]
      #Extract positions from query sequence
      query = namesa[0]
      query_seqnr = muscle_names.index(query)
      query_seq = muscle_seqs[query_seqnr]
      seq = ""
      for j in poslist:
        aa = query_seq[j]
        k = j
        l = j
        if aa == "-":
          k += 1
          l = l - 1
          if l not in poslist:
            aa = query_seq[l]
          elif (j + 1) not in poslist:
            aa = query_seq[k]
        seq = seq + aa
      querysig34codes.append(seq)

    #Load reference NRPS signatures
    infile3 = utils.get_full_path(__file__, "knowncodes.fasta")
    signaturesdict = fastadict(infile3)
    signaturenames = fastanames(infile3)
    signatureseqs = fastaseqs(signaturenames,signaturesdict)

    outfile1 = open(out_file1,"w")
    outfile2 = open(out_file2,"w")
    #Compare NRPS signature with database of signatures and write output to txt file
    for k in querysignames:
      querysigseq = querysigseqs[querysignames.index(k)]
      scoredict = {}
      for i in signaturenames:
        sigseq = signatureseqs[signaturenames.index(i)]
        positions  = range(len(querysigseq))
        score = 0
        for j in positions:
          if querysigseq[j] == sigseq[j]:
            score += 1
        score = ((float(score) / 10) * 100)
        scoredict[i] = score
      sortedhits = sortdictkeysbyvalues(scoredict)
      sortedhits = sortedhits
      sortedscores = []
      sortedhits2 = []
      for i in sortedhits:
        score = scoredict[i]
        if score > 40:
          score = "%.0f"%(score)
          sortedscores.append(score)
          sortedhits2.append(i)
      allsortedhits = sortedhits
      sortedhits = sortedhits2
      #Find all other scores after best score
      nextbesthitsdict = {}
      nextbesthits = []
      for hit in allsortedhits:
        aa = hit.split("__")[-1]
        score = scoredict[hit]
        if not nextbesthitsdict.has_key(aa) and aa != "xxx":
          nextbesthitsdict[aa] = score
          nextbesthits.append(aa)
      #Write output to txt file
      outfile1.write(querysig34codes[querysignames.index(k)] + "\t" + k + "\n")
      if len(sortedhits) > 0:
        outfile2.write(k + "\t" + sortedhits[0].split("__")[-1] + "\t")
      else:
        outfile2.write(k + "\t" + "N/A" + "\t")
      outfile2.write(";".join(["%s(%s)" % (aa, nextbesthitsdict[aa]) for aa in nextbesthits]) + "\n")
    outfile1.close()
    outfile2.close()

    return
