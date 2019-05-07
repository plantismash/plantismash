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


def run_pkssignature_analysis(infile2, outfile):
    ##Core script
    #Extract PKS signature from AT domains
    infile = utils.get_full_path(__file__, "AT_domains_muscle.fasta")
    muscle_file = "muscle.fasta"
    dict2 = fastadict(infile2)
    namesb = fastanames(infile2)
    seqsb = fastaseqs(namesb,dict2)
    startpos = 7
    querysignames = []
    querysigseqs = []
    for i in namesb:
      seq = seqsb[namesb.index(i)]
      querysignames.append(i)
      writefasta([i],[seq],"infile.fasta")
      infile2 = "infile.fasta"
      refsequence = "P0AAI9_AT1"
      namesa = [i]
      #Run muscle and collect sequence positions from file
      utils.execute(["muscle", "-profile", "-quiet", "-in1", infile, "-in2", infile2, "-out", "muscle.fasta"])
      file = open(utils.get_full_path(__file__, "ATpositions.txt"), "r")
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
        seq = seq + aa
      querysigseqs.append(seq)

    #Load reference PKS signatures
    infile3 = utils.get_full_path(__file__, "pks_signatures.fasta")
    signaturesdict = fastadict(infile3)
    signaturenames = fastanames(infile3)
    signatureseqs = fastaseqs(signaturenames,signaturesdict)

    out_file = open(outfile,"w")
    #Compare PKS signature with database of signatures and write output to txt file
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
        score = ((float(score) / 24) * 100)
        scoredict[i] = score
      sortedhits = sortdictkeysbyvalues(scoredict)
      sortedhits = sortedhits[:10]
      sortedscores = []
      sortedhits2 = []
      for i in sortedhits:
        score = scoredict[i]
        if score > 50:
          score = "%.0f"%(score)
          sortedscores.append(score)
          sortedhits2.append(i)
      sortedhits = sortedhits2
      #Write output to txt file
      out_file.write("//\n" + k + "\t" + querysigseq + "\n")
      a = 0
      for i in sortedhits:
        out_file.write(i + "\t" + signatureseqs[signaturenames.index(i)] + "\t" + sortedscores[a] + "\n")
        a += 1
      out_file.write("\n\n")
