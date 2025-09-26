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

import sys
import os
import shutil
import subprocess
import time
import logging
import signal
from antismash import utils

def smcog_analysis(inputgenes, inputnr, seq_record, smcogdict, smcogsoutputfolder):
    "run smCOG search on all gene cluster CDS features"
    for feature in inputgenes:
        k = utils.get_gene_id(feature)
        tag = k
        seq = str(utils.get_aa_sequence(feature))
        #create input.fasta file with single query sequence to be used as input for MSA
        utils.writefasta([tag], [seq], "input" + str(inputnr) + ".fasta")
        if k in smcogdict and len(smcogdict[k]) > 0:
            smcog = (smcogdict[k][0][0]).split(":")[0]
            alignsmcogs(smcog, inputnr)
            #Generate trimmed alignment
            trimalignment(inputnr)
            #Draw phylogenetic tree
            drawtree(inputnr)
            #Convert tree to draw PNG image
            converttree(inputnr, smcogsoutputfolder, tag)

def alignsmcogs(smcog, inputnr):
     #Align to multiple sequence alignment, output as fasta file
     infile1 = utils.get_full_path(__file__, "%s_muscle.fasta" % str(smcog).lower())
     if sys.platform == ('linux2') or sys.platform == ('win32'):
         musclecommand = ["muscle", "-quiet", "-profile", "-in1", infile1, "-in2", "input" + str(inputnr) + ".fasta", "-out", "muscle" + str(inputnr) + ".fasta"]
     elif sys.platform == ('darwin'):
         musclecommand = ["muscle", "-quiet", "-profile", "-in1", infile1, "-in2", "input" + str(inputnr) + ".fasta", "-out", "muscle" + str(inputnr) + ".fasta"]
     utils.execute(musclecommand)

def trimalignment(inputnr):
    #Trim alignment
    #edit muscle fasta file: remove all positions before the first and after the last position shared by >33% of all sequences
    musclefile = open("muscle" + str(inputnr) + ".fasta","r")
    filetext = musclefile.read()
    filetext = filetext.replace("\r","\n")
    lines = filetext.split("\n")
    ##Combine all sequence lines into single lines
    lines2 = []
    seq = ""
    nrlines = len(lines)
    a = 0
    lines = lines[:-1]
    for i in lines:
        if a == (nrlines - 2):
            seq = seq + i
            lines2.append(seq)
        if i[0] == ">":
            lines2.append(seq)
            seq = ""
            lines2.append(i)
        else:
            seq = seq + i
        a += 1
    lines = lines2[1:]
    #Retrieve names and seqs from muscle fasta lines
    seqs = []
    names = []
    for i in lines:
        if len(i) > 0 and i[0] == ">":
            name = i[1:]
            names.append(name)
        else:
            seq = i
            seqs.append(seq)
    #Find first and last amino acids shared conserved >33%
    #Create list system to store conservation of residues
    conservationlist = []
    lenseqs = len(seqs[0])
    nrseqs = len(seqs)
    for i in range(lenseqs):
        conservationlist.append({"A":0,"B":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"J":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"U":0,"V":0,"W":0,"X":0,"Y":0,"Z":0,"-":0})
    a = 0
    for i in seqs:
        aa = list(i)
        for i in aa:
            conservationlist[a][i] += 1
            a += 1
        a = 0
    firstsharedaa = 0
    lastsharedaa = lenseqs
    #Find first amino acid shared
    first = "yes"
    nr = 0
    for i in conservationlist:
        aa = utils.sortdictkeysbyvaluesrev(i)
        if aa[0] != "-" and i[aa[1]] > (nrseqs / 3) and first == "yes":
            firstsharedaa = nr
            first = "no"
        nr += 1
    #Find last amino acid shared
    conservationlist.reverse()
    first = "yes"
    nr = 0
    for i in conservationlist:
        aa = utils.sortdictkeysbyvaluesrev(i)
        if aa[0] != "-" and i[aa[1]] > (nrseqs / 3) and first == "yes":
            lastsharedaa = lenseqs - nr
            first = "no"
        nr += 1
    #Shorten sequences to detected conserved regions
    seqs2 = []
    for i in seqs:
        seq = i[firstsharedaa:lastsharedaa]
        seqs2.append(seq)
    seqs = seqs2
    seedfastaname = "trimmed_alignment" + str(inputnr) + ".fasta"
    utils.writefasta(names, seqs, seedfastaname)

def drawtree(inputnr):
    #Draw phylogenetic tree with fasttree 2.1.1
    nwkfile = "tree" + str(inputnr) + ".nwk"
    fasttreecommand = "fasttree -quiet -fastest -noml trimmed_alignment" + str(inputnr) + ".fasta > " + nwkfile
    os.system(fasttreecommand)

def converttree(inputnr, smcogsoutputfolder, tag):
     #Convert tree to XTG and draw PNG image using TreeGraph
     command = ['java', '-Djava.awt.headless=true', '-jar', utils.get_full_path(__file__, 'TreeGraph.jar'),
                '-convert', 'tree%s.nwk'% inputnr, '-xtg', 'tree%s.xtg' % inputnr ]
     p = subprocess.Popen(command, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
     processes_starttime = time.time()
     while True:
         if (time.time() - processes_starttime) > 1200:
             if sys.platform == ('linux2') or sys.platform == ('darwin'):
                 os.kill(p.pid,signal.SIGKILL)
                 logging.info("Now in " + os.getcwd() + " TreeGraph -convert on tree" + str(inputnr) + " ran out out of time")
                 break
             elif sys.platform == ('win32'):
                 subprocess.Popen("taskkill /F /T /PID %i"%p.pid , shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                 logging.info("Now in " + os.getcwd() + " TreeGraph -convert on tree " + str(inputnr) + " ran out out of time")
                 break
         if p.poll() == 0:
             break
         time.sleep(2)
     out, err = p.communicate()
     output = out
     if "exception" not in output and "Exception" not in output:
         command = ['java', '-Djava.awt.headless=true', '-jar', utils.get_full_path(__file__, 'TreeGraph.jar'),
                    '-image', 'tree%s.xtg'% inputnr, "%s.png" % tag.split('.')[0] ]
         p = subprocess.Popen(command, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
         processes_starttime = time.time()
         while True:
             if (time.time() - processes_starttime) > 1200:
                 if sys.platform == ('linux2') or sys.platform == ('darwin'):
                     os.kill(p.pid,signal.SIGKILL)
                     logging.info("Now in " + os.getcwd() + " TreeGraph -image on tree " + str(inputnr) + " ran out out of time")
                     break
                 elif sys.platform == ('win32'):
                     subprocess.Popen("taskkill /F /T /PID %i"%p.pid , shell=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                     logging.info("Now in " + os.getcwd() + " TreeGraph -image on tree " + str(inputnr) + " ran out out of time")
                     break
             if p.poll() == 0:
                 break
             time.sleep(2)
         out, err = p.communicate()
         output = out
         if "exception" not in output and "Exception" not in output:
             shutil.copy(tag.split(".")[0] + '.png', smcogsoutputfolder)
             os.remove(tag.split(".")[0] + ".png")
             os.remove("tree" + str(inputnr) + ".xtg")
             os.remove("trimmed_alignment" + str(inputnr) + ".fasta")
