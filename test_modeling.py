#!/usr/bin/env python
# Small script for testing the active site finder library without the whole antiSMASH Overhead

import antismash.generic_genome_modules.metabolicmodel.automodel
import antismash.output_modules.metabolicmodel
import antismash.config
import logging
import sys
from Bio import SeqIO 
from argparse import Namespace
from antismash import utils
import os
from Bio.Alphabet import generic_dna


def main():

    GenbankFile = ""
    logging.basicConfig(format='%(levelname)s: %(message)s', level="DEBUG")
    
    if len(sys.argv)>1:
        GenbankFile = sys.argv[1]
    else:
        GenbankFile = "NC_021985EC.gb"
    #GenbankFile = "/Users/karl/Documents/test/NC_018750.1.final.gbk"
    #OutFile = "/Users/karl/Documents/test/test.gbk"
    


    options = Namespace()
    antismash.config.load_config(options)
    options.statusfile = "test.status"
    options.outputfoldername = "modeling"
    options.cpus = 8
    options.metabolicmodeldir = options.outputfoldername
    options.automodel = Namespace()
    options.automodel.solver = "glpk"
    
    if len(sys.argv)>2:
        options.modeling = sys.argv[2]
    else:
        options.modeling = "sco"
    if not os.path.exists(options.metabolicmodeldir):
        os.mkdir(options.metabolicmodeldir)
    fh = open(GenbankFile, "r")
    seq_records = list(SeqIO.parse(fh, "genbank", generic_dna))
    out_recs = []
    logging.debug("There are %s records in file:" % len(seq_records))
  
    antismash.generic_genome_modules.metabolicmodel.automodel.run_automodel(seq_records, options)
    antismash.output_modules.metabolicmodel.write(seq_records, options)
       
    
    
    

if __name__ == "__main__":
    main()