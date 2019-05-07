#!/usr/bin/env python
# Small script for testing the active site finder library without the whole antiSMASH Overhead

import helperlibs.bio.seqio
import antismash.generic_modules.ecpredictor
import antismash.config
import logging
import sys
from Bio import SeqIO 
from argparse import Namespace
from antismash import utils


def main():

    GenbankFile = ""
    logging.basicConfig(format='%(levelname)s: %(message)s', level="DEBUG")
    
    if len(sys.argv)>1:
        GenbankFile = sys.argv[1]
    else:
        GenbankFile = "NC_021985.gb"
    #GenbankFile = "/Users/karl/Documents/test/NC_018750.1.final.gbk"
    #OutFile = "/Users/karl/Documents/test/test.gbk"
    OutFile = "test.gbk"


    options = Namespace()
    antismash.config.load_config(options)
    options.statusfile = "test.status"
    options.outputfoldername = "."
    options.cpus = 8
    
    if len(sys.argv)>2:
        options.ecpred = sys.argv[2]
    else:
        options.ecpred = "eficaz"

    seq_records = helperlibs.bio.seqio.parse(GenbankFile)
    out_recs = []
    for seq_record in seq_records:
        antismash.generic_modules.ecpredictor.run(seq_record, options)
        out_recs.append(seq_record)

    logging.debug("Writing output file")
    
    SeqIO.write(out_recs, OutFile, "genbank")
    

if __name__ == "__main__":
    main()