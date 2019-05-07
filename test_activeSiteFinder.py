#!/usr/bin/env python
# Small script for testing the active site finder library without the whole antiSMASH Overhead

import helperlibs.bio.seqio
import antismash.generic_modules.active_site_finder
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
        if utils.locate_file(sys.argv[1]):
            GenbankFile = sys.argv[1]
    else:
        GenbankFile = "/Users/karl/test/test.gb"
    #GenbankFile = "/Users/karl/Documents/test/NC_018750.1.final.gbk"
    #OutFile = "/Users/karl/Documents/test/test.gbk"
    OutFile = "/Users/karl/test-asf.gbk"

    options = Namespace()
    antismash.config.load_config(options)

    seq_records = helperlibs.bio.seqio.parse(GenbankFile)
    out_recs = []
    for seq_record in seq_records:
        ASFObj = antismash.generic_modules.active_site_finder.active_site_finder(seq_record, options)
        
        status = ASFObj.execute()
        
        if status == True:
            logging.debug("preparation/execution successful")
            out_recs.append(seq_record)
        else:
            logging.debug("Error in preparation")
            
    logging.debug("Writing output file")
    
    SeqIO.write(out_recs, OutFile, "genbank")
    

if __name__ == "__main__":
    main()