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

"""GenBank output format module

"""
import logging
from helperlibs.bio import seqio
from os import path

name = "embl"
short_description = "EMBL output"
priority = 1

def write(seq_records, options):
    basename = seq_records[0].id
    output_name = path.join(options.outputfoldername, "%s.final.embl" % basename)
    logging.debug("Writing seq_records to %r" % output_name)
    if options.input_type == 'nucl':
        seqio.write(seq_records, output_name, 'embl')