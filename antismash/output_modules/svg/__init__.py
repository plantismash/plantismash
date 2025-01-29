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

"""SVG output format module

"""
import logging
import os
from os import path
import shutil
from antismash import utils
from .data_loading import prepare_visualization
from .svg_drawer import create_svgs

name = "svg"
short_description = "SVG output"
priority = 1

def write(seq_records, options):
    basename = options.outputfoldername
    options.svgdir = path.join(basename, "svg")
    logging.debug("Writing seq_records SVGs to %r" % options.svgdir)
    if not path.exists(options.svgdir):
        os.mkdir(options.svgdir)
    for seq_record in seq_records:
        if len(utils.get_cluster_features(seq_record)) > 0:
            #Parse clusterblast output to prepare visualization
            prepare_visualization(options, seq_record)
            create_svgs(options, seq_record)


