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

"""Gene finding using Glimmer

"""

import logging
from os import path
from antismash import utils
from helperlibs.wrappers.io import TemporaryDirectory
from helperlibs.bio import seqio
from antismash.utils import execute
from Bio.SeqFeature import SeqFeature, FeatureLocation

def run_glimmer(seq_record, options):
    "Run glimmer3 to annotate prokaryotic sequences"
    basedir = utils.get_genefinding_basedir(options)
    with TemporaryDirectory(change=True):
        utils.fix_record_name_id(seq_record, options)
        name = seq_record.id
        while len(name) > 0 and name[0] == '-':
            name = name[1:]
        if name == "":
            name = "unknown"
        fasta_file = '%s.fasta' % name
        longorfs_file = '%s.longorfs' % name
        icm_file = '%s.icm' % name
        result_file = '%s.predict' % name

        # run long-orfs
        with open(fasta_file, 'w') as handle:
            seqio.write([seq_record], handle, 'fasta')
        long_orfs = [path.join(basedir, 'long-orfs')]
        long_orfs.extend(['-l', '-n', '-t', '1.15',
                          '--trans_table', '11',
                          fasta_file,
                          longorfs_file
                         ])
        out, err, retcode = execute(long_orfs)
        if err.find('ERROR') > -1:
            logging.error("Locating long orfs failed: %r" % err)
            return

        # run extract
        extract = [path.join(basedir, 'extract'), '-t', fasta_file,
                   longorfs_file]
        out, err, retcode = execute(extract)
        if out == '':
            logging.error("Failed to extract genes from model, aborting: %r" % err)
            return

        build_icm = [path.join(basedir, 'build-icm'), '-r', icm_file]
        out, err, retcode = execute(build_icm, input=out)
        if err != '':
            logging.error("Failed to build gene model: %r" % err)
            return

        # run glimmer3
        glimmer = [path.join(basedir, 'glimmer3')]
        glimmer.extend(['-l', '-o', '50', '-g', '90', '-q', '3000', '-t', '30',
                        '--trans_table', '11', fasta_file, icm_file, name ])

        out, err, retcode = execute(glimmer)
        if err.find('ERROR') > -1:
            logging.error("Failed to run glimmer3: %r" % err)
            return
        for line in open(result_file, 'r'):
            # skip first line
            if line.startswith('>'):
                continue

            name, start, end, strand, score = line.split()

            try:
                start = int(start)
                end = int(end)
                strand = int(strand)
            except ValueError:
                logging.error('Malformatted glimmer output line %r' % line.rstrip())

            if start > end:
                bpy_strand = -1
                tmp = start
                start = end
                end = tmp
            else:
                bpy_strand = 1

            loc = FeatureLocation(start-1, end, strand=bpy_strand)
            feature = SeqFeature(location=loc, id=name, type="CDS",
                        qualifiers={'locus_tag': ['ctg%s_%s' % (options.record_idx, name)],
                                    'note': ['Glimmer score: %s' %score]})
            seq_record.features.append(feature)
