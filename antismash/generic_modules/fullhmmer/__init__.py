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

"""Full genome PFAM anotation"""

import logging
from os import path
from antismash import utils
from collections import defaultdict
from Bio.SeqFeature import SeqFeature, FeatureLocation
from .name2pfamid import name_to_pfamid

name = "fullhmmer"
short_description = name.capitalize()
priority = 10000


# Tuple is ( binary_name, optional)
_required_binaries = [
    ('hmmscan', False)
]

_required_files = [
    ('Pfam-A.hmm', False),
    ('Pfam-A.hmm.h3f', False),
    ('Pfam-A.hmm.h3i', False),
    ('Pfam-A.hmm.h3m', False),
    ('Pfam-A.hmm.h3p', False)
]

def check_prereqs(options):
    "Check if all required applications are around"
    
    if 'pfamdir' not in options:
        options.pfamdir = utils.get_full_path(__file__, '')

    failure_messages = []
    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)

    for file_name, optional in _required_files:
        if utils.locate_file(path.join(options.pfamdir, file_name)) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % file_name)

    return failure_messages

def run(seq_record, options):
    "run hmmsearch against PFAM for all CDS features"
    if 'pfamdir' not in options:
        options.pfamdir = utils.get_full_path(__file__, '')

    query_sequence = utils.get_multifasta(seq_record)

    target_hmmfile = path.join(options.pfamdir, 'Pfam-A.hmm')

    logging.info('Running whole-genome pfam search')
    results = utils.run_hmmscan(target_hmmfile, query_sequence)

    _annotate(seq_record, options, results)

def _annotate(seq_record, options, results):
    "Annotate seq_record with CDS_motifs for the result"
    logging.debug("generating feature objects for PFAM hits")
    min_score = _min_score(options)
    max_evalue = _max_evalue(options)

    feature_by_id = utils.get_feature_dict(seq_record)
    
    for r in results:
        i = 1
        for hsp in r.hsps:
            if hsp.bitscore <= min_score or hsp.evalue >= max_evalue:
                continue

            if hsp.query_id not in feature_by_id:
                continue

            feature = feature_by_id[hsp.query_id]

            start, end = _calculate_start_end(feature, hsp)
            loc = FeatureLocation(start, end, strand=feature.strand)
            
            newFeature = SeqFeature(location=loc, type=options.FeatureTags.fullhmmer_tag)
            
            quals = defaultdict(list)
            
            quals['label'].append(r.id)
            if 'locus_tag' in feature.qualifiers:       
                quals['locus_tag'] = feature.qualifiers['locus_tag']
            else:
                quals['locus_tag'] = [hsp.query_id]
            quals['domain'] = [hsp.hit_id]
            quals['asDomain_id'] = ['fullhmmer_'+'_'.join(quals['locus_tag'])+'_'+'{:04d}'.format(i)]
            i += 1
            
            quals['evalue'] = [str("{:.2E}".format(float(hsp.evalue)))]
            quals['score'] = [str(hsp.bitscore)]
            quals['aSTool'] = ["fullhmmer"]
            quals['detection'] = ["hmmscan"]
            quals['database'] = [path.basename(r.target)]
            if 'transl_table' in feature.qualifiers:
                [transl_table] = feature.qualifiers['transl_table']
            else:
                transl_table = 1
            quals['translation'] = [str(newFeature.extract(seq_record.seq).translate(table=transl_table))]

            quals['note'].append("%s-Hit: %s. Score: %s. E-value: %s. Domain range: %s..%s." % \
                    (path.basename(r.target), hsp.hit_id, hsp.bitscore, hsp.evalue,
                     hsp.hit_start, hsp.hit_end))

            quals['description'] = [hsp.hit_description]

            try:
                pfamid = name_to_pfamid[hsp.hit_id]
                if 'db_xref' in quals:
                    quals['db_xref'].append("PFAM: %s" % pfamid)
                else:
                    quals['db_xref'] = ["PFAM: %s" % pfamid]    
            except KeyError:
                pass
            
            newFeature.qualifiers=quals
            seq_record.features.append(newFeature)


def _calculate_start_end(feature, result):
    "Calculate start and end of a result"
    if feature.strand == 1:
        start = feature.location.start + (3 * result.query_start)
        end = feature.location.start + (3 * result.query_end)
    else:
        end = feature.location.end - (3 * result.query_start)
        start = feature.location.end - (3 * result.query_end)

    return start, end

def _min_score(options):
    "Load the minimal score from the configuration or set a default"
    try:
        score = float(options.fullhmmer.score)
    except:
        score = 0

    return score

def _max_evalue(options):
    "Load the maximal evalue from the configuration or set a default"
    try:
        evalue = float(options.fullhmmer.evalue)
    except:
        evalue = 0.01

    return evalue
