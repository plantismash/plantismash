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
import warnings
from helperlibs.bio import seqio
from os import path
from antismash import utils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation

name = "genbank"
short_description = "GenBank output"
priority = 1

def seq_record_convert_nucl_to_prot(seq_records, options):
    seq_record = seq_records[0]
    cdsfeatures = utils.get_cds_features(seq_record)
    cdsmotifs = utils.get_all_features_of_type(seq_record, ["CDS_motif"])
    #Find corresponding cdsmotifs for each cdsfeature
    cdsmotifdict = {}
    for cdsfeature in cdsfeatures:
        for cdsmotif in cdsmotifs:
            if cdsfeature.location.start <= cdsmotif.location.start <= cdsfeature.location.end:
                if not cdsmotifdict.has_key(cdsfeature.qualifiers['product'][0]):
                    cdsmotifdict[cdsfeature.qualifiers['product'][0]] = [cdsmotif]
                else:
                    cdsmotifdict[cdsfeature.qualifiers['product'][0]].append(cdsmotif)
    #For each cdsfeature, write a protein SeqRecord with CDS_motif features (abMotifs AND sec_met)
    prot_seq_records = []
    for cdsfeature in cdsfeatures:
        cds_domains = []
        #Extract sec_met info from feature
        if 'sec_met' in cdsfeature.qualifiers:
            if len([qual for qual in cdsfeature.qualifiers['sec_met'] if "NRPS/PKS subtype: " in qual]) > 0:
                cds_description = [qual for qual in cdsfeature.qualifiers['sec_met'] if "NRPS/PKS subtype: " in qual][0].partition("NRPS/PKS subtype: ")[2]
            else:
                cds_description = "Unknown protein"
            cds_domains = [qual for qual in cdsfeature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qual]
        else:
            cds_description = "Unknown protein"
        #Create protein seq_record
        prot_seq_record = SeqRecord(Seq(cdsfeature.qualifiers['translation'][0], IUPAC.protein),
                                    id=cdsfeature.qualifiers['product'][0], name=cdsfeature.qualifiers['product'][0],
                                    description=cds_description)
        utils.fix_record_name_id(prot_seq_record, options)
        #Add CDS_motif features based on NRPS/PKS domains
        cdsmotif_features = []
        for cds_domain in cds_domains:
            domainstart, domainend = cds_domain.partition(" (")[2].partition("). ")[0].split("-")
            domainlocation = FeatureLocation(int(domainstart), int(domainend))
            domain_feature = SeqFeature(domainlocation, type="CDS_motif")
            domain_feature.qualifiers['note'] = [cds_domain]
            cdsmotif_features.append(domain_feature)
        #Add CDS_motif features based on NRPS/PKS abMotifs
        if cdsmotifdict.has_key(cdsfeature.qualifiers['product'][0]):
            for cdsmotif in cdsmotifdict[cdsfeature.qualifiers['product'][0]]:
                oldstart, oldend = cdsmotif.location.start, cdsmotif.location.end
                newstart = (oldstart - cdsfeature.location.start) / 3
                newend = (oldend - cdsfeature.location.start) / 3
                newlocation = FeatureLocation(newstart, newend)
                cdsmotif.location = newlocation
                cdsmotif_features.append(cdsmotif)
        prot_seq_record.features.extend(cdsmotif_features)
        prot_seq_records.append(prot_seq_record)
    return prot_seq_records

def write(seq_records, options):
    basename = seq_records[0].id
    if options.input_type == 'nucl':
        output_name = path.join(options.outputfoldername, "%s.final.gbk" % basename)
        logging.debug("Writing %s seq_records to %r" % (len(seq_records), output_name))
        seqio.write(seq_records, output_name, 'genbank')
        i=1
        for rec in seq_records:
            # For compatibility with the database importer, we have to check whether we are dealing
            # with a seq_record obtained from a file (then its class will be SeqRecord) or from a#
            # database (then its class will be DBSeqRecord)
            #
            # running the cluster extraction on a DBSeqRecord will throw an exception, as splitting the object is not supported
            if rec.__class__.__name__ == 'SeqRecord':
                for cluster in utils.get_cluster_features(rec):
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        cluster_rec = rec[cluster.location.start:cluster.location.end]
                    output_name = path.join(options.outputfoldername,
                                            "%s.cluster%03d.gbk" % (basename, i))
                    seqio.write([cluster_rec], output_name, 'genbank')
                    i += 1
    else:
        seq_records = seq_record_convert_nucl_to_prot(seq_records, options)
        output_name = path.join(options.outputfoldername, "%s.final.gp" % basename)
        logging.debug("Writing seq_records to %r" % output_name)
        seqio.write(seq_records, output_name, 'genbank')
