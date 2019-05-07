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

"""General HMM detection module

"""
import logging
import pprint
import argparse
from os import path
from antismash import utils
from antismash import config
import os
import sys
from Bio.SeqFeature import SeqFeature, FeatureLocation
from cluster_predict import do_HMM as predict_clusters
from cluster_predict import A as trans_prob
from cluster_predict import B as emit_prob
from cluster_predict import start_p, index_dict
from numpy import mean

name = "clusterfinder"
short_description = name.capitalize()
priority = 2000

# The tuple is the name of the binary and whether it is an optional requirement
_required_binaries = [
    ]

class ClusterFinderHit(object):
    "A putative cluster identified by ClusterFinder"
    def __init__(self, positions, probability):
        self.location = FeatureLocation(positions[0], positions[1])
        self.probability = probability

def check_prereqs():
    """Don't check for prerequisites, we don't have any"""
    return []

def run_cluster_finder(seq_record, options):
    """Load the data and run the cluster_predict tool"""
    logging.info('Running ClusterFinder HMM to detect gene clusters')
    pfam_features = utils.get_pfam_features(seq_record)
    if len(pfam_features) == 0:
        return 0
    pfam_ids = [[xref for xref in feature.qualifiers['db_xref'] if "PFAM: " in xref][0].partition("PFAM: ")[2] for feature in pfam_features]
    if len(pfam_ids) == 0:
        return 0
    pfam_features = [feature for feature in pfam_features if len([xref for xref in feature.qualifiers['db_xref'] if "PFAM: " in xref]) > 0]

    options.clusterfinderfolder = options.full_outputfolder_path + os.sep + "clusterfinder"
    if not os.path.exists(options.clusterfinderfolder):
        os.mkdir(options.clusterfinderfolder)

    predictions = predict_clusters(pfam_ids, trans_prob, emit_prob,
                                   index_dict, start_p, graph=True,
                                   name='ctg_' + str(options.record_idx), FilterRepeats=0,
                                   path=options.clusterfinderfolder, outs=pfam_ids)

    if len(pfam_ids) != len(predictions):
        print "Warning: ClusterFinder input length %s != output length %s" % (len(pfam_ids), len(predictions))

    #Save ClusterFinder probabilities within CDS_motif features
    for idx in range(len(pfam_features)):
        feature = pfam_features[idx]
        prediction = str(predictions[idx])
        feature.qualifiers['note'].append("ClusterFinder probability: %s" % prediction)
    annotate_geneclusters(seq_record, options)
    return 0

def find_nr_cds(clusterpositions, seq_record):
    #Find the number of CDSs in candidate cluster and adjust the cluster starts and ends to match the CDS starts and ends
    cdsfeatures = utils.get_cds_features(seq_record)
    withinclustercdsfeatures = []
    for cds in cdsfeatures:
         if clusterpositions[0] <= int(cds.location.start) <= clusterpositions[1] or \
            clusterpositions[0] <= int(cds.location.end) <= clusterpositions[1] or \
            int(cds.location.start) <= clusterpositions[0] <= int(cds.location.end) or \
            int(cds.location.start) <= clusterpositions[1] <= int(cds.location.end):
            withinclustercdsfeatures.append(cds)
    if len(withinclustercdsfeatures) == 0:
        return clusterpositions, 0
    newclusterstart = min([int(cds.location.start) for cds in withinclustercdsfeatures])
    newclusterend = max([int(cds.location.end) for cds in withinclustercdsfeatures])
    newclusterpositions = [newclusterstart, newclusterend]
    return newclusterpositions, len(withinclustercdsfeatures)

def find_cf_clusters(pfam_features, seq_record, options):
    #Find clusters based on ClusterFinder probabilities of CDS_motif features
    cf_clusters = []
    state = "seed"
    clusterpositions = [0,0]
    pfam_ids = []
    loop_index = 1
    probabilities = [0]
    for feature in pfam_features:
        featurepositions = int(feature.location.start), int(feature.location.end)
        try:
            cf_probability = get_cf_probability(feature)
        except:
            loop_index += 1
            continue
        if cf_probability >= 0.3:
            if state == "seed":
                state = "extend"
                probabilities = [cf_probability]
                clusterpositions = [min(featurepositions), max(featurepositions)]
                pfam_ids = []
                pfam_ids.append([xref for xref in feature.qualifiers['db_xref'] if "PFAM: " in xref][0].partition("PFAM: ")[2])
            else:
                probabilities.append(cf_probability)
                if max(featurepositions) > clusterpositions[1]:
                    clusterpositions[1] = max(featurepositions)
                pfam_ids.append([xref for xref in feature.qualifiers['db_xref'] if "PFAM: " in xref][0].partition("PFAM: ")[2])
        else:
            if state == "extend":
                state = "seed"
                clusterpositions, cdsnr = find_nr_cds(clusterpositions, seq_record)
                if is_good_cluster_hit(cdsnr, probabilities, pfam_ids, options):
                    logging.info('Adding cluster ' + str(len(cf_clusters)+1) + \
                                 ' at position ' + str(clusterpositions[0]) + \
                                 ' to ' + str(clusterpositions[1]))
                    cf_clusters.append(ClusterFinderHit(clusterpositions, mean(probabilities)))
                clusterpositions = []
                pfam_ids = []
        if loop_index == len(pfam_features):
            if len(clusterpositions) > 0:
                clusterpositions, cdsnr = find_nr_cds(clusterpositions, seq_record)
                if is_good_cluster_hit(cdsnr, probabilities, pfam_ids, options):
                    logging.info('Adding cluster ' + str(len(cf_clusters)+1) + \
                                 ' at position ' + str(clusterpositions[0]) + \
                                 ' to ' + str(clusterpositions[1]))
                    cf_clusters.append(ClusterFinderHit(clusterpositions, mean(probabilities)))
            clusterpositions = []
            pfam_ids = []
        loop_index += 1
    return cf_clusters

def annotate_geneclusters(seq_record, options):
    """Re-annotate gene clusters in the seq_record"""
    pfam_features = utils.get_pfam_features(seq_record)
    cf_clusters = find_cf_clusters(pfam_features, seq_record, options)
    #Integrate ClusterFinder clusters with existing cluster features
    newclusters = []
    cluster_features = utils.get_cluster_features(seq_record)
    for cf_cluster in cf_clusters:
        overlaps = False
        cf_type = "cf_putative"
        for cluster in cluster_features:
            cluster_sig_genes = [gene for gene in utils.get_secmet_cds_features(seq_record) if gene in utils.get_cluster_cds_features(cluster, seq_record)]
            if utils.features_overlap(cf_cluster, cluster):
                overlaps = True
                if options.borderpredict: #Predict gene cluster borders using ClusterFinder
                    if ((cluster.location.end + cluster.location.start) / 2) in cf_cluster.location:
                        cluster.location = cf_cluster.location
                        for sig_gene in cluster_sig_genes:
                            startpoint = min([sig_gene.location.start, sig_gene.location.end])
                            endpoint = max([sig_gene.location.start, sig_gene.location.end])
                            if cluster.location.start > startpoint:
                                cluster.location = FeatureLocation(startpoint, cluster.location.end)
                            if cluster.location.end < endpoint:
                                cluster.location = FeatureLocation(cluster.location.start, endpoint)
                elif cf_cluster.location.start < cluster.location.start and cf_cluster.location.end > cluster.location.end:
                    cluster.location = cf_cluster.location
                elif cf_cluster.location.start < cluster.location.start:
                    cluster.location = FeatureLocation(cf_cluster.location.start, cluster.location.end)
                elif cf_cluster.location.end > cluster.location.end:
                    cluster.location = FeatureLocation(cluster.location.start, cf_cluster.location.end)
                cluster.qualifiers['probability'] = [ "%01.4f" % cf_cluster.probability ]
        if not overlaps:
            cf_cluster_CDSs = utils.get_cluster_cds_features(cf_cluster, seq_record)
            for CDS in cf_cluster_CDSs:
                if 'sec_met' in CDS.qualifiers:
                    type_sec_met_qualifiers = [feat for feat in CDS.qualifiers['sec_met'] if "Type: " in feat]
                    for qualifier in type_sec_met_qualifiers:
                        if "cf_fatty_acid" in qualifier:
                            if cf_type == "cf_putative":
                                cf_type = "cf_fatty_acid"
                            elif cf_type == "cf_saccharide":
                                cf_type = "cf_fatty_acid-saccharide"
                        if "cf_saccharide" in qualifier:
                            if cf_type == "cf_putative":
                                cf_type = "cf_saccharide"
                            elif cf_type == "cf_fatty_acid":
                                cf_type = "cf_fatty_acid-saccharide"
            new_cluster = SeqFeature(cf_cluster.location, type="cluster")
            new_cluster.qualifiers['product'] = [cf_type]
            new_cluster.qualifiers['probability'] = [ "%01.4f" % cf_cluster.probability ]
            newclusters.append(new_cluster)
    seq_record.features.extend(newclusters)
    #Re-number clusters
    clusters = utils.get_cluster_features(seq_record)
    clusters.sort(compare_feature_locations)
    clusternr = options.clusternr_offset
    for cluster in clusters:
        cluster.qualifiers['note'] = ["Cluster number: %s" % clusternr]
        clusternr += 1
    options.next_clusternr = clusternr

def compare_feature_locations(feature1, feature2):
    return cmp(feature1.location.start, feature2.location.start)

def get_cf_probability(feature):
    "Get ClusterFinder probability from feature"
    cf_probability_string = [note for note in feature.qualifiers['note'] if "ClusterFinder probability: " in note][0]
    return float(cf_probability_string.partition("ClusterFinder probability: ")[2])

def is_good_cluster_hit(cdsnr, probabilities, pfam_ids, options):
    "Check if the current cluster is a good hit"
    biosynthetic_pfams = ["PF00109", "PF02801", "PF08659", "PF00378", "PF08541", "PF08545", "PF02803", "PF00108",
    "PF02706", "PF03364", "PF08990", "PF00501", "PF00668", "PF08415", "PF00975", "PF03061", "PF00432", "PF00494", "PF03936",
    "PF01397", "PF00432", "PF04275", "PF00348", "PF02401", "PF04551", "PF00368", "PF00534", "PF00535", "PF02922", "PF01041",
    "PF00128", "PF00908", "PF02719", "PF04321", "PF01943", "PF02806", "PF02350", "PF02397", "PF04932", "PF01075", "PF00953",
    "PF01050", "PF03033", "PF01501", "PF05159", "PF04101", "PF02563", "PF08437", "PF02585", "PF01721", "PF02052", "PF02674",
    "PF03515", "PF04369", "PF08109", "PF08129", "PF09221", "PF09683", "PF10439", "PF11420", "PF11632", "PF11758", "PF12173",
    "PF04738", "PF04737", "PF04604", "PF05147", "PF08109", "PF08129", "PF08130", "PF00155", "PF00202", "PF00702", "PF06339",
    "PF04183", "PF10331", "PF03756", "PF00106", "PF01370", "PF00107", "PF08240", "PF00441", "PF02770", "PF02771", "PF08028",
    "PF01408", "PF02894", "PF00984", "PF00725", "PF03720", "PF03721", "PF07993", "PF02737", "PF00903", "PF00037", "PF04055",
    "PF00171", "PF00067", "PF01266", "PF01118", "PF02668", "PF00248", "PF01494", "PF01593", "PF03992", "PF00355", "PF01243",
    "PF00384", "PF01488", "PF00857", "PF04879", "PF08241", "PF08242", "PF00698", "PF00483", "PF00561", "PF00583", "PF01636",
    "PF01039", "PF00288", "PF00289", "PF02786", "PF01757", "PF02785", "PF02409", "PF01553", "PF02348", "PF00891", "PF01596",
    "PF04820", "PF02522", "PF08484", "PF08421"]
    npfams = len(set(biosynthetic_pfams) & set(pfam_ids))
    rc = cdsnr >= options.cdsnr and \
         mean(probabilities) >= options.cf_prob_thres and \
         npfams >= options.cf_npfams
    logging.info('rc = ' + str(rc) + \
                 ' cdsnr = ' + str(cdsnr) + \
                 ' mean(prob) = ' + str(mean(probabilities)) + \
                 ' npfams = ' + str(npfams))
    return rc;

__all__ = [ check_prereqs, run_cluster_finder]
