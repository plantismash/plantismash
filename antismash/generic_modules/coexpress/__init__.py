# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2016 Satria A Kautsar, Hernando G Suarez
# Wageningen University & Research, NL
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Compare detected clusters against prepared (GEO) Expression data"""


import logging
import json
import re
import numpy as np
from os import path
from antismash import utils
from antismash.config import get_config
from Bio import SeqIO
from Bio.Seq import translate
from helperlibs.wrappers.io import TemporaryDirectory

import networkx as nx
import community

name = "coexpress"
short_description = name.capitalize()
priority = 10000

_required_binaries = [
]


def check_prereqs(options):
    "Check if all required applications are around"
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate file: %r" % binary_name)
    return failure_messages


def run_coexpress(seq_record, all_gene_expressions, geo):
    options = get_config()
    cl_count = 1
    cl_list = utils.get_cluster_features(seq_record)

    gene_expressions = all_gene_expressions[seq_record.id]

    logging.info('Running CoExpress analysis on the clusters..')
    for cluster in cl_list:
        logging.debug('Running CoExpress analysis on record "%s".. (Cluster %s of %s)' % (geo["info"]["id"], cl_count, len(cl_list)))
        features = utils.get_cluster_cds_features(cluster, seq_record)
        cl_count += 1
        cluster_genes = {}

        for feature in features:
            gene_id =  utils.get_gene_id(feature)
            if gene_id in gene_expressions:
                cluster_genes[gene_id] = gene_expressions[gene_id]

        #calculate correlation value between genes
        for gene_1 in cluster_genes:
            if "cor" not in cluster_genes[gene_1]:
                cluster_genes[gene_1]["cor"] = {}
            if "exp" not in cluster_genes[gene_1]:
                continue
            for gene_2 in cluster_genes:
                if "cor" not in cluster_genes[gene_2]:
                    cluster_genes[gene_2]["cor"] = {}
                if gene_2 == gene_1:
                    continue
                if "exp" not in cluster_genes[gene_2]:
                    continue
                if gene_1 in cluster_genes[gene_2]["cor"]:
                    continue
                cor_val = calc_correlation_value(cluster_genes[gene_1], cluster_genes[gene_2])
                cluster_genes[gene_1]["cor"][gene_2] = cor_val
                cluster_genes[gene_2]["cor"][gene_1] = cor_val

        #calculate distance value for building dendogram
        for gene_1 in cluster_genes:
            if "dist" not in cluster_genes[gene_1]:
                cluster_genes[gene_1]["dist"] = {}
            for gene_2 in cluster_genes:
                if "dist" not in cluster_genes[gene_2]:
                    cluster_genes[gene_2]["dist"] = {}
                dist = 100.0
                if "cor" in cluster_genes[gene_1] and gene_2 in cluster_genes[gene_1]["cor"]:
                    cor_val = min(1.00, cluster_genes[gene_1]["cor"][gene_2])
                    dist = 100.0 * (1.0 - cor_val)
                cluster_genes[gene_1]["dist"][gene_2] = dist
                cluster_genes[gene_2]["dist"][gene_1] = dist

        # check for remote genes, add if correlation value >= 0.9
        for gene_1 in cluster_genes:
            for seqid in all_gene_expressions:
                prefix = "%s:" % seqid.replace(":", "_")
                for gene_2 in all_gene_expressions[seqid]:
                    if (prefix + gene_2) not in options.hmm_results: # only add biosynthetic remote genes
                        continue
                    if gene_2 == gene_1:
                        continue
                    if gene_2 in cluster_genes:
                        continue
                    cor_val = min(1.00, calc_correlation_value(cluster_genes[gene_1], all_gene_expressions[seqid][gene_2]))
                    if 1.00 > cor_val >= 0.9:
                        cluster_genes[gene_1]["dist"][gene_2] = 100.0 * (1.0 - cor_val)

        # review the remote genes, discard genes with less than 2 edges
        if True:
            edges_count = {}
            for gene_1 in cluster_genes:
                for gene_2 in cluster_genes[gene_1]["dist"]:
                    if gene_2 not in cluster_genes:
                        if gene_2 not in edges_count:
                            edges_count[gene_2] = 0
                        edges_count[gene_2] += 1
            for gene_1 in cluster_genes:
                new_dists = {}
                for gene_2 in cluster_genes[gene_1]["dist"]:
                    if (gene_2 in cluster_genes) or (edges_count[gene_2] >= 2):
                        new_dists[gene_2] = cluster_genes[gene_1]["dist"][gene_2]
                cluster_genes[gene_1]["dist"] = new_dists

        # review the remote genes, discard genes without any connection to cluster's biosynthetic genes
        if True:
            have_connections = []
            prefix = "%s:" % seq_record.id.replace(":", "_")
            for gene_1 in cluster_genes:
                if (prefix + gene_1) in options.hmm_results:
                    for gene_2 in cluster_genes[gene_1]["dist"]:
                        if (gene_2 not in cluster_genes) and (gene_2 not in have_connections):
                            have_connections.append(gene_2)
            for gene_1 in cluster_genes:
                new_dists = {}
                for gene_2 in cluster_genes[gene_1]["dist"]:
                    if (gene_2 in cluster_genes) or (gene_2 in have_connections):
                        new_dists[gene_2] = cluster_genes[gene_1]["dist"][gene_2]
                cluster_genes[gene_1]["dist"] = new_dists

        #update seq_record
        update_features(features, cluster_genes, geo)

    if False: #This feature is temporarily disabled, saved for next version #options.coexpress_signal_cluster_size < len(overlaps):
        logging.info('Running expression signal analysis on seq_record..')
        signals = [];
        n = options.coexpress_signal_cluster_size - 1
        #build list of cluster locations (for annotating signal regions)
        clrefs = []
        for cluster in cl_list:
            clrefs.append(((cluster.location.start, cluster.location.end), utils.get_cluster_number(cluster)))
        clrefs = sorted(clrefs, key=lambda cl: cl[0][0])
        #build signals
        for i in range(0, len(overlaps) - n):
            genes = []
            for overlap in overlaps[i:i+n]:
                gene = overlap[0]
                for feature in overlap:
                    if utils.get_gene_id(feature) in gene_expressions:
                        gene = feature
                        break
                genes.append(gene)
            cors = []
            checked = []
            hits = []
            for x in range(0, len(genes)):
                gene_x = utils.get_gene_id(genes[x])
                if prefix + gene_x in options.hmm_results:
                    hits.append(options.hmm_results[prefix + gene_x][0].query_id)
                for y in range(0, len(genes)):
                    if ((x,y) in checked) or ((y,x) in checked):
                        continue
                    cor_val = 0
                    gene_y = utils.get_gene_id(genes[y])
                    if (gene_x in gene_expressions) and (gene_y in gene_expressions):
                        cor_val = calc_correlation_value(gene_expressions[gene_x], gene_expressions[gene_y])
                    cors.append(cor_val)
                    checked.append((x, y))
            sloc = (genes[0].location.start + genes[-1].location.end) / 2
            cor_val = 0
            if len(cors) > 0 and len(list(set(hits))) > 1:
                cor_val = np.median(cors)
            cl_idx = -1
            for clref in clrefs:
                if sloc < clref[0][0]:
                    continue
                if sloc <= clref[0][1]:
                    cl_idx = clref[1]
                    break
            signals.append((sloc, cor_val, cl_idx))
        if "coexpress_signal" not in options:
            options.coexpress_signal = {}
        if geo["info"]["id"] not in options.coexpress_signal:
            options.coexpress_signal[geo["info"]["id"]] = {}
        options.coexpress_signal[geo["info"]["id"]][seq_record.id] = signals


def update_dist_between_clusters(seq_records, all_gene_expressions, geo):
    """Check and add remote genes that have > 0.9 PCC and in a cluster"""
    cluster_genes = {}
    for seq_record in seq_records:
        gene_expressions = all_gene_expressions[seq_record.id]
        for feature in utils.get_withincluster_cds_features(seq_record):
            gene_id =  utils.get_gene_id(feature)
            if gene_id in gene_expressions:
                cluster_genes[gene_id] = gene_expressions[gene_id]

    for gene_1 in cluster_genes:
        for gene_2 in cluster_genes:
            if gene_2 == gene_1:
                continue
            if (gene_1 in cluster_genes[gene_2]["dist"]) or (gene_2 in cluster_genes[gene_1]["dist"]):
                if (gene_1 not in cluster_genes[gene_2]["dist"]):
                    cluster_genes[gene_2]["dist"][gene_1] = cluster_genes[gene_1]["dist"][gene_2]
                if (gene_2 not in cluster_genes[gene_1]["dist"]):
                    cluster_genes[gene_1]["dist"][gene_2] = cluster_genes[gene_2]["dist"][gene_1]
                continue
            cor_val = min(1.00, calc_correlation_value(cluster_genes[gene_1], cluster_genes[gene_2]))
            if 1.00 > cor_val >= 0.9:
                cluster_genes[gene_1]["dist"][gene_2] = 100.0 * (1.0 - cor_val)
                cluster_genes[gene_2]["dist"][gene_1] = 100.0 * (1.0 - cor_val)

    for seq_record in seq_records:
        update_features(utils.get_withincluster_cds_features(seq_record), cluster_genes, geo)


def update_features(features, cluster_genes, geo):
    for feature in features:
        gene_id = utils.get_gene_id(feature)
        if gene_id in cluster_genes and "ref" in cluster_genes[gene_id]:
            cg = cluster_genes[gene_id]
            removed_lines = []
            if "geo" in feature.qualifiers:
                for line in feature.qualifiers["geo"]:
                    if line.startswith("Record: %s;" % geo["info"]["id"]):
                        removed_lines.append(line)
                for line in removed_lines:
                    feature.qualifiers["geo"].remove(line)
            else:
                feature.qualifiers["geo"] = []
            geo_str = "Record: %s;%s(%s);%s;%s" % (
                    geo["info"]["id"],
                    cg["ref"],
                    cg["evalue"],
                    ",".join(["%s=%s(%.2f)" % (sample, 0, cg["exp"][sample]) for sample in cg["exp"]]),
                    ",".join(["%s=%.2f" % (gene_2, cg["dist"][gene_2]) for gene_2 in cg["dist"]])
                    )
            feature.qualifiers["geo"].append(geo_str)


def parse_geofiles(geo_paths):
    result = []
    for geo_path in geo_paths:
        result.append(parse_geofile(geo_path))
    return result


def parse_geofile(geo_path):
    temp_data = {}
    dataset_info = {"id": "", "title": "", "samples": [], "desc": ""}
    dataset_data = {}

    # parse geo file to get dataset_info and temp_data
    with open(geo_path) as geo_file:
        pstate = ""
        sstate = None
        temp_table_header = []
        for line in iter(geo_file.readline, ''):
            line = line.lstrip().rstrip()
            if line.startswith("^"):
                # is a header line
                ss = line.split("=")
                if (len(ss) < 2):
                    ss.append("")
                left = ss[0][1:].rstrip().upper()
                right = ss[1].lstrip()
                pstate = left
                sstate = None
                temp_table_header = []
                if left in ["DATASET", "SERIES"]:
                    dataset_info["id"] = right
                    if left == "DATASET":
                        dataset_info["type"] = "GDS"
                    if left == "SERIES":
                        dataset_info["type"] = "GSE"
                elif left == "SAMPLE":
                    sample_id = get_safe_sample_name(right.upper())
                    if sample_id not in dataset_info["samples"]:
                        dataset_info["samples"].append(sample_id)
                    sstate = sample_id
            elif line.startswith("!"):
                # is an attribute line
                ss = line.split("=")
                if (len(ss) < 2):
                    ss.append("")
                left = ss[0][1:].rstrip().upper()
                right = ss[1].lstrip()
                temp_table_header = []
                if (pstate == "DATASET"):
                    if left == "DATASET_TITLE":
                        dataset_info["title"] = right
                    elif left == "DATASET_DESCRIPTION":
                        dataset_info["desc"] = right
                    elif left == "DATASET_TABLE_BEGIN":
                        pstate = "parsing_table"
                elif (pstate == "SERIES"):
                    if left == "SERIES_TITLE":
                        dataset_info["title"] = right
                    elif left == "SERIES_SUMMARY":
                        dataset_info["desc"] = right
                elif (pstate == "PLATFORM"):
                    if left == "PLATFORM_TABLE_BEGIN":
                        pstate = "parsing_table"
                elif (pstate == "SAMPLE"):
                    if left == "SAMPLE_TABLE_BEGIN":
                        pstate = "parsing_table"
                elif (pstate == "SUBSET"):
                    if left == "SUBSET_SAMPLE_ID":
                        for sample_id in right.split(","):
                            sample_id = get_safe_sample_name(sample_id.upper())
                            if sample_id not in dataset_info["samples"]:
                                dataset_info["samples"].append(sample_id)
            elif line.startswith("#"):
                continue
            else:
                if pstate == "parsing_table":
                    if len(temp_table_header) < 1:
                        # assign table header
                        temp_table_header = [ss.upper() for ss in line.split("\t")]
                    else:
                        # parse row
                        row = {}
                        i = 0
                        cols = [ss for ss in line.split("\t")]
                        for col in temp_table_header:
                            if i == len(cols):
                                cols.append(None)
                            if sstate != None and col == "VALUE":
                                row[sstate] = cols[i]
                            else:
                                row[col] = cols[i]
                            i += 1
                        if "ID_REF" in row:
                            if row["ID_REF"] not in temp_data:
                                temp_data[row["ID_REF"]] = {}
                            for col in row:
                                if col != "ID_REF":
                                    temp_data[row["ID_REF"]][col] = row[col]
                        elif "ID" in row:
                            if row["ID"] not in temp_data:
                                temp_data[row["ID"]] = {}
                            for col in row:
                                if col != "ID":
                                    temp_data[row["ID"]][col] = row[col]

    # process temp_data to dataset_data
    for id_ref in [id_ref for id_ref in temp_data]:
        identifiers = []
        values = {}
        for col in temp_data[id_ref]:
            val = temp_data[id_ref][col]
            if val != None:
                if get_safe_sample_name(col) in dataset_info["samples"]:
                    values[get_safe_sample_name(col)] = float(val)
                elif col in ["GB_ACC", "PT_ACC", "SP_ACC", "ORF", "GI", "GENBANK ACCESSION", "PLATFORM_CLONEID", "PLATFORM_ORF", "PLATFORM_SPOTID"]:
                    val = val.upper()
                    if len(val) > 0 and val not in identifiers:
                        identifiers.append(val)
        dataset_data[id_ref] = (identifiers, values)
        del temp_data[id_ref]

    return {"info" : dataset_info, "data" : dataset_data}


def calc_sample_ranges(geo_dataset):
    # calc ranges (10, 50, 90 percentile) of each sample sets
    geo_dataset["info"]["ranges"] = {}
    for sample in geo_dataset["info"]["samples"]:
        expressions = []
        for id_ref in geo_dataset["data"]:
            if sample in list(geo_dataset["data"][id_ref][1].keys()):
                expressions.append(geo_dataset["data"][id_ref][1][sample])
        expressions = np.array(expressions)
        geo_dataset["info"]["ranges"][sample] = (np.percentile(expressions, 10), np.percentile(expressions, 50), np.percentile(expressions, 90))


def parse_csvfiles(csv_paths):
    result = []
    for csv_path in csv_paths:
        result.append(parse_csv_file(csv_path))
    return result


def parse_csv_file(csv_path):
    header=[]
    dataset_data={}
    dataset_info = {"id": get_safe_sample_name(path.basename(csv_path).strip(".csv")), "title": "", "samples": [], "desc": "", "type": "CSV"}
    with open(csv_path) as f:
        skip_cols = []
        for line in f:
            if line.startswith('#'):
                if line.startswith('#title:'):
                    dataset_info['title']=line[8:].rstrip()
                if line.startswith('#desc:'):
                    dataset_info['desc']=line[7:].rstrip()
            else:
                Bline=re.split(',',line.rstrip())
                if not header:
                    header=Bline
                else:
                    cur_id=Bline[0]
                    dataset_data[cur_id]=[]
                    dataset_data[cur_id].append([cur_id])
                    dataset_data[cur_id].append({})
                    for i,exp in enumerate(Bline):
                        if (i != 0) and (i not in skip_cols):
                            if exp=='': exp=0
                            if exp=='NA': exp=float('nan')
                            else:
                                try:
                                    exp=float(exp)
                                except ValueError:
                                    logging.warning("Incorrectly formatted CSV file (%s)! column '%d: %s' contains non-numeric values and will be skipped." % (dataset_info["id"], i, header[i]))
                                    skip_cols.append(i)
                                    continue
                            dataset_data[cur_id][1][get_safe_sample_name(header[i])]=exp

        for i in skip_cols:
            skipped_col = get_safe_sample_name(header[i])
            for cur_id in dataset_data:
                if skipped_col in dataset_data[cur_id][1]:
                    del dataset_data[cur_id][1][skipped_col]
        temp_header = []
        for i in range(0, len(header)):
            if i not in skip_cols:
                temp_header.append(header[i])
        header = temp_header
    if len(header) > 1:
        dataset_info["samples"] = [get_safe_sample_name(h) for h in header[1:]]
    return {"info" : dataset_info, "data" : dataset_data}


def match_exp_to_genes(features, geo_dataset):
    cluster_genes = {}
    geo_info = geo_dataset["info"]
    geo_data = geo_dataset["data"]

    # get gene_id to ref_id table
    gene_to_ref = {}
    col_gene_id = geo_info["col_id"]
    if col_gene_id < 0:
        return {} # gene_id columns not found
    for id_ref, data in list(geo_data.items()):
        gene_to_ref[data[0][col_gene_id].upper()] = id_ref

    # fill cluster_genes
    for feature in features:
        gene_id = utils.get_gene_id(feature)
        if gene_id.upper() in gene_to_ref:
            cluster_genes[gene_id] = {}
            cluster_genes[gene_id]["ref"] = gene_to_ref[gene_id.upper()]
            cluster_genes[gene_id]["evalue"] = float(-1)

    #calculate scaled value for each hits
    for gene_id in cluster_genes:
        cg = cluster_genes[gene_id]
        if "ref" in cg:
            cg["exp"] = {}
            for sample, value in list(geo_data[cg["ref"]][1].items()):
                cg["exp"][sample] = value

    return cluster_genes


def calc_correlation_value(gene_1, gene_2):
    options = get_config()
    set_sample = []
    for s in gene_1["exp"]:
        if s not in set_sample:
            set_sample.append(s)
    for s in gene_2["exp"]:
        if s not in set_sample:
            set_sample.append(s)
    if len(set_sample) < 1:
        return 0
    #using Pearson correlation coefficient
    Xi = []
    Yi = []
    for s in set_sample:
        val_1 = None
        val_2 = None
        if s in gene_1["exp"]:
            val_1 = gene_1["exp"][s]
        if s in gene_2["exp"]:
            val_2 = gene_2["exp"][s]
        if (val_1 != None and val_2 != None):
            Xi.append(val_1)
            Yi.append(val_2)

    # MAD filter
    if options.coexpress_min_MAD == '': # Default case, geo_min_MAD is an empty string => filter out genes with 0 MAD.
        if (calc_MAD(Xi) == 0) or (calc_MAD(Yi) == 0):
            return 0
    else: # Filter out genes with MAD below specified min.
        if (calc_MAD(Xi) < float(options.coexpress_min_MAD)) or (calc_MAD(Yi) < float(options.coexpress_min_MAD)):
            return 0

    cor_val = (((len(set_sample) * sum([x**2 for x in Xi])) - (sum(Xi)**2)) ** 0.5)
    cor_val *= (((len(set_sample) * sum([y**2 for y in Yi])) - (sum(Yi)**2)) ** 0.5)
    if cor_val != 0:
        cor_val = ((len(set_sample) * sum([Xi[i]*Yi[i] for i in range(0, len(set_sample))])) - (sum(Xi) * sum(Yi))) / cor_val
    return cor_val


def calc_MAD(x):
    """ MAD = median( | Xi - median(X) |  ) """
    med=np.median(x)
    x2=np.absolute(x-med)
    MAD=np.median(x2)
    return MAD


def find_col_id(geo_dataset, seq_records):
    if geo_dataset["info"]["type"] == "CSV":
        geo_dataset["info"]["col_id"] = 0
        return geo_dataset
    for id_ref, data in list(geo_dataset["data"].items()):
        for i in range(0, len(data[0])):
            for seq_record in seq_records:
                for feature in utils.get_cds_features(seq_record):
                    gene_id = utils.get_gene_id(feature)
                    if gene_id.upper() == data[0][i].upper():
                        geo_dataset["info"]["col_id"] = i
                        return geo_dataset
    geo_dataset["info"]["col_id"] = -1
    return geo_dataset


def update_g(cur_gene1, interactions, distances, full_g):
    for cur_gene2 in interactions:
        # dist = 100 - (PCC * 100)
        # PCC = (100 - dist) / 100
        score = distances[cur_gene2]
        score = (100 - score) / 100

        #For good scores, add an edge to full_g
        if score >= 0.9:
            full_g.add_edge(cur_gene1, cur_gene2, weight = score)


def get_inter_cluster_relation(seq_records, geo_id):
    logging.debug('Calculating inter cluster relations on geo_record "%s"..' % (geo_id))
    data=[]
    full_g=nx.Graph()
    cluster_genes={}
    bio_genes=set()
    cur_cluster1=0
    # First, inspect all cluster to get cluster_genes
    for record in seq_records:
        for cluster in utils.get_cluster_features(record):
            cur_cluster1+=1
            cluster_genes[cur_cluster1]=set()

            for cluster_gene in utils.get_cluster_cds_features(cluster,record):
                # We only care about cluster_genes that have a geo match
                for cluster_gene_geo in utils.parse_geo_feature(cluster_gene):
                    # We only care about data from the current geo_id
                    if cluster_gene_geo['rec_id']==geo_id:
                        cur_gene1=utils.get_gene_id(cluster_gene)
                        cur_gene1_distances=cluster_gene_geo['dist']
                        cur_gene1_neighbors=set(cur_gene1_distances)

                        # Add each gene to cluster_genes, and to the full_g(raph) and to bio_genes
                        cluster_genes[cur_cluster1].add(cur_gene1)
                        full_g.add_node(cur_gene1)
                        if 'sec_met' in cluster_gene.qualifiers:
                            bio_genes.add(cur_gene1)

                        # Get intra-cluster edges
                        interactions=cur_gene1_neighbors.intersection(cluster_genes[cur_cluster1])
                        update_g(cur_gene1,interactions,cur_gene1_distances,full_g)

                        # From the second cluster onwards, we'll add inter-cluster edges backwards, i.e.: 2-1, 3-1, 3-2, 4-1, 4-2, etc...
                        if cur_cluster1 != 1:
                            for cur_cluster2 in cluster_genes:
                                if cur_cluster1 is not cur_cluster2:
                                    interactions=cur_gene1_neighbors.intersection(cluster_genes[cur_cluster2])
                                    update_g(cur_gene1,interactions,cur_gene1_distances,full_g)

    # Remove single nodes
    for node in full_g.nodes():
        if full_g.degree(node)==0:
            full_g.remove_node(node)

    # Get communities
    community_dict  =community.best_partition(full_g)

    number_of_clusters  =len(cluster_genes)

    # Now check inter-cluster interactions
    for i in range(1,number_of_clusters+1):
        cluster1    =cluster_genes[i]

        for j in range(i+1,number_of_clusters+1):
            cluster2    =cluster_genes[j]
            cluster3    =cluster1.union(cluster2)

            cluster_pair_g =full_g.subgraph(cluster3)

            communities_present =np.unique([community_dict[n] for n in cluster3 if n in community_dict])

            # CRITERIA 1 = only intra-community edges
            for cur_community in communities_present:
                cur_community_nodes =[n for n in cluster3 if n in community_dict and community_dict[n]==cur_community]
                cur_community_g =cluster_pair_g.subgraph(cur_community_nodes)

                decomposed_g=list(nx.connected_component_subgraphs(cur_community_g))
                for cur_g in decomposed_g:
                    # CRITERIA 2 = no isolates. anything with a clustering_coefficient=0 will be pruned out.
                    clustering_coefficient  =nx.clustering(cur_g)

                    pred_nodes  =[n for n in clustering_coefficient if clustering_coefficient[n]>0]
                    pred_g      =cur_g.subgraph(pred_nodes)
                    pred_edges  =pred_g.edges()

                    prediction  =set(pred_g.nodes())
                    prediction_cluster1 =prediction.intersection(cluster1)
                    prediction_cluster2 =prediction.intersection(cluster2)

                    bio_prediction  =prediction.intersection(bio_genes)
                    bio_prediction_cluster1 =prediction_cluster1.intersection(bio_genes)
                    bio_prediction_cluster2 =prediction_cluster2.intersection(bio_genes)

                    #CRITERIA 3 = at least 2 genes per cluster
                    #CRITERIA 5 = at least 1 bio per cluster
                    #CRITERIA 4 = at least 3 bio
                    if (    len(prediction_cluster1)        >=2 and
                            len(prediction_cluster2)        >=2 and
                            len(bio_prediction_cluster1)    >=1 and
                            len(bio_prediction_cluster2)    >=1 and
                            len(bio_prediction)             >=3
                            ):

                            pred_edges1    =[n for n in pred_edges if n[0] in cluster1 and n[1] in cluster1]
                            pred_edges2    =[n for n in pred_edges if n[0] in cluster2 and n[1] in cluster2]

                            pred_edges12    =[n for n in pred_edges if n[0] in cluster1 and n[1] in cluster2]
                            pred_edges21    =[n for n in pred_edges if n[0] in cluster2 and n[1] in cluster1]
                            inter_cluster_edges =pred_edges12 + pred_edges21

                            data.append({})
                            data[-1]['source']={}
                            data[-1]['source']['id']    =i
                            data[-1]['source']['links'] =pred_edges1

                            data[-1]['target']={}
                            data[-1]['target']['id']    =j
                            data[-1]['target']['links'] =pred_edges2

                            data[-1]['links']=inter_cluster_edges
    return data


def get_safe_sample_name(sample_name):
    safe_sample_name = re.sub("[^A-Za-z0-9]", "_", sample_name)
    return safe_sample_name
