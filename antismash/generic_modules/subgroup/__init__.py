# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2024 Ziqiang Luo
# Wageningen University & Research, NL
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""collect families proteins sequence from the clusters and run hmmsearch to identify the subgroups using set models"""

import os
from helperlibs.wrappers.io import TemporaryDirectory
from antismash import utils
from antismash.config import get_config
import logging
import re
import subprocess

defult_folder = os.path.dirname(os.path.abspath(__file__))


def load_subgroup_model(input_dir):
    """check whether needed subgroup models from the Subgroup_Model.csv file are available in family folders;
    if is, return a dictionary with family names as keys and the paths to the subgroup models as values."""

    subgroup_model_dic = {}
    missing_subgroup_model = False

    csv_path = os.path.join(input_dir, 'Subgroup_Model.csv')
    if not os.path.exists(csv_path):
        utils.log_status("Cannot find the required Subgroup_Model.csv file. Skip the subgroup identification.")
        logging.debug("Cannot find the required Subgroup_Model.csv file. Skip the subgroup identification.")
        return False

    with open(csv_path, 'r') as csv:
        next(csv)  # Skip the header line
        for line in csv:
            columns = line.strip().split('\t')
            if columns[0] == 'Y':  # Check if "Enable" is "Y"

                # subgroup models should be in their family folders indicated in the Subgroup_Model.csv
                family_model = columns[2][:-4]  # Remove the last 4 characters (.hmm)
                subgroup_model = columns[3]
                subgroup_model_path = os.path.join(os.path.join(input_dir, family_model), subgroup_model)
                if os.path.exists(subgroup_model_path) is False:
                    utils.log_status("Cannot find the required subgroup_model %s" % subgroup_model_path)
                    logging.debug("Cannot find the required subgroup_model %s" % subgroup_model_path)
                    missing_subgroup_model = True
                elif missing_subgroup_model is False:
                    if family_model in subgroup_model_dic:
                        subgroup_model_dic[family_model].append(subgroup_model_path)
                    else:
                        subgroup_model_dic[family_model] = [subgroup_model_path]

    if subgroup_model_dic == {} or missing_subgroup_model:
        utils.log_status("No enabled (Y) in Subgroup_Model.csv or missing the required missing subgroup model(s). "
                         "Skip the subgroup identification.")
        logging.debug("No enabled (Y) in Subgroup_Model.csv or missing the required missing subgroup model(s). "
                        "Skip the subgroup identification.")
        return False
    else:
        utils.log_status("Subgroup_Model.csv required Subgroup Models are loaded")
        logging.debug("Subgroup_Model.csv required Subgroup Models are loaded")
        return subgroup_model_dic


def fam_aa_fas(seq_records, subgroup_model_dic, output_dir):
    """Extract protein sequences of the genes in clusters that also belong to the families of the required Subgroup
    Models; Output the sequences in fasta files for each family in the output_dir. Return the family names."""

    fam_sub_model_set = set(list(subgroup_model_dic.keys()))
    fam_id_fas_dic = {}

    for seq_record in seq_records:
        features = utils.get_withincluster_cds_features(seq_record)

        # each feature contains a gene CDS information same as it in a gbk file.
        for feature in features:
            domain_list = []
            if 'sec_met' in feature.qualifiers:
                if feature.qualifiers['sec_met'][1].startswith("Domains detected:"):
                    item_without_prefix = feature.qualifiers['sec_met'][1].replace("Domains detected: ", "", 1)
                    raw_domain_list = item_without_prefix.split(";")
                    for raw_domain in raw_domain_list:
                        domain = raw_domain.split(" (E-value")

                        # because of filterhmmdetails.txt, in results only shows UDPGT_2 when also found UDPGT.
                        if domain[0] == 'UDPGT_2':
                            domain_list.append('UDPGT')
                        else:
                            domain_list.append(domain[0])

            # get the common families between the domain_list of a gene and subgroup_model_dic.keys()
            common_fam = fam_sub_model_set & set(domain_list)
            if common_fam:
                id_fas = (utils.get_gene_id(feature), feature.qualifiers['translation'][0])
                for fam_name in common_fam:
                    # add the gene id and protein sequence to the corresponding family in fam_id_fas_dic
                    if fam_name in fam_id_fas_dic:
                        fam_id_fas_dic[fam_name].append(id_fas)
                    else:
                        fam_id_fas_dic[fam_name] = [id_fas]

    # use fam_id_fas_dic to output collected protein sequences in fasta files for each family
    if fam_id_fas_dic:
        fam_names = list(fam_id_fas_dic.keys())
        utils.log_status(str(fam_names) + " domain(s) found in clusters and will be used for subgroup identification.")
        logging.debug(str(fam_names) + " domain(s) found in clusters and will be used for subgroup identification.")
        fam_fas_paths = {}
        for fam_name in fam_names:
            fam_fas_path = os.path.join(output_dir, "{}(within_clusters).fasta".format(fam_name))
            with open(fam_fas_path, "w") as fam_fas:
                fam_id_fas = fam_id_fas_dic[fam_name]
                for id_fas in fam_id_fas:
                    fam_fas.write('>' + id_fas[0] + '\n')
                    fam_fas.write(id_fas[1] + '\n')
            fam_fas_paths[fam_name] = fam_fas_path
        return fam_names, fam_fas_paths
    else:
        utils.log_status("No required family domains found in any cluster. Skip the subgroup identification.")
        logging.debug("No required family domains found in any cluster. Skip the subgroup identification.")
        return False


def merge_subgroup_model(subgroup_model_dic, fam_name, temp_dir):
    """ In temp_dir, marge the subgroup models of a given family into one .hmm file and hmmpress it.
    Return its path"""

    subgroup_model_paths = subgroup_model_dic[fam_name]
    merged_hmm_path = os.path.join(temp_dir, "{}_merged.hmm".format(fam_name))

    # in web version, the merged hmm file not needed to remove
    with open(merged_hmm_path, 'w') as merged_hmm:
        subprocess.call(['cat'] + subgroup_model_paths, stdout=merged_hmm)
    utils.run_hmmpress(merged_hmm_path)
    utils.log_status("required Subgroup Model for {} family is merged and hmm pressed".format(fam_name))
    logging.debug("required Subgroup Model for {} family is merged and hmm pressed".format(fam_name))
    return merged_hmm_path


def _update_feature(seq_records, hmmscan_result_dict):
    """add subgroup feature in seq_records for genes used in hmmscan for subgroup_identification."""

    for seq_record in seq_records:
        features = utils.get_withincluster_cds_features(seq_record)
        for feature in features:
            if utils.get_gene_id(feature) in hmmscan_result_dict:
                hits = hmmscan_result_dict[utils.get_gene_id(feature)]
                result = "; ".join(["%s (E-value: %s, bitscore: %s)" % (hit[1], hit[2], hit[3]) for hit in hits])
                feature.qualifiers['subgroup'] = [result]


def subgroup_identification(seq_records, output_dir, input_dir):
    """ Identify subgroups within the clusters using the set models in Subgroup_Model.csv."""

    if input_dir:
        utils.log_status("use {} as the input subgroup model directory".format(input_dir))
        logging.debug("use {} as the input subgroup model directory".format(input_dir))
    else:
        input_dir = defult_folder

    # create a subgroup result folder in the output_dir
    output_dir = os.path.join(output_dir, "subgroup")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # load the subgroup models from Subgroup_Model.csv
    subgroup_model_dic = load_subgroup_model(input_dir)
    if subgroup_model_dic:
        hmmscan_result_dict = {}
        fam_names, fam_fas_paths = fam_aa_fas(seq_records, subgroup_model_dic, output_dir)

        # for each family found in clusters, merge its subgroup models and run hmmscan to identify subgroups
        for fam_name in fam_names:
            result_lines = []
            utils.log_status("detecting {} family".format(fam_name))
            logging.debug("detecting {} family".format(fam_name))

            with TemporaryDirectory() as temp_dir:
                merged_hmm_path = merge_subgroup_model(subgroup_model_dic, fam_name, temp_dir)
                hmm_result = utils.run_hmmscan(merged_hmm_path, fam_fas_paths[fam_name], query_sequence_path=True)

            # parse the hmmscan result
            for gene_hits in hmm_result:
                gene_hits_unsort = []
                if gene_hits.hsps:
                    for hsp in gene_hits.hsps:
                        result_lines.append(str(hsp))   # for saving hmmscan raw results in a file

                        # only keep 4 attributes of a hit in gene_hits_unsort
                        gene_hits_unsort.append((hsp.query_id, hsp.hit_id, hsp.evalue, hsp.bitscore))

                    # Only keep the top 3 bitscore of hits in hmmscan_result_dict for updating the seq_records
                    gene_hits_sorted = sorted(gene_hits_unsort, key=lambda x: x[3], reverse=True)
                    hmmscan_result_dict[gene_hits_sorted[0][0]] = gene_hits_sorted[:3]

            # save hmmscan raw results in a file
            if result_lines:
                result_path = os.path.join(output_dir, "{}_subgroup_hmm_result(within_clusters).txt".format(fam_name))
                with open(result_path, "w") as file:
                    for hsp in result_lines:
                        file.write(hsp + '\n')

        # update the seq_records with the subgroup information
        if hmmscan_result_dict:
            utils.log_status("subfamilies detected, adding results in the gbk and output")
            logging.debug("subfamilies detected, adding results in the gbk and output")
            _update_feature(seq_records, hmmscan_result_dict)
        else:
            utils.log_status("no subfamily models does not match genes in the cluster")
            logging.debug("no subfamily models does not match genes in the cluster")