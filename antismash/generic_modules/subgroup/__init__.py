# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2024 Ziqiang Luo
# Wageningen University & Research, NL
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""collect families proteins sequence from the clusters and run hmmsearch and pplacer to identify the subgroups using
set models"""

import os
import logging
import subprocess
import shutil
from helperlibs.wrappers.io import TemporaryDirectory
from antismash import utils
from run_pplacer_tree import run_subtree_subgroup
from make_tree_svg import run_graphlan

default_folder = os.path.dirname(os.path.abspath(__file__))


def load_subgroup_dict(input_dir, file_name, output_dir):
    """check whether needed subgroup models from the Subgroup_Model.txt file are available in family folders;
    if is, return a dictionary with family names as keys and the paths to the subgroup models as values."""

    subgroup_dict = {}
    subgroup_substrate_dict = {}
    folder_family_dict = {}
    missing_subgroup_model = False

    txt_path = os.path.join(input_dir, file_name)
    if not os.path.exists(txt_path):
        logging.debug("Cannot find the required {} file. Skip the subgroup identification.".format(file_name))
        return False
    with open(txt_path, 'r') as txt:
        next(txt)  # Skip the header line
        for line in txt:
            columns = line.strip().split('\t')
            if columns[0] == 'Y':  # Check if "Enable" is "Y"

                # subgroup models should be in their family folders indicated in the Subgroup_Model.txt
                folder_name = columns[2]  # in other words, the family name in short
                family_models =  [model.strip() for model in columns[3].split("/")] # get the family model name
                subgroup_model = columns[4]
                subgroup_model_path = os.path.join(os.path.join(input_dir, folder_name), subgroup_model)
                if os.path.exists(subgroup_model_path) is False:
                    logging.debug("Cannot find the required %s" % subgroup_model_path)
                    missing_subgroup_model = True
                elif missing_subgroup_model is False:
                    if folder_name in subgroup_dict:
                        subgroup_dict[folder_name].append(subgroup_model_path)
                    else:
                        subgroup_dict[folder_name] = [subgroup_model_path]
                    # this allows a family corresponding to multiple subgroup models and vice versa
                    # todo: consider the case that a family member must have multiple models in the same time
                    for family_model in family_models:
                        if family_model not in folder_family_dict:
                            folder_family_dict[family_model] = [folder_name]
                        elif folder_name not in folder_family_dict[family_model]:
                            folder_family_dict[family_model].append(folder_name)

            # get the Match type of the subgroups
            if len(columns) == 7:
                subgroup_substrate_dict[columns[4][:-4]] = [substrate.strip() for substrate in columns[6].split(",")]

    if subgroup_dict == {} or missing_subgroup_model:
        logging.debug("No enabled (Y) in {} or missing the required missing subgroup model/tree package(s). "
                      "Skip the subgroup identification.".format(file_name))
        return [], {}
    else:
        logging.debug("{} required Subgroup Models/Tree package are loaded".format(file_name))
        shutil.copy(txt_path, os.path.join(output_dir, "0_"+file_name))
        return [folder_family_dict, subgroup_dict], subgroup_substrate_dict


def fam_aa_fas(seq_records, subgroup_model_dic, subgroup_tree_dic, output_dir, options):
    """Extract protein sequences of the genes in clusters that also belong to the families of the required Subgroup
    Models; Output the sequences in fasta files for each family in the output_dir. Return the family names."""

    fam_model = list(subgroup_model_dic[0].keys()) if subgroup_model_dic else []
    fam_tree = list(subgroup_tree_dic[0].keys()) if subgroup_tree_dic else []

    fam_sub_model_set = set(fam_model + fam_tree)
    fam_id_fas_dic = {}

    for seq_record in seq_records:
        features = utils.get_withincluster_cds_features(seq_record)
        # get genes features within clusters
        for feature in features:
            domain_list = []
            id = utils.get_gene_id(feature)

            # get the domain list of a gene from the hmm_results; so ignored the filterhmmdetails.txt
            # prefix_query = "%s:" % seq_record.id.replace(":", "_") + id
            # if prefix_query in options.hmm_results:
            #     for hsp in options.hmm_results[prefix_query]:
            #         domain_list.append(hsp.query_id.split("/")[1])

            if 'sec_met' in feature.qualifiers:
                if feature.qualifiers['sec_met'][1].startswith("Domains detected:"):
                    item_without_prefix = feature.qualifiers['sec_met'][1].replace("Domains detected: ", "", 1)
                    raw_domain_list = item_without_prefix.split(";")
                    for raw_domain in raw_domain_list:
                        domain = raw_domain.split(" (E-value")
                        domain_list.append(domain[0])

            # get the common families between the domain_list of a gene and subgroup_model_dic.keys()
            common_model = fam_sub_model_set & set(domain_list)

            if common_model:
                id_fas = (id, feature.qualifiers['translation'][0])
                for model_name in common_model:
                    for fam_name in subgroup_model_dic[0][model_name]:
                        # add the gene id and protein sequence to the corresponding family in fam_id_fas_dic
                        if fam_name not in fam_id_fas_dic:
                            fam_id_fas_dic[fam_name] = [id_fas]
                        elif id_fas not in fam_id_fas_dic[fam_name]:
                            fam_id_fas_dic[fam_name].append(id_fas)
    # use fam_id_fas_dic to output collected protein sequences in fasta files for each family
    if fam_id_fas_dic:
        fam_names = list(fam_id_fas_dic.keys())
        logging.debug(str(fam_names) + " domain(s) found in clusters and will be used for subgroup identification.")
        fam_fas_paths = {}
        for fam_name in fam_names:
            fam_fas_path = os.path.join(output_dir, "{}_within_clusters.fasta".format(fam_name))
            with open(fam_fas_path, "w") as fam_fas:
                fam_id_fas = fam_id_fas_dic[fam_name]
                for id_fas in fam_id_fas:
                    fam_fas.write('>plantismash_' + id_fas[0] + '\n') # add prefix to avoid repeated id in tree making
                    fam_fas.write(id_fas[1] + '\n')
            fam_fas_paths[fam_name] = fam_fas_path
        return fam_names, fam_fas_paths
    else:
        logging.debug("No required family domains found in any cluster. Skip the subgroup identification.")
        return [], []


def merge_subgroup_model(subgroup_model_dic, fam_name, temp_dir):
    """ In temp_dir, marge the subgroup models of a given family into one .hmm file and hmmpress it.
    Return its path"""

    if fam_name in subgroup_model_dic:
        subgroup_model_paths = subgroup_model_dic[fam_name]
    else:
        return None
    merged_hmm_path = os.path.join(temp_dir, "{}_merged.hmm".format(fam_name))

    # in web version, the merged hmm file not needed to remove
    with open(merged_hmm_path, 'w') as merged_hmm:
        subprocess.call(['cat'] + subgroup_model_paths, stdout=merged_hmm)
    utils.run_hmmpress(merged_hmm_path)
    logging.debug("required Subgroup Model for {} family is merged and hmm pressed".format(fam_name))
    return merged_hmm_path


def _update_feature(seq_records, hmmscan_result_dicts, tree_query_subgroups, output_dir, fam_names, subgroup_substrate_dict):
    """add subgroup feature in seq_records for genes used in hmmscan for subgroup_identification."""
    unsure_subgroups = []
    for seq_record in seq_records:
        features = utils.get_withincluster_cds_features(seq_record)
        for feature in features:
            query = utils.get_gene_id(feature)

            for fam_name in fam_names:

                hmmscan_result_dict = hmmscan_result_dicts[fam_name] if fam_name in hmmscan_result_dicts else None
                tree_query_subgroup = tree_query_subgroups[fam_name] if fam_name in tree_query_subgroups else None

                if hmmscan_result_dict and query in hmmscan_result_dict:
                    hits = hmmscan_result_dict[query]
                    result_hmm = ";".join(["%s (E=%s, bitscore=%s)" % (hit[1], hit[2], hit[3]) for hit in hits])
                    # considered the case that a gene has multiple domains from more one required families
                    if "subgroup_hmm" in feature.qualifiers:
                        feature.qualifiers['subgroup_hmm'].append(result_hmm)
                    else:
                        feature.qualifiers['subgroup_hmm'] = [result_hmm]
                else:
                    hits = result_hmm = None

                if tree_query_subgroup and query in tree_query_subgroup:
                    subgroup_likelihood = tree_query_subgroup[query]
                    subgroup_max = max(subgroup_likelihood.items(), key=lambda item: item[1])[0]  # get a tuple
                    len_likelihood = len(subgroup_likelihood)
                    len_subgroup_max = len(subgroup_max)
                    likelihood_max = subgroup_likelihood[subgroup_max]
                    if "subgroup_tree" in feature.qualifiers:
                        feature.qualifiers['subgroup_tree'].append(str(subgroup_likelihood))
                    else:
                        feature.qualifiers['subgroup_tree'] = [str(subgroup_likelihood)]
                else:
                    likelihood_max = len_likelihood = len_subgroup_max = subgroup_likelihood = subgroup_max = None

                # Only if the grouping results of the tree is same to hmm, the subgroup is sure
                if 'subgroup' not in feature.qualifiers and (hits or likelihood_max):
                    feature.qualifiers['subgroup'] = []

                if len_likelihood == 1 and len_subgroup_max == 1 and hits and hits[0][1] == subgroup_max[0]:
                    feature.qualifiers['subgroup'].append(hits[0][1])
                    if hits[0][1] in subgroup_substrate_dict:
                        # todo: Eliminate self-reference (e.g. A-Glycosidic_branch_elongating in saccharide) and
                        #  consider more conbinations
                        substrates_list = subgroup_substrate_dict[hits[0][1]]
                        if 'substrates' not in feature.qualifiers:
                            feature.qualifiers['substrates'] = substrates_list
                        else:
                            feature.qualifiers['substrates'] += substrates_list
                elif hits or likelihood_max:
                    unsure_subgroups.append([query, subgroup_likelihood, result_hmm, fam_name])
                    if hits:
                        feature.qualifiers['subgroup'].append(hits[0][1] + " (E: %s)" % hits[0][2])
                    elif likelihood_max > 0.8 and len_subgroup_max == 1:
                        feature.qualifiers['subgroup'].append(";".join(subgroup_max) + " (tree>0.8)")
                    elif subgroup_max:
                        feature.qualifiers['subgroup'].append("unsure tree")

    write_unsure_subgroup(unsure_subgroups, output_dir) if unsure_subgroups else None


def write_unsure_subgroup(unsure_subgroups, output_dir):
    """write the unsure subgroups results (trees and Hmm gave different results) in a file"""

    output_file = os.path.join(output_dir, "0_unsure_subgroup.txt")
    with open(output_file, "w") as file:
        for unsure_subgroup in unsure_subgroups:
            query, subgroup_likelihood, result_hmm, fam_name = unsure_subgroup
            file.write(query + " for " + fam_name + '\n')
            if subgroup_likelihood:
                subgroup_likelihood_list = sorted(subgroup_likelihood.items(), key=lambda x: x[1], reverse=True)
                file.write("subgroup_likelihood \t" + str(subgroup_likelihood_list) + '\n')
            file.write("result_hmm \t" + result_hmm + '\n') if result_hmm else None


def subgroup_identification(seq_records, output_dir, options):
    """ Identify subgroups within the clusters using the set models in Subgroup_Model.txt."""

    if options.subgroup_inputpath:
        input_dir = options.subgroup_inputpath
        logging.debug("use {} as the input subgroup model directory".format(input_dir))
    else:
        input_dir = default_folder

    # create a subgroup result folder in the output_dir
    output_dir = os.path.join(output_dir, "subgroup")
    os.mkdir(output_dir) if not os.path.exists(output_dir) else None

    # load the subgroup models from Subgroup_Model.txt
    subgroup_model_dic, subgroup_substrate_dict = load_subgroup_dict(input_dir, "Subgroup_Model.txt", output_dir)
    subgroup_tree_dic = load_subgroup_dict(input_dir, "Subgroup_Tree.txt", output_dir)[0]
    if subgroup_model_dic or subgroup_tree_dic:
        fam_names, fam_fas_paths = fam_aa_fas(seq_records, subgroup_model_dic, subgroup_tree_dic, output_dir, options)

        # test the module on QS UGT sequences which needed to put in  you/outputfolder/subgroup ahead
        # fam_fas_paths["UDPGT"] = os.path.join(output_dir, "sequence.fasta")
        # # Read the contents of the file
        # with open(fam_fas_paths["UDPGT"], 'r') as file:
        #     lines = file.readlines()

        # Open the file in write mode to overwrite its contents
        # with open(fam_fas_paths["UDPGT"], 'w') as file:
        #     for line in lines:
        #         # Replace '>' with '>plantismash_'
        #         modified_line = line.replace(">", ">plantismash_")
        #         file.write(modified_line)
        # fam_names = ["UDPGT"]

        hmmscan_result_dicts = {}
        tree_query_subgroups = {}

        # for each family found in clusters, merge its subgroup models and run hmmscan to identify subgroups
        for fam_name in fam_names:
            with TemporaryDirectory() as temp_dir:
                merged_hmm_path = merge_subgroup_model(subgroup_model_dic[1], fam_name, temp_dir)
                if not merged_hmm_path:
                    continue
                logging.debug("detecting {} family".format(fam_name))
                hmm_result = utils.run_hmmscan(merged_hmm_path, fam_fas_paths[fam_name], query_sequence_path=True)

            # parse the hmmscan result
            hmmscan_result_dict = {}
            result_lines = []
            for gene_hits in hmm_result:
                gene_hits_unsort = []
                if gene_hits.hsps:
                    for hsp in gene_hits.hsps:
                        # result_lines.append(str(hsp))  # for saving hmmscan raw results in a file
                        # only keep 4 attributes of a hit in gene_hits_unsort
                        # remove plantismash_ prefix in the query_id
                        gene_hits_unsort.append((hsp.query_id[12:], hsp.hit_id, hsp.evalue, hsp.bitscore))
                    # Only keep the top 3 bitscore of hits in hmmscan_result_dict for updating the seq_records
                    gene_hits_sorted = sorted(gene_hits_unsort, key=lambda x: x[3], reverse=True)
                    hmmscan_result_dict[gene_hits_sorted[0][0]] = gene_hits_sorted[:3]
                    result_lines.append(gene_hits_sorted[:3])
            hmmscan_result_dicts[fam_name] = hmmscan_result_dict if hmmscan_result_dict else None

            ## save hmmscan raw results in a file
            if result_lines:
                result_path = os.path.join(output_dir, "{}_top3_hmm.txt".format(fam_name))
                with open(result_path, "w") as file:
                    for hsp in result_lines:
                        file.write(str(hsp) + '\n')

        for fam_name in fam_names:
            # use pplacer to place the query sequences on the reference tree and get the subgroups and likelihoods
            if subgroup_tree_dic and fam_name in subgroup_tree_dic[1]:
                refpkg_dir = subgroup_tree_dic[1].get(fam_name)[0]
                if os.path.exists(refpkg_dir):
                    # todo fam_fas_paths[fam_name] should give the DNA sequences files if datatype is DNA in
                    #  phylo_model.json of refpkg_dir.
                    # todo the same fam_name with more than one refpkg_dir wanted to be used in the future

                    logging.debug("use pplacer to place the query protein sequences on {} family tree".format(fam_name))
                    tree_query_subgroup, prefix_singletrees, ref_group_file = run_subtree_subgroup(fam_name, refpkg_dir,
                                                                                    output_dir, fam_fas_paths[fam_name])

                    if not options.disable_treesvg:
                        tree_svg_dir = os.path.join(output_dir, "tree_svg")
                        os.mkdir(tree_svg_dir) if not os.path.exists(tree_svg_dir) else None
                        logging.debug("run graphlan to make tree.svg for each query")
                        for prefix_singletree in prefix_singletrees:
                            prefix, singletrees = prefix_singletree
                            for singletree in singletrees:
                                run_graphlan(ref_group_file, tree_svg_dir, prefix, singletree, fam_name)

                else:
                    logging.debug(
                        "no refpkg for {} family provide in input directory {}, skip".format(fam_name, refpkg_dir))
                    tree_query_subgroup = {}
                tree_query_subgroups[fam_name] = tree_query_subgroup

        # update the seq_records with the subgroup information
        if hmmscan_result_dicts:
            logging.debug("subfamilies detected, adding results in the gbk and output")
            _update_feature(seq_records, hmmscan_result_dicts, tree_query_subgroups,
                            output_dir, fam_names, subgroup_substrate_dict)
        else:
            logging.debug("no subfamily models does not match genes in the cluster")
