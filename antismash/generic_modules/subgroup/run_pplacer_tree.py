# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2024 Ziqiang Luo
# Wageningen University & Research, NL
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""run pplacer to place query sequence on the reference tree and get the subgroup of the subtree of query sequence"""

import subprocess
import json
from ete3 import Tree
import os


def run_pplacer(fam_fasta, refpkg_dir, fam_name, output_dir):
    """run pplacer to place query protein sequences on the reference tree in refpkg_dir"""

    # get reference alignment and ref_subgroup_txt file path from refpkg_dir
    contents_file = os.path.join(refpkg_dir, "CONTENTS.json")
    with open(contents_file, 'r') as f:
        contents_dict = json.load(f)
        ref_alignment_path = os.path.join(refpkg_dir, contents_dict['files']["aln_fasta"])
        ref_group_file = os.path.join(refpkg_dir, contents_dict['files']["seq_info"])

    # output file paths for hmmbuild, hmmalign, pplacer and guppy
    hmm_path = os.path.join(output_dir, "{}.hmm".format(fam_name))
    align_sto = os.path.join(output_dir, "{}_query_ref_aligned.sto".format(fam_name))
    jplace_file = os.path.join(output_dir, "{}.jplace".format(fam_name))
    trees_file = os.path.join(output_dir, "{}_total.tre".format(fam_name))

    # run hmmalign to align query sequences to reference alignment hmm made by hmmbuild
    hmmbuild_command = ["hmmbuild", hmm_path, ref_alignment_path]
    print(hmmbuild_command)
    subprocess.call(' '.join(hmmbuild_command), shell=True)
    hmmalign_command = ["hmmalign", "-o", align_sto, "--mapali", ref_alignment_path, hmm_path, fam_fasta]
    subprocess.call(' '.join(hmmalign_command), shell=True)

    # run pplacer
    pplacer_command = ["pplacer", "-o ", jplace_file, "-c", refpkg_dir, align_sto]
    subprocess.call(' '.join(pplacer_command), shell=True)

    # run guppy
    subprocess.call("export LC_ALL=C", shell=True)
    guppy_command = ["guppy sing", "-o", trees_file, jplace_file]
    subprocess.call(' '.join(guppy_command), shell=True)

    return jplace_file, trees_file, ref_group_file


def make_ref_group_dict(ref_group_file):
    """make a dictionary of reference node names as key and its subgroup as value"""

    ref_group_dict = {}
    with open(ref_group_file, 'r') as f:
        for line in f:
            line = line.rstrip().split('\t')
            ref_group_dict[line[0]] = line[1]
    return ref_group_dict


def get_prefix_list(jplace_file):
    """get the query sequence ID from jplace file"""

    prefix_list = []
    with open(jplace_file, 'r') as f:
        for placement in json.load(f)['placements']:
            prefix_list.append([ID[0] for ID in placement.get('nm', [])])
    return prefix_list


def get_single_tree(trees_file, prefix_list, output_dir, fam_name):
    """get single tree for each query sequence from the total tree file grabbed by guppy"""

    output_dir = os.path.join(output_dir, "{}_single_tree".format(fam_name))
    os.mkdir(output_dir) if not os.path.exists(output_dir) else None
    prefix_singletrees = []
    with open(trees_file, 'r') as f:
        i = 0
        for line in f:
            prefix = prefix_list[i]
            i += 1
            singletrees = []
            for ID in prefix:
                # there may be multiple IDs with the same sequence, the guppy will merge them in one tree
                # the merge ID starts with the first ID in the list
                # but the oder of the trees for different sequence is the same as in the jplace file
                # remove plantismash_ prefix in the query_id
                singletree = os.path.join(output_dir, "{}.tre".format(ID[12:]))
                with open(singletree, "w") as file:
                    file.write(line)
                singletrees.append(singletree)
            prefix_singletrees.append((prefix,singletrees))
    return prefix_singletrees


def get_subtree_subgroup(trees_file, prefix_list, ref_group_dict):
    """get the subgroups of the subtrees of query sequences; query_subgroups is a dictionary with
    dictionary {(subgroup1, subgroup2, ...): likelihood; (subgroup1, ...): likelihood} as value
    of each query sequence ID as key"""

    query_subgroups = {}
    with open(trees_file, 'r') as f:
        i = 0
        for line in f:
            prefix = prefix_list[i]
            i += 1
            tree = Tree(line, format=1)
            # get the leave names of the same query sequence ID in different places on the tree
            target_leaves = [leaf for leaf in tree.get_leaf_names() if leaf.startswith(prefix[0])]

            subgroup_likelihood = {}
            for target_leaf in target_leaves:
                likelihood = float(target_leaf.split('=')[-1])

                neighbor_subgroup_set = set()
                # get the parent node of the target leaf
                parent_node = tree.search_nodes(name=target_leaf)[0].up
                # get the subgroups of the all leaves under parent node
                for leaf in parent_node.iter_leaves():
                    subgroup = ref_group_dict.get(leaf.name)
                    neighbor_subgroup_set.add(subgroup) if subgroup else None
                subgroup = tuple(neighbor_subgroup_set)

                # if the subgroups result is the same to it in other places case, add the likelihood
                subgroup_likelihood[subgroup] = subgroup_likelihood.get(subgroup, 0) + likelihood
            # if the query sequence only has one subgroups result on the tree, set the likelihood to 1
            if len(subgroup_likelihood) == 1:
                subgroup_likelihood[subgroup_likelihood.keys()[0]] = 1
            for ID in prefix:
                # remove plantismash_ prefix in the query_id
                query_subgroups[ID[12:]] = subgroup_likelihood
    return query_subgroups


# def write_subtree_subgroup(result_path, query_subgroups):
#     with open(result_path, "w") as file:
#         file.write("# query" + '\n')
#         file.write("# subgroup1,subgroup2" + '\t likelihood' + '\n')
#         file.write("# subgroup3" + '\t likelihood' + '\n')
#         for query, subgroup_likelihood in query_subgroups.items():
#             file.write(query + '\n')
#             for subgroup, likelihood in subgroup_likelihood.items():
#                 file.write(', '.join(subgroup) + '\t' + str(likelihood) + '\n')


def run_subtree_subgroup(fam_name, refpkg_dir, output_dir, fam_fasta):
    """run pplacer to place query on the reference tree and get the subgroup of the subtree of query"""

    jplace_file, trees_file, ref_group_file = run_pplacer(fam_fasta, refpkg_dir, fam_name, output_dir)

    prefix_list = get_prefix_list(jplace_file)

    ref_group_dict = make_ref_group_dict(ref_group_file)

    prefix_singletrees = get_single_tree(trees_file, prefix_list, output_dir, fam_name)

    tree_query_subgroup = get_subtree_subgroup(trees_file, prefix_list, ref_group_dict)

    return tree_query_subgroup, prefix_singletrees, ref_group_file
