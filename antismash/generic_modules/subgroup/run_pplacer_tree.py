# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2024 Ziqiang Luo
# Wageningen University & Research, NL
# Bioinformatics Group, Department of Plant Sciences
#
# Copyright (C) 2025 Elena Del Pup
# Wageningen University & Research, NL
# Bioinformatics Group, Department of Plant Sciences
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Run pplacer to place query sequence on the reference tree and get the subgroup of the subtree of query sequence"""
import subprocess
import json
from ete3 import Tree
import os

def run_pplacer(fam_fasta, refpkg_dir, fam_name, output_dir):
    """Run pplacer to place query protein sequences on the reference tree in refpkg_dir"""

    # Get reference alignment and group info paths from refpkg_dir
    contents_file = os.path.join(refpkg_dir, "CONTENTS.json")
    if not os.path.exists(contents_file):
        raise FileNotFoundError(f"Missing required file: {contents_file}")

    with open(contents_file, 'r') as f:
        contents_dict = json.load(f)
        ref_alignment_path = os.path.join(refpkg_dir, contents_dict['files']["aln_fasta"])
        ref_group_file = os.path.join(refpkg_dir, contents_dict['files']["seq_info"])

    if not os.path.exists(ref_alignment_path):
        raise FileNotFoundError(f"Missing alignment: {ref_alignment_path}")

    # Output paths
    hmm_path = os.path.join(output_dir, f"{fam_name}.hmm")
    align_sto = os.path.join(output_dir, f"{fam_name}_query_ref_aligned.sto")
    jplace_file = os.path.join(output_dir, f"{fam_name}.jplace")
    trees_file = os.path.join(output_dir, f"{fam_name}_total.tre")

    # Run commands safely with error checking
    subprocess.run(["hmmbuild", hmm_path, ref_alignment_path], check=True)
    subprocess.run(["hmmalign", "-o", align_sto, "--mapali", ref_alignment_path, hmm_path, fam_fasta], check=True)
    subprocess.run(["pplacer", "-o", jplace_file, "-c", refpkg_dir, align_sto], check=True)
    subprocess.run(["guppy", "sing", "-o", trees_file, jplace_file], check=True)

    return jplace_file, trees_file, ref_group_file

def make_ref_group_dict(ref_group_file):
    ref_group_dict = {}
    with open(ref_group_file, 'r') as f:
        for line in f:
            line = line.rstrip().split('\t')
            ref_group_dict[line[0]] = line[1]
    return ref_group_dict

def get_prefix_list(jplace_file):
    prefix_list = []
    with open(jplace_file, 'r') as f:
        for placement in json.load(f)['placements']:
            prefix_list.append([ID[0] for ID in placement.get('nm', [])])
    return prefix_list

def get_single_tree(trees_file, prefix_list, output_dir, fam_name):
    output_dir = os.path.join(output_dir, f"{fam_name}_single_tree")
    os.makedirs(output_dir, exist_ok=True)

    prefix_singletrees = []
    with open(trees_file, 'r') as f:
        for i, line in enumerate(f):
            prefix = prefix_list[i]
            singletrees = []
            for ID in prefix:
                singletree = os.path.join(output_dir, f"{ID[12:]}.tre")
                with open(singletree, "w") as file:
                    file.write(line)
                singletrees.append(singletree)
            prefix_singletrees.append((prefix, singletrees))
    return prefix_singletrees

def get_subtree_subgroup(trees_file, prefix_list, ref_group_dict):
    query_subgroups = {}
    with open(trees_file, 'r') as f:
        for i, line in enumerate(f):
            prefix = prefix_list[i]
            tree = Tree(line, format=1)
            target_leaves = [leaf for leaf in tree.get_leaf_names() if leaf.startswith(prefix[0])]

            subgroup_likelihood = {}
            for target_leaf in target_leaves:
                likelihood = float(target_leaf.split('=')[-1])
                parent_node = tree.search_nodes(name=target_leaf)[0].up
                neighbor_subgroup_set = {
                    ref_group_dict.get(leaf.name)
                    for leaf in parent_node.iter_leaves()
                    if ref_group_dict.get(leaf.name)
                }
                subgroup = tuple(neighbor_subgroup_set)
                subgroup_likelihood[subgroup] = subgroup_likelihood.get(subgroup, 0) + likelihood

            if len(subgroup_likelihood) == 1:
                subgroup_likelihood[list(subgroup_likelihood.keys())[0]] = 1

            for ID in prefix:
                query_subgroups[ID[12:]] = subgroup_likelihood

    return query_subgroups

def run_subtree_subgroup(fam_name, refpkg_dir, output_dir, fam_fasta):
    jplace_file, trees_file, ref_group_file = run_pplacer(fam_fasta, refpkg_dir, fam_name, output_dir)
    prefix_list = get_prefix_list(jplace_file)
    ref_group_dict = make_ref_group_dict(ref_group_file)
    prefix_singletrees = get_single_tree(trees_file, prefix_list, output_dir, fam_name)
    tree_query_subgroup = get_subtree_subgroup(trees_file, prefix_list, ref_group_dict)
    return tree_query_subgroup, prefix_singletrees, ref_group_file
