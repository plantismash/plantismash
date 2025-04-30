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

import os
import re
import tempfile
import subprocess
import argparse
from ete3 import Tree
import colorsys

def generate_color(hue):
    rgb = colorsys.hls_to_rgb(hue, 0.5, 1.0)
    return '#{:02x}{:02x}{:02x}'.format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))

def replace_special_characters(name):
    return re.sub(r'[^\w]', '_', name)

def modify_leaf_names(tree_file, modified_prefix, verbose=False):
    tree = Tree(tree_file, format=1)
    modified_names = {}
    for leaf in tree.iter_leaves():
        original_name = leaf.name
        modified_name = replace_special_characters(original_name)
        leaf.name = modified_name
        modified_names[modified_name] = original_name

    target_leaves = [leaf for leaf in tree.get_leaf_names() if leaf.startswith(modified_prefix)]
    neighbor_leave_set = set()
    for target_leaf in target_leaves:
        parent_node = tree.search_nodes(name=target_leaf)[0].up
        for leaf in parent_node.iter_leaves():
            neighbor_leave_set.add(leaf.name)

    neighbor_names_dic = {name: modified_names[name] for name in neighbor_leave_set}

    if verbose:
        print(f"[DEBUG] Modified {len(modified_names)} leaf names. {len(neighbor_leave_set)} neighbors identified.")

    return tree, neighbor_names_dic, neighbor_leave_set

def make_group_color(ref_group_file, neighbor_leave_set, verbose=False):
    with open(ref_group_file, 'r') as infile:
        groups = {}
        next(infile)
        for row in infile:
            row = row.strip().split('\t')
            leaf_name = replace_special_characters(row[0])
            group_name = row[1]
            if group_name not in groups:
                groups[group_name] = []
            groups[group_name].append(leaf_name)

    make_color = {g: l for g, l in groups.items() if neighbor_leave_set & set(l)}

    group_colors = {}
    total_groups = len(make_color)
    hue_step = 1.0 / total_groups if total_groups else 1.0
    for i, group_name in enumerate(make_color):
        hue = i * hue_step
        color = generate_color(hue)
        group_colors[group_name] = color

    if verbose:
        print(f"[DEBUG] Assigned colors to {len(group_colors)} subgroups.")

    return group_colors, make_color

def write_annotation_file(annotations, prefix, ann_file, query_id, fam_name, group_colors, make_color):
    with open(ann_file, 'w') as f:
        f.write('# all special characters changed to underscore in leaf names for Graphlan compatibility\n')
        f.write("clade_marker_size\t0\n")
        f.write(f"{fam_name} family tree\tclade_marker_font_size\t20\n")
        f.write("shown leaves under queries` parent_nodes\tclade_marker_font_size\t20\n")

        query = f"{query_id}_query_node"
        f.write(f"{query}\tclade_marker_color\tk\n")
        f.write(f"{query}\tclade_marker_size\t30\n")
        f.write(f"{query}\tclade_marker_font_size\t20\n")
        f.write(f"{query}\tclade_marker_shape\t*\n")

        for group_name, color in group_colors.items():
            f.write(f"{group_name}\tclade_marker_color\t{color}\n")
            f.write(f"{group_name}\tclade_marker_size\t30\n")
            f.write(f"{group_name}\tclade_marker_font_size\t20\n")
            f.write(f"{group_name}\tclade_marker_shape\ts\n")
            for leaf_name in make_color[group_name]:
                f.write(f"{leaf_name}\tannotation_background_color\t{color}\n")

        for modified_name, original_name in annotations.items():
            f.write(f"{modified_name}\tannotation_rotation\t90\n")
            f.write(f"{modified_name}\tannotation_font_size\t4\n")
            if original_name.startswith(prefix[0]):
                f.write(f"{modified_name}\tannotation_background_color\twhite\n")
                f.write(f"{modified_name}\tclade_marker_shape\t*\n")
                f.write(f"{modified_name}\tclade_marker_size\t30\n")
                color = 'r' if (prefix[-1] + "_#0_") in original_name else 'k'
                f.write(f"{modified_name}\tclade_marker_color\t{color}\n")
                label = query_id + "#" + original_name.split('#')[-1]
                f.write(f"{modified_name}\tannotation\t{label}\n")
            else:
                f.write(f"{modified_name}\tannotation\t{original_name}\n")

def run_graphlan(ref_group_file, tree_svg_dir, prefix, singletree, fam_name, verbose=False):
    os.makedirs(tree_svg_dir, exist_ok=True)

    temp_dir = tempfile.mkdtemp()
    modified_prefix = replace_special_characters(prefix[0])
    query_id = os.path.basename(singletree)[:-4]

    tree_file_modified = os.path.join(temp_dir, 'tree_modified.tre')
    ann_file = os.path.join(temp_dir, 'ann.txt')
    tree_xml_file = os.path.join(temp_dir, 'tree.xml')
    tree_svg_file = os.path.join(tree_svg_dir, f'{query_id}.svg')

    if verbose:
        print(f"[INFO] Generating Graphlan SVG for {query_id} -> {tree_svg_file}")

    tree, modified_names, neighbor_leave_set = modify_leaf_names(singletree, modified_prefix, verbose=verbose)
    group_colors, make_color = make_group_color(ref_group_file, neighbor_leave_set, verbose=verbose)
    tree.write(outfile=tree_file_modified)

    write_annotation_file(modified_names, prefix, ann_file, query_id, fam_name, group_colors, make_color)

    try:
        subprocess.run(["graphlan_annotate.py", "--annot", ann_file, tree_file_modified, tree_xml_file], check=True)
        subprocess.run(["graphlan.py", tree_xml_file, tree_svg_file, "--external_legends"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Graphlan failed for {query_id}: {e}")

    if not os.path.exists(tree_svg_file):
        print(f"[WARNING] SVG not created for {query_id}. Check Graphlan installation and input files.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Graphlan to generate tree SVGs from subgroup analysis")
    parser.add_argument("--ref-group-file", required=True, help="Reference group annotation file")
    parser.add_argument("--tree-dir", required=True, help="Directory to output SVG trees")
    parser.add_argument("--prefix", nargs='+', required=True, help="Prefix used to identify query sequences")
    parser.add_argument("--tree", required=True, help="Input single tree file")
    parser.add_argument("--family", required=True, help="Name of the protein/domain family")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    run_graphlan(
        ref_group_file=args.ref_group_file,
        tree_svg_dir=args.tree_dir,
        prefix=args.prefix,
        singletree=args.tree,
        fam_name=args.family,
        verbose=args.verbose
    )
