# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2024 Ziqiang Luo
# Wageningen University & Research, NL
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os
import re
import tempfile
from ete3 import Tree
import colorsys


def generate_color(hue):
    rgb = colorsys.hls_to_rgb(hue, 0.5, 1.0)
    return '#{:02x}{:02x}{:02x}'.format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))


def replace_special_characters(name):
    return re.sub(r'[^\w]', '_', name)


def modify_leaf_names(tree_file, modified_prefix):
    tree = Tree(tree_file, format=1)

    # change special characters to underscore in leaf names for Graphlan compatibility
    modified_names = {}
    for leaf in tree.iter_leaves():
        original_name = leaf.name
        modified_name = replace_special_characters(original_name)
        leaf.name = modified_name
        modified_names[modified_name] = original_name

    # get the leave names of the same query sequence ID in different places on the tree
    target_leaves = [leaf for leaf in tree.get_leaf_names() if leaf.startswith(modified_prefix)]
    neighbor_leave_set = set()
    for target_leaf in target_leaves:
        # get the parent node of the target leaf
        parent_node = tree.search_nodes(name=target_leaf)[0].up
        for leaf in parent_node.iter_leaves():
            neighbor_leave_set.add(leaf.name)

    neighbor_names_dic = {}
    for neighbor_leave in neighbor_leave_set:
        neighbor_names_dic[neighbor_leave] = modified_names[neighbor_leave]

    return tree, neighbor_names_dic, neighbor_leave_set


def make_group_color(ref_group_file, neighbor_leave_set):
    with open(ref_group_file, 'r') as infile:
        groups = {}
        next(infile)
        for row in infile:
            row = row.strip().split('\t')
            # change special characters to underscore in leaf names for Graphlan compatibility
            leaf_name = replace_special_characters(row[0])
            group_name = row[1]
            if group_name not in groups:
                groups[group_name] = []
            groups[group_name].append(leaf_name)

    make_color = {}
    for group_name, leaves in groups.items():
        if neighbor_leave_set & set(leaves):
            make_color[group_name] = leaves

    # generate colors for each group
    group_colors = {}
    total_groups = len(make_color)
    hue_step = 1.0 / total_groups

    for i, group_name in enumerate(make_color):
        hue = i * hue_step
        color = generate_color(hue)
        group_colors[group_name] = color

    return group_colors, make_color


def write_annotation_file(annotations, prefix, ann_file, query_id, fam_name, group_colors, make_color):
    with open(ann_file, 'w') as f:
        f.write('# all special characters changed to underscore in leaf names for Graphlan compatibility\n')
        # f.write("title\t{} in the subtree(s) of {} family\n".format(query_id, fam_name))
        f.write("clade_marker_size\t0\n")
        f.write('{}\tclade_marker_font_size\t20\n'.format(fam_name + " family tree"))
        f.write('shown leaves under queries` parent_nodes\tclade_marker_font_size\t20\n')

        query = query_id + "_#No._M=likelyhood"
        f.write('{}\tclade_marker_color\tk\n'.format(query))
        f.write('{}\tclade_marker_size\t30\n'.format(query))
        f.write('{}\tclade_marker_font_size\t20\n'.format(query))
        f.write('{}\tclade_marker_shape\t*\n'.format(query))
        for group_name, color in group_colors.items():
            f.write('{}\tclade_marker_color\t{}\n'.format(group_name, color))
            f.write('{}\tclade_marker_size\t30\n'.format(group_name))
            f.write('{}\tclade_marker_font_size\t20\n'.format(group_name))
            f.write('{}\tclade_marker_shape\ts\n'.format(group_name))
            for leaf_name in make_color[group_name]:
                f.write('{}\tannotation_background_color\t{}\n'.format(leaf_name, color))

        for modified_name, original_name in annotations.items():
            f.write('{}\tannotation_rotation\t90\n'.format(modified_name))
            f.write("{}\tannotation_font_size\t4\n".format(modified_name))
            if original_name.startswith(prefix[0]):
                f.write('{}\tannotation_background_color\twhite\n'.format(modified_name))
                f.write('{}\tclade_marker_shape\t*\n'.format(modified_name))
                f.write('{}\tclade_marker_size\t30\n'.format(modified_name))
                if (prefix[-1] + "_#0_") in original_name:
                    f.write('{}\tclade_marker_color\tr\n'.format(modified_name))
                else:
                    f.write('{}\tclade_marker_color\tk\n'.format(modified_name))
                f.write('{}\tannotation\t{}\n'.format(modified_name, query_id +"#"+ original_name.split('#')[-1]))
            else:
                f.write('{}\tannotation\t{}\n'.format(modified_name, original_name))


def run_graphlan(ref_group_file, tree_svg_dir, prefix, singletree, fam_name):

    temp_dir = tempfile.mkdtemp()
    modified_prefix = replace_special_characters(prefix[0])
    query_id = os.path.basename(singletree)[:-4]

    tree_file_modified = os.path.join(temp_dir, 'tree_modified.tre')
    ann_file = os.path.join(temp_dir, 'ann.txt')
    tree_xml_file = os.path.join(temp_dir, 'tree.xml')

    tree_svg_file = os.path.join(tree_svg_dir, '{}.svg'.format(query_id))

    # Modify the leaf node name and save the modified tree file
    tree, modified_names, neighbor_leave_set = modify_leaf_names(singletree, modified_prefix)
    group_colors, make_color = make_group_color(ref_group_file, neighbor_leave_set)
    tree.write(outfile=tree_file_modified)

    write_annotation_file(modified_names, prefix, ann_file, query_id, fam_name, group_colors, make_color)

    # run graphlan
    os.system('graphlan_annotate.py --annot {} {} {}'.format(ann_file, tree_file_modified, tree_xml_file))
    os.system('graphlan.py {} {}  --external_legends'.format(tree_xml_file, tree_svg_file))
