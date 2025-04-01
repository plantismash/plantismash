#!/usr/bin/env python3

"""
Author: Hannah E. Augustijn
University: Wageningen University and Research & Leiden University
Department: Department of Bioinformatics & Institute of Biology Leiden
Date: 20/09/2024
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from Bio import SeqIO
import subprocess
import sys
import matplotlib.gridspec as gridspec


def get_sequence_motifs(fasta_s60, fasta_s200, gene_name, outdir):
    """ Wrapper function for creating sequence motifs of DAP-seq peaks """
    out_path = os.path.join(outdir, "meme_out")
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    meme_out = run_meme(fasta_s60, out_path, 5)
    con_motif, motifs = parse_meme(meme_out)
    if con_motif:
        pfm = create_pfm(motifs)
        shan_ic = get_ic(pfm)
        create_logo(shan_ic, gene_name, out_path)
    else:
        meme_out = run_meme(fasta_s200, out_path, 5)
        con_motif, motifs = parse_meme(meme_out)
        if con_motif:
            pfm = create_pfm(motifs)
            shan_ic = get_ic(pfm)
            create_logo(shan_ic, gene_name, out_path)
        else:
            print("No motif found")
            pfm = False
            shan_ic = False
    return pfm, shan_ic


def run_meme(fasta_file, outdir, min_width):
    """ Run MEME to identify motifs in binding sequences """
    reg_name = fasta_file.split("/")[-1].split(".")[0]
    out_file = f"{outdir}/{reg_name}.meme"
    cmd_meme = f"meme {fasta_file} -o {outdir} -text -dna -minw " \
               f"{min_width} -maxw 30 -nmotifs 3 -evt 0.05 -revcomp > {out_file}"
    try:
        if not os.path.exists(out_file):
            subprocess.check_output(cmd_meme, shell=True,
                                    stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print(f"Unable to run meme with command: {cmd_meme}")
        sys.exit()
    return out_file


def parse_motif_coordinates(meme_out):
    """ Extract the coordinates and consensus of the motif from MEME results """
    consensus_motif = ""
    start_motif = 0
    end_motif = 0
    for num, line in enumerate(meme_out, 1):
        if "E-value = " in line:
            consensus_motif = line.split()[1]
        elif line.startswith("BL   MOTIF"):
            start_motif = num
        elif line.startswith("//"):
            end_motif = num
    return consensus_motif, start_motif, end_motif


def parse_meme(meme_results):
    """ parses the meme output to lists of motifs and consensus sequences """
    motif_list = []
    with open(meme_results, "r") as meme_out:
        consensus_motif, start_motif, end_motif = parse_motif_coordinates(meme_out)
        if not consensus_motif:
            print(f"Could not parse results from the MEME output. Please check the MEME output")
            return False, False

        else:
            # seek the index in the file and extract the motifs
            meme_out.seek(0)
            txt = meme_out.readlines()
            motifs = txt[start_motif: end_motif - 1]
            for motif in motifs:
                motif_list.append(motif.split(" ")[-4])

            return consensus_motif, motif_list


def create_pfm(motifs):
    """ constructs a position frequency matrix from a list of binding sites

     :param motifs: list of binding sites
     :returns: pfm_dict: position frequency matrix as dictionary
     """
    pfm_dict = {"A": [], "C": [], "G": [], "T": []}

    for j in range(len(motifs[0])):
        pfm_dict["A"].append(([i[j] for i in motifs].count("A")))
        pfm_dict["C"].append(([i[j] for i in motifs].count("C")))
        pfm_dict["G"].append(([i[j] for i in motifs].count("G")))
        pfm_dict["T"].append(([i[j] for i in motifs].count("T")))
    return pfm_dict


def get_ic(pfm):
    """ Create information content profiles for Shannon entropy

     :param pfm: position frequency matrix as dictionary
     :returns: shan_ent: information content profiles for Shannon entropy as dataframe
     """
    shan_ent = {"A": [], "C": [], "G": [], "T": []}

    df = pd.DataFrame.from_dict(pfm)
    df["sum"] = df.sum(axis=1)
    ppm = df.loc[:, "A":"T"].div(df["sum"], axis=0)

    for index, row in ppm.iterrows():
        # SeqLogo information content profile (Shannon uncertainty measure)
        entropy = 2 + (row["A"] * (0 if (row["A"]) == 0 else np.log2(row["A"]))) + \
                  (row["C"] * (0 if (row["C"]) == 0 else np.log2(row["C"]))) + \
                  (row["G"] * (0 if (row["G"]) == 0 else np.log2(row["G"]))) + \
                  (row["T"] * (0 if (row["T"]) == 0 else np.log2(row["T"])))

        shan_ent["A"].append(row["A"] * entropy)
        shan_ent["C"].append(row["C"] * entropy)
        shan_ent["G"].append(row["G"] * entropy)
        shan_ent["T"].append(row["T"] * entropy)

    return pd.DataFrame.from_dict(shan_ent)


def create_logo(ic, reg_name, out_dir):
    """ Creates sequence logo in pdf format

     :param ic: information content profiles for Shannon entropy
     :param reg_name: name of the transcription regulator
     :param out_dir: output directory
     :returns: None
     """
    pdf_file = os.path.join(out_dir, f"{reg_name}_logo.pdf")
    png_file = os.path.join(out_dir, f"{reg_name}_logo.png")
    if not os.path.exists(pdf_file):
        plt.figure(figsize=(90, 35))

        logo = logomaker.Logo(ic, alpha=.75)
        logo.style_spines(visible=False)
        logo.style_spines(spines=('left', 'bottom'), visible=True)
        logo.style_xticks(rotation=0, fmt='%d', anchor=0)

        logo.ax.set_ylabel("Information content", labelpad=-1)
        logo.ax.xaxis.set_ticks_position('none')
        logo.ax.xaxis.set_tick_params(pad=-1)
        logo.ax.set_xticklabels([x + 1 for x in ic.index.values])
        logo.ax.set_ylim([0, 2])

        plt.savefig(pdf_file, format="pdf")
        plt.savefig(png_file, format="png")
    return