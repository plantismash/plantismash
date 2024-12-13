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
# Copyright (C) 2024 Ziqiang Luo
# Wageningen University & Research, NL
#
# Copyright (C) 2024 Elena Del Pup 
# Wageningen University & Research, NL
# Bioinformatics Group, Department of Plant Sciences 
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""ClusterBlast comparative gene cluster analysis"""

import logging
import subprocess
from os import path
from antismash import utils
from clusterblast import load_clusterblast_database, internal_homology_blast, perform_clusterblast, filter_overlap
from data_loading import prepare_data, generate_Storage_for_cb
import os 
name = "clusterblast"
short_description = name.capitalize()
priority = 10000


# Tuple is ( binary_name, optional)
_required_binaries = [
    ('blastp', False),
    ('makeblastdb', False),
    ('diamond', False),
]

_required_files = [
    # ('geneclusterprots.dmnd', False),
    # ('geneclusterprots.fasta', False),
    # ('geneclusters.txt', False),
    ('plantgeneclusterprots.dmnd', False),
    ('plantgeneclusterprots.fasta', False),
    ('plantgeneclusters.txt', False),
]

def merge_files_in_directory(directory, merged_fasta, merged_txt):
    """Merge all .fasta and .txt files in a directory into the provided output files."""
    fasta_files = [path.join(directory, f) for f in os.listdir(directory) if f.endswith('.fasta')]
    txt_files = [path.join(directory, f) for f in os.listdir(directory) if f.endswith('.txt')]

    # Handle case with only one FASTA and one TXT file
    if len(fasta_files) == 1 and len(txt_files) == 1:
        logging.info("Only one FASTA and one TXT file found.")
        logging.info("Renaming {} to {}".format(fasta_files[0], merged_fasta))
        logging.info("Renaming {} to {}".format(txt_files[0], merged_txt))
        os.rename(fasta_files[0], merged_fasta)
        os.rename(txt_files[0], merged_txt)
    else:
        # Merge FASTA files
        with open(merged_fasta, 'w') as fasta_out:
            for fasta_file in fasta_files:
                logging.info("Merging FASTA file: {}".format(fasta_file))
                with open(fasta_file, 'r') as f:
                    fasta_out.write(f.read())
        
        # Merge TXT files
        with open(merged_txt, 'w') as txt_out:
            for txt_file in txt_files:
                logging.info("Merging TXT file: {}".format(txt_file))
                with open(txt_file, 'r') as f:
                    txt_out.write(f.read())


def generate_diamond_database(fasta_file, output_dmnd):
    """Generate a Diamond database from the given FASTA file."""
    logging.info("Generating Diamond database from FASTA file")
    diamond_command = [
        "diamond", "makedb",
        "--in", fasta_file,
        "-d", output_dmnd
    ]
    try:
        subprocess.check_call(diamond_command)
        logging.info("Diamond database created: {}".format(output_dmnd))
    except subprocess.CalledProcessError as e:
        logging.error("Error generating Diamond database: {}".format(e))
        raise


def ensure_required_files():
    """Ensure all required files are present or generated."""
    clusterblast_dir = path.dirname(path.abspath(__file__))  # Directory where this script is located
    merged_fasta = path.join(clusterblast_dir, "plantgeneclusterprots.fasta")
    merged_txt = path.join(clusterblast_dir, "plantgeneclusters.txt")
    dmnd_file = path.join(clusterblast_dir, "plantgeneclusterprots")

    # Merge .fasta and .txt files into consolidated files if missing
    if not path.exists(merged_fasta) or not path.exists(merged_txt):
        logging.info("Merging all FASTA and TXT files in the clusterblast directory")
        merge_files_in_directory(clusterblast_dir, merged_fasta, merged_txt)

    # Generate Diamond database if not present
    if not path.exists("{}.dmnd".format(dmnd_file)):
        logging.info("Generating the Diamond database file")
        generate_diamond_database(merged_fasta, dmnd_file)


def check_prereqs(options):
    """Check if all required applications and files are available or generate them."""
    if options.clusterblastdir == "":
        options.clusterblastdir = utils.get_full_path(__file__, '')

    # Check binaries
    logging.info("Checking prerequisites for ClusterBlast")
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate binary: %r" % binary_name)

    # Ensure required files are available or generated
    clusterblast_dir = options.clusterblastdir
    merged_fasta = path.join(clusterblast_dir, "plantgeneclusterprots.fasta")
    merged_txt = path.join(clusterblast_dir, "plantgeneclusters.txt")
    dmnd_file = path.join(clusterblast_dir, "plantgeneclusterprots")

    if not path.exists(merged_fasta) or not path.exists(merged_txt):
        logging.info("Required FASTA/TXT files not found. Merging files in ClusterBlast directory.")
        merge_files_in_directory(clusterblast_dir, merged_fasta, merged_txt)

    if not path.exists("{}.dmnd".format(dmnd_file)):
        logging.info("Diamond database not found. Generating it from the FASTA file.")
        generate_diamond_database(merged_fasta, dmnd_file)

    # Check required files explicitly
    for file_name, optional in _required_files:
        file_path = path.join(clusterblast_dir, file_name)
        if not path.exists(file_path) and not optional:
            failure_messages.append("Failed to locate required file: %r" % file_name)

    if failure_messages:
        for message in failure_messages:
            logging.error(message)
    return failure_messages


def run_clusterblast(seq_record, options):
    logging.info('Running ClusterBlast')
    clusters, proteinlocations, proteinstrands, proteinannotations, proteintags = load_clusterblast_database(seq_record)
    seq_record.internalhomologygroupsdict = internal_homology_blast(seq_record)
    perform_clusterblast(options, seq_record, clusters, proteinlocations, proteinstrands, proteinannotations, proteintags)
    prepare_data(seq_record, options, searchtype="general")
    generate_Storage_for_cb(options, seq_record)



def make_geneclusterprots(seq_records, options, output_filename="plantgeneclusterprots.fasta"):
    """make gene cluster proteins fasta file for clusterblast"""
    # Check if seq_records is empty
    if not seq_records:
        logging.warning("No sequence records provided to make_geneclusterprots.")
        return
    
    names = []
    seqs = []
    logging.info("Received {} sequence records for processing.".format(len(seq_records)))

    for seq_record in seq_records:
        geneclusters = utils.get_sorted_cluster_features(seq_record)
        if not geneclusters:
            logging.warning("No gene clusters found for sequence record {}.".format(seq_record.id))
            continue

        for genecluster in geneclusters:
            queryclusternames = []
            queryclusterseqs = []
            strand_start_ends = []
            queryclusterprots = filter_overlap(utils.get_cluster_cds_features(genecluster, seq_record))

            if not queryclusterprots:
                logging.warning("No CDS features found for cluster {} in {}".format(genecluster, seq_record.id))
                continue

            # logging.info("Total CDS features for cluster {} in {}: {}".format(genecluster, seq_record.id, len(queryclusterprots)))

            for cds in queryclusterprots:
                # logging.debug("Processing CDS: {}".format(cds))
                if cds.strand == 1:
                    strand = "+"
                else:
                    strand = "-"
                start = str(cds.location.start).replace(">", "").replace("<", "")
                end = str(cds.location.end).replace(">", "").replace("<", "")
                strand_start_end = (strand, start, end)
                
                logging.debug("Strand-start-end for CDS: {}".format(strand_start_end))
                logging.debug("Current strand_start_ends: {}".format(strand_start_ends))

                if strand_start_end not in strand_start_ends:
                    # todo: Incompletely overlapping splicing should be treated as one gene if their sequences are 50% similar？
                    strand_start_ends.append(strand_start_end)
                    annotation = utils.get_gene_annotation(cds)

                    fullname = "|".join([seq_record.id, "c" + str(utils.get_cluster_number(genecluster)),
                                         start + "-" + end,
                                         strand, utils.get_gene_id(cds), annotation.replace(' ', '_') , utils.get_gene_acc(cds)])
                    # scaffold|cluster number|location|strand|locustag|annotation|protein or accession
                    queryclusternames.append(fullname)
                    queryclusterseqs.append(str(utils.get_aa_sequence(cds)))

            for i in range(len(queryclusternames)):
                names.append(queryclusternames[i])
                seqs.append(queryclusterseqs[i])
    
    # Check if no sequences were collected
    if not names or not seqs:
        logging.warning("No sequences found in make_geneclusterprots!")
        return

    if not path.exists(options.clusterblastdir):
        try:
            os.makedirs(options.clusterblastdir)
            logging.info("Created directory for clusterblast output: {}".format(options.clusterblastdir))
        except Exception as e:
            logging.error("Failed to create directory {}: {}".format(options.clusterblastdir, e))
            return

    outputname = path.join(options.clusterblastdir, output_filename)
    
    logging.info("Writing the following sequences to FASTA: ")
    for i in range(len(names)):
        logging.info("> {}\n{}".format(names[i], seqs[i]))
        # how many sequences are written to the file

    try:
        utils.writefasta(names, seqs, outputname)
        logging.info("FASTA file {} with {} sequences written to {}".format(output_filename, len(names), outputname))
    except Exception as e:
        logging.error("Failed to write FASTA file {}: {}".format(outputname, e))

def where_is_clusterblast():
    return utils.get_full_path(__file__, '')