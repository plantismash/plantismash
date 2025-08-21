#!/usr/bin/env python
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
# Copyright (C) 2024 Elena Del Pup 
# Wageningen University & Research, NL
# Bioinformatics Group, Department of Plant Sciences 
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""Run the antiSMASH pipeline"""

import sys
import os
import random

if sys.platform ==  ('win32') or sys.platform == ('darwin'):
    os.environ['EXEC'] = os.getcwd() + os.sep + "exec"
    os.environ['PYTHON'] = os.getcwd() + os.sep + "python"
    sys.path.append(os.sep.join([os.getcwd(), "python", "Lib", "site-packages"]))
    os.environ['PATH'] = os.pathsep + os.environ['PYTHON'] + os.pathsep + os.environ['PATH']
    os.environ['PATH'] = os.pathsep + os.environ['EXEC'] + os.pathsep + os.environ['PATH']

import logging
import argparse
from os import path
import multiprocessing
from straight.plugin import load
from helperlibs.bio import seqio
from antismash.config import load_config, set_config
from antismash import utils
from antismash.output_modules import xls
from antismash.generic_modules import check_prereqs as generic_check_prereqs
from antismash.generic_genome_modules import check_prereqs as ggm_check_prereqs
from antismash.generic_modules import (
    hmm_detection,
    genefinding,
    fullhmmer,
    smcogs,
    clusterblast,
    subclusterblast,
    knownclusterblast,
    active_site_finder,
    coexpress,
    gff_parser,
    subgroup, 
    tfbs_finder 
)
try:
    from antismash.db.biosql import get_record
    from antismash.db.extradata import getExtradata
    USE_BIOSQL = True
except ImportError: 
    USE_BIOSQL = False
from numpy import array_split, array
import urllib.request, urllib.error, urllib.parse
from urllib.error import URLError
import http.client
import time
from collections import defaultdict
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
try:
    from io import StringIO
except ImportError:
    from io import StringIO
from datetime import datetime
import re


def ValidateDetectionTypes(detection_models):
    class Validator(argparse.Action):
        def __call__(self, parser, args, values, option_string = None):
            if values == "all":
                values = ",".join(detection_models).replace(";",",").split(",")
            else:
                values = values.replace(";",",").split(",")
                try:
                    for value in values:
                        if value not in detection_models:
                            raise ValueError('invalid detection models {s!r}'.format(s = value))
                except ValueError as e:
                    print("\nInput error:", e, "\n")
                    print("Please choose from the following list:\n", "\n".join(detection_models))
                    sys.exit(1)
            setattr(args, self.dest, values)
    return Validator


def ValidateClusterTypes(clustertypes):
    #TODO : list according to selected detection models
    class Validator(argparse.Action):
        def __call__(self, parser, args, values, option_string = None):
            values = values.replace(";",",").split(",")
            try:
                for value in values:
                    if value not in clustertypes:
                        raise ValueError('invalid clustertype {s!r}'.format(s = value))
            except ValueError as e:
                print("\nInput error:", e, "\n")
                print("Please choose from the following list:\n", "\n".join(clustertypes), "\n\nExample: --enable t1pks,nrps,other")
                sys.exit(1)
            setattr(args, self.dest, values)
    return Validator


def main():
    "Actually run the pipeline"
    multiprocessing.freeze_support()

 # We print the logging debug message as soon as the logging object is initialized after parsing the cmdline
    start_time = datetime.now()

    # First load the output plugins so we can present appropriate options
    output_plugins = load_output_plugins()

    detection_models = hmm_detection.get_supported_detection_models()
    clustertypes = hmm_detection.get_supported_cluster_types()

    parser = utils.getArgParser()

    #positional arguments
    parser.add_argument('sequences',
                        metavar='sequence',
                        nargs="*",
                        help="GenBank/EMBL/FASTA file(s) containing DNA.")

    #optional non-grouped arguments
    parser.add_argument('-h', '--help',
                        dest='help',
                        action='store_true',
                        default=False,
                        help="Show this help text.")
    parser.add_argument('--help-showall',
                        dest='help_showall',
                        action='store_true',
                        default=False,
                        help="Show full lists of arguments on this help text.")
    parser.add_argument('-c', '--cpus',
                        dest='cpus',
                        type=int,
                        default=multiprocessing.cpu_count(),
                        help="How many CPUs to use in parallel. (default: %(default)s)")

    ## grouped arguments

    group = parser.add_argument_group('Basic analysis options', '', basic=True)
    group.add_argument('--input-type',
                       dest='input_type',
                       default='nucl',
                       choices=['nucl', 'prot'],
                       help="Determine input type: amino acid sequence(s) or nucleotide sequence(s). (default: %(default)s)")
    group.add_argument('--taxon',
                        dest='taxon',
                        default='bacteria',
                        choices=['bacteria', 'fungi', 'plants'],
                        help="Determine the taxon from which the sequence(s) came from. (default: plants)")
    group.add_argument('--hmmsearch-chunk',
                       dest='hmmsearch_chunk',
                       default=10000,
                       type=int,
                       help="Number of sequences per hmm_search run (set lower for CPU with less RAM). (default: %(default)s)")

    group = parser.add_argument_group('Additional analysis',
                                      'Use -h along with the corresponding options to show '
                                      'additional analysis specific parameters. (e.g.: '
                                      '"{prog} -h --clusterblast --smcogs")'.format(prog=parser.prog),
                                      basic=True)
    group.add_argument('--cdhit',
                       dest='enable_cdhit',
                       action='store_true',
                       default=False,
                       help="Use CD-HIT to filter tandem arrays on hmm_detection.")
    group.add_argument('--clusterblast',
                       dest='clusterblast',
                       action='store_true',
                       default=False,
                       help="Compare identified clusters against a database of antiSMASH-predicted clusters.")
    group.add_argument('--update_clusterblast',
                       dest='update_clusterblast',
                       action='store_true',
                       default=False,
                       help="update the database of plantiSMASH-predicted clusters (plantgeneclusters.txt, plantgeneclusterprots.fasta, plantgeneclusterprots.dmnd) using this time found clusters.")
    group.add_argument('--clusterblastdir',
                       dest='clusterblastdir',
                       type=str,
                       default="",
                       help="the clusterblastdir contain the database of antiSMASH-predicted clusters (plantgeneclusters.txt, plantgeneclusterprots.fasta).")
    group.add_argument('--subclusterblast',
                       dest='subclusterblast',
                       action='store_true',
                       default=False,
                       help="Compare identified clusters against known subclusters responsible for synthesising precursors.")
    group.add_argument('--knownclusterblast',
                       dest='knownclusterblast',
                       action='store_true',
                       default=False,
                       help="Compare identified clusters against known gene clusters from the MIBiG database.")
    group.add_argument('--coexpress',
                       dest='coexpress',
                       action='store_true',
                       default=False,
                       help="Compare identified clusters against a prepared Expression data.")
    group.add_argument('--tfbs',
                   dest='tfbs_detection',
                   action='store_true',
                   default=False,
                   help="Run transcription factor binding site (TFBS) prediction.")
    group.add_argument('--tfbs-pvalue', dest='tfbs_pvalue', type=float, default= 1e-4,
                   help="TFBS p-value cutoff (default from config)")
    group.add_argument('--tfbs-range', dest='tfbs_range', type=int, default=500,
                    help="Range around clusters to scan for TFBS (default from config)")
    group.add_argument('--disable_subgroup',
                       dest='disable_subgroup',
                       action='store_true',
                       default=False,
                       help="disable identifying the subgroup.")
    group.add_argument('--disable_treesvg',
                       dest='disable_treesvg',
                       action='store_true',
                       default=False,
                       help="disable making svg pictures of subgroupping tree to save time.")
    group.add_argument('--subgroup_inputpath',
                       dest="subgroup_inputpath",
                       type=str,
                       default="",
                       help="give path of folder with the structure same to subgroup folder in antismash/generic_modules.")
    group.add_argument('--disable_specific_modules',
                       dest='disable_specific_modules',
                       action='store_true',
                       default=False,
                       help="disable specific modules.")
    group.add_argument('--smcogs',
                       dest='smcogs',
                       action='store_true',
                       default=False,
                       help="Look for sec. met. clusters of orthologous groups.")
    group.add_argument('--inclusive',
                       dest='inclusive',
                       action='store_true',
                       default=False,
                       help="Use inclusive ClusterFinder algorithm for additional cluster detection.")
    group.add_argument('--borderpredict',
                       dest='borderpredict',
                       action='store_true',
                       default=False,
                       help="Use ClusterFinder algorithm to predict gene cluster borders.")
    group.add_argument('--full-hmmer',
                       dest='full_hmmer',
                       action='store_true',
                       default=False,
                       help="Run a whole-genome HMMer analysis.")
    group.add_argument('--asf',
                       dest='run_asf',
                       action='store_true',
                       default=False,
                       help='Run active site finder module.')
    group.add_argument('--ecpred',
                       dest='ecpred',
                       default='none',
                       choices=['none', 'kegg', 'eficaz'],
                       help="Run EC prediction. kegg = query KEGG REST API (Internet connection required!) "
                            "eficaz = use EFICAZ2.5a EC predictor. (default: %(default)s)")
    group.add_argument('--modeling',
                       dest="modeling",
                       default='none',
                       choices=['none', 'sco', 'eco'],
                       help="Run homology based metabolic modeling pipeline against template model "
                            "eco or sco. (default: %(default)s)")

    group = parser.add_argument_group('Advanced options')
    group.add_argument('--limit',
                       dest="limit",
                       type=int,
                       default=1000,
                       help="Only process the first <limit> records (default: %(default)s). -1 to disable")
    group.add_argument('--from',
                       dest='start',
                       type=int,
                       default=-1,
                       help="Start analysis at nucleotide specified.")
    group.add_argument('--to',
                       dest='end',
                       type=int,
                       default=-1,
                       help="End analysis at nucleotide specified")
    group.add_argument('--pfamdir',
                       dest='pfamdir',
                       default=argparse.SUPPRESS,
                       help="Directory the Pfam-A.hmm file is located in.")
    group.add_argument('--models',
                        metavar="TYPES",
                        dest='enabled_detection_models',
                        action=ValidateDetectionTypes(detection_models),
                        default=["default"],
                        help="Select user-defined hmm detection models to be used.")
    group.add_argument('--enable',
                       metavar="TYPES",
                       dest='enabled_cluster_types',
                       action=ValidateClusterTypes(clustertypes),
                       default=clustertypes,
                       help="Select sec. met. cluster types to search for.")
    group.add_argument('--fix-id-line',
                       dest='fix_id_line',
                       action='store_true',
                       default=False,
                       help="Try to fix invalid sequence file ID lines.")
    group.add_argument('--min-domain-number',
                       dest='min_domain_number',
                       type=int,
                       default=0,
                       help="Minimum unique domains required for minimum_rules hmm detection.")
    group.add_argument('--dynamic-cutoff',
                       dest='enable_dynamic_cutoff',
                       action='store_true',
                       default=False,
                       help="Use dynamic cutoff multiplier based on local intergenic distances.")
    group.add_argument('--cutoff-multiplier',
                       dest='cutoff_multiplier',
                       type=float,
                       default=1.00,
                       help="Multiplier for tuning hmm_detection cutoff globally.")
    group.add_argument('--gene-num-cutoff',
                       dest='gene_num_cutoff',
                       type=int,
                       default=0,
                       help="Accept core genes within kb distance or having separated "
                            "by equal to or less than this number of intervening genes.")
    group.add_argument('--gene-num-cutoff-only',
                        dest='gene_num_cutoff_only',
                        action='store_true',
                        default=False,
                        help="Only use number of intervening genes for the cutoff.")
    group.add_argument('--ignore-short-aa',
                       dest='ignore_short_aa',
                       action='store_true',
                       default=True,
                       help="In hmm detection, ignore short AA without significant pfam hits (affect dynamic cutoffs).")
    group.add_argument('--require-internal-cyclopeptide-repeats',
                        dest='require_internal_cyclopeptide_repeats',
                        action='store_true',
                        default=False,
                        help="Only accept cyclopeptide clusters where repeats occur inside BURP features (default: repeats can be anywhere if BURP is present)")
    group = parser.add_argument_group('Gene finding options (ignored when ORFs are annotated)')
    group.add_argument('--genefinding',
                       dest='genefinding',
                       default='none',
                       choices=['glimmer', 'prodigal', 'prodigal-m', 'none'],
                       help="Specify algorithm used for gene finding: Glimmer, "
                            "Prodigal, or Prodigal Metagenomic/Anonymous mode. (default: %(default)s).")
    group.add_argument('--eukaryotic',
                       dest='eukaryotic',
                       action='store_true',
                       default=False,
                       help="DNA is of eukaryotic origin.")
    group.add_argument('--all-orfs',
                       dest='all_orfs',
                       action='store_true',
                       default=False,
                       help="Use all ORFs > 60 nt instead of running gene finding.")
    group.add_argument('--gff3',
                       dest='gff3',
                       default=False,
                       help="Specify GFF3 file to extract features from.")
    group.add_argument('--use_phase',
                       dest='use_phase',
                       default=False,
                       action='store_true',
                       help="This flag forces CDS coordinates in the GFF3 to be modified by phase before translation. Otherwise this field of the GFF3 file will be ignored. (Phytozome GFF3s require applying phase before translation, whereas NCBI, augustus and glimmerhmm GFF3s do not.)")
    group.add_argument('--glimmerhmm-train_folder',
                       dest='glimmerhmm_train_folder',
                       default='crypto',
                       help="Training models for the glimmerHMM (see antismash/generic_modules/genefinding/).")


    group = parser.add_argument_group('CD-HIT specific options', '',
                                      param=["--cdhit"])
    group.add_argument('--cdh-memory',
                       dest='cdh_memory',
                       type=str,
                       default="2000",
                       help="memory of using CD-HIT, default 2000M.")
    group.add_argument('--cdh-cutoff',
                       dest='cdh_cutoff',
                       default=0.5,
                       help="Cutoff for CD-HIT filtering, lower means more strict filtering. (0.2 - 1.0)")
    group.add_argument('--cdh-display-cutoff',
                       dest='cdh_display_cutoff',
                       default=0.2,
                       help="Display percentage identities for cluster of core genes passing this cutoff value. (0.2-1.0)")

    group = parser.add_argument_group('ClusterBlast specific options', '',
                                      param=["--clusterblast", "--subclusterblast", "--knownclusterblast"])
    group.add_argument('--nclusters',
                       dest='nclusters',
                       type=int,
                       default=10,
                       help="Number of clusters from ClusterBlast to display.")
    group.add_argument('--seed',
                       dest='seed',
                       type=int,
                       default=0,
                       help="Random number seed for ClusterBlast coloring.")
    group.add_argument('--homologyscalelimit',
                       dest='homologyscalelimit',
                       type=float,
                       default=0.0,
                       help="If positive float number greater than 1, limit horizontal shrinkage "
                            "of the graphical display of the query BGC in ClusterBlast results to "
                            "this ratio. Warning: some homologous genes may no longer be visible!")
    group.add_argument('--dbgclusterblast',
                       dest='dbgclusterblast',
                       default="",
                       help='Retrieve the resuts from previous antiSMASH run stored in the specified directory.')

    group = parser.add_argument_group('CoExpress specific options', '',
                                      param=["--coexpress"])
    group.add_argument('--coexpress-soft_file',
                       dest='coexpress_soft_file',
                       default="",
                       help="GEO soft formatted file (.soft) from either GEO Series (GSE) or GEO Dataset (GDS)."
                            "Use ',' to separate multiple input files.")
    group.add_argument('--coexpress-csv_file',
                      dest='coexpress_csv_file',
                      default="",
                      help="CSV formatted gene expression file."
                           "Use ',' to separate multiple input files.")
    group.add_argument('--coexpress-min_MAD',
                       dest='coexpress_min_MAD',
                       default="",
                       help="By default, genes with a Mean Absolute Deviation (MAD) of 0 get filtered out of any correlation calculation. "
                            "A minimum MAD value can be defined to filter out genes with values below the specified one. "
                            "If option is set to 0, the feature is turned off entirely as MAD is never below 0.")

    group = parser.add_argument_group('ClusterFinder specific options', '',
                                      param=["--inclusive", "--borderpredict"])
    group.add_argument('--cf_cdsnr',
                       dest='cdsnr',
                       type=int,
                       default=5,
                       help="Minimum size of a ClusterFinder cluster, in number of CDS.")
    group.add_argument('--cf_threshold',
                       dest='cf_prob_thres',
                       type=float,
                       default=0.6,
                       help="ClusterFinder mean probability threshold.")
    group.add_argument('--cf_npfams',
                       dest='cf_npfams',
                       type=int,
                       default=5,
                       help="Minimum number of biosynthetic Pfam domains in a ClusterFinder cluster.")

    group = parser.add_argument_group('Output options', '', basic=True)
    group.add_argument('--outputfolder',
                       dest='outputfoldername',
                       default=argparse.SUPPRESS,
                       help="Directory to write results to.")
    group.add_argument('--db',
                       dest='use_db',
                       action='store_true',
                       default=False,
                       help="write results to BioSQL database.")
    group.add_argument('--dboverwrite',
                       dest='db_overwrite',
                       action="store_true",
                       default=False,
                       help="Re-run antiSMASH analysis and overwrite results with same accession number in BioSQL database.")
    for plugin in output_plugins:
        group.add_argument('--disable-%s' % plugin.name,
                           dest=plugin.name,
                           action='store_false',
                           default=argparse.SUPPRESS,
                           help="Disable %s" % plugin.short_description)

    group = parser.add_argument_group("Debugging & Logging options", '', basic=True)
    group.add_argument('-v', '--verbose',
                       dest='verbose',
                       action='store_true',
                       default=False,
                       help="Print verbose status information to stderr.")
    group.add_argument('-d', '--debug',
                       dest='debug',
                       action='store_true',
                       default=False,
                       help="Print debugging information to stderr.")
    group.add_argument('--logfile',
                       dest='logfile',
                       default=argparse.SUPPRESS,
                       help="Also write logging output to a file.")
    group.add_argument('--statusfile',
                       dest='statusfile',
                       default=argparse.SUPPRESS,
                       help="Write the current status to a file.")
    group.add_argument('--list-plugins',
                       dest='list_plugins',
                       action='store_true',
                       default=False,
                       help="List all available sec. met. detection modules.")
    group.add_argument('--check-prereqs',
                       dest='check_prereqs_only',
                       action='store_true',
                       default=False,
                       help="Just check if all prerequisites are met.")
    group.add_argument('-V', '--version',
                       dest='version',
                       action='store_true',
                       default=False,
                       help="Display the version number and exit.")

    ## endof grouped arguments

    #if --help, show help texts and exit
    if (list(set(["-h", "--help", "--help-showall"]) & set(sys.argv))):
        parser.print_help(None, "--help-showall" in sys.argv)
        sys.exit(0)

    #Parse arguments, removing hyphens from the beginning of file names to avoid conflicts with argparse
    infile_extensions = ('.fasta', '.fas', '.fa', '.gb', '.gbk', '.emb', '.embl')
    sys.argv = [arg.replace("-","< > HYPHEN < >") if (arg.endswith(infile_extensions) and arg[0] == "-") else arg for arg in sys.argv]
    options = parser.parse_args()
    options.sequences = [filename.replace("< > HYPHEN < >","-") for filename in options.sequences]

    #Remove ClusterFinder cluster types when ClusterFinder is not turned on
    if not options.inclusive:
        options.enabled_cluster_types = [cl_type for cl_type in options.enabled_cluster_types if not cl_type.startswith("cf_")]

    # Logging is useful for all the following code, so make sure that is set up
    # right after parsing the arguments.
    setup_logging(options)

    #preset options according to taxon
    apply_taxon_preset(options)

    #Remove cluster types not included in selected detection models
    temp_ct = []
    for cl_type in options.enabled_cluster_types:
        if len(cl_type.split("/")) > 1:
            if cl_type.split("/")[0] in options.enabled_detection_models:
                temp_ct.append(cl_type)
        elif "default" in options.enabled_detection_models :
            temp_ct.append(cl_type)
    options.enabled_cluster_types = temp_ct

    #if -V, show version texts and exit
    if options.version:
        print("antiSMASH %s" % utils.get_version())
        sys.exit(0)

    logging.info("starting plantiSMASH {version}, called with {cmdline}".format(version=utils.get_version(), cmdline=" ".join(sys.argv)))
    logging.debug("plantismash analysis started at %s", str(start_time))
    options.run_info = {}
    options.run_info["ver"] = utils.get_version()
    options.run_info["param"] = " ".join(sys.argv[1:-1])
    options.run_info["start"] = str(start_time)


    if options.nclusters > 50:
        logging.info("Number of clusters (%d) is too large. Reducing to 50.", options.nclusters)
        options.nclusters = 50
    logging.info("Number of clusters = %d", options.nclusters)
    if options.seed != 0:
        random.seed(options.seed)

    if options.enable_cdhit:
        cdh_cutoff = float(options.cdh_cutoff)
        if not 0.2 <= cdh_cutoff <= 1.0:
            logging.error('cd-hit cutoff is on wrong value (use 0.2 - 1.0)')
            sys.exit(1)
        cdh_display_cutoff = float(options.cdh_display_cutoff)
        if not 0.2 <= cdh_display_cutoff <= 1.0:
            logging.error('cd-hit display cutoff is on wrong value (use 0.2 - 1.0)')
            sys.exit(1)

    if options.input_type == 'prot' and (
            options.clusterblast or
            options.knownclusterblast or
            options.subclusterblast or
            options.inclusive or
            options.full_hmmer):
        logging.error("Protein input option is not compatible with --clusterblast, --subclusterblast, " \
                      "--knownclusterblast, --inclusive, and --full-hmmer")
        sys.exit(2)
    if options.input_type == 'prot' and (options.start != -1 or options.end != -1):
        logging.error("Protein input option is not compatible with --start and --end.")
        sys.exit(2)

    if options.coexpress:
        # check if soft file exist
        if (options.coexpress_soft_file == "" and options.coexpress_csv_file == ""):
            logging.error("Required .soft or .csv not found (set with --coexpress-soft_file PATH or --coexpress-csv_file PATH).")
            sys.exit(2)
        options.coexpress_soft_files = []
        options.coexpress_csv_files = []
        if options.coexpress_soft_file:
            options.coexpress_soft_file = path.abspath(options.coexpress_soft_file)
            soft_files = options.coexpress_soft_file.split(",")
            for fname in soft_files:
                if fname == "" or not path.isfile(fname):
                    logging.error("Required .soft not found (%s)." % fname)
                    sys.exit(2)
                else:
                    options.coexpress_soft_files.append(path.abspath(fname))
        if options.coexpress_csv_file:
            options.coexpress_csv_file = path.abspath(options.coexpress_csv_file)
            csv_files = options.coexpress_csv_file.split(",")
            for fname in csv_files:
                if fname == "" or not path.isfile(fname):
                    logging.error("Required .csv not found (%s)." % fname)
                    sys.exit(2)
                else:
                    options.coexpress_csv_files.append(path.abspath(fname))

    if options.coexpress_min_MAD != '':
        try:
            float(options.coexpress_min_MAD)
        except ValueError:
            logging.error("coexpress_min_MAD value must be a float.")
            sys.exit(2)

    if not USE_BIOSQL:
        logging.info("BioSQL not loaded. Turning off DB access")
        options.use_db = False

    #Load configuration data from config file
    load_config(options)
    set_config(options)

    # Unpack TFBS config section into top-level options
    options.tfbs_pvalue = getattr(options, "tfbs_pvalue", 1e-4)
    options.tfbs_range = getattr(options, "tfbs_range", 1000)

    # set up standard DB namespace
    options.dbnamespace = options.BioSQLconfig.dbgenomenamespace

    #Load and filter plugins
    utils.log_status("Loading detection plugins")
    if not options.disable_specific_modules:
        plugins = load_specific_modules()
    else:
        plugins = []
    if options.list_plugins:
        list_available_plugins(plugins, output_plugins)
        sys.exit(0)
    filter_plugins(plugins, options, clustertypes)

    filter_outputs(output_plugins, options)

    # Check prerequisites
    if check_prereqs(plugins, options) > 0:
        logging.error("Not all prerequisites met")
        sys.exit(1)
    if options.check_prereqs_only:
        logging.info("All prerequisites are met")
        sys.exit(0)
    if not options.sequences:
        parser.error("Please specify at least one sequence file")
    if options.gff3 and not os.path.exists(options.gff3):
        logging.error('No file found at %r', options.gff3)
        sys.exit(1)

    if 'outputfoldername' not in options:
        options.outputfoldername = path.splitext(path.basename(options.sequences[0]))[0]
    if not os.path.exists(options.outputfoldername):
        os.mkdir(options.outputfoldername)
    options.full_outputfolder_path = path.abspath(options.outputfoldername)

    #Remove old clusterblast files
    clusterblast_files = ["subclusterblastoutput.txt", "clusterblastoutput.txt", "knownclusterblastoutput.txt"]
    for clusterblast_file in clusterblast_files:
        clusterblastoutput_path = os.path.join(options.full_outputfolder_path, clusterblast_file)
        if os.path.exists(clusterblastoutput_path):
            os.remove(clusterblastoutput_path)

    if options.debug and os.path.exists(options.dbgclusterblast):
        logging.debug("Using %s instead of computing Clusterblasts and variantes!", options.dbgclusterblast)

    if ("train_%s" % options.glimmerhmm_train_folder) not in genefinding.get_supported_training_models():
        logging.error("Specified glimmerHMM training models not found (%s)!", options.glimmerhmm_train_folder)
        sys.exit(1)

    # Structure for options.extrarecord
    #
    # options.extrarecord[seq_record_id]=Namespace()
    # options.extrarecord[seq_record_id].extradata[extradataID] can contain the different (deserialized) objects,
    # e.g. ...extradata['SubClusterBlastData'] will contain the storage object for subclsuterblast results
    options.extrarecord = {}

    #Parse input sequence
    try:
        utils.log_status("Parsing the input sequence(s)")
        seq_records = parse_input_sequences(options)
        if options.input_type == 'nucl':
            seq_records = [record for record in seq_records if len(record.seq) > 1000]
            if len(seq_records) == 0:
                logging.error("Input does not contain contigs larger than minimum size of 1000 bp.")
                sys.exit(1)
    except IOError as e:
        logging.error(str(e))
        sys.exit(1)

    if len(seq_records) < 1:
        logging.error("Sequence file is incorrectly formatted or contains no sequences of sufficient quality.")
        sys.exit(1)

    # parse coexpress input file
    if options.coexpress:
        logging.info("Parsing CoExpress input file...")
        options.geo_dataset = coexpress.parse_geofiles(options.coexpress_soft_files)
        options.geo_dataset.extend(coexpress.parse_csvfiles(options.coexpress_csv_files))
        options.gene_expressions = []
        geo_ids = []
        for geo in options.geo_dataset:
            options.gene_expressions.append({})
            for seq_record in seq_records:
                options.gene_expressions[-1][seq_record.id] = {}
            geo_id = geo["info"]["id"]
            i = 0
            while geo_id in geo_ids:
                i += 1
                geo_id = "%s_%d" % (geo["info"]["id"], i)
            if i > 0:
                logging.warning("Duplicated GEO rec_id found, renaming to %s..." % geo_id)
                geo["info"]["id"] = geo_id
            geo_ids.append(geo_id)

        # calculate samples ranges
        logging.info("Finding CoExpress locus tag columns...")
        for geo in options.geo_dataset:
            coexpress.find_col_id(geo, seq_records)

        # calculate samples ranges
        logging.info("Calculating CoExpress sample ranges...")
        for geo in options.geo_dataset:
            coexpress.calc_sample_ranges(geo)

        features_to_match = {}
        logging.info("Doing whole genome expression matching...")
        for seq_record in seq_records:
            features_to_match[seq_record.id] = utils.get_cds_features(seq_record)
        for i in range(0, len(options.geo_dataset)):
            logging.debug("Matching expression for dataset (%d/%d).." % (i + 1, len(options.geo_dataset)))
            for seq_record in seq_records:
                logging.debug("Seq_record %s.." % (seq_record.id))
                options.gene_expressions[i][seq_record.id] = coexpress.match_exp_to_genes(features_to_match[seq_record.id], options.geo_dataset[i])
                # filter out gene_expressions in case of physical overlaps,
                # choose 1) the one with pfam hits (assuming has been filtered in hmm_detection filtering)
                # and/or 2) the one with highest avg expression values
                overlaps = utils.get_overlaps_table(seq_record)[0]
                prefix = "%s:" % seq_record.id.replace(":", "_")
                gene_expressions = {}
                all_gene_expressions = options.gene_expressions[i][seq_record.id]
                for overlap in overlaps:
                    chosen_gene = -1
                    best_exp = 0
                    for j in range(0, len(overlap)):
                        feature = overlap[j]
                        gene_id = utils.get_gene_id(feature)
                        if gene_id in all_gene_expressions:
                            if (prefix + gene_id) in options.hmm_results:
                                chosen_gene = j
                                break # this assumes one overlap cluster only have 1 gene with hits
                            exps = [all_gene_expressions[gene_id]["exp"][sample] for sample in all_gene_expressions[gene_id]["exp"]]
                            avg_exp = sum(exps) / len(exps)
                            if chosen_gene < 0 or best_exp < avg_exp:
                                chosen_gene = j
                    if chosen_gene >= 0:
                        gene_id = utils.get_gene_id(overlap[chosen_gene])
                        gene_expressions[gene_id] = all_gene_expressions[gene_id]
                options.gene_expressions[i][seq_record.id] = gene_expressions
    # The current implementation of the modeling pipeline works on seq_record level; to support multi-contig sequences, we have to modify it
    # to work with seq_records instead and be called outside the "for seq_record in seq_records" loop

#    if (not options.modeling == "none") and (len(seq_records) > 1):
#        logging.error('Metabolic modeling is currently only available for finished genome sequences, i.e. sequences subitted as 1 contig')
#        options.modeling = "none"



    options.record_idx = 1
    options.orig_record_idx =1


    temp_seq_records = []
    for seq_record in seq_records:
        old_seq_id = seq_record.id
        utils.log_status("Analyzing record %d" % options.record_idx)
        logging.info("Analyzing record %d", options.record_idx)
        utils.sort_features(seq_record)
        strip_record(seq_record)
        utils.fix_record_name_id(seq_record, options)
        logging.debug("record name/id: %s", seq_record.id)

        if (seq_record.id != old_seq_id):
            # update seq_record id in dictionaries to match new seq_record id
            # this assume utils.fix_record_name_id retained the unique ids naming
            new_hmm_results = {}
            for composite_id in options.hmm_results:
                new_id = composite_id
                if (composite_id.split(":")[0] == old_seq_id):
                    new_id = ":".join([seq_record.id.replace(":", "_"), composite_id.split(":")[-1]])
                new_hmm_results[new_id] = options.hmm_results[composite_id]
            options.hmm_results = new_hmm_results
            if options.coexpress:
                for i in range(0, len(options.gene_expressions)):
                    if old_seq_id in options.gene_expressions[i]:
                        options.gene_expressions[i][seq_record.id] = options.gene_expressions[i][old_seq_id]
                        del options.gene_expressions[i][old_seq_id]

        if options.use_db and not options.db_overwrite:
            if utils.check_if_dbrecord_exists(seq_record.name, options):
                logging.warn("A record with accession number %s exits in antiSMASH-DB; " \
                             "using this record instead of re-calculating everything\n" \
                             "If a recalculation is desired, please add --dboverwrite parameter to command line", seq_record.name)
                try:
                    utils.log_status("retrieving record")
                    #TODO: This new association is not!!! Transferred to seq_records!!!
                    seq_record = get_record(seq_record.name, options)
                except Exception as e:
                    logging.exception("Uncaptured error %s when reading entries from antiSMASH-DB. This should not have happened :-(", e)
                    sys.exit(1)
                # set from_database flag to prevent updating the record in the output module
                temp_seq_records.append(seq_record)
                options.from_database = True

                # and change options for database export
                options.clusterblast = None
                options.subclusterblast = None
                options.knownclusterblast = None
                options.modeling = "none"
                options.smcogs="TRUE"
                options.enabled_cluster_types = ValidateClusterTypes(clustertypes)

                options.extrarecord[seq_record.id] = argparse.Namespace()

                # Get extradata from DB
                extradataHash = getExtradata(options, seq_record.id)

                options.extrarecord[seq_record.id].extradata = extradataHash

                if 'ClusterBlastData' in options.extrarecord[seq_record.id].extradata:
                    options.clusterblast = True
                if 'SubClusterBlastData' in options.extrarecord[seq_record.id].extradata:
                    options.subclusterblast = True
                if 'KnownClusterBlastData' in options.extrarecord[seq_record.id].extradata:
                    options.knownclusterblast = None
                if 'MetabolicModelDataObj' in options.extrarecord[seq_record.id].extradata:
                    options.modeling = "db"
                    options.metabolicmodeldir = options.outputfoldername + os.sep + "metabolicModel"

            # There must be a nicer solution...
            else:
                options.from_database = False
                run_analyses(seq_record, options, plugins)

        else:
            options.from_database = False
            run_analyses(seq_record, options, plugins)
            
        # Add seq_record to temp_seq_records 
        temp_seq_records.append(seq_record)

        options.record_idx += 1
        options.orig_record_idx += 1
        logging.debug("The record %s is originating from db %s", seq_record.name, options.from_database)

    # update coexpress data to calculate all possible high corelation betwen genes in clusters
    if options.coexpress:
        for i in range(0, len(options.geo_dataset)):
            coexpress.update_dist_between_clusters(seq_records, options.gene_expressions[i], options.geo_dataset[i])

    #check CoExpressinter-cluster relation
    if options.coexpress:
        options.coexpress_inter_clusters = {}
        for geo_id in utils.get_all_geo_rec_ids(seq_records):
            options.coexpress_inter_clusters[geo_id] = coexpress.get_inter_cluster_relation(seq_records,geo_id)

    #overwrite seq_records with newly assembled temp_seq_records
    seq_records = temp_seq_records
    del temp_seq_records

    # execute whole genome-scope modules
    run_generic_genome_modules(seq_records, options)

    options.run_info["end"] = str(datetime.now())

    #before writing to output, remove all hmm_detection's subdir prefixes from clustertype
    for seq_record in seq_records:
        for cluster in utils.get_cluster_features(seq_record):
            prod_names = []
            for prod in cluster.qualifiers['product']:
                prod_name = []
                for name in prod.split('-'):
                    prod_name.append(name.split('/')[-1])
                prod_names.append("-".join(prod_name))
            cluster.qualifiers['product'] = prod_names
        for cds in utils.get_cds_features(seq_record):
            if 'sec_met' in cds.qualifiers:
                temp_qual = []
                for row in cds.qualifiers['sec_met']:
                    if row.startswith('Type: '):
                        clustertypes = [(ct.split('/')[-1]) for ct in row.split('Type: ')[-1].split('-')]
                        temp_qual.append('Type: ' + "-".join(clustertypes))
                    elif row.startswith('Domains detected: '):
                        cluster_results = []
                        for cluster_result in row.split('Domains detected: ')[-1].split(';'):
                            cluster_results.append(cluster_result.split(' (E-value')[0].split('/')[-1] + ' (E-value' + cluster_result.split(' (E-value')[-1])
                        temp_qual.append('Domains detected: ' + ";".join(cluster_results))
                    else:
                        temp_qual.append(row)
                cds.qualifiers['sec_met'] = temp_qual

    #on plants, remove plant clustertype from hybrid types, and replace single
    #plant clustertype with "putative"
    if options.taxon == "plants":
        for seq_record in seq_records:
            for cluster in utils.get_cluster_features(seq_record):
                prod_names = []
                for prod in cluster.qualifiers['product']:
                    prod_name = list(set(prod.split('-')))
                    if (len(prod_name) > 1) and ("plant" in prod_name):
                        prod_name.remove("plant")
                    elif prod_name == ["plant"]:
                        prod_name = ["putative"]
                    prod_names.append("-".join(prod_name))
                cluster.qualifiers['product'] = prod_names
            for cds in utils.get_cds_features(seq_record):
                if 'sec_met' in cds.qualifiers:
                    temp_qual = []
                    for row in cds.qualifiers['sec_met']:
                        if row.startswith('Type: '):
                            clustertypes = list(set(row.split('Type: ')[-1].split('-')))
                            if (len(clustertypes) > 1) and ("plant" in clustertypes):
                                clustertypes.remove("plant")
                            elif clustertypes == ["plant"]:
                                clustertypes = ["putative"]
                            temp_qual.append('Type: ' + "-".join(clustertypes))
                        else:
                            temp_qual.append(row)
                    cds.qualifiers['sec_met'] = temp_qual

    # run subgroup identification
    # output family sequences fasta files and hmmscan results txt in the output folder, and update seq_records
    if not options.disable_subgroup:
        logging.info("identificating subgroup")
        subgroup.subgroup_identification(seq_records,path.abspath(options.outputfoldername), options)
    else:
        logging.info("subgroup identification is disabled")

    if options.update_clusterblast:
        # Updated function to generate unique output files using the input sequence name.
        logging.info("Updating ClusterBlast database for input: {0}".format(options.sequences))

        # Make sure the clusterblastdir directory exists
        if options.clusterblastdir == "":
            options.clusterblastdir = clusterblast.where_is_clusterblast()
        if not os.path.exists(options.clusterblastdir):
            os.makedirs(options.clusterblastdir)

        # Generate unique identifiers based on the input sequence
        try:
            # Split the input path into parts
            input_basename = os.path.normpath(options.sequences[0]).split(os.sep)
            # Get the last two elements of the filepath name 
            last_two_parts = "_".join(input_basename[-2:])  # Join with an underscore
            # Sanitize the combined parts to ensure filename safety
            sanitized_path = re.sub(r'[^\w\-_\.]', '_', last_two_parts)
            # Add more unique identifiers like timestamp if needed
            unique_name = sanitized_path
        except IndexError:
            raise ValueError("No sequences provided in options.sequences!")

        # File paths
        clusteblast_txt = os.path.join(options.clusterblastdir, "plantgeneclusters_{0}.txt".format(unique_name))
        clusteblast_fasta = os.path.join(options.clusterblastdir, "plantgeneclusterprots_{0}.fasta".format(unique_name))

         # Generate and write the plantgeneclusterprots.fasta file
        try:
            logging.debug("Calling make_geneclusterprots to generate FASTA file...")
            clusterblast.make_geneclusterprots(seq_records, options, "plantgeneclusterprots_{0}.fasta".format(unique_name))
            generated_fasta = os.path.join(options.clusterblastdir, "plantgeneclusterprots_{0}.fasta".format(unique_name))
            logging.debug("Expected generated FASTA file path: {0}".format(generated_fasta))

            if os.path.exists(generated_fasta) and os.path.getsize(generated_fasta) > 0:
                logging.debug("FASTA file {} exists and is non-empty.".format(generated_fasta))
            else:
                logging.error("FASTA file {} is missing or empty!".format(generated_fasta))
                logging.debug("ClusterBlast directory: {0}".format(options.clusterblastdir))
                raise IOError("FASTA file generation failed.")
        except Exception as e:
            logging.error("Error generating ClusterBlast FASTA file: {0}".format(e))
            raise

        # Generate and write the plantgeneclusters.txt file
        try:
            # generate the file in the result directory 
            xls.write(seq_records, options)
            print("plantgeneclusters.txt saved in the results directory")

             # Write the TXT file explicitly for the clusterblast directory
            with open(clusteblast_txt, "w") as clusteblast_file:
                for seq_record in seq_records:
                    clusters = utils.get_sorted_cluster_features(seq_record)
                    for cluster in clusters:
                        clustertype = utils.get_cluster_type(cluster)
                        clusternr = utils.get_cluster_number(cluster)
                        clustergenes = [utils.get_gene_id(cds) for cds in utils.get_cluster_cds_features(cluster, seq_record)]
                        accessions = [utils.get_gene_acc(cds) for cds in utils.get_cluster_cds_features(cluster, seq_record)]
                        # Write cluster data to TXT
                        clusteblast_file.write(
                            "\t".join(
                                [
                                    seq_record.id,
                                    seq_record.description,
                                    "c{0}".format(clusternr),
                                    clustertype,
                                    ";".join(clustergenes),
                                    ";".join(accessions),
                                ]
                            )
                            + "\n"
                        )
            print(("plantgeneclusters_{0}.txt saved in the clusterblast directory".format(unique_name)))
            
            # Validation: Ensure the file was created
            if not os.path.exists(clusteblast_txt):
                raise IOError("TXT file generation failed for clusterblast directory.")
            logging.info("TXT file successfully created at {0}".format(clusteblast_txt))

            # Count non-empty lines in the generated TXT file
            try:
                non_empty_line_count = 0
                with open(clusteblast_txt, "r") as txtfile:
                    non_empty_line_count = sum(1 for line in txtfile if line.strip())
                logging.info(
                    "Number of clusters for {0}: {1}, {2}".format(
                        unique_name, options.sequences, non_empty_line_count
                    )
                )
            except Exception as e:
                logging.error("Error counting non-empty lines in TXT file: {0}".format(e))
                raise

        except Exception as e:
            logging.error("Error generating ClusterBlast TXT file: {0}".format(e))
            raise

    # Write results and complete the process
    try:
        options.plugins = plugins
        utils.log_status("Writing the output files")
        logging.debug("Writing output for {0} sequence records".format(len(seq_records)))
        write_results(output_plugins, seq_records, options)
        zip_results(seq_records, options)

        # Log runtime
        end_time = datetime.now()
        running_time = end_time - start_time
        logging.debug(
            "antiSMASH calculation finished at {0}; runtime: {1}".format(
                str(end_time), str(running_time)
            )
        )
        utils.log_status("antiSMASH status: SUCCESS")
        logging.debug("antiSMASH status: SUCCESS")
    except Exception as e:
        logging.error("Error during results writing or finalizing: {0}".format(e))
        raise


def strip_record(seq_record):
    features = utils.get_cds_features(seq_record)
    for feature in features:
        if 'sec_met' in feature.qualifiers:
            del feature.qualifiers['sec_met']

        # remove aStorage note qualifiers
        if 'note' in feature.qualifiers:
            note_list=[]
            for notetext in feature.qualifiers['note']:
                if "aSstorage" not in notetext:
                    note_list.append(notetext)
            feature.qualifiers['note']=note_list

    clusters = utils.get_cluster_features(seq_record)
    deletefeatures = []
    for f in seq_record.features:
        if f in clusters:
            deletefeatures.append(seq_record.features.index(f))
    deletefeatures.reverse()
    for featurenr in deletefeatures:
        del seq_record.features[featurenr]


def apply_taxon_preset(options):
    "Apply presets according to specific --taxon"

    if options.taxon == "bacteria":
        logging.debug("Applying preset for %s", options.taxon)

    if options.taxon == "fungi":
        logging.debug("Applying preset for %s", options.taxon)
        options.eukaryotic = True

    # force the default preset for plants
    if not options.taxon:
        options.taxon = "plants"

    if options.taxon == "plants":
        logging.debug("Applying preset for %s", options.taxon)
        options.eukaryotic = True
        options.enable_cdhit = True
        if not "--limit" in sys.argv:
            options.limit = 9999
        if not "--cdh-cutoff" in sys.argv:
            options.cdh_cutoff = 0.5
        if not "--cdh-display-cutoff" in sys.argv:
            options.cdh_display_cutoff = 0.2
        if not "--min-domain-number" in sys.argv:
            options.min_domain_number = 2
        options.enable_dynamic_cutoff = True
        if not "--glimmerhmm-train_folder" in sys.argv:
            options.glimmerhmm_train_folder = "arabidopsis"
        if "default" in options.enabled_detection_models:
            options.enabled_detection_models.remove("default")
            if not "plants" in options.enabled_detection_models:
                options.enabled_detection_models.append("plants")

def run_analyses(seq_record, options, plugins):
    """Run antiSMASH analyses for a single SeqRecord"""

    if 'next_clusternr' not in options:
        options.next_clusternr = 1

    options.clusternr_offset = options.next_clusternr

    # Detect gene clusters
    detect_geneclusters(seq_record, options)

    for f in utils.get_cluster_features(seq_record):
        logging.debug(f)

    logging.debug(f"[DEBUG] run_analyses: disable_specific_modules = {options.disable_specific_modules}")
    logging.debug(f"[DEBUG] run_analyses: plugins passed in: {[p.name for p in plugins]}")

    if len(utils.get_cluster_features(seq_record)) > 0:
        # Run specific analyses first
        logging.debug("Running specific analyses for %s", seq_record.id)
        run_specific_analyses(seq_record, options, plugins)

        # Renumber the clusters to maintain contiguous numbering
        renumber_clusters(seq_record, options)

        # Run general analyses
        logging.debug("Running general analyses for %s", seq_record.id)
        run_general_analyses(seq_record, options)


def run_specific_analyses(seq_record, options, plugins):
    """Run specific cluster related analyses"""
    if options.disable_specific_modules:
        logging.debug("Skipping run_specific_analyses because --disable_specific_modules is set.")
        return
    cluster_specific_analysis(plugins, seq_record, options)



def run_general_analyses(seq_record, options):
    """Run general analyses that are independent of specific cluster types"""
    
    # Run smCOG analysis
    if options.smcogs:
        utils.log_status("Detecting smCOGs for contig #%d" % options.record_idx)
        smcogs.run_smcog_analysis(seq_record, options)

    # Run ClusterBlast
    if options.clusterblast:
        utils.log_status("ClusterBlast analysis for contig #%d" % options.record_idx)
        clusterblast.run_clusterblast(seq_record, options)
    
    # Run SubClusterBlast
    if options.subclusterblast:
        utils.log_status("SubclusterBlast analysis for contig #%d" % options.record_idx)
        subclusterblast.run_subclusterblast(seq_record, options)
    
    # Run KnownClusterBlast
    if options.knownclusterblast:
        utils.log_status("KnownclusterBlast analysis for contig #%d" % options.record_idx)
        knownclusterblast.run_knownclusterblast(seq_record, options)
    
    # Run CoExpress
    if options.coexpress:
        utils.log_status("Coexpression analysis for contig #%d" % options.record_idx)
        for i in range(0, len(options.geo_dataset)):
            coexpress.run_coexpress(seq_record, options.gene_expressions[i], options.geo_dataset[i])

    # Run Transcription Factor Binding Site (TFBS) analysis
    if options.tfbs_detection:
        utils.log_status(
            f"TFBS analysis for contig #{options.record_idx} (p-value={options.tfbs_pvalue}, region ±{options.tfbs_range} bp)")
        tfbs_finder.run_tfbs_finder_for_record(seq_record, options)
    
    # Run Active Site Finder
    if options.run_asf:
        ASFObj = active_site_finder.active_site_finder(seq_record, options)
        status = ASFObj.execute()
        if status:
            logging.debug("Active site finder execution successful")
        else:
            logging.error("Error in active site finder module!")



def renumber_clusters(seq_record, options):
    """Renumber clusters in the SeqRecord to be contiguous across all chromosomes."""
    if 'global_cluster_counter' not in options:
        options.global_cluster_counter = 1  
    for cluster in utils.get_cluster_features(seq_record):
        cluster.qualifiers['note'] = [note if not note.startswith("Cluster number") else "Cluster number: %d" % options.global_cluster_counter for note in cluster.qualifiers.get('note', [])]
        options.global_cluster_counter += 1



def list_available_plugins(plugins, output_plugins):
    print("Support for detecting the following secondary metabolites:")
    for plugin in plugins:
        print((" * %s" % plugin.short_description))

    print("\nSupport for the following modules acting in whole-genome scope")
    for plugin in load_generic_genome_plugins():
        print((" * %s: %s" % (plugin.name, plugin.short_description)))
    print("\nSupport for the following output formats:")
    for plugin in output_plugins:
        print((" * %s" % plugin.short_description))


def filter_plugins(plugins, options, clustertypes):
    if options.enabled_cluster_types is None or options.enabled_cluster_types == clustertypes:
        return

    for plugin in plugins:
        if plugin.name in clustertypes and plugin.name not in options.enabled_cluster_types:
            plugins.remove(plugin)

    if not options.disable_specific_modules and plugins == []:
        print("No plugins enabled, use --list-plugins to show available plugins")
        sys.exit(1)


def filter_outputs(plugins, options):
    for plugin in plugins:
        if plugin.name in options:
            logging.debug("Removing plugin %r", plugin.name)
            plugins.remove(plugin)

    if plugins == []:
        print("No plugins enabled, use --list-plugins to show available plugins")
        sys.exit(1)


def write_results(plugins, seq_records, options):
    for plugin in plugins:
        plugin.write(seq_records, options)


def zip_results(seq_records, options):
    "Create a zip archive with all the results generated so far"
    zip_name = '%s.zip' % seq_records[0].id
    utils.zip_path(path.abspath(options.outputfoldername), zip_name)


def setup_logging(options):
    "Set up the logging output"
    if options.debug:
        log_level = logging.DEBUG
    elif options.verbose:
        log_level = logging.INFO
    else:
        log_level = logging.WARNING

    logging.basicConfig(format='%(levelname)s: %(message)s',
                        level=log_level)
    if 'logfile' in options:
        if not (path.dirname(options.logfile) == "" or os.path.exists(path.dirname(options.logfile))):
            os.mkdir(path.dirname(options.logfile))
        fh = logging.FileHandler(options.logfile)
        fh.setLevel(logging.INFO)
        fh.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        logging.getLogger('').addHandler(fh)


def load_specific_modules():
    "Load available secondary metabolite detection modules"
    logging.info('Loading specific modules')
    detection_plugins = list(load('antismash.specific_modules'))

    logging.info("The following modules were loaded:%s "%(detection_plugins))
    #TODO remove this logging block

    # Sort after priority to ensure correct order of execution
    #detection_plugins.sort(cmp=lambda x, y: cmp(x.priority, y.priority))

    return detection_plugins


def load_output_plugins():
    "Load available output formats"
    plugins = list(load('antismash.output_modules'))
    plugins.sort(key=lambda x: x.priority)
    return plugins

def load_generic_genome_plugins():
    "Load available output formats"
    plugins = list(load('antismash.generic_genome_modules'))
    plugins.sort(key=lambda x: x.priority)

    return plugins

def fetch_entries_from_ncbi(efetch_url):
    urltry = "n"
    nrtries = 0
    output = ""
    while urltry == "n" and nrtries < 4:
        try:
            nrtries += 1
            time.sleep(3)
            req = urllib.request.Request(efetch_url)
            response = urllib.request.urlopen(req)
            output = response.read()
            if len(output) > 5:
                urltry = "y"
        except (IOError,http.client.BadStatusLine,URLError,http.client.HTTPException):
            logging.error("Entry fetching from NCBI failed. Waiting for connection...")
            time.sleep(5)
    return output


def fix_wgs_master_record(seq_record):
    updated_seq_records = []
    #If seq_record is a WGS master record, parse out contig accession numbers and download these as separate seq_records
    if 'wgs_scafld' in seq_record.annotations:
        contigranges = seq_record.annotations['wgs_scafld']
    else:
        contigranges = [seq_record.annotations['wgs']]
    allcontigs = []
    for contigrange in contigranges:
        if len(contigrange) == 1:
            allcontigs.extend(contigrange)
            continue
        startnumber, endnumber = '', ''
        alpha_tag = ''
        for char in contigrange[0].partition(".")[0]:
            if char.isdigit():
                startnumber += char
            else:
                alpha_tag += char
        for char in contigrange[1].partition(".")[0]:
            if char.isdigit():
                endnumber += char
        nrzeros = 0
        for char in startnumber:
            if char == "0":
                nrzeros += 1
            else:
                break
        contigrange = [alpha_tag + nrzeros * "0" + str(number) for number in range(int(startnumber), int(endnumber))]
        allcontigs.extend(contigrange)
    #Create contig groups of 50 (reasonable download size per download)
    nr_groups = len(allcontigs) / 50 + 1
    contig_groups = [list(np_array) for np_array in array_split(array(allcontigs), nr_groups)]
    #Return unchanged if no contigs supplied
    if len(contig_groups[0]) == 0 or (len(contig_groups[0]) == 1 and contig_groups[0][0] == ""):
        return [seq_record]
    #Download contigs and parse into seq_record objects
    for contig_group in contig_groups:
        efetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
        efetch_url = efetch_url + ",".join(contig_group) + '&rettype=gbwithparts&retmode=text'
        output = fetch_entries_from_ncbi(efetch_url)
        if not len(output) > 5:
            break
        if "Resource temporarily unavailable" in output[:200] or "<h1>Server Error</h1>" in output[:500] or "NCBI - WWW Error" in output[:500]:
            logging.error('ERROR: NCBI server temporarily unavailable: downloading contigs failed.')
            sys.exit(1)
        try:
            handle = StringIO(output)
            updated_seq_records.extend(list(SeqIO.parse(handle, 'genbank')))
        except ValueError as e:
            logging.error('Parsing %r failed: %s', "temporary contig file", e)
    return updated_seq_records


def fix_supercontig_record(seq_record):
    updated_seq_records = []
    #If seq_record is a supercontig record, reconstruct sequence and replace CONTIG feature by ORIGIN feature
    contig_info = seq_record.annotations['contig'].replace("join(","")
    if contig_info[-1] == ")":
        contig_info = contig_info[:-1]
    allcontigparts = contig_info.split(",")
    accessions = []
    for part in contig_info.split(","):
        if "gap(" not in part:
            accessions.append(part.partition(":")[0].partition(".")[0].rpartition("complement(")[2].replace(")",""))
    #Create contig groups of 50 (reasonable download size per download)
    nr_groups = len(accessions) / 50 + 1
    contig_groups = [list(np_array) for np_array in array_split(array(accessions), nr_groups)]
    #Return unchanged if no contigs supplied
    if len(contig_groups[0]) == 0 or (len(contig_groups[0]) == 1 and contig_groups[0][0] == ""):
        return [seq_record]
    #Download contig sequences based on their accessions
    contigseqdict = {}
    for contig_group in contig_groups:
        efetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
        efetch_url = efetch_url + ",".join(contig_group) + '&rettype=fasta&retmode=text'
        output = fetch_entries_from_ncbi(efetch_url)
        if not len(output) > 5:
            break
        if "Resource temporarily unavailable" in output[:200] or "<h1>Server Error</h1>" in output[:500]:
            logging.error('ERROR: NCBI server temporarily unavailable: downloading contigs failed.')
            sys.exit(1)
        sequences = [seq for seq in output.split(">") if len(seq) > 5]
        for sequence in sequences:
            for contig_acc in contig_group:
                if contig_acc in sequence.partition("\n")[0]:
                    contigseqdict[contig_acc] = sequence.partition("\n")[2].replace("\n","")
    #Reconstruct supercontig sequence based on contig sequences and gaps
    fullsequence = ''
    for part in allcontigparts:
        if "gap(" in part:
            candidate_gap_int = part.partition('gap(')[2][:-1]
            if "unk" in candidate_gap_int:
                candidate_gap_int = candidate_gap_int.partition("unk")[2]
            if candidate_gap_int.isdigit():
                fullsequence += int(candidate_gap_int) * 'N'
            else:
                logging.error('Parsing supercontig file failed: faulty gap identifier' + part)
                sys.exit(1)
        else:
            accession = part.partition(":")[0].partition(".")[0].rpartition("complement(")[2]
            sequence = contigseqdict[accession]
            if "complement(" in part:
                sequence = str(Seq(sequence, generic_dna).reverse_complement())
            if ":" in part and ".." in part:
                seqstart, seqend = part.partition(":")[2].replace(")","").replace("(","").split("..")
                if int(seqstart) > 0:
                    seqstart = int(seqstart) - 1
                sequence = sequence[int(seqstart) : int(seqend)]
            fullsequence += sequence
    #Add compiled sequence to seq_record
    SeqObject = Seq(fullsequence, generic_dna)
    seq_record.seq = SeqObject
    updated_seq_records.append(seq_record)
    return updated_seq_records


def process_wgs_master_scaffolds(seq_records):
    updated_seq_records = []
    for seq_record in seq_records:
        #Check if seq_record is a WGS master record or a supercontig record
        if not ('wgs_scafld' in seq_record.annotations or 'wgs' in seq_record.annotations or 'contig' in seq_record.annotations) or ('contig' in seq_record.annotations and len(str(seq_record.seq).replace("N","")) > 0):
            updated_seq_records.append(seq_record)
        else:
            if 'wgs_scafld' in seq_record.annotations or 'wgs' in seq_record.annotations:
                updated_seq_records.extend(fix_wgs_master_record(seq_record))
            elif 'contig' in seq_record.annotations and len(str(seq_record.seq).replace("N","")) == 0:
                updated_seq_records.extend(fix_supercontig_record(seq_record))
    return updated_seq_records


def is_nucl_seq(sequence):
    if len(str(sequence).lower().replace("a","").replace("c","").replace("g","").replace("t","").replace("n","")) < 0.2 * len(sequence):
        return True
    else:
        return False


def generate_nucl_seq_record(sequences):
    "Generate nucleotide seq_record"
    if len(sequences) == 0:
        return []
    seq_record = SeqRecord(Seq(""),id="Protein_Input", name="ProteinInput",
                   description="antiSMASH protein input")
    position = 0
    cds_features = []
    cdsnames = []
    for sequence in sequences:
        startpos = position
        endpos = position + len(sequence) * 3
        position += len(sequence) * 3 + 1000
        location = FeatureLocation(startpos, endpos)
        cdsfeature = SeqFeature(location, type="CDS")
        cdsfeature.strand = 1
        sequence_id = sequence.id[:15].replace(" ","_")
        if sequence_id not in cdsnames:
            cdsfeature.qualifiers['product'] = [sequence_id]
            cdsfeature.qualifiers['locus_tag'] = [sequence_id]
            cdsnames.append(sequence_id)
        else:
            x = 1
            while sequence_id[:8] + "_" + str(x) in cdsnames:
                x += 1
            cdsfeature.qualifiers['product'] = [sequence_id[:8] + "_" + str(x)]
            cdsfeature.qualifiers['locus_tag'] = [sequence_id[:8] + "_" + str(x)]
            cdsnames.append(sequence_id[:8] + "_" + str(x))
        cdsfeature.qualifiers['translation'] = [str(sequence.seq).replace('.', 'X')]
        cds_features.append(cdsfeature)
    seq_record.features.extend(cds_features)
    return [seq_record]


def add_translations(seq_records):
    "Add a translation qualifier to all CDS features"
    for seq_record in seq_records:
        cdsfeatures = utils.get_cds_features(seq_record)
        for cdsfeature in cdsfeatures:
            if 'translation' not in cdsfeature.qualifiers or len(cdsfeature.qualifiers['translation']) == 0:
                if len(seq_record.seq) == 0:
                    logging.error('No amino acid sequence in input entry for CDS %r, ' \
                            'and no nucleotide sequence provided to translate it from.', cdsfeature.id)
                    sys.exit(1)
                else:
                    import Bio.Data.CodonTable
                    try:
                        translation = str(utils.get_aa_translation(seq_record, cdsfeature))
                    except Bio.Data.CodonTable.TranslationError as e:
                        logging.error('Getting amino acid sequences from %s, CDS %r failed: %s',
                                seq_record.name, cdsfeature.id, e)
                        sys.exit(1)
                    cdsfeature.qualifiers['translation'] = [translation]


def add_seq_record_seq(seq_records):
    for seq_record in seq_records:
        if len(seq_record.seq) == 0:
            seqmax = max([cds.location.start for cds in utils.get_cds_features(seq_record)] + [cds.location.end for cds in utils.get_cds_features(seq_record)])
            seq_record.seq = Seq(seqmax * "n")


def check_duplicate_gene_ids(sequences):
    "Fix duplicate locus tags so that they are different"
    NO_TAG = "no_tag_found"
    high_water_mark = 0
    all_ids = defaultdict(lambda: False)
    for sequence in sequences:
        seq_ids = utils.get_cds_features(sequence)
        for cdsfeature in seq_ids:
            gene_id = utils.get_gene_id(cdsfeature)
            if not all_ids[gene_id]:
                all_ids[gene_id] = True
            else:
                if gene_id == NO_TAG:
                    x = high_water_mark + 1
                else:
                    x = 1
                id_str = "%s_%s" % ( gene_id[:8], x)
                while all_ids[id_str]:
                    x += 1
                    id_str = "%s_%s" % ( gene_id[:8], x)
                logging.debug("generated id %r", id_str)
                cdsfeature.qualifiers['product'] = [id_str]
                cdsfeature.qualifiers['locus_tag'] = [id_str]
                all_ids[id_str] = True
                if gene_id == NO_TAG:
                    high_water_mark = x


def fix_id_lines(options, filename):

    with open(filename, 'r') as fh:
        file_content = fh.read()

    while file_content[:2] == "\n\n":
        file_content = file_content[1:]

    try:
        filetype = seqio._get_seqtype_from_ext(filename)
    except ValueError as e:
        logging.debug("Got ValueError %s trying to guess the file format", e)
        return filename

    needed_fixing = False

    if filetype == "fasta":
        if len([line for line in file_content.split("\n") if (len(line) > 0 and line[0] == ">" and " " in line)]) > 0:
            lines = []
            for line in file_content.split("\n"):
                if (len(line) > 0 and line[0] == ">" and " " in line):
                    lines.append(line.replace(" ", ""))
                else:
                    lines.append(line)
            file_content = "\n".join(lines)
            needed_fixing = True
    elif filetype == "genbank":
        if "LOCUS       " not in file_content.partition("\n")[0]:
            file_content = "LOCUS       A01                    0 bp        DNA              BCT 01-JAN-2000\n" + file_content
            needed_fixing = True
    elif filetype == "embl":
        if "ID   " not in file_content.partition("\n")[0]:
            file_content = "ID   A01; SV 1; linear; unassigned DNA; STD; PRO; 0 BP.\nXX\n" + file_content
            needed_fixing = True

    if needed_fixing:
        relative_name = path.basename(filename)
        name, ext = path.splitext(relative_name)
        new_name = path.join(options.full_outputfolder_path, "{}_fixed{}".format(name, ext))
        filename = new_name
        with open(filename, 'w') as fh:
            fh.write(file_content)

    return filename


def parse_input_sequences(options):
    "Parse the input sequences from given filename"
    filenames = options.sequences
    logging.info('Parsing input sequences %r', filenames)

    sequences = []
    for filename in filenames:

        if not path.exists(filename):
            logging.error('No sequence file found at %r', filename)
            sys.exit(1)

        infile = open(filename, "r")
        file_content = infile.read()
        infile.close()
        if "Resource temporarily unavailable" in file_content[:200] or \
                "<h1>Server Error</h1>" in file_content[:500] or \
                "NCBI - WWW Error" in file_content[:500]:
            logging.error('ERROR: NCBI server temporarily unavailable: downloading %s failed.', os.path.basename(filename))
            sys.exit(1)

        if options.fix_id_line:
            filename = fix_id_lines(options, filename)

        try:
            record_list = list(seqio.parse(filename))
            if len(record_list) == 0:
                logging.error('No sequence in file %r', filename)
            sequences.extend(record_list)
        except ValueError as e:
            logging.error('Parsing %r failed: %s', filename, e)
            sys.exit(1)
        except AssertionError as e:
            logging.error('Parsing %r failed: %s', filename, e)
            sys.exit(1)
        except Exception as e:
            logging.error('Parsing %r failed with unhandled exception: %s',
                          filename, e)
            sys.exit(1)
    #Check if seq_records have appropriate content
    i = 0
    new_seqs = []
    for sequence in sequences:
        cleanseq = str(sequence.seq).replace("-","")
        cleanseq = cleanseq.replace(":","")
        sequence.seq = Seq(cleanseq)
        #Check if seq_record has either a sequence or has at least 80% of CDS features with 'translation' qualifier
        cdsfeatures = utils.get_cds_features(sequence)
        cdsfeatures_with_translations = [cdsfeature for cdsfeature in cdsfeatures if 'translation' in cdsfeature.qualifiers]
        if len(sequence.seq) == 0 or (
                options.input_type == 'nucl' and \
                len(str(sequence.seq).replace("N","")) == 0 and \
                len(cdsfeatures_with_translations) < 0.8 * len(cdsfeatures)):
            logging.error("Record %s has no sequence, skipping.", sequence.id)
            continue

        if options.input_type == 'prot':
            if is_nucl_seq(sequence.seq):
                logging.error("Record %s is a nucleotide record, skipping.", sequence.id)
                continue
        elif options.input_type == 'nucl':
            if not isinstance(sequence.seq.alphabet, NucleotideAlphabet) and not is_nucl_seq(sequence.seq):
                logging.error("Record %s is a protein record, skipping.", sequence.id)
                continue
            if sequence.seq.alphabet != generic_dna:
                sequence.seq.alphabet = generic_dna

        new_seqs.append(sequence)

    sequences = new_seqs

    #If protein input, convert all protein seq_records to one nucleotide seq_record
    if options.input_type == 'prot':
        sequences = generate_nucl_seq_record(sequences)

    #Handle WGS master or supercontig entries
    sequences = process_wgs_master_scaffolds(sequences)

    # Make sure we don't waste weeks of runtime on huge records, unless requested by the user
    old_len = len(sequences)
    if options.limit > -1:
        # sort seq_record arrays before applying --limit, prioritize larger sequences, re-sort by acc# afterwards
        sequences = sorted(sequences, key = lambda sequence: len(sequence.seq), reverse = True)
        sequences = sequences[:options.limit]
        sequences = sorted(sequences, key = lambda sequence: sequence.id)

    new_len = len(sequences)
    if new_len < old_len:
        options.triggered_limit = True
        logging.warning("Only analysing the first %d records (increase via --limit)" % options.limit)

    #Check if no duplicate locus tags / gene IDs are found
    check_duplicate_gene_ids(sequences)

    #If no CDS entries in records, run gene finding
    options.record_idx = 1


    #Store IDs for all entries
    options.all_record_ids = [seq.id for seq in sequences]

    # For retaining the correct contig numbers, a second counter is required also including removed sequences without genes
    options.orig_record_idx =1

    # Check GFF suitability
    gff_parser.check_gff_suitability(options,sequences)

    # Remove contigs < 100 bp without annotations
    new_seqs = []
    for sequence in sequences:
        if len(utils.get_cds_features(sequence)) < 1 and len(sequence.seq) < 100:
            continue
        new_seqs.append(sequence)
    sequences = new_seqs

    # Concatenate contigs < 100kbp that are having no annotations and run genefinding at once
    if options.genefinding != 'none':  # changed
        new_seqs = []
        concated_seq_loc = []
        seq_type = generic_dna
        num_concated = 0
        if options.input_type == 'prot':
            seq_type = generic_protein
        concat_seq_text = ""
        for sequence in sequences:
            concated_seq_loc.append((-1, -1))
            if len(utils.get_cds_features(sequence)) < 1:
                if options.gff3:
                    if sequence.id in options.gff_ids: # sequence is covered in provided gff, skip genefinding
                        continue
                if len(concat_seq_text) > 1:
                    concat_seq_text += "TGA-TGA--TGA" # 3-frames stop codon to prevent cross-contig gene prediction
                    for ix in range(0, 30): # create gaps to separate between contigs
                        concat_seq_text += "-"
                s_pos = len(concat_seq_text)
                concat_seq_text += "%s" % sequence.seq
                e_pos = len(concat_seq_text) - 1
                concated_seq_loc[-1] =  (s_pos, e_pos) # take note of the contig's start and end position
                num_concated += 1
        if num_concated > 0:
            concat_seq = SeqRecord(Seq(concat_seq_text, seq_type))
            logging.info("Running genefinding on %s Contigs without annotations.." % num_concated)
            genefinding.find_genes(concat_seq, options)
            for det_gene in utils.get_cds_features(concat_seq): # pass the annotations back into the sequences
                gene_start = det_gene.location.start
                gene_end = det_gene.location.end
                for ci in range(0, len(concated_seq_loc)):
                    s_pos, e_pos = concated_seq_loc[ci]
                    if (e_pos >= gene_start >= s_pos) and (e_pos >= gene_end >= s_pos):
                        det_gene.location = FeatureLocation((gene_start - s_pos), (e_pos - gene_end))
                        sequences[ci].features.append(det_gene)
                        break
            del concat_seq
        del concated_seq_loc

    # re-filter the resulting sequences, discard checked contigs without genes
    new_seqs = []
    for sequence in sequences:
        if len(utils.get_cds_features(sequence)) < 1 and ((not options.gff3) or (sequence.id not in options.gff_ids)):
            continue
        new_seqs.append(sequence)
    if ((len(sequences) > len(new_seqs))):
        logging.info("Discarding %s contigs without annotations after gene finding.." % (len(sequences) - len(new_seqs)))
    sequences = new_seqs

    i = 0
    new_seqs = []
    for sequence in sequences:
        #Fix sequence name (will be ID) if it contains illegal chars
        illegal_chars  = '''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ '''
        for char in sequence.name:
            if char in illegal_chars:
                sequence.name = sequence.name.replace(char, "_")
        #Iterate through sequence objects
        if len(utils.get_cds_features(sequence)) < 1:
            if options.gff3:
                logging.info("No CDS features found in record %r but GFF3 file provided, running GFF parser.", sequence.id)
                gff_parser.run(sequence, options)
                check_duplicate_gene_ids(sequences)
            if len(utils.get_cds_features(sequence)) < 1:
                logging.info("No genes found, skipping record")
                options.orig_record_idx += 1
                continue
        #Fix locus tags
        utils.fix_locus_tags(sequence, options)
        options.record_idx += 1
        options.orig_record_idx += 1
        new_seqs.append(sequence)
    sequences = new_seqs

    #Fix sequence record IDs to be unique
    ids_used = []
    for sequence in sequences:
        seq_id = sequence.id
        if seq_id not in ids_used:
            ids_used.append(seq_id)
        else:
            x = 0
            #Make sure the length of the ID does not exceed 16
            if len(seq_id) <= 12:
                while "%s_%i" % (seq_id, x) in ids_used:
                    x += 1
                sequence.id = "%s_%i" % (seq_id, x)
            else:
                while "%s_%i" % (seq_id[:-4], x) in ids_used:
                    x += 1
                sequence.id = "%s_%i" % (seq_id[:-4], x)
            ids_used.append(sequence.id)
            options.all_record_ids.append(sequence.id) #Update all_record_ids with new record

    #Make sure that all CDS entries in all seq_records have translation tags, otherwise add them
    add_translations(sequences)

    #Do hmm_search in advance, since it will be used by multiple modules beside hmm_detection
    logging.info("Running hmmsearch..")
    schunk_size = options.hmmsearch_chunk # will do hmmsearch per-n sequences
    cur_count = 0
    cur_total = 0
    total_count = 0
    #hmmdetails = [line.split("\t") for line in open(utils.get_full_path(__file__, path.join("antismash", "generic_modules", "hmm_detection", "hmmdetails.txt")),"r").read().split("\n") if line.count("\t") == 3]
    _signature_profiles = [(sig.name, sig.description, sig.cutoff, sig.path) for sig in hmm_detection.get_sig_profiles()]
    full_fasta = ""
    results_by_id = {}
    for seq_record in sequences:
        total_count += len(utils.get_cds_features(seq_record))
    for seq_record in sequences:
        farray = utils.get_cds_features(seq_record)
        for feature in farray:
            prefix = "%s:" % seq_record.id.replace(":", "_")
            gene_id = prefix + utils.get_gene_id(feature)
            fasta_seq = feature.qualifiers['translation'][0]
            if "-" in str(fasta_seq):
                fasta_seq = Seq(str(fasta_seq).replace("-",""), generic_protein)
            if len(fasta_seq) == 0:
                continue
            if cur_count < schunk_size:
                full_fasta += "\n>%s\n%s" % (gene_id, fasta_seq)
            else:
                logging.debug("Running hmmsearch for %s sequences X %s models (%s/%s)" % (cur_count, len(_signature_profiles), cur_total, total_count))
                for sig in _signature_profiles:
                    runresults = utils.run_hmmsearch(utils.get_full_path(__file__, path.join("antismash", "generic_modules", "hmm_detection", sig[3])), full_fasta, sig[2])
                    for runresult in runresults:
                        #Store result if it is above cut-off
                        for hsp in runresult.hsps:
                            if hsp.bitscore > sig[2]:
                                if len(sig[0].split("/")) > 1:
                                    hsp.query_id = sig[0].split("/")[0] + "/" + hsp.query_id
                                if hsp.hit_id not in results_by_id:
                                    results_by_id[hsp.hit_id] = [hsp]
                                else:
                                    results_by_id[hsp.hit_id].append(hsp)
                full_fasta = ">%s\n%s" % (gene_id, fasta_seq)
                cur_count = 0
            cur_count += 1
            cur_total += 1
    if len(full_fasta) > 0:
        logging.debug("Running hmmsearch for %s sequences X %s models (%s/%s)" % (cur_count, len(_signature_profiles), cur_total, total_count))
        for sig in _signature_profiles:
            runresults = utils.run_hmmsearch(utils.get_full_path(__file__, path.join("antismash", "generic_modules", "hmm_detection", sig[3])), full_fasta, sig[2])
            for runresult in runresults:
                #Store result if it is above cut-off
                for hsp in runresult.hsps:
                    if hsp.bitscore > sig[2]:
                        if len(sig[0].split("/")) > 1:
                            hsp.query_id = sig[0].split("/")[0] + "/" + hsp.query_id
                        if hsp.hit_id not in results_by_id:
                            results_by_id[hsp.hit_id] = [hsp]
                        else:
                            results_by_id[hsp.hit_id].append(hsp)
    for hit_id in results_by_id:
        for hsp in results_by_id[hit_id]:
            hsp.hit_id = hsp.hit_id.split(":", 1)[1]
    options.hmm_results = results_by_id
    del results_by_id

    #Remove records without potential hits, checked by HMMer - saves time for inputs with many contigs
    contigs_hits = check_signature_gene_presence(sequences, options)
    new_sequences = []
    for idx in range(len(sequences)):
        if len(contigs_hits[idx]) > 0:
            new_sequences.append(sequences[idx])
    del contigs_hits
    logging.info("Removing %i out of %i contigs without any pHMM hits" % (len(sequences) - len(new_sequences), len(sequences)))
    sequences = new_sequences

    #Make sure that all seq_records have a sequence
    add_seq_record_seq(sequences)

    if len(sequences) > 1:
        options.start = -1
        options.end = -1
        logging.info("Discarding --from and --to options, as multiple entries are used.")

    i = 0
    while i < len(sequences):
        sequence = sequences[i]

        if options.start > 1:
            if options.start > len(sequence):
                logging.error('Specified analysis start point is at %r, which is larger ' \
                              'than record size %r', options.start, len(sequence))
                sys.exit(1)
            sequence = sequence[options.start-1:]
            # new sequence is shorter, so fix the end calculation
            options.end -= options.start
            sequences[i] = sequence

        if options.end > 0:
            if options.end > len(sequence):
                logging.error('Specified analysis end point is at %r, which is larger ' \
                              'than record size %r', options.end, len(sequence))
                sys.exit(1)
            sequence = sequence[:options.end]
            sequences[i] = sequence

        # Some programs write gaps as - not N, but translate() hates that
        if sequence.seq.find('-') > -1:
            sequence.seq = Seq(str(sequence.seq).replace('-', 'N'),
                               alphabet=sequence.seq.alphabet)

        # Some programs like to write gaps as X, translate() hates that
        if sequence.seq.find('X') > -1:
            sequence.seq = Seq(str(sequence.seq).replace('X', 'N'),
                               alphabet=sequence.seq.alphabet)
        if sequence.seq.find('x') > -1:
            sequence.seq = Seq(str(sequence.seq).replace('x', 'N'),
                               alphabet=sequence.seq.alphabet)

        i += 1

    return sequences


def check_prereqs(plugins, options):
    failure_messages = []
    failure_messages.extend(generic_check_prereqs(options))
    failure_messages.extend(ggm_check_prereqs(options))
    for plugin in plugins:
        if 'check_prereqs' in dir(plugin):
            failure_messages.extend(plugin.check_prereqs())

    for msg in failure_messages:
        logging.error(msg)

    return len(failure_messages)


def detect_signature_genes(seq_record, clustertypes, options):
    "Detect different secondary metabolite clusters based on HMM signatures"
    logging.info('Looking for secondary metabolite cluster signatures')
    hmm_detection.detect_signature_genes(seq_record, clustertypes, options)


def detect_geneclusters(seq_record, options):
    if options.input_type == 'nucl':
        utils.log_status("Detecting secondary metabolite signature genes " \
                         "for contig #%d" % options.record_idx)
    else:
        utils.log_status("Detecting secondary metabolite signature genes " \
                         "for supplied amino acid sequences")
    detect_signature_genes(seq_record, options.enabled_cluster_types, options)
    if options.inclusive:
        utils.log_status("Detecting secondary metabolite clusters using "\
                         "inclusive ClusterFinder algorithm for contig #%d" % options.record_idx)
        fullhmmer.run(seq_record, options)
        clusterfinder.run_cluster_finder(seq_record, options)
        options.full_hmmer = False


def cluster_specific_analysis(plugins, seq_record, options):
    if options.disable_specific_modules:
        logging.debug("Skipping cluster-specific analyses because --disable_specific_modules is set.")
        return

    logging.info('Running cluster-specific analyses')

    for plugin in plugins:
        if options.taxon == "plants" and not plugin.name.startswith("plant_"):
            continue
        if hasattr(plugin, "specific_analysis"):
            logging.debug('Running analyses specific to %s clusters', plugin.short_description)
            plugin.specific_analysis(seq_record, options)
        else:
            logging.debug('No specific analyses implemented for %s clusters', plugin.short_description)


def run_generic_genome_modules(seq_records, options):
    "Run genome wide analysis modules"
    logging.info('Running genome-wide analysis modules')

    plugins = load_generic_genome_plugins()

    for plugin in plugins:
        logging.debug("Executing genome-wide plugin: %s", plugin.name)
        plugin.run_analyses(seq_records, options)

def unspecific_analysis(seq_record, options):
    "Run analyses independent of specific clusters"
    logging.info('Running analyses independent of specific cluster types')

    if options.full_hmmer:
        utils.log_status("Running full-genome PFAM analysis for contig #%d" % \
                         options.record_idx)
        fullhmmer.run(seq_record, options)

    if not options.ecpred == 'none':
        utils.log_status("Running EC prediction")
        ecpredictor.run(seq_record, options)
    # TODO: offer fullblast here


def check_signature_gene_presence(seq_records, options):
    "For each contig, check if contig contains signature genes, in order to discard empty ones"
    logging.info('Checking gene clusters for presence of profile hits using HMM library')
    #Identify contigs with HMM hits
    contigs_hits = []
    for idx in range(len(seq_records)):
        contigs_hits.append([])
        seq_record = seq_records[idx]
        for feature in utils.get_cds_features(seq_record):
            prefix = "%s:" % seq_record.id.replace(":", "_")
            gene_id = utils.get_gene_id(feature)
            if (prefix + gene_id) in options.hmm_results:
                for hsp in options.hmm_results[prefix + gene_id]:
                    if hsp.query_id not in contigs_hits[-1]:
                        contigs_hits[-1].append(hsp.query_id)
    return contigs_hits


if __name__ == "__main__":
    main()
