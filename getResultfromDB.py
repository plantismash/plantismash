#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2014 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011-2014 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
#
# Copyright (C) 2014 Tilmann Weber
# The Novo Nordisk Foundation Center for Biosustainability
# Technical University of Denmark
# Section: Metabolic Engineering for Natural Compounds / New Bioactive Compounds
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""Get entry from database and generate antiSMASH output files from entry"""

import sys
import os

if sys.platform ==  ('win32') or sys.platform == ('darwin'):
    os.environ['EXEC'] = os.getcwd() + "\exec"
    os.environ['PYTHON'] = os.getcwd() + "\python"
    sys.path.append(os.getcwd() + "\python\Lib\site-packages")
    os.environ['PATH'] = os.pathsep + os.environ['PYTHON'] + os.pathsep + os.environ['PATH']
    os.environ['PATH'] = os.pathsep + os.environ['EXEC'] + os.pathsep + os.environ['PATH']

import random
import logging
import argparse
from os import path
import straight.plugin
from antismash.config import load_config, set_config
from antismash import utils
from antismash.db.biosql import get_records
from antismash.db.extradata import getExtradata
from antismash.generic_modules import hmm_detection
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO



def main():
    "Retrieve antiSMASH entry from database"


    # First load the output plugins so we can present appropriate options
    output_plugins = load_output_plugins()


    parser = argparse.ArgumentParser(description='Retrieve entry from database')
    parser.add_argument('seq_ids', metavar='seq_ids', nargs="*",
                        help="accession numbers of antiSMASH-DB entries")

    parser.add_argument('-d', '--debug', dest='debug',
                        action='store_true', default=False,
                        help="Print debugging information to stderr")
    parser.add_argument('--list-plugins', dest='list_plugins',
                        action='store_true', default=False,
                        help="List all available sec. met. detection modules")
    parser.add_argument('-v', '--verbose', dest='verbose',
                        action='store_true', default=False,
                        help="Print verbose status information to stderr")
    parser.add_argument('--logfile', dest='logfile',
                        default=argparse.SUPPRESS,
                        help="Also write logging output to a file")
    parser.add_argument('--statusfile', dest='statusfile',
                        default=argparse.SUPPRESS,
                        help="Write the current status to a file")

    group = parser.add_argument_group('Output options')
    for plugin in output_plugins:
        group.add_argument('--disable-%s' % plugin.name, dest=plugin.name,
                           action='store_false', default=argparse.SUPPRESS,
                           help="Disable %s" % plugin.short_description)

    group = parser.add_argument_group('Settings')
    group.add_argument('--outputfolder', dest='outputfoldername',
                        default=argparse.SUPPRESS,
                        help="Directory to write results to")
    group.add_argument('--dbnamespace', dest='dbnamespace',
                       help="Define BioSQL namespace to search")
    
    group.add_argument('--nclusters', dest='nclusters',
                        default=10, type=int,
                        help="Number of clusters from ClusterBlast to display")
    group.add_argument('--seed', dest='seed', default=0, type=int,
                        help="Random number seed for ClusterBlast coloring")
    options = parser.parse_args()
    
    # Logging is useful for all the following code, so make sure that is set up
    # right after parsing the arguments.
    setup_logging(options)

    if options.nclusters > 50:
        logging.info("Number of clusters (" + str(options.nclusters) + ") is too large. Reducing to 50.")
        options.nclusters = 50
    logging.debug("Number of clusters to show in clusterblast = " + str(options.nclusters))
    if options.seed != 0:
        random.seed(options.seed)
        
    # Load list of clutertypes
    clustertypes = hmm_detection.get_supported_cluster_types()
    
    # Manually set some opions that are required for working with the same codebase as run_antismash.ph
    options.input_type="nucl"
    
    # Note: the clusterblast/subclusterblast options are automatically activated if aSstorage object with data is found
    options.clusterblast=None
    options.subclusterblast=None
    options.knownclusterblast=None
    options.smcogs="TRUE"
    options.modeling = "none"
    options.enabled_cluster_types = ValidateClusterTypes(clustertypes)


    
    #Load configuration data from config file
    load_config(options)
    set_config(options)

    
    
    #Load and filter plugins
    utils.log_status("Loading detection plugins")
    plugins = load_detection_plugins()
    filter_plugins(plugins, options, clustertypes)
    filter_outputs(output_plugins, options)
    
    options.plugins = plugins
    
    if options.list_plugins:
        list_available_plugins(output_plugins)
        sys.exit(0)

    filter_outputs(output_plugins, options)

    #Check prerequisites
    if not options.seq_ids:
        parser.error("Please specify at least one antiSMASH-DB accession number")

    if not 'outputfoldername' in options:
        options.outputfoldername = path.splitext(path.basename(options.sequences[0]))[0]
    if not os.path.exists(options.outputfoldername):
        os.mkdir(options.outputfoldername)
    options.full_outputfolder_path = path.abspath(options.outputfoldername)

    if not options.dbnamespace in [options.BioSQLconfig.dbgenomenamespace, options.BioSQLconfig.dbclusternamespace]:

        logging.warn("DBnamespace %s not defined in default.cfg, switching to standard namespace %s." % (options.dbnamespace, options.BioSQLconfig.dbgenomenamespace))
        options.dbnamespace = options.BioSQLconfig.dbgenomenamespace

    #Parse input sequence
    try:
        utils.log_status("retrieving record")
        seq_records = get_records(options)
    except:
        logging.exception("Uncaptured error when reading entries from antiSMASH-DB. This should not have happened :-(")
        sys.exit(1)
    
    options.extrarecord = {}
    
    for seq_record in seq_records:
        options.extrarecord[seq_record.id] = argparse.Namespace()
        
        logging.debug("DB retrieval: trying to find extra data for %s" % seq_record.id)
        extradataHash = getExtradata(options, seq_record.id)
        logging.debug("Keys of extradataHash: %s" % ", ".join(extradataHash.keys()))
        options.extrarecord[seq_record.id].extradata = extradataHash
        
        if options.extrarecord[seq_record.id].extradata.has_key('ClusterBlastData'):
            logging.debug("DB retrieval: Found extra data for ClusterBlast")
            options.clusterblast = True
        if options.extrarecord[seq_record.id].extradata.has_key('SubClusterBlastData'):
            logging.debug("DB retrieval: Found extra data for SubClusterBlast")
            options.subclusterblast = True
        if options.extrarecord[seq_record.id].extradata.has_key('KnownClusterBlastData'):
            logging.debug("DB retrieval: Found extra data for KnownClusterBlast")
            options.knownclusterblast = True
        if options.extrarecord[seq_record.id].extradata.has_key('MetabolicModel'):
            logging.debug("DB retrieval: Found extra data for Modeling")
            options.modeling = "db"

    #Write results
    utils.log_status("Writing the output files")
    write_results(output_plugins, seq_records, options)
    zip_results(seq_records, options)




def list_available_plugins(output_plugins):
    print("Support for the following output formats:")
    for plugin in output_plugins:
        print(" * %s" % plugin.short_description)

def filter_plugins(plugins, options, clustertypes):
    if options.enabled_cluster_types is None or options.enabled_cluster_types == clustertypes:
        return

    for plugin in plugins:
        if plugin.name in clustertypes and plugin.name not in options.enabled_cluster_types:
            plugins.remove(plugin)

    if plugins == []:
        print("No plugins enabled, use --list-plugins to show available plugins")
        sys.exit(1)
        
def load_detection_plugins():
    "Load available secondary metabolite detection modules"
    logging.info('Loading detection modules')
    return straight.plugin.load('antismash.specific_modules')


def filter_outputs(plugins, options):
    for plugin in plugins:
        if plugin.name in options:
            logging.debug("Removing plugin %r" % plugin.name)
            plugins.remove(plugin)

    if plugins == []:
        print("No plugins enabled, use --list-plugins to show available plugins")
        sys.exit(1)

def ValidateClusterTypes(clustertypes):
    class Validator(argparse.Action):
        def __call__(self, parser, args, values, option_string = None):
            values = values.replace(";",",").split(",")
            try:
                for value in values:
                    if value not in clustertypes:
                        raise ValueError('invalid clustertype {s!r}'.format(s = value))
            except ValueError as e:
                print "\nInput error:", e, "\n"
                print "Please choose from the following list:\n", "\n".join(clustertypes), "\n\nExample: --enable t1pks,nrps,other"
                sys.exit(1)
            setattr(args, self.dest, values)
    return Validator


def write_results(plugins, seq_records, options):
    for plugin in plugins:
        plugin.write(seq_records, options)


def zip_results(seq_records, options):
    "Create a zip archive with all the results generated so far"
    zip_name = '%s.zip' % seq_records[0].id
    utils.zip(path.abspath(options.outputfoldername), zip_name)


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
        if not os.path.exists(path.dirname(options.logfile)):
            os.mkdir(path.dirname(options.logfile))
        fh = logging.FileHandler(options.logfile)
        fh.setLevel(logging.INFO)
        fh.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        logging.getLogger('').addHandler(fh)


def load_output_plugins():
    "Load available output formats"
    plugins = list(straight.plugin.load('antismash.output_modules'))
    plugins.sort(cmp=lambda x, y: cmp(x.priority, y.priority))
    
    # We have to remove the BioSQL exporter from the output plugins-List
    for plugin in plugins:
        if plugin.name == 'BioSQL':
            plugins.remove(plugin)
    return plugins



if __name__ == "__main__":
    main()
