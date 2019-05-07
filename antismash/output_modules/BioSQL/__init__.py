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
# Copyright (C) 2014 Tilmann Weber
# The Novo Nordisk Foundation Center for Biosustainability
# Technical University of Denmark
# Section: Metabolic Engineering for Natural Compounds / New Bioactive Compounds
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.


"""
BioSQL export module
"""

import logging
import warnings
import sys
import re
try:
    from psycopg2 import IntegrityError
    from psycopg2.extensions import TransactionRollbackError
    have_psycopg2 = True
except ImportError:
    have_psycopg2 = False
from antismash import utils
from antismash.db import biosql
from antismash.db import extradata
from time import strftime, gmtime, sleep

name = "BioSQL"
short_description = "export to BioSQL database"
priority = 1


# This plugin requires a functional BioSQL setup
# Configuration in default.cfg, section BioSQLconfig
#
# Demo config:

# [BioSQLconfig]
# name = BioSQL
# DBdriver = psycopg2
# DBuser = biosql
# DBpass = biosql
# DBhost = localhost
# DBport = 5432
# DBdb = biosql
# DBnamespace = antismash


def write(seq_records, options):
    if not options.use_db:
        logging.debug('--db option not set; skipping BioSQL export')
        return 1
    elif not have_psycopg2:
        logging.debug('failed to import psycopg2; skipping BioSQL export')
        return 1
    else:
        logging.debug("starting DB export...")

    if options.input_type == 'nucl':



        if not "BioSQLconfig" in options:
            logging.warning("Parameters for database access not defined in default.cfg. Skipping database export")

        # Set up database object
        myDB = biosql.aSDB(options)

        # connect to namespace for full genomes (dbgenomesnamespace)
        try:
            myDB.connect(namespace=options.BioSQLconfig.dbgenomenamespace)
        except Exception as e:
            logging.warning("Could not assign BioSQL namespace/database %s for storing genome seq_records, probably namespace does not exist.\nError message: %s" , options.BioSQLconfig.dbgenomenamespace, e)

            try:
                myDB.generate_namespace(namespace = options.BioSQLconfig.dbgenomenamespace,
                                        description = "namespace for storing antiSMASH genome results")
                myDB.commit()
            except Exception as e:
                logging.exception("Could not generate namespace/db %s for storing genome seq_records in database %s\nError message: %s", options.BioSQLconfig.dbgenomenamespace, options.BioSQLconfig.dbdb, e)
                sys.exit(1)

        for seq_record in seq_records:
            # if entry was extracted from database we can simply continue without updating
            if options.from_database:
                logging.warn("Seq record %s was retrieved from database; skipping loading to BioSQL", seq_record.name)
                return
            # As the biopyhton BioSQL-adapter doesn't allow updates of records we first have to delete records with the same accession_number
            # first delete genome entry
            # myDB.switch_namespace(options.BioSQLconfig.dbgenomenamespace)
            entryid = myDB.fetch_entryid_by_name(name=seq_record.name)
            if entryid:
                logging.debug("deleting existing database entry with ID: %s", entryid)
                myDB.delete(entryid)
            if 'extrarecord' in options:
                if seq_record.id in options.extrarecord:
                    if'extradata' in options.extrarecord[seq_record.id]:
                        logging.debug("Uploading %s rows to extraTable", len(options.extrarecord[seq_record.id].extradata.keys()))
                        loadExtraData(seq_record.id, options.extrarecord[seq_record.id].extradata, options)
            else:
                logging.debug("No extradata available, skipping")

        # Store genome seq_record
        # myDB.switch_namespace(options.BioSQLconfig.dbgenomenamespace)

        dbretries = 1

        # try to load entries to BioSQL database
        while dbretries <= int(options.BioSQLconfig.dbretries):
            try:
                logging.debug("Try {0} of {1}: Writing genome seq_records to database {2}".format(dbretries, options.BioSQLconfig.dbretries, options.BioSQLconfig.dbdb))
                count = myDB.load_records(seq_records)
                break
            except (IntegrityError, TransactionRollbackError) as err:
                logging.error('{0}: {1}'.format(err.diag.severity, err.diag.message_primary))
                logging.error(err.diag.message_detail)
                logging.error(err.diag.message_hint)
                logging.error('Error on try {0} of {1}: sleeping for 30 sec, rolling back transaction and trying again'.format(dbretries, options.BioSQLconfig.dbretries))
                sleep(30)
                myDB.rollback()
                dbretries += 1

        if dbretries > int(options.BioSQLconfig.dbretries):
            logging.exception('ERROR: Failed to write genome record to database!')
            sys.exit(1)

        logging.debug("commiting %d changes to database", count)
        myDB.commit()


    else:
        logging.error("Protein data import to database not implemented")
        sys.exit(1)


def loadExtraData(seqrecord_id, extradataHash, options):
    "upload extra data"
    logging.debug("writing extra data to database")

    myExtraDB = extradata.extradata(options)
    
    #before uploading the new data, all existing data for genbank_id should be deleted
    if myExtraDB.check_SeqID(seqrecord_id, "%"):
        logging.debug("Deleting extradata for %s", seqrecord_id)
        myExtraDB.delete(seqrecord_id, "%")
    
    dbresult = myExtraDB.putExtradataHash(seqrecord_id,extradataHash)
    myExtraDB.close()
    return dbresult
