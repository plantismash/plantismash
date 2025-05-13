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
Module for handling database access

This library contains the generic antiSMASH database class aSDB which handles all communication with the DBI
Basically, it is a wrapper around the Biopython BioSeqDatabase object
"""
import logging
import sys
from BioSQL import BioSeqDatabase

# Although the database import currently is not implemented as a plugin, we define the required parameters
name = "BioSQLImporter"
short_description = "BioSQL importer to retrieve seq_records from BioSQL database"
priority = 1


class aSDB(object):

    def __init__(self, options):
        """
        constructor method to setup object and connect to the database

         This object requires a functional BioSQL setup
         Configuration in default.cfg, section BioSQLconfig

         Demo config:

         [BioSQLconfig]
         name = BioSQL
         DBdriver = psycopg2
         DBuser = biosql
         DBpass = biosql
         DBhost = localhost
         DBport = 5432
         DBdb = biosql
        """
        # Database namespace must be provided via options.dbnamespace

        if "BioSQLconfig" not in options:
            logging.exception("Parameters for database access not defined in default.cfg. Bailing out...")
            sys.exit(1) # FIXME: Library functions should never call sys.exit()
        try:
            self.server = BioSeqDatabase.open_database(driver=options.BioSQLconfig.dbdriver,
                                                      user=options.BioSQLconfig.dbuser,
                                                      passwd=options.BioSQLconfig.dbpass,
                                                      host=options.BioSQLconfig.dbhost,
                                                      port=options.BioSQLconfig.dbport,
                                                      db=options.BioSQLconfig.dbdb)
        except Exception as e:
            logging.exception("Error setting up the database object; Check parameters in default.cfg")
            raise e
        self.db = {}


    def connect(self, namespace):
        """This method connects to a BioSQL namespace

        For convenience, the method returns the database adaptor object. But as this
        obect is also stored within the aSDB object, normally shouldn't be used directly.
        """

        if namespace not in self.db:
            try:
                self.db[namespace] = self.server[namespace]
            except Exception as e:
                logging.exception("Could not assign BioSQL namespace/database %s, probably namespace does not exist", namespace)
                raise e
        self.namespace = namespace
        return(self.db[namespace])


    def close(self):
        "Close connection to the database server"
        self.server.close()
        del self


    def commit(self):
        "Commit changes to SQL server"
        if 'namespace' in self.db:
            self.server[self.namespace].commit()
        else:
            self.server.commit()

    def rollback(self):
        "Rollback changes to SQL database"
        if 'namespace' in self.db:
            self.server[self.namespace].rollback()
        else:
            self.server.rollback()

    def delete(self, entryid):
        "Delete record from BioSQL database using its entry_id"

        try:
            del self.db[self.namespace][entryid]
            logging.debug("Delete record %s from database", entryid)
        except Exception as e:
            logging.warn("Could not delete record %s from database/namespace %s", entryid, self.namespace)
            raise e


    def fetch_entryid_by_name(self, name):
        "get database id from display_id"

        try:
            namespace_id = self.db[self.namespace].adaptor.fetch_dbid_by_dbname(self.namespace)
            dbid=self.db[self.namespace].adaptor.fetch_seqid_by_display_id(dbid=namespace_id, name=name)
        except Exception as e:
            logging.warn("%s: Could not find database id by searching for %s",  e, name)
            dbid = None
        return dbid


    def generate_namespace(self, namespace, description):
        try:
            self.db[namespace] = self.server.new_database(db_name=namespace, description=description)
            self.namespace = namespace
            # self.commit()
            logging.debug("namespace %s generated for storing genome seq_records", namespace)
        except Exception as e:
            logging.exception("Could not generate namespace/db %s", namespace)
            raise e


    def get_current_namespace(self):
        """"Get the name of last connected namespace.

        Note: aSDB object can simultaneously connect to multiple namespaces, this method only returns the name of the last connected namespace
        """
        return self.namespace


    def get_record_by_name(self, name):
        "Get a database DBSeqRecord object from database by selected by name"

        try:
            seq_record = self.db[self.namespace].lookup(name=name)
        except Exception as e:
            logging.warn("%s: Database entry with name %s does not exist", e, name)
            seq_record = None
        return(seq_record)


    def load_records(self, seq_records):
        "Load sequence records into current namespace of database"

        count = self.db[self.namespace].load(seq_records)
        return count


    def switch_namespace(self, new_namespace):
        if new_namespace not in self.db:
            logging.exception("cannot switch to namespace %s; namespace in not connected", new_namespace)
            sys.exit(1)
        self.namespace = new_namespace

# Maybe we should move this function to a own module or the main script as it is not directly related to the aSDB object
def get_record(seq_id, options):
    myDB = aSDB(options)
    myDB.connect(namespace=options.dbnamespace)
    # OK, now let's retrieve the sequence record
    seq_record = myDB.get_record_by_name(seq_id)
    if not seq_record:
        logging.exception("Database entry with id %s does not exist", seq_id)
        sys.exit(1)
    logging.debug("Retrieved entry %s from database", seq_record.name)
    return seq_record

def get_records(options):
    "get list of DBSeqRecords from database"

#    myDB = aSDB(options)
#    myDB.connect(namespace=options.dbnamespace)
    seq_records=[]

    for seq_id in options.seq_ids:#
        seq_record = get_record(seq_id, options)
        # OK, now let's retrieve the sequence record
#        seq_record = myDB.get_record_by_name(seq_id)
#        if not seq_record:
#            logging.exception("Database entry with id %s does not exist", seq_id)
#            sys.exit(1)
        seq_records.append(seq_record)

    return seq_records

