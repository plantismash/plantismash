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
database access module to read/write extra data, e.g. clusterblastresults, models,...
"""

import logging

try:
    import psycopg2
    got_psycopg2 = True
except ImportError:
    got_psycopg2 = False

    
import sys

try:
    import pickle as pickle
except ImportError:
    import pickle


name = 'hmmdb'
short_description ="database access module to read/write extra data, e.g. clusterblastresults, models,..."
priority = 1
extradataTableName = "extradata"

class extradata(object):

    def __init__(self, options):
        # Just cause the ImportError to be raised again
        
        if not got_psycopg2:
            raise ImportError("psycopg2 not imported")
       
        
        if "BioSQLconfig" not in options:
            logging.exception("Parameters for database access not defined in default.cfg. Bailing out...")
            sys.exit(1)

        try:
            dbconn = psycopg2.connect(database=options.BioSQLconfig.dbdb,
                                      user=options.BioSQLconfig.dbuser,
                                      password=options.BioSQLconfig.dbpass,
                                      host=options.BioSQLconfig.dbhost,
                                      port=options.BioSQLconfig.dbport)
        except Exception as e:
            logging.exception("failed to connect to database!")
            raise e

        self.dbconn = dbconn
        self.extraschema = options.BioSQLconfig.extraschema
        self.extradataTablename = extradataTableName
        self.options=options



    def close(self):
        try:
            del self.dbconn
        except Exception as e:
            logging.warning("could not close database connection: %s!", e)

    def check_SeqID(self, genbank_id, secondary_id):
        "check if there are entries for genbank_id"

        db_cur = self.dbconn.cursor()

        SQL = """SELECT
                    1
                 FROM
                    {schema}.{table}
                 WHERE
                    {schema}.{table}.genbank_id = %s AND
                    {schema}.{table}.secondary_id LIKE %s
                 LIMIT 1;""".format(schema=self.extraschema, table=self.extradataTablename)

        db_cur.execute(SQL, (genbank_id, secondary_id,))
        # I have no idea why this query results in a (1,) tupel and not in [(1)]
        if db_cur.fetchall() == [(1,)]:
            return(True)
        else:
            return(False)

    def delete(self, genbank_id, secondary_id):
        "delete entry(row) in table"

        db_cur = self.dbconn.cursor()

        SQL = """DELETE FROM
                    {schema}.{table}
                 WHERE
                    {schema}.{table}.genbank_id = %s AND
                    {schema}.{table}.secondary_id LIKE %s;""".format(schema=self.extraschema, table=self.extradataTablename)

        try:
            db_cur.execute(SQL, (genbank_id, secondary_id,))
        except psycopg2.Error as e:
            logging.warning("Could not delete record with genbank_id, secondary_id %s, %s from database\n" \
                            "psycopg error is: %s\n", genbank_id, secondary_id, e)

        if not db_cur.rowcount == 1:
            logging.debug("psycopg2 delete operation returned %s deleted rows", db_cur.rowcount)
        self.commit()


    def commit(self):
        "Commit changes to SQL server"

        self.dbconn.commit()

    def rollback(self):
        "Rollback changes to SQL database"

        self.dbconn.rollback()


    def insert(self, genbank_id, secondary_id, string):
        "Insert row into table"

        db_cur = self.dbconn.cursor()
        status = 'failed'

        SQL = """INSERT INTO
                    {schema}.{table} (genbank_id, secondary_id, asdata)
                 VALUES
                    (%s, %s, %s);""".format(schema=self.extraschema, table=self.extradataTablename)
        try:
            db_cur.execute(SQL, (genbank_id, secondary_id, string,))
            status = 'success'
        except psycopg2.Error as e:
            logging.warning("Could not insert record with genbank_id %s\npsycopg error is: %s\n", genbank_id, e)

        try:
            self.commit()
        except psycopg2.Error as e:
            logging.warning("Error commiting insert into db")
            status = 'failed'

        return status

    def update(self, genbank_id, secondary_id, string):
        "Insert row into table"

        db_cur = self.dbconn.cursor()
        status='failed'

        SQL = """UPDATE
                    {schema}.{table}
                SET
                    genbank_id = %s, secondary_id = %s, asdata = %s
                WHERE
                    {schema}.{table}.genbank_id = %s AND
                    {schema}.{table}.secondary_id = %s;""".format(schema=self.extraschema, table=self.extradataTablename)
        try:
            db_cur.execute(SQL, (genbank_id, secondary_id, string, genbank_id, secondary_id,))
            status = "success"
        except psycopg2.Error as e:
            logging.warning("Could not insert record with genbank_id %s\npsycopg error is: %s\n", genbank_id, e)


        try:
            self.commit()
        except psycopg2.Error as e:
            logging.warning("Error commiting insert into db")
            status = 'failed'

        return status

    def insertOrUpdate(self, genbank_id, secondary_id, string):
        "Checks if record exists and then either updates record (if it exists) or inserts a new row"

        status = 'failed'
        if self.check_SeqID(genbank_id, secondary_id):
            status = self.update(genbank_id, secondary_id, string)
        else:
            status = self.insert(genbank_id, secondary_id, string)

        return status

    def fetchOneRow(self, genbank_id, secondary_id):
        "Fetch entry from db"

        db_cur = self.dbconn.cursor()

        SQL = """SELECT
                    genbank_id, secondary_id, asdata
                 FROM
                    {schema}.{table}
                 WHERE
                    {schema}.{table}.genbank_id = %s AND
                    {schema}.{table}.secondary_id = %s;""".format(schema=self.extraschema, table=self.extradataTablename)


        try:
            db_cur.execute(SQL, (genbank_id, secondary_id, ))

        except psycopg2.Error as e:
            logging.warning("Could not fetch row with genbank_id, secondary_id %s, %s from database\n" \
                            "psycopg error is: %s\n" , genbank_id, secondary_id, e)

        if not db_cur.rowcount == 1:
            logging.warning("psycopg2 SELECT operation returned %s rows", db_cur.rowcount)

        # the data is in field [1] in table
        result = db_cur.fetchone()[1]

        return result

    def getExtradataHash(self, genbank_id):
        "Fetch entry from db"

        db_cur = self.dbconn.cursor()

        SQL = """SELECT
                    genbank_id, secondary_id, asdata
                 FROM
                    {schema}.{table}
                 WHERE
                    {schema}.{table}.genbank_id = %s ;""".format(schema=self.extraschema, table=self.extradataTablename)


        try:
            db_cur.execute(SQL, (genbank_id, ))

        except psycopg2.Error as e:
            logging.warning("Could not get ExtradataHash record for genbank_id, %s from database\n" \
                            "psycopg error is: %s\n", genbank_id, e)



        # the data is in field [1] in table
        lines = db_cur.fetchall();
        logging.debug("getExtradataHas retrieved %s lines from DB", len(lines))

        results = {}

        if len(lines) > 0:
            for line in lines:
                results[line[1]] = pickle.loads((line[2]))


        return results

    def putExtradataHash(self, genbank_id, extradataHash):

        status = "failed"
        for key in extradataHash:
            serializedData = pickle.dumps(extradataHash[key])
            # logging.debug("DB-extra: length of serialized data for %s: %s" % (len(serializedData), key))
            if len(serializedData) > 0:
                status = self.insertOrUpdate(genbank_id, key, serializedData)
            self.dbconn.commit()
        return status

def getExtradata(options, genbank_id):
    myExtraDB = extradata(options)
    extradataHash = myExtraDB.getExtradataHash(genbank_id)
    myExtraDB.close()
    return extradataHash
