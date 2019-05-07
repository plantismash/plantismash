"""
Test suite for antismash.db module
"""

try:
    import unittest2
except ImportError:
    import unittest as unittest2

from os import path
from argparse import Namespace
from minimock import mock, restore, TraceTracker, assert_same_trace
from Bio import SeqIO
from BioSQL import BioSeqDatabase
skip_db_tests = False
try:
    from antismash.db import biosql as db
except ImportError:
    skip_db_tests = True

@unittest2.skipIf(skip_db_tests, "Skipping DB-related tests")
class TestBasicDBconnection(unittest2.TestCase):
    "Test connection to database"
    def setUp(self):
        "set up database test framework and test for db connection"
        self.options = Namespace()
        self.options.BioSQLconfig = Namespace()
        self.options.BioSQLconfig.dbdriver = 'psycopg2'
        self.options.BioSQLconfig.dbuser = 'biosql'
        self.options.BioSQLconfig.dbpass = 'biosql'
        self.options.BioSQLconfig.dbhost = 'localhost'
        self.options.BioSQLconfig.dbport = '5432'
        self.options.BioSQLconfig.dbdb = 'antiSMASHnosetest'
        self.options.BioSQLnamespace = 'test'

        self.trace_tracker = TraceTracker()

        # test correct method call for db.asDB)()
        mock('BioSeqDatabase.open_database', tracker=self.trace_tracker, returns=['mock result'])

        expected = """Called BioSeqDatabase.open_database(
        db='antiSMASHnosetest',
        driver='psycopg2',
        host='localhost',
        passwd='biosql',
        port='5432',
        user='biosql')
        """

        mydb = db.aSDB(self.options)
        assert_same_trace(self.trace_tracker, expected)

        # restore mocked objects that other tests can run with the original code
        restore()
        self.trace_tracker = TraceTracker()

        #test initializing database object db.__init__"

        try:
            mydb = db.aSDB(self.options)
        except Exception:
            self.skipTest('database connection could not be established, skipping tests')

        self.assertIsInstance(mydb, db.aSDB)


        self.mydb = mydb

    def test_generate_namespace(self):
        "test generate_namespace() method"

        mydb = db.aSDB(self.options)

        mock('BioSeqDatabase.DBServer.new_database', tracker=self.trace_tracker, returns=["mock result"])
        expected = """Called BioSeqDatabase.DBServer.new_database(
        db_name='test',
        description='test namespace for nosetests')
        """

        mydb.generate_namespace(self.options.BioSQLnamespace, 'test namespace for nosetests')
        assert_same_trace(self.trace_tracker, expected)


class TestDBAccess(unittest2.TestCase):
    "Test access to database"

    def setUp(self):
        # restore all mocked methods
        restore()

        #Open testfile and parse it with Biopython
        testfile = path.join(path.dirname(__file__), "AAA26493.embl")
        testfile_fh = open(testfile, "rU")
        self.seq_record = SeqIO.read(testfile_fh, 'embl')

        self.options = Namespace()
        self.options.BioSQLconfig = Namespace()
        self.options.BioSQLconfig.dbdriver = 'psycopg2'
        self.options.BioSQLconfig.dbuser = 'biosql'
        self.options.BioSQLconfig.dbpass = 'biosql'
        self.options.BioSQLconfig.dbhost = 'localhost'
        self.options.BioSQLconfig.dbport = '5432'
        self.options.BioSQLconfig.dbdb = 'antiSMASHnosetest'
        self.options.BioSQLnamespace = 'test'

        self.trace_tracker = TraceTracker()

        try:
            mydb = db.aSDB(self.options)
        except Exception:
            self.skipTest('database connection could not be established, skipping tests')

        try:
            del mydb.server[self.options.BioSQLnamespace]
            mydb.commit()
        except Exception:
            pass

        mydb.generate_namespace(self.options.BioSQLnamespace, 'test namespace for nosetests')
        mydb.commit()

        mydb.connect(self.options.BioSQLnamespace)
        self.mydb = mydb

    def tearDown(self):
        self.mydb.close()

    def test_get_current_namespace(self):
        "test get_current_namespace() method"

        mydb = self.mydb

        mynamespace = mydb.get_current_namespace()
        self.assertEqual(mynamespace, self.options.BioSQLnamespace)

    def test_connect(self):
        "test database connection with connect() and get_current_namespace() method"

        mydb = self.mydb
        # get list of namespaces in DB and check for our test namespace
        dblist = mydb.db[self.options.BioSQLnamespace].adaptor.list_biodatabase_names()
        self.assertTrue(self.options.BioSQLnamespace in dblist)

        self.assertEqual(mydb.get_current_namespace(), self.options.BioSQLnamespace)



    def test_commit(self):
        "test commit method"

        mydb = self.mydb

        mock('BioSeqDatabase.DBServer.commit', tracker=self.trace_tracker, returns=["mock result"])

        mydb.commit()
        expected = "Called BioSeqDatabase.DBServer.commit()"
        assert_same_trace(self.trace_tracker, expected)
        
        
    def test_rollback(self):
        "test rollback method"
        
        mydb = self.mydb
        
        mock('BioSeqDatabase.DBServer.rollback', tracker=self.trace_tracker, returns=["mock result"])
        
        mydb.rollback()
        expected = 'Called BioSeqDatabase.DBServer.rollback()'
        assert_same_trace(self.trace_tracker, expected)


    def test_database_operations(self):
        "Test for methods: load_records(), fetch_entryid_by_name(), get_record_by_name(); delete()"

        # needs to be monolithic to ensure that record is loaded to DB

        mydb = self.mydb

        mock('BioSeqDatabase.BioSeqDatabase.load', tracker=self.trace_tracker, returns=["mock result"])

        rec_no = mydb.load_records([self.seq_record])
        self.assertListEqual(rec_no, ['mock result'])

        expected = """Called BioSeqDatabase.BioSeqDatabase.load(
        [SeqRecord(seq=Seq('GTGTCGGGGCCGCGCTCGCGCACGACGAGCAGGCGGACGCCGGTCCGCATCGGC...TGA', IUPACAmbiguousDNA()), id='AAA26493.1', name='AAA26493', description='Saccharopolyspora erythraea EryA', dbxrefs=[])])
        """
        assert_same_trace(self.trace_tracker, expected)
        restore()
        self.trace_tracker = TraceTracker()

        # now do the actual loading...
        rec_no = mydb.load_records([self.seq_record])
        self.assertEqual(rec_no, 1)

        commit_res = mydb.commit()
        self.assertIsNone(commit_res)

        mock('BioSeqDatabase.Adaptor.fetch_dbid_by_dbname', tracker=self.trace_tracker, returns=1234)
        mock('BioSeqDatabase.Adaptor.fetch_seqid_by_display_id', tracker=self.trace_tracker, returns=4321)

        db_seq_rec_id = mydb.fetch_entryid_by_name('testname')


        expected = """Called BioSeqDatabase.Adaptor.fetch_dbid_by_dbname('test')
Called BioSeqDatabase.Adaptor.fetch_seqid_by_display_id(
    dbid=1234,
        name='testname')
        """

        assert_same_trace(self.trace_tracker, expected)
        restore()
        self.trace_tracker = TraceTracker()

        # database record does not exist, should deliver None
        db_seq_rec_id = mydb.fetch_entryid_by_name('test')
        self.assertIsNone(db_seq_rec_id)

        db_seq_rec_id = mydb.fetch_entryid_by_name('AAA26493')
        self.assertIsInstance(db_seq_rec_id, int)

        # retrieve entry from database

        mock('BioSeqDatabase.BioSeqDatabase.lookup', tracker=self.trace_tracker, returns=['mock results'])
        myrecord = mydb.get_record_by_name('AAA26493')

        self.assertListEqual(myrecord, ['mock results'])

        expected = "Called BioSeqDatabase.BioSeqDatabase.lookup(name='AAA26493')"
        assert_same_trace(self.trace_tracker, expected)
        restore()
        self.trace_tracker = TraceTracker()

        myrecord = mydb.get_record_by_name('AAA26493')
        myrecord_class = myrecord.__class__.__name__

        self.assertEqual(myrecord_class, "DBSeqRecord")

        self.assertEqual(str(myrecord.seq), str(self.seq_record.seq))
        self.assertEqual(myrecord.id, self.seq_record.id)
        self.assertEqual(myrecord.name, self.seq_record.name)

        myrecord = mydb.get_record_by_name('Idonotexist')
        self.assertIsNone(myrecord)

        # test delete method
        mock('BioSeqDatabase.BioSeqDatabase.__delitem__', tracker=self.trace_tracker, returns=['mock results'])

        mydb.delete(db_seq_rec_id)
        expected = "Called BioSeqDatabase.BioSeqDatabase.__delitem__(...)"
        assert_same_trace(self.trace_tracker, expected)
        restore()

        mydb.delete(db_seq_rec_id)
        db_seq_rec_id = mydb.fetch_entryid_by_name('AAA26493')
        self.assertIsNone(db_seq_rec_id, 'Entry %s still present in database' % db_seq_rec_id)

class Test_get_functions(unittest2.TestCase):
    "Test get_record and get_records functions of antismash.db module"

    def setUp(self):
        "set up database connection"
        # restore all mocked methods
        restore()

        #Open testfile and parse it with Biopython
        testfile = path.join(path.dirname(__file__), "AAA26493.embl")
        testfile_fh = open(testfile, "rU")
        self.seq_record = SeqIO.read(testfile_fh, 'embl')

        self.options = Namespace()
        self.options.BioSQLconfig = Namespace()
        self.options.BioSQLconfig.dbdriver = 'psycopg2'
        self.options.BioSQLconfig.dbuser = 'biosql'
        self.options.BioSQLconfig.dbpass = 'biosql'
        self.options.BioSQLconfig.dbhost = 'localhost'
        self.options.BioSQLconfig.dbport = '5432'
        self.options.BioSQLconfig.dbdb = 'antiSMASHnosetest'
        self.options.BioSQLnamespace = 'test'
        self.options.dbnamespace = self.options.BioSQLnamespace
        self.trace_tracker = TraceTracker()

        try:
            mydb = db.aSDB(self.options)
        except Exception:
            self.skipTest('database connection could not be established, skipping tests')

        if not mydb:
            self.skipTest('database connection impossible; skipping tests')

        try:
            del mydb.server[self.options.BioSQLnamespace]
            mydb.commit()
        except Exception:
            pass

        mydb.generate_namespace(self.options.BioSQLnamespace, 'test namespace for nosetests')

        mydb.connect(self.options.BioSQLnamespace)

        rec_no = mydb.load_records([self.seq_record])
        mydb.commit()
        self.mydb = mydb

    def tearDown(self):
        "lose database connection"
        self.mydb.close()

    def test_get_record(self):
        "Test get_record() function"

        Fake_seq_rec = Namespace()
        Fake_seq_rec.name = 'fakedName'

        mock('db.aSDB.get_record_by_name', tracker=self.trace_tracker, returns=Fake_seq_rec)

        db_seq_rec = db.get_record('AAA26493', self.options)
        expected = "Called db.aSDB.get_record_by_name('AAA26493')"
        assert_same_trace(self.trace_tracker, expected)
        restore()
        self.trace_tracker = TraceTracker()

        db_seq_rec = db.get_record('AAA26493', self.options)
        myrecord_class = db_seq_rec.__class__.__name__

        self.assertEqual(myrecord_class, "DBSeqRecord")
        self.assertEqual(str(db_seq_rec.seq), str(self.seq_record.seq))
        self.assertEqual(db_seq_rec.id, self.seq_record.id)
        self.assertEqual(db_seq_rec.name, self.seq_record.name)

    def test_get_records(self):
        "Test get_records() function"

        self.options.seq_ids = ['AAA26493', 'AAA26493']

        mydbrecords = db.get_records(self.options)
        self.assertEqual(len(mydbrecords), 2)
        self.assertEqual(mydbrecords[0].name, self.seq_record.name)
