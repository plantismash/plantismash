"""Integration test running on the balhimycin cluster"""

import sys
import unittest
from os import path
from helperlibs.wrappers.io import TemporaryDirectory
from helperlibs.bio import seqio
from minimock import mock, restore, TraceTracker, assert_same_trace
from run_antismash import main

class IntTestExit(Exception):
    pass

class IntegrationBalh(unittest.TestCase):
    def setUp(self):
        self.trace_tracker = TraceTracker()
        def fake_exit(code):
            raise IntTestExit(code)
        mock("sys.exit", tt=self.trace_tracker, returns_func=fake_exit)

    def tearDown(self):
        restore()

    def integration_run(self):
        """Run a sanity check on the Balhimycin cluster"""
        infile = path.join(path.dirname(__file__), 'Y16952.gbk')
        with TemporaryDirectory() as tempdir:
            argv = [
                'run_antismash.py',
                '--outputfolder', tempdir,
                infile,
            ]
            sys.argv = argv
            main()
            assert_same_trace(self.trace_tracker, "")
            outfile = path.join(tempdir, 'Y16952.3.final.gbk')
            self.assertTrue(path.exists(outfile), "Failed to create ouput")
            rec = seqio.read(outfile)
            self.assertIsNotNone(rec, "Failed to parse output")
