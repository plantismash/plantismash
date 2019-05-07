try:
    import unittest2
except ImportError:
    import unittest as unittest2
from minimock import mock, restore, TraceTracker, assert_same_trace
from argparse import Namespace
import antismash.generic_modules

class GenericModulesTest(unittest2.TestCase):
    def setUp(self):
        self.tt = TraceTracker()
        self.options = Namespace(full_hmmer=True, clusterblast=True,
                                 subclusterblast=True, knownclusterblast=True,
                                 smcogs=True)
        mock('antismash.generic_modules.fullhmmer.check_prereqs',
             returns=['fullhmmer'], tracker=self.tt)
        mock('antismash.generic_modules.genefinding.check_prereqs',
             returns=['genefinding'], tracker=self.tt)
        mock('antismash.generic_modules.hmm_detection.check_prereqs',
             returns=['hmm_detection'], tracker=self.tt)
        mock('antismash.generic_modules.clusterblast.check_prereqs',
             returns=['clusterblast'], tracker=self.tt)
        mock('antismash.generic_modules.subclusterblast.check_prereqs',
             returns=['subclusterblast'], tracker=self.tt)
        mock('antismash.generic_modules.knownclusterblast.check_prereqs',
             returns=['knownclusterblast'], tracker=self.tt)
        mock('antismash.generic_modules.smcogs.check_prereqs',
             returns=['smcogs'], tracker=self.tt)

    def tearDown(self):
        restore()

    def test_check_prereqs(self):
        "Test antismash.generic_modules.check_prereqs()"
        trace = """    Called antismash.generic_modules.fullhmmer.check_prereqs(Namespace(...))
    Called antismash.generic_modules.genefinding.check_prereqs(Namespace(...))
    Called antismash.generic_modules.hmm_detection.check_prereqs()
    Called antismash.generic_modules.smcogs.check_prereqs(Namespace(...))
    Called antismash.generic_modules.clusterblast.check_prereqs(Namespace(...))
    Called antismash.generic_modules.subclusterblast.check_prereqs(Namespace(...))
    Called antismash.generic_modules.knownclusterblast.check_prereqs(Namespace(...))"""
        expected = ['fullhmmer', 'genefinding', 'hmm_detection',
                    'smcogs', 'clusterblast', 'subclusterblast', 'knownclusterblast']
        returned = antismash.generic_modules.check_prereqs(self.options)

        self.assertListEqual(expected, returned)
        assert_same_trace(self.tt, trace)
