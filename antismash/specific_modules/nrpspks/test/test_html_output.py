#set fileencoding: utf-8
try: import unittest2
except ImportError: import unittest as unittest2

from minimock import mock, restore, TraceTracker, assert_same_trace
from antismash import utils
from antismash.specific_modules.nrpspks import html_output as ho

class TestHtml(unittest2.TestCase):
    def setUp(self):
        self.tt = TraceTracker()
        mock('utils.get_aa_sequence', tracker=self.tt, returns="FAKESEQ")


    def tearDown(self):
        restore()


    def test__parse_substrate_predictions(self):
        """Test _parse_substrate_predictions()"""
        domain = "NRPS/PKS Domain: PKS_AT (565-854). E-value: " \
                 "2.2e-103. Score: 336.9; Substrate specificity predictions: " \
                 "emal (PKS signature), mmal (Minowa), xmal (consensus);"
        expected = [('PKS signature', 'emal'), ('Minowa', 'mmal'),
                    ('consensus', 'xmal')]
        ret = ho._parse_substrate_predictions(domain)
        self.assertEqual(expected, ret)


    def test__parse_substrate_predictions_no_specificity(self):
        """Test _parse_substrate_predictions() without specificity string"""
        domain = "Nothing to see here, move along"
        expected = []
        ret = ho._parse_substrate_predictions(domain)
        self.assertEqual(expected, ret)


    def test__get_monomer_prediction(self):
        """Test _get_monomer_prediction()"""
        class FakeFeature(object): pass
        feature = FakeFeature()
        feature.qualifiers = dict(note=['Nothing to see here',
                                  'Monomers prediction: foo bar baz'])
        expected = "foo bar baz"
        ret = ho._get_monomer_prediction(feature)
        self.assertEqual(expected, ret)


    def test__get_monomer_prediction_no_monomers(self):
        """Test _get_monomer_prediction() without monomers"""
        class FakeFeature(object): pass
        feature = FakeFeature()
        feature.qualifiers = dict(note=['Nothing to see here'])
        expected = "N/A"
        ret = ho._get_monomer_prediction(feature)
        self.assertEqual(expected, ret)
