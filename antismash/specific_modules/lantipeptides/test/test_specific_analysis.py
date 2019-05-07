try:
    import unittest2 as unittest
except ImportError:
    import unittest
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqFeature import SeqFeature, FeatureLocation
from minimock import mock, restore, TraceTracker
from antismash import utils
from antismash.specific_modules.lantipeptides.specific_analysis import (
    Lantipeptide,
    predict_cleavage_site,
    result_vec_to_features,
)

class TestLantipeptide(unittest.TestCase):
    def test_init(self):
        "Test Lantipeptide instantiation"
        lant = Lantipeptide(23, 42, 17, 'Class-I')
        self.assertTrue(isinstance(lant, Lantipeptide))
        self.assertEqual(23, lant.start)
        self.assertEqual(42, lant.end)
        self.assertEqual(17, lant.score)
        self.assertEqual('Class-I', lant.lantype)
        self.assertEqual('', lant.core)
        with self.assertRaises(ValueError):
            lant.molecular_weight

    def test_repr(self):
        "Test Lantipeptide representation"
        lant = Lantipeptide(23, 42, 17, 'Class-I')
        expected = "Lantipeptide(23..42, 17, 'Class-I', '', -1, -1(-1))"
        self.assertEqual(expected, repr(lant))

    def test_core(self):
        "Test Lantipeptide.core"
        lant = Lantipeptide(23, 42, 17, 'Class-I')
        self.assertEqual('', lant.core)
        self.assertFalse(hasattr(lant, 'core_analysis'))
        lant.core = "MAGICHAT"
        self.assertEqual('MAGICHAT', lant.core)
        self.assertTrue(hasattr(lant, 'core_analysis'))

    def test_core_ignore_invalid(self):
        "Test Lantipeptide.core ignores invalid amino acids"
        lant = Lantipeptide(23, 42, 17, 'Class-I')
        self.assertEqual('', lant.core)
        self.assertFalse(hasattr(lant, 'core_analysis'))
        lant.core = "MAGICXHAT"
        self.assertEqual('MAGICHAT', lant.core)
        self.assertTrue(hasattr(lant, 'core_analysis'))

    def test_number_of_lan_bridges(self):
        "Test Lantipeptide.number_of_lan_bridges"
        lant = Lantipeptide(23, 42, 17, 'Class-I')
        lant.core = "MAGICHAT"
        self.assertEqual(-1, lant._lan_bridges)
        self.assertEqual(1, lant.number_of_lan_bridges)
        lant._lan_bridges = 2
        self.assertEqual(2, lant._lan_bridges)
        self.assertEqual(2, lant.number_of_lan_bridges)

    def test_monoisotopic_mass(self):
        "Test Lantipeptide.monoisotopic_mass"
        lant = Lantipeptide(23, 42, 17, 'Class-I')
        lant.core = "MAGICHAT"
        analysis = ProteinAnalysis("MAGICHAT", monoisotopic=True)
        mw = analysis.molecular_weight()
        # Thr is assumed to be dehydrated
        mw -= 18
        self.assertAlmostEqual(mw, lant.monoisotopic_mass)
        self.assertAlmostEqual(mw, lant._monoisotopic_weight)

        lant._weight = 42
        self.assertEqual(42, lant.molecular_weight)

    def test_molecular_weight(self):
        "Test Lantipeptide.molecular_weight"
        lant = Lantipeptide(23, 42, 17, 'Class-I')
        lant.core = "MAGICHAT"
        analysis = ProteinAnalysis("MAGICHAT", monoisotopic=False)
        mw = analysis.molecular_weight()
        # Thr is assumed to be dehydrated
        mw -= 18.02
        self.assertAlmostEqual(mw, lant.molecular_weight)
        self.assertAlmostEqual(mw, lant._weight)

        lant._weight = 42
        self.assertEqual(42, lant.molecular_weight)

    def test_alternative_weights(self):
        "Test Lantipeptide.alt_weights"
        lant = Lantipeptide(23, 42, 17, 'Class-I')
        lant.core = "MAGICHATS"
        analysis = ProteinAnalysis("MAGICHATS", monoisotopic=False)
        mw = analysis.molecular_weight()
        # One Ser/Thr is assumed to be dehydrated, but not the other
        mw -= 18.02
        self.assertEqual([mw], lant.alternative_weights)


class TestSpecificAnalysis(unittest.TestCase):
    class FakeHit(object):
        class FakeHsp(object):
            def __init__(self, start, end, score):
                self.query_start = start
                self.query_end = end
                self.bitscore = score

        def __init__(self, start, end, score, desc):
            self.hsps = [ self.FakeHsp(start, end, score) ]
            self.description = desc

        def __iter__(self):
            return iter(self.hsps)

    def setUp(self):
        self.tt = TraceTracker()
        self.hmmpfam_return_vals = []
        mock('utils.run_hmmpfam2', tracker=self.tt, returns=self.hmmpfam_return_vals)

    def tearDown(self):
        restore()

    def test_predict_cleavage_site(self):
        "Test lantipeptides.predict_cleavage_site()"
        resvec = predict_cleavage_site('foo', 'bar')
        self.assertEqual(None, resvec)
        fake_hit = self.FakeHit(24, 42, 17, 'fake')
        self.hmmpfam_return_vals.append([fake_hit])

        res = predict_cleavage_site('foo', 'bar')

        self.assertEqual(23, res.start)
        self.assertEqual(42, res.end)
        self.assertEqual(17, res.score)
        self.assertEqual('fake', res.lantype)

    def test_result_vec_to_features(self):
        "Test lantipeptides.result_vec_to_features()"
        loc = FeatureLocation(0, 165, strand=1)
        orig_feature = SeqFeature(loc, 'CDS')
        orig_feature.qualifiers['locus_tag'] = ['FAKE0001']
        vec = Lantipeptide(17, 23, 42, 'Class-I')
        seq = "TAILTAILTAILTAILTAILTAILTAILTAILTAILCC"
        vec.core = seq
        vec.leader = "HEADHEADHEAD"
        new_features = result_vec_to_features(orig_feature, vec)
        self.assertEqual(2, len(new_features))
        leader, core = new_features

        self.assertEqual(loc.start, leader.location.start)
        self.assertEqual(loc.start + (23 * 3), leader.location.end)
        self.assertEqual(loc.strand, leader.location.strand)
        self.assertEqual('CDS_motif', leader.type)
        self.assertEqual(orig_feature.qualifiers['locus_tag'],
                         leader.qualifiers['locus_tag'])
        self.assertEqual(['leader peptide',
                          'predicted leader seq: HEADHEADHEAD'], leader.qualifiers['note'])

        self.assertEqual(leader.location.end, core.location.start)
        self.assertEqual(loc.end, core.location.end)
        self.assertEqual(loc.strand, core.location.strand)
        self.assertEqual('CDS_motif', core.type)
        expected = ['core peptide', 'monoisotopic mass: 3646.3',
                    'molecular weight: 3648.6',
                    'alternative weights: 3666.6; 3684.6; 3702.7; 3720.7; 3738.7; 3756.7; 3774.7',
                    'number of bridges: 2',
                    'predicted core seq: TAILTAILTAILTAILTAILTAILTAILTAILTAILCC',
                    'predicted class: Class-I',
                    'score: 42.00'
                   ]
        self.assertEqual(expected, core.qualifiers['note'])
        self.assertEqual(orig_feature.qualifiers['locus_tag'],
                         core.qualifiers['locus_tag'])
