try:
    import unittest2
except ImportError:
    import unittest as unittest2
from minimock import mock, restore, Mock, TraceTracker, assert_same_trace
from argparse import Namespace
import antismash.utils
from antismash.generic_modules import clusterfinder
from Bio.SeqFeature import SeqFeature, FeatureLocation

def _make_fake_feature(start, end, probability=None, pfam_id=None, type_=None):
    location = FeatureLocation(start, end)
    feature = SeqFeature(location)
    feature.qualifiers = {'note': [], 'db_xref': []}
    if probability is not None:
        feature.qualifiers['note'].append('ClusterFinder probability: %02.4f' % probability)
    if pfam_id is not None:
        feature.qualifiers['db_xref'].append('PFAM: %s' % pfam_id)
    if type_ is not None:
        feature.type = type_
    return feature

class ClusterFinderTest(unittest2.TestCase):
    def setUp(self):
        self.config = Namespace()
        self.tt = TraceTracker()
        self.fake_features = [
            _make_fake_feature(10, 20, 0.1, 'FAKE007'),
            _make_fake_feature(30, 40, 0.3, 'PF00106'),
            _make_fake_feature(50, 60, 0.4, 'PF00107'),
            _make_fake_feature(60, 70, 0.7, 'PF00109'),
            _make_fake_feature(70, 80, 0.98, 'PF08484'),
            _make_fake_feature(90, 100, 0.8, 'PF02401'),
            _make_fake_feature(100, 110, 0.32, 'PF04369'),
            _make_fake_feature(110, 120, 1.0, 'PF00128'),
            _make_fake_feature(130, 140, 0.2, 'FAKE234'),
            _make_fake_feature(500, 505, pfam_id='FAKE505'),
            _make_fake_feature(1010, 1020, 0.1, 'FAKE007'),
            _make_fake_feature(1030, 1040, 0.3, 'PF00106'),
            _make_fake_feature(1050, 1060, 0.4, 'PF00107'),
            _make_fake_feature(1060, 1070, 0.7, 'PF00109'),
            _make_fake_feature(1070, 1080, 0.98, 'PF08484'),
            _make_fake_feature(1090, 1100, 0.8, 'PF02401'),
            _make_fake_feature(1100, 1110, 0.32, 'PF04369'),
            _make_fake_feature(1110, 1120, 1.0, 'PF00128'),
        ]
        mock('antismash.utils.get_cds_features', tracker=self.tt, returns=self.fake_features)
        self.config.cdsnr = 5
        self.config.cf_prob_thres = 0.6
        self.config.cf_npfams = 5
        self.config.borderpredict = False

    def tearDown(self):
        restore()

    def test_check_prereqs(self):
        "Test clusterfinder.check_prereqs()"
        self.assertEqual([], clusterfinder.check_prereqs())

    def test_compare_feature_locations(self):
        "Test clusterfinder.compare_feature_locations()"
        feature1 = _make_fake_feature(24, 25)
        feature2 = _make_fake_feature(42, 45)

        self.assertEqual(0, clusterfinder.compare_feature_locations(feature1, feature1))
        self.assertEqual(-1, clusterfinder.compare_feature_locations(feature1, feature2))
        self.assertEqual(1, clusterfinder.compare_feature_locations(feature2, feature1))

    def test_find_nr_cds(self):
        "Test clusterfinder.find_nr_cds"
        left = (0, 5)
        newpos, num = clusterfinder.find_nr_cds(left, None)
        self.assertEqual(left, newpos)
        self.assertEqual(0, num)

        right = (150, 160)
        newpos, num = clusterfinder.find_nr_cds(right, None)
        self.assertEqual(right, newpos)
        self.assertEqual(0, num)

        middle = (35, 115)
        newpos, num = clusterfinder.find_nr_cds(middle, None)
        self.assertEqual([30, 120], newpos)
        self.assertEqual(7, num)

        small = (501, 504)
        newpos, num = clusterfinder.find_nr_cds(small, None)
        self.assertEqual([500, 505], newpos)
        self.assertEqual(1, num)

    def test_find_cf_clusters(self):
        "Test clusterfinder.find_cf_clusters()"
        ret = clusterfinder.find_cf_clusters(self.fake_features, None, self.config)
        self.assertEqual(30, ret[0].location.start)
        self.assertEqual(120, ret[0].location.end)
        self.assertEqual(1030, ret[1].location.start)
        self.assertEqual(1120, ret[1].location.end)

    def test_annotate_geneclusters_no_overlaps(self):
        "Test clusterfinder.annotate_geneclusters() without overlaps"
        mock("antismash.utils.get_pfam_features", tracker=self.tt,
             returns=self.fake_features)
        self.config.clusternr_offset = 1
        self.config.next_clusternr = -5
        seq_record = Mock('seq_record', tracker=self.tt)
        seq_record.features = []
        clusterfinder.annotate_geneclusters(seq_record, self.config)

        self.assertEqual(3, self.config.next_clusternr)
        self.assertEqual(2, len(seq_record.features))
        cluster1, cluster2 = seq_record.features

        self.assertIn('probability', cluster1.qualifiers.keys())
        self.assertAlmostEqual(0.6429, float(cluster1.qualifiers['probability'][0]))

        self.assertIn('probability', cluster2.qualifiers.keys())
        self.assertAlmostEqual(0.6429, float(cluster2.qualifiers['probability'][0]))

    def test_annotate_geneclusters_overlaps(self):
        "Test clusterfinder.annotate_geneclusters() with overlaps"
        mock("antismash.utils.get_pfam_features", tracker=self.tt,
             returns=self.fake_features)
        self.config.clusternr_offset = 1
        self.config.next_clusternr = -5
        seq_record = Mock('seq_record', tracker=self.tt)
        seq_record.features = [
            _make_fake_feature(10, 40, type_='cluster'),
            _make_fake_feature(1040, 1050, type_='cluster'),
            _make_fake_feature(110, 400, type_='cluster'),
        ]
        clusterfinder.annotate_geneclusters(seq_record, self.config)

        self.assertEqual(4, self.config.next_clusternr)
        self.assertEqual(3, len(seq_record.features))

        self.assertEqual(10, int(seq_record.features[0].location.start))
        self.assertEqual(120, int(seq_record.features[0].location.end))
        self.assertAlmostEqual(0.6429, float(seq_record.features[0].qualifiers['probability'][0]))
        self.assertEqual(1030, int(seq_record.features[1].location.start))
        self.assertEqual(1120, int(seq_record.features[1].location.end))
        self.assertAlmostEqual(0.6429, float(seq_record.features[1].qualifiers['probability'][0]))
        self.assertEqual(30, int(seq_record.features[2].location.start))
        self.assertEqual(400, int(seq_record.features[2].location.end))
        self.assertAlmostEqual(0.6429, float(seq_record.features[2].qualifiers['probability'][0]))
