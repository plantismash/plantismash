try:
    import unittest2
except ImportError:
    import unittest as unittest2
from minimock import mock, restore, TraceTracker, assert_same_trace
from argparse import Namespace
import antismash.utils
from antismash.generic_modules import fullhmmer
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

def _create_dummy_record():
    seq = Seq('GTGGAGCGGTACTAAATGTACTCCACTATCTGCTGATTGGAAACCACGGAGCGCTCTTAG',
              generic_dna)

    rec = SeqRecord(seq, id="FAKE")

    locations = (FeatureLocation(0,15, strand=1), FeatureLocation(15,36, strand=1),
                 FeatureLocation(36,60, strand=1))
    idx = 1
    for loc in locations:
        f = SeqFeature(loc, type='CDS')
        f.qualifiers['locus_tag'] = ['orf%04d' % idx]
        rec.features.append(f)
        idx += 1

    return rec


class FakeHSP(object):
    def __init__(self, start, end):
        self.query_start = start
        self.query_end = end


def _create_dummy_results():
    pass



class TestFullhmmer(unittest2.TestCase):
    def setUp(self):
        f_config = Namespace(score='1', evalue='0.02')
        self.config = Namespace(fullhmmer=f_config)
        self.tt = TraceTracker()
        mock('antismash.utils.locate_executable', returns='hmmsearch',
             tracker=self.tt)
        self.file_list = ['Pfam-A.hmm', 'Pfam-A.hmm.h3f', 'Pfam-A.hmm.h3i',
                          'Pfam-A.hmm.h3m', 'Pfam-A.hmm.h3p']
        mock('antismash.utils.locate_file', returns_iter=self.file_list,
             tracker=self.tt)
        mock('antismash.utils.run_hmmscan', returns=[])

    def tearDown(self):
        restore()

    def test_check_prereqs(self):
        "Test fullhmmer.check_prereqs()"
        config = Namespace()
        trace = """    Called antismash.utils.locate_executable('hmmscan')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3f')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3i')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3m')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3p')"""
        expected = []
        returned = fullhmmer.check_prereqs(config)

        self.assertListEqual(expected, returned)
        assert_same_trace(self.tt, trace)


    def test_check_prereqs_missing_exe(self):
        "Test fullhmmer.check_prereqs() with a missing executable"
        config = Namespace()
        antismash.utils.locate_executable.mock_returns = None
        trace = """    Called antismash.utils.locate_executable('hmmscan')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3f')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3i')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3m')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3p')"""
        expected = ["Failed to locate file: 'hmmscan'"]
        returned = fullhmmer.check_prereqs(config)

        self.assertListEqual(expected, returned)
        assert_same_trace(self.tt, trace)


    def test_check_prereqs_missing_file(self):
        "Test fullhmmer.check_prereqs() with a missing file"
        config = Namespace()
        self.file_list[0] = None
        antismash.utils.locate_file.mock_returns_iter = self.file_list
        trace = """    Called antismash.utils.locate_executable('hmmscan')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3f')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3i')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3m')
    Called antismash.utils.locate_file(
        '.../antismash/generic_modules/fullhmmer/Pfam-A.hmm.h3p')"""
        expected = ["Failed to locate file: 'Pfam-A.hmm'"]
        returned = fullhmmer.check_prereqs(config)

        self.assertListEqual(expected, returned)
        assert_same_trace(self.tt, trace)


    def test_min_score(self):
        "Test fullhmmer._min_score()"
        self.assertAlmostEqual(1, fullhmmer._min_score(self.config))
        self.assertAlmostEqual(0, fullhmmer._min_score(Namespace()))

    def test_max_evalue(self):
        "Test fullhmmer._max_evalue()"
        self.assertAlmostEqual(0.02, fullhmmer._max_evalue(self.config))
        self.assertAlmostEqual(0.01, fullhmmer._max_evalue(Namespace()))


    def test_calculate_start_end(self):
        "Test fullhmmer._calculate_start_end()"
        hsp = FakeHSP(1, 5)
        rec = _create_dummy_record()
        feature = rec.features[1]
        start, end = fullhmmer._calculate_start_end(feature, hsp)
        self.assertEqual((18, 30), (start, end))
        self.assertEqual("YSTI", str(rec.seq[start:end].translate()))

        rev = rec.reverse_complement()
        feature = rev.features[1]
        start, end = fullhmmer._calculate_start_end(feature, hsp)
        self.assertEqual((30, 42), (start, end))
        self.assertEqual("YSTI", str(rev.seq[start:end].reverse_complement().translate()))
