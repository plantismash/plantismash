try:
    import unittest2
except ImportError:
    import unittest as unittest2
import helperlibs.bio.seqio
from minimock import mock, restore, TraceTracker, assert_same_trace, Mock
import antismash.output_modules.genbank as gb

class TestGenbankWrite(unittest2.TestCase):
    def setUp(self):
        self.tt = TraceTracker()
        mock('helperlibs.bio.seqio.write', tracker=self.tt)
        self.options = Mock('options', tracker=self.tt, outputfoldername='test', input_type='nucl')
        record1 = Mock('record1', tracker=self.tt, id='record1', features=[])
        record2 = Mock('record2', tracker=self.tt, id='record2', features=[])
        self.records = [record1, record2]

    def tearDown(self):
        restore()

    def test_write(self):
        "Test output_modules.genbank.write()"
        gb.write(self.records, self.options)
        expected = "Called helperlibs.bio.seqio.write([<Mock ... record1>, <Mock ... record2>], 'test/record1.final.gbk', 'genbank')"
        assert_same_trace(self.tt, expected)

    def test_write_abspath(self):
        "Test output_modules.genbank.write() with absolute path output folder name"
        self.options.outputfoldername = "/test"
        gb.write(self.records, self.options)
        expected = "Called helperlibs.bio.seqio.write([<Mock ... record1>, <Mock ... record2>], '/test/record1.final.gbk', 'genbank')"
        assert_same_trace(self.tt, expected)
