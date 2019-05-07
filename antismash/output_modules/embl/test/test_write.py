try:
    import unittest2
except ImportError:
    import unittest as unittest2
import helperlibs.bio.seqio
from minimock import mock, restore, TraceTracker, assert_same_trace, Mock
import antismash.output_modules.embl as embl

class TestEmblWrite(unittest2.TestCase):
    def setUp(self):
        self.tt = TraceTracker()
        mock('helperlibs.bio.seqio.write', tracker=self.tt)
        self.options = Mock('options', tracker=self.tt, outputfoldername='test', input_type="nucl")
        record1 = Mock('record1', tracker=self.tt, id='record1')
        record2 = Mock('record2', tracker=self.tt, id='record2')
        self.records = [record1, record2]

    def tearDown(self):
        restore()

    def test_write(self):
        "Test output_modules.embl.write()"
        embl.write(self.records, self.options)
        expected = "Called helperlibs.bio.seqio.write([<Mock ... record1>, <Mock ... record2>], 'test/record1.final.embl', 'embl')"
        assert_same_trace(self.tt, expected)

    def test_write_abspath(self):
        "Test output_modules.embl.write() with absolute path outputfoldername"
        self.options.outputfoldername = "/long/test"
        embl.write(self.records, self.options)
        expected = "Called helperlibs.bio.seqio.write([<Mock ... record1>, <Mock ... record2>], '/long/test/record1.final.embl', 'embl')"
        assert_same_trace(self.tt, expected)
