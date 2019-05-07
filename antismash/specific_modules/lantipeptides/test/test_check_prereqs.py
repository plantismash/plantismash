try:
    import unittest2
except ImportError:
    import unittest as unittest2

from minimock import Mock, mock, restore, TraceTracker, assert_same_trace
import antismash.specific_modules.lantipeptides
from antismash.specific_modules.lantipeptides import check_prereqs

class TestCheckPrereqs(unittest2.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.tt = TraceTracker()
        self.locate_exe = Mock('antismash.utils.locate_executable',
                               tracker=self.tt, returns="/fake/path/to/binary")
        mock('antismash.utils.locate_executable',
             mock_obj=self.locate_exe, tracker=self.tt)

    def tearDown(self):
        restore()

    def test_check_prereqs(self):
        "Test lantipeptides.check_prereqs()"
        ret = check_prereqs()
        self.assertEqual(ret, [])
        expected = "    Called antismash.utils.locate_executable('hmmpfam2')"
        assert_same_trace(self.tt, expected)

    def test_check_binary_prereqs_failing(self):
        "Test lantipeptides.check_prereqs() returns 'missing binary' error"
        self.locate_exe.mock_returns = None
        ret = check_prereqs()
        self.assertEqual(len(ret), 1)
        self.assertIn("Failed to locate executable for 'hmmpfam2'", ret)
