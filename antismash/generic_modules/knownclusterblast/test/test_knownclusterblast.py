import codecs
import unittest
from os import path
import antismash.generic_modules.knownclusterblast as kcb

class KnownClusterBlastTest(unittest.TestCase):
    def test_no_unicode_in_knownclusters(self):
        """Test knownclusters.txt doesn't contain unicode codepoints"""
        knownclusters_path = path.join(path.dirname(kcb.__file__), "knownclusters.txt")
        self.assertTrue(path.exists(knownclusters_path))
        try:
            fh = codecs.open(knownclusters_path, 'r', 'ascii')
            fh.read()
            fh.close()
        except ValueError:
            with codecs.open(knownclusters_path, 'r', 'utf-8') as fh:
                lines = fh.readlines()
            counter = 1
            for line in lines:
                try:
                    line.encode('ascii')
                    counter += 1
                except UnicodeEncodeError:
                    assert False, "Non-ASCII char in line %s: %r" % (counter, line)
