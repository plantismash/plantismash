from nose.tools import assert_is_instance
try:
    import unittest2
except ImportError:
    import unittest as unittest2
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
from minimock import mock, restore, TraceTracker, assert_same_trace
from argparse import Namespace
from Bio import SearchIO
from helperlibs.bio import seqio
from antismash import utils
import antismash.generic_modules.active_site_finder


class TestASF(unittest2.TestCase):
    "Test Active Site Finder class"
    
    def setUp(self):
        "set-up required variables, etc; load demo sequence and parse XML"
        self.options = Namespace()
        QualifierTags = Namespace()
        QualifierTags.ASF_scaffold = 'aSASF_scaffold'
        QualifierTags.ASF_choice = 'aSASF_choice'
        QualifierTags.ASF_prediction = 'aSASF_prediction'
        
        result = Namespace()
        result.id = "fullhmmer_oxyB_0001"
        result.hsps = [Namespace]
        result.hsps[0] = Namespace()
        result.hsps[0].query_start = 1
        result.hsps[0].query_end = 99
        result.hsps[0].hit_start = 322
        result.hsps[0].hit_end = 428
        result.hsps[0].aln=[Namespace(), Namespace()]
        result.hsps[0].aln[0].seq = "RAVDELIRYLTVPYGPTPRIAKQDVTVGDQVIKAGESVICSLPAANRDPALVPDADRLDVTR--------DPVPHVAFGHGIHHCLGAALARLELRTVFTALWRRF"
        result.hsps[0].aln[1].seq = "avikEtLRlhpvvplllpRevtkdvvirgylipkGtevivnlyalhrdpevfpnPeeFdpeRFldekgsrksfaflPFGaGpRnCiGerlArmelklflatlLqnF"
        
        self.result = result
        
        self.options.QualifierTags = QualifierTags
        
        self.seq_record = seqio.read(utils.get_full_path(__file__, 'Y16952.3.final.gbk'))
        
        self.assertEqual(361, len(self.seq_record.features))
        
        xmltext ="""<?xml version="1.0" encoding="UTF-8"?>

<resource xmlns:xsd="http://www.w3.org/2001/XMLSchema"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
    <analysis name='ASP_P450Oxy' type="active_site"> 
        <Prerequisite>
            <primary_tag_type>PFAM_domain</primary_tag_type>
            <tag>domain</tag>
            <tag_value>p450</tag_value>
        </Prerequisite>
        <Execute program="hmmscan" CaptureConsole="TRUE">
            <!-- Currently, the location of the hmmpfam2 binary and the database location 
                is inferred from the antismash configuration file! -->

            <parameters>
                <!-- If prefixes for parameters are required they can be added as attribute 
                    prefix -->
                <parameter name="evalue" prefix="--domE">0.1</parameter>
                <parameter name="cpus" prefix="--cpu">1</parameter>
            </parameters>
            <database>p450.hmm3</database>
            <db_source>PFAM24</db_source>
            <BioPythonParser>hmmer3-text</BioPythonParser>
        </Execute>
        <Alignment>
            <scaffold>
                <scaffoldOffset>327,330,400,403,409</scaffoldOffset>
                <scaffoldValue>E,R,F,G,G</scaffoldValue>
            </scaffold>
            <choice result="active site cystein present">
                <offset>407</offset>
                <value>C</value>
                <comment>Cytochrome P450 oxygenase active site cystein; coordinates heme Fe ligand</comment>
            </choice>
        </Alignment>
        <description>Pediction of cytochrome P450 active site cystein</description>
        <referenceList>
            <reference>Del Vecchio, F., H. Petkovic, S. G. Kendrew, L. Low, B.
                Wilkinson, R. Lill, J. Cortes, B. A. Rudd, J. Staunton, and P. F.
                Leadlay. 2003. Active-site residue, domain and module swaps in
                modular polyketide synthases. J Ind. Microbiol Biotechnol
                30:489-494.</reference>
        </referenceList>
    </analysis>
</resource>
"""
        
        ETObj = ET.fromstring(xmltext)
        test = ETObj.find('./analysis')
        self.ETObj = ETObj
        
        # now acutally generate ASF object
        myASF = antismash.generic_modules.active_site_finder.active_site_finder(self.seq_record, self.options)
        
        assert_is_instance(myASF, antismash.generic_modules.active_site_finder.active_site_finder)
        
        self.my_ASF = myASF
        
        
    def test__init(self):
        "Test active_site_finder.__init__ method"
        
        self.tt = TraceTracker()
        
        class dummyET:
            def __init__(self):
                pass
            def getroot(self):
                return dummyET()
            def findall (self, a):
                return "ET-Subtree"
        
        
        mock('ET.parse', returns=dummyET(), tracker=self.tt)
        #mock('ET.getroot', returns="XMLTree", tracker=self.tt)
        #mock('ET.findall', returns="XMLSubtree", tracker=self.tt)
        
        mock('antismash.generic_modules.active_site_finder.active_site_finder.check_prereqs', returns=[], tracker=self.tt)
        
        expected = """    Called ET.parse(
        '.../antismash/generic_modules/active_site_finder/config/SignatureResources.xml')
    Called antismash.generic_modules.active_site_finder.active_site_finder.check_prereqs(
        )"""
        
        ASF = antismash.generic_modules.active_site_finder.active_site_finder(self.seq_record, self.options)
        
        assert_same_trace(self.tt, expected)

        # restore mocked objects that other tests can run with the original code
        restore()
        self.tt = TraceTracker()
        
        # now acutally generate ASF object
        myASF = antismash.generic_modules.active_site_finder.active_site_finder(self.seq_record, self.options)
        
        assert_is_instance(myASF, antismash.generic_modules.active_site_finder.active_site_finder)
        
        self.my_ASF = myASF


    def test__get_scaffold_annotation(self):
        "Test active_site_finder._get_scaffold_annotation method"
        
        my_ScaffXML = self.ETObj.find('./analysis/Alignment/scaffold')
        
        resultline = self.my_ASF._get_scaffold_annotation(self.result, my_ScaffXML)
        
        expectation = "Scaffold coordinates: (327,330,400,403,409); scaffold residues: (E,R,F,G,G); expected: (E,R,F,G,G); matchArray: (True,True,True,True,True); emission probability array (n.d.,n.d.,n.d.,n.d.,n.d.); overall match: TRUE"
        
        self.assertEqual(resultline, expectation, "scaffold line mismatch")
        
        
    def test__get_prediction_annotation(self):
        "Test active_site_finder._get_prediction_annotation method"
        my_ChoicesXML = self.ETObj.findall('./analysis/Alignment/choice')
        
        (choiceList, predictionList) = self.my_ASF._get_prediction_annotation(self.result, my_ChoicesXML)
        
        expected = ["Description: Cytochrome P450 oxygenase active site cystein; coordinates heme Fe ligand, choice result: active site cystein present, choice coordinates: (407); residues: (C); expected for choice: (C); matchArray: (True); emission probability array (n.d.); overall match: TRUE"]
        self.assertListEqual(choiceList, expected, "prediction choices mismatch")
        
        expected = ["Full match for prediction: active site cystein present"]
        self.assertListEqual(predictionList, expected, "prediction string mismatch")
        
        
    def test__execute_tool(self):
        "Test active_site_finder._execute_tool method"
        self.tt=TraceTracker()
        self.tt2=TraceTracker()
        mock('utils.execute', returns=["shell output", 0, 0], tracker=self.tt)
        
        mock('SearchIO.parse', returns = ["SearchIO object"], tracker=self.tt2)
        
        result = self.my_ASF._execute_tool(self.ETObj.find('./analysis'), fileName="testTempfile")
        
        self.assertListEqual(result, ["SearchIO object"])
        
        expected = """        Called utils.execute(
            ['hmmscan', '--domE', '0.1', '--cpu', '1', '.../antismash/generic_modules/active_site_finder/hmm/p450.hmm3', 'testTempfile'])
      """
        assert_same_trace(self.tt, expected)
        
        expected = "   Called SearchIO.parse(<cStringIO.StringI object at ...>, 'hmmer3-text')"
        assert_same_trace(self.tt2, expected)
        
        restore()
        
        self.tt=TraceTracker()
        self.tt2=TraceTracker()
        mock('utils.execute', returns=["shell output", 0, 0], tracker=self.tt)
        
        mock('SearchIO.parse', returns = ["SearchIO object"], tracker=self.tt2)
        
        result = self.my_ASF._execute_tool(self.ETObj.find('./analysis'), stdin_data="fasta sequence from stdin")
        
        self.assertListEqual(result, ["SearchIO object"])
        
        expected = """        Called utils.execute(
            ['hmmscan', '--domE', '0.1', '--cpu', '1', '.../antismash/generic_modules/active_site_finder/hmm/p450.hmm3'],
            input='fasta sequence from stdin')
      """
        assert_same_trace(self.tt, expected)
        
        expected = "   Called SearchIO.parse(<cStringIO.StringI object at ...>, 'hmmer3-text')"
        assert_same_trace(self.tt2, expected)
        
        restore()
        self.tt=TraceTracker()
        del self.tt2
        
        
    def test__run_external_tool(self):
        "Test active_site_finder._run_external_tool method"

        self.tt = TraceTracker()
        
        SeqFeatureList = utils.get_all_features_of_type_with_query(self.seq_record, "PFAM_domain", "domain", "p450")
        self.assertEqual(6, len(SeqFeatureList))
       
        mock('antismash.generic_modules.active_site_finder.active_site_finder._execute_tool', returns=["external program call successful"], tracker=self.tt)
        
        result = self.my_ASF._run_external_tool(self.ETObj.find('./analysis'), SeqFeatureList)
        expected = ["external program call successful"]
        self.assertListEqual(result, expected)
        
        
    def test__fix_coordinates(self):
        "Test active_site_finder._fix_coordinates method"
        
        seq = "lsgpeavkevlikkgeefs..grgdeallatsrkafkgkgvlfangekwkklRrfltptltsf.klsleelveeeaedlv"
        
        coord = self.my_ASF._fix_coordinates(10, seq)
        self.assertEqual(coord, 10, "assignment without gap")
        
        coord = self.my_ASF._fix_coordinates(25, seq)
        self.assertEqual(coord, 27, "assignment with 2 gaps (at pos 20,21 of original string)")
        
        coord = self.my_ASF._fix_coordinates(20, seq)
        self.assertEqual(coord, 22, "assignment with 2 gaps (at end position 20 and 21)")
        
        coord = self.my_ASF._fix_coordinates(76, seq)
        self.assertEqual(coord, 79, "assignment with 3 gaps (at end position 20,21 and 64)")
                
        with self.assertRaises(ValueError):
            coord = self.my_ASF._fix_coordinates(77, seq)
        
    def test_check_prereqs(self):
        "Test active_site_finder.check_prereqs method"
        # THIS TEST HAS TO BE UPDATED WHEN NEW PROFILES ARE ADDED TO THE MODULE!
        
        self.tt = TraceTracker()
        
        mock('utils.locate_executable', returns="/my/path/to/executable", tracker=self.tt)
        mock('utils.locate_file', returns="/my/path/to/file", tracker=self.tt)
        
        result = self.my_ASF.check_prereqs()
        
        self.assertListEqual(result, [], "return empty list if executables/files are found")
        
        expected = """    Called utils.locate_executable('blastp')
    Called utils.locate_executable('hmmpfam2')
    Called utils.locate_executable('hmmscan')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/PKSI-KR.hmm2')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/PKSI-KS_N.hmm2')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/PKSI-KS_C.hmm2')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/PKSI-AT.hmm2')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/PKSI-ACP.hmm2')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/PKSI-DH.hmm2')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/PKSI-KR.hmm2')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/Thioesterase.hmm2')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/PKSI-ER.hmm2')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/PKSI-AT.hmm2')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/aa-activating.aroundLys.hmm2')
    Called utils.locate_file(
        '.../antismash/generic_modules/active_site_finder/hmm/p450.hmm3')"""
        assert_same_trace(self.tt, expected)
        
        restore()
        self.tt = TraceTracker()
        
        mock('utils.locate_executable', returns=None, tracker=self.tt)
        mock('utils.locate_file', returns="/my/path/to/file", tracker=self.tt)
        
        result = self.my_ASF.check_prereqs()
        expected = ["Failed to locate file: 'blastp'", "Failed to locate file: 'hmmpfam2'", "Failed to locate file: 'hmmscan'"]
        self.assertListEqual(result, expected, "test result if file not found")
        restore()