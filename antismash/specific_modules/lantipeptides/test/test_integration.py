try: import unittest2 as unittest
except ImportError: import unittest
from helperlibs.bio import seqio
from antismash import utils
from antismash.config import set_config
from antismash.specific_modules.lantipeptides import specific_analysis
from antismash.specific_modules.lantipeptides import html_output as h

class TestIntegration(unittest.TestCase):
    def setUp(self):
        from argparse import Namespace
        conf = Namespace()
        conf.cpus = 1
        set_config(conf)

    def tearDown(self):
        set_config(None)

    def test_nisin(self):
        "Test lantipeptide prediction for nisin A"
        rec = seqio.read(utils.get_full_path(__file__, 'nisin.gbk'))
        self.assertEqual(38, len(rec.features))

        specific_analysis(rec, None)
        self.assertEqual(40, len(rec.features))
        prepeptides = h._find_core_peptides(utils.get_cluster_by_nr(rec, 1), rec)
        self.assertEqual(1, len(prepeptides))
        prepeptide = prepeptides[0]
        leaders = h._find_leader_peptides(utils.get_cluster_by_nr(rec, 1), rec)
        self.assertEqual(1, len(leaders))
        leader = leaders[0]
        # real monoisotopic mass is 3351.51, but we overpredict a Dha
        self.assertAlmostEqual(3333.6, h._get_monoisotopic_mass(prepeptide))
        # real mw is 3354.5, see above
        self.assertAlmostEqual(3336.0, h._get_molecular_weight(prepeptide))
        self.assertEqual([3354.0, 3372.1, 3390.1, 3408.1], h._get_alternative_weights(prepeptide))
        self.assertEqual(5, h._get_number_bridges(prepeptide))
        self.assertEqual("MSTKDFNLDLVSVSKKDSGASPR", h._get_leader_peptide_sequence(leader))
        self.assertEqual("ITSISLCTPGCKTGALMGCNMKTATCHCSIHVSK", h._get_core_peptide_sequence(prepeptide))
        self.assertEqual('Class I', h._get_core_peptide_class(prepeptide))


    def test_epidermin(self):
        "Test lantipeptide prediction for epidermin"
        rec = seqio.read(utils.get_full_path(__file__, 'epidermin.gbk'))
        self.assertEqual(18, len(rec.features))

        specific_analysis(rec, None)
        self.assertEqual(20, len(rec.features))
        prepeptides = h._find_core_peptides(utils.get_cluster_by_nr(rec, 1), rec)
        self.assertEqual(1, len(prepeptides))
        prepeptide = prepeptides[0]
        leaders = h._find_leader_peptides(utils.get_cluster_by_nr(rec, 1), rec)
        self.assertEqual(1, len(leaders))
        leader = leaders[0]
        self.assertAlmostEqual(2164, h._get_monoisotopic_mass(prepeptide))
        self.assertAlmostEqual(2165.6, h._get_molecular_weight(prepeptide))
        self.assertEqual(3, h._get_number_bridges(prepeptide))
        self.assertEqual("MEAVKEKNDLFNLDVKVNAKESNDSGAEPR", h._get_leader_peptide_sequence(leader))
        self.assertEqual("IASKFICTPGCAKTGSFNSYCC", h._get_core_peptide_sequence(prepeptide))
        self.assertEqual('Class I', h._get_core_peptide_class(prepeptide))
        self.assertEqual(['AviCys'], h._get_core_peptide_extra_modifications(prepeptide))


    def test_microbisporicin(self):
        "Test lantipeptide prediction for microbisporicin"
        rec = seqio.read(utils.get_full_path(__file__, 'microbisporicin.gbk'))
        self.assertEqual(56, len(rec.features))

        specific_analysis(rec, None)
        self.assertEqual(58, len(rec.features))
        prepeptides = h._find_core_peptides(utils.get_cluster_by_nr(rec, 1), rec)
        self.assertEqual(1, len(prepeptides))
        prepeptide = prepeptides[0]
        leaders = h._find_leader_peptides(utils.get_cluster_by_nr(rec, 1), rec)
        self.assertEqual(1, len(leaders))
        leader = leaders[0]
        # NOTE: this is not the correct weight for microbisporicin
        # there are some additional modifications we do not predict yet
        self.assertAlmostEqual(2212.9, h._get_monoisotopic_mass(prepeptide))
        self.assertAlmostEqual(2214.5, h._get_molecular_weight(prepeptide))
        self.assertEqual(4, h._get_number_bridges(prepeptide))
        self.assertEqual("MPADILETRTSETEDLLDLDLSIGVEEITAGPA", h._get_leader_peptide_sequence(leader))
        self.assertEqual("VTSWSLCTPGCTSPGGGSNCSFCC", h._get_core_peptide_sequence(prepeptide))
        self.assertEqual('Class I', h._get_core_peptide_class(prepeptide))
        self.assertEqual(['AviCys', 'Cl', 'OH'], h._get_core_peptide_extra_modifications(prepeptide))


    def test_epicidin(self):
        "Test lantipeptide prediction for epicidin 280"
        rec = seqio.read(utils.get_full_path(__file__, 'epicidin_280.gbk'))
        self.assertEqual(21, len(rec.features))

        specific_analysis(rec, None)
        self.assertEqual(23, len(rec.features))
        prepeptides = h._find_core_peptides(utils.get_cluster_by_nr(rec, 1), rec)
        self.assertEqual(1, len(prepeptides))
        prepeptide = prepeptides[0]
        leaders = h._find_leader_peptides(utils.get_cluster_by_nr(rec, 1), rec)
        self.assertEqual(1, len(leaders))
        leader = leaders[0]
        self.assertAlmostEqual(3115.7, h._get_monoisotopic_mass(prepeptide))
        self.assertAlmostEqual(3117.7, h._get_molecular_weight(prepeptide))
        self.assertEqual([3135.7, 3153.7, 3171.7], h._get_alternative_weights(prepeptide))
        self.assertEqual(3, h._get_number_bridges(prepeptide))
        self.assertEqual("MENKKDLFDLEIKKDNMENNNELEAQ", h._get_leader_peptide_sequence(leader))
        self.assertEqual("SLGPAIKATRQVCPKATRFVTVSCKKSDCQ", h._get_core_peptide_sequence(prepeptide))
        self.assertEqual('Class I', h._get_core_peptide_class(prepeptide))
        self.assertEqual(['Lac'], h._get_core_peptide_extra_modifications(prepeptide))


    def test_labyrinthopeptin(self):
        "Test lantipeptide prediction for labyrinthopeptin"
        rec = seqio.read(utils.get_full_path(__file__, 'labyrinthopeptin.gbk'))
        self.assertEqual(7, len(rec.features))

        specific_analysis(rec, None)
        self.assertEqual(11, len(rec.features))


    def test_sco_cluster3(self):
        "Test lantipeptide prediction for SCO cluster #3"
        rec = seqio.read(utils.get_full_path(__file__, 'sco_cluster3.gbk'))
        self.assertEqual(69, len(rec.features))

        specific_analysis(rec, None)
        self.assertEqual(71, len(rec.features))
        prepeptides = h._find_core_peptides(utils.get_cluster_by_nr(rec, 1), rec)
        self.assertEqual(1, len(prepeptides))
        prepeptide = prepeptides[0]
        self.assertEqual('Class I', h._get_core_peptide_class(prepeptide))
