import unittest

from glycopeptidepy.structure import sequence, modification, residue


class TestModifcationTarget(unittest.TestCase):
    def test_n_term_special_flag(self):
        mt = modification.ModificationTable()
        rule = mt["Pyro-glu from Q"]
        peptide = sequence.parse("QVPQLQSQ")
        sites = rule.find_valid_sites(peptide)
        self.assertTrue(sites[0] is modification.SequenceLocation.n_term)
        deamidated = mt["Deamidated"]
        sites = deamidated.find_valid_sites(peptide)
        self.assertEqual(sites, [0, 3, 5, 7])


if __name__ == '__main__':
    unittest.main()
