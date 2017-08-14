import unittest

from glycopeptidepy.structure import sequence, modification, residue


n_glycan_rule = '#:glycoct:RES 1b:b-dglc-HEX-1:5 2s:n-acetyl 3b:b-dglc-HEX-1:5 \
4s:n-acetyl 5b:b-dman-HEX-1:5 6b:a-dman-HEX-1:5 7b:a-dman-HEX-1:5 LIN \
1:1d(2+1)2n 2:1o(4+1)3d 3:3d(2+1)4n 4:3o(4+1)5d 5:5o(3+1)6d 6:5o(6+1)7d'


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

    def test_parse_target(self):
        spaced = 'Q @ N-term'
        part_spaced = 'Q @N-Term'
        part_spaced2 = 'Q@ n-term'
        no_space = 'Q@N-term'
        spaced_parsed = modification.extract_targets_from_string(spaced)
        self.assertEqual(spaced_parsed, modification.extract_targets_from_string(part_spaced))
        self.assertEqual(spaced_parsed, modification.extract_targets_from_string(part_spaced2))
        self.assertEqual(spaced_parsed, modification.extract_targets_from_string(no_space))


class TestGlycosylationRule(unittest.TestCase):
    def test_parse(self):
        rule = modification.Glycosylation.try_parse(n_glycan_rule)
        base_mass = modification.glycoct.loads("\
            RES 1b:b-dglc-HEX-1:5 2s:n-acetyl \
            3b:b-dglc-HEX-1:5 4s:n-acetyl \
            5b:b-dman-HEX-1:5 6b:a-dman-HEX-1:5 \
            7b:a-dman-HEX-1:5 LIN 1:1d(2+1)2n \
            2:1o(4+1)3d 3:3d(2+1)4n 4:3o(4+1)5d \
            5:5o(3+1)6d 6:5o(6+1)7d").mass()
        dehydrated_mass = base_mass - modification.Composition("H2O").mass
        self.assertAlmostEqual(rule.mass, dehydrated_mass, 3)


if __name__ == '__main__':
    unittest.main()
