import unittest

from glycopeptidepy.structure import sequence, modification, residue
from glycopeptidepy.structure import Composition

from glypy.io import glycoct

n_glycan_rule = '#:glycoct:RES 1b:b-dglc-HEX-1:5 2s:n-acetyl 3b:b-dglc-HEX-1:5 \
4s:n-acetyl 5b:b-dman-HEX-1:5 6b:a-dman-HEX-1:5 7b:a-dman-HEX-1:5 LIN \
1:1d(2+1)2n 2:1o(4+1)3d 3:3d(2+1)4n 4:3o(4+1)5d 5:5o(3+1)6d 6:5o(6+1)7d'


class TestModifcationTarget(unittest.TestCase):
    def test_n_term_special_flag(self):
        mt = modification.ModificationTable()
        rule = mt["Acetylation"]
        peptide = sequence.parse("QVPQLQSQ")
        sites = rule.find_valid_sites(peptide)
        self.assertTrue(sites[0] is modification.SequenceLocation.n_term)

        rule = mt['Gln->pyro-Glu']
        sites = rule.find_valid_sites(peptide)
        self.assertEqual(sites[0], 0)

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

    def test_serialize(self):
        spaced = 'Q @ N-term'
        spaced_parsed = modification.extract_targets_from_string(spaced)
        assert spaced == spaced_parsed.serialize()

    def test__len__(self):
        spaced = 'Q @ N-term'
        spaced_parsed = modification.extract_targets_from_string(spaced)
        assert len(spaced_parsed) == 1


class TestGlycosylationRule(unittest.TestCase):
    def test_parse(self):
        rule = modification.Glycosylation.try_parse(n_glycan_rule)
        base_mass = glycoct.loads("\
            RES 1b:b-dglc-HEX-1:5 2s:n-acetyl \
            3b:b-dglc-HEX-1:5 4s:n-acetyl \
            5b:b-dman-HEX-1:5 6b:a-dman-HEX-1:5 \
            7b:a-dman-HEX-1:5 LIN 1:1d(2+1)2n \
            2:1o(4+1)3d 3:3d(2+1)4n 4:3o(4+1)5d \
            5:5o(3+1)6d 6:5o(6+1)7d").mass()
        dehydrated_mass = base_mass - Composition("H2O").mass
        self.assertAlmostEqual(rule.mass, dehydrated_mass, 3)


class TestModificationTable(unittest.TestCase):
    def test_restricted_table(self):
        table = modification.RestrictedModificationTable(
            constant_modifications=['Carbamidomethyl (C)'],
            variable_modifications=['Deamidated (N)'])

        total = len(table.rules())
        given = len(table.other_modifications)
        assert total - given == 2

    def test_by_unimod(self):
        mod = modification.Modification("UNIMOD:4")
        assert mod.name == "Carbamidomethyl"


class TestModificationRule(unittest.TestCase):
    def test_targets(self):
        mod = modification.Modification("Deamidation")
        assert len(mod.rule.targets) == 5
        assert len(mod.rule.n_term_targets) == 1
        assert len(mod.rule.c_term_targets) == 0

    def test_subtraction(self):
        mod = modification.Modification("Deamidation")
        assert len((mod.rule - mod.rule).targets) == 0

    def test_serialize(self):
        mod = modification.Modification("Deamidation")
        assert mod.rule.serialize() == mod.rule.name


class TestAminoAcidSubstitution(unittest.TestCase):
    def test_parse(self):
        rule_name = "@Asn->Pro"
        rule = modification.AminoAcidSubstitution.try_parse(rule_name)
        assert rule is not None
        assert rule.serialize() == rule_name


class AnonymousModificationRule(unittest.TestCase):
    def test_parse(self):
        rule_name = "@Name-1.0"
        rule = modification.AnonymousModificationRule.try_parse(rule_name)
        assert rule.name == "Name"
        assert rule.mass == 1.0
        assert rule.serialize() == rule_name
        assert not rule.is_standard
        assert rule == rule.clone()


if __name__ == '__main__':
    unittest.main()
