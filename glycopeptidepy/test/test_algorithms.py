import unittest

from glycopeptidepy.structure import sequence
from glycopeptidepy.algorithm import reverse_preserve_sequon, edit_distance

p2 = "YPVLN(N-Glycosylation)VTMPN(Deamidation)NGKFDK{Hex:9; HexNAc:2}"


heparanase = ''.join("""MLLRSKPALPPPLMLLLLGPLGPLSPGALPRPAQAQDVVDLDFFT
QEPLHLVSPSFLSVTIDANLATDPRFLILLGSPKLRTLARGLSPA
YLRFGGTKTDFLIFDPKKESTFEERSYWQSQVNQDICKYGSIPPD
VEEKLRLEWPYQEQLLLREHYQKKFKNSTYSRSSVDVLYTFANCS
GLDLIFGLNALLRTADLQWNSSNAQLLLDYCSSKGYNISWELGNE
PNSFLKKADIFINGSQLGEDFIQLHKLLRKSTFKNAKLYGPDVGQ
PRRKTAKMLKSFLKAGGEVIDSVTWHHYYLNGRTATKEDFLNPDV
LDIFISSVQKVFQVVESTRPGKKVWLGETSSAYGGGAPLLSDTFA
AGFMWLDKLGLSARMGIEVVMRQVFFGAGNYHLVDENFDPLPDYW
LSLLFKKLVGTKVLMASVQGSKRRKLRVYLHCTNTDNPRYKEGDL
TLYAINLHNVTKYLRLPYPFSNKQVDKYLLRPLGPHGLLSKSVQL
NGLTLKMVDDQTLPPLMEKPLRPGSSLGLPAFSYSFFVIRNAKVA
ACI""".splitlines())


class TestReversePreserveSequonTarget(unittest.TestCase):
    def test_with_glycan(self):
        ref = sequence.PeptideSequence(p2)
        mass = ref.total_mass
        composition = ref.total_composition()

        perm = reverse_preserve_sequon(ref)
        self.assertAlmostEqual(mass, perm.total_mass, 4)
        self.assertAlmostEqual(ref.total_mass, perm.total_mass, 4)
        self.assertEqual(composition, perm.total_composition())
        self.assertEqual(ref.total_composition(), perm.total_composition())

    def test_detect_similarity(self):
        fwd = sequence.PeptideSequence("PEPEPR")
        rev = reverse_preserve_sequon(fwd)
        assert edit_distance(fwd, rev) > 0


class TestEditDistance(unittest.TestCase):
    def test_distance(self):
        assert edit_distance("PEPTIDE", "PIPTIDE") == 1
        assert edit_distance(sequence.PeptideSequence("PEPTIDE"),
                             sequence.PeptideSequence("PIPTIDE")) == 1

    def test_reverse_distance(self):
        fwd = sequence.PeptideSequence(heparanase)
        rev = reverse_preserve_sequon(fwd)
        assert edit_distance(fwd, rev) > 0


if __name__ == '__main__':
    unittest.main()
