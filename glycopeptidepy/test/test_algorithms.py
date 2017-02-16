import unittest

from glycopeptidepy.structure import sequence
from glycopeptidepy.algorithm import reverse_preserve_sequon

p2 = "YPVLN(N-Glycosylation)VTMPN(Deamidation)NGKFDK{Hex:9; HexNAc:2}"


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


if __name__ == '__main__':
    unittest.main()
