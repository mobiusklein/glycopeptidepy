import unittest

from glycopeptidepy.structure import residue


class AminoAcidResidueTests(unittest.TestCase):
    def test_residue_creation(self):
        leucine = residue.AminoAcidResidue(symbol='L')
        isoleucine = residue.AminoAcidResidue(name='Ile')
        self.assertEqual(leucine.composition, isoleucine.composition)

    def test_degeneracy(self):
        xle = residue.AminoAcidResidue("J")
        leucine = residue.AminoAcidResidue(symbol='L')
        self.assertTrue(xle.is_degenerate)
        self.assertFalse(leucine.is_degenerate)
        self.assertTrue(leucine in residue.degeneracy_index[xle.name])

        self.assertNotIn(leucine, residue.get_all_sequencing_residues())
