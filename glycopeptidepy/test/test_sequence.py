import unittest
import pickle

from glycopeptidepy.structure import sequence, modification, residue, composition
from glypy import GlycanComposition, Glycan, MonosaccharideResidue


HYDROGEN = composition.Composition("H").mass
R = residue.Residue


p1 = "PEPTIDE"
p2 = "YPVLN(N-Glycosylation)VTMPN(Deamidation)NGKFDK{Hex:9; HexNAc:2}"
p3 = "NEEYN(N-Glycosylation)K{Hex:5; HexNAc:4; NeuAc:2}"
hexnac_mass = MonosaccharideResidue.from_iupac_lite("HexNAc").mass()
hexose_mass = MonosaccharideResidue.from_iupac_lite("Hex").mass()


class TestSequenceParser(unittest.TestCase):
    def test_parser(self):
        chunks, mods, glycan, n_term, c_term = sequence.sequence_tokenizer(p1)
        self.assertEqual(len(mods), 0)
        self.assertEqual(len(chunks), len(p1))
        self.assertEqual(glycan, "")

        chunks, mods, glycan, n_term, c_term = sequence.sequence_tokenizer(p2)
        self.assertEqual(GlycanComposition.parse("{Hex:9; HexNAc:2}"), glycan)
        self.assertEqual(len(mods), 2)
        self.assertEqual(len(chunks), 16)


class PeptideSequenceSuiteBase(object):
    def parse_sequence(self, seqstr):
        raise NotImplementedError()

    def test_fragmentation(self):
        seq = self.parse_sequence(p1)
        mapping = {
            "b3": 324.1554 - HYDROGEN,
            "b6": 653.3141 - HYDROGEN,
            "y2": 263.0874 - HYDROGEN,
            "y5": 574.2719 - HYDROGEN,
            "a4": 397.2082 - HYDROGEN,
            "a5": 510.2922 - HYDROGEN,
            "z3": 360.1527 - HYDROGEN * 2,
        }
        for fragment, mass in mapping.items():
            fmass = seq.fragment(fragment).mass
            self.assertAlmostEqual(mass, fmass, 2, (mass, fmass, fragment))
            self.assertAlmostEqual(mass, seq.fragment(fragment).total_composition().mass, 2, (mass, fmass, fragment))

    def test_stub_ions(self):
        peptide = self.parse_sequence(p3)
        composition = peptide.total_composition()
        stubs = sorted({f.mass for f in peptide.stub_fragments()})
        self.assertAlmostEqual(stubs[0], 795.3399, 3)
        self.assertAlmostEqual(stubs[1], 998.4193, 3)
        self.assertAlmostEqual(stubs[2], 1201.4986, 3)
        self.assertAlmostEqual(stubs[3], 1363.5515, 3)
        self.assertAlmostEqual(stubs[4], 1525.6043, 3)
        self.assertAlmostEqual(stubs[5], 1687.6571, 3)
        self.assertEqual(composition, peptide.total_composition())

    def test_glycan_fragments_stubs(self):
        peptide = self.parse_sequence(p3)
        stubs = {f.name: f.mass for f in peptide.glycan_fragments(all_series=True, allow_ambiguous=True)}
        self.assertAlmostEqual(stubs["peptide+{Hex:3; HexNAc:2}"], 1687.6571, 3)
        self.assertAlmostEqual(stubs["peptide+{Hex:3; HexNAc:3}"],
                               1687.6571 + hexnac_mass, 3)
        self.assertAlmostEqual(stubs["peptide+{Hex:4; HexNAc:3}"],
                               1687.6571 + hexnac_mass + hexose_mass, 3)

    def test_clone(self):
        peptide = self.parse_sequence(p3)
        self.assertEqual(peptide, peptide.clone())
        self.assertEqual(peptide.total_mass, peptide.clone().total_mass)
        self.assertEqual(peptide.total_composition(), peptide.clone().total_composition())
        self.assertTrue(peptide.full_structure_equality(peptide.clone()))

    def test_add_remove_modification(self):
        p = self.parse_sequence("ENGTISR")
        ref_mass = p.mass
        p.add_modification(1, "Deamidated")
        self.assertAlmostEqual(self.parse_sequence(str(p)).mass, p.mass, 5)

        p.drop_modification(1, "Deamidated")

        self.assertAlmostEqual(p.mass, ref_mass, 5)

    def test_picklability(self):
        original = self.parse_sequence(p3)
        duplicate = pickle.loads(pickle.dumps(original))
        self.assertEqual(
            sequence._total_composition(original),
            sequence._total_composition(duplicate))
        self.assertEqual(original, duplicate)
        self.assertEqual(original.total_mass, duplicate.total_mass)
        self.assertEqual(original.total_composition(), duplicate.total_composition())
        self.assertTrue(original.full_structure_equality(duplicate))

    def test_total_mass(self):
        seq = self.parse_sequence(p3)
        self.assertAlmostEqual(seq.total_mass, 3000.1123374719496, 3)


class TestPeptideSequence(PeptideSequenceSuiteBase, unittest.TestCase):
    def parse_sequence(self, seqstr):
        return sequence.PeptideSequence(seqstr)


class TestNamedSequence(PeptideSequenceSuiteBase, unittest.TestCase):
    def parse_sequence(self, seqstr):
        return sequence.NamedSequence("spam", seqstr)


if __name__ == '__main__':
    unittest.main()
