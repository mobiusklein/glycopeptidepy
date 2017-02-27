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


class TestPeptideSequence(unittest.TestCase):
    def test_parser(self):
        chunks, mods, glycan, n_term, c_term = sequence.sequence_tokenizer(p1)
        self.assertEqual(len(mods), 0)
        self.assertEqual(len(chunks), len(p1))
        self.assertEqual(glycan, "")

        chunks, mods, glycan, n_term, c_term = sequence.sequence_tokenizer(p2)
        self.assertEqual(GlycanComposition.parse("{Hex:9; HexNAc:2}"), glycan)
        self.assertEqual(len(mods), 2)
        self.assertEqual(len(chunks), 16)

    def test_fragmentation(self):
        seq = sequence.parse(p1)
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
        peptide = sequence.parse(p3)
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
        peptide = sequence.parse(p3)
        stubs = {f.name: f.mass for f in peptide.glycan_fragments(all_series=True, allow_ambiguous=True)}
        self.assertAlmostEqual(stubs["peptide+{Hex:3; HexNAc:2}"], 1687.6571, 3)
        self.assertAlmostEqual(stubs["peptide+{Hex:3; HexNAc:3}"],
                               1687.6571 + hexnac_mass, 3)
        self.assertAlmostEqual(stubs["peptide+{Hex:4; HexNAc:3}"],
                               1687.6571 + hexnac_mass + hexose_mass, 3)

    def test_clone(self):
        peptide = sequence.parse(p3)
        self.assertEqual(peptide, peptide.clone())

    def test_add_remove_modification(self):
        p = sequence.parse("ENGTISR")
        ref_mass = p.mass
        p.add_modification(1, "Deamidated")
        self.assertAlmostEqual(sequence.parse(str(p)).mass, p.mass, 5)

        p.drop_modification(1, "Deamidated")

        self.assertAlmostEqual(p.mass, ref_mass, 5)

    def test_picklability(self):
        original = sequence.parse(p3)
        duplicate = pickle.loads(pickle.dumps(original))
        self.assertEqual(
            sequence._total_composition(original),
            sequence._total_composition(duplicate))


if __name__ == '__main__':
    unittest.main()
