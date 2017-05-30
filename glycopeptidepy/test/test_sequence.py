import unittest
import pickle

from glycopeptidepy.structure import sequence, modification, residue, composition
from glypy import GlycanComposition, Glycan, MonosaccharideResidue


HYDROGEN = composition.Composition("H").mass
R = residue.Residue


p1 = "PEPTIDE"
p2 = "YPVLN(N-Glycosylation)VTMPN(Deamidation)NGKFDK{Hex:9; HexNAc:2}"
p3 = "NEEYN(N-Glycosylation)K{Hex:5; HexNAc:4; NeuAc:2}"
p4 = "TVDGT(O-Glycosylation)AR{Fuc:1; Hex:1; HexNAc:1; Neu5Ac:1}"
p5 = "(Carbamidomethyl)-FFYFTPNK"
p6 = "NEEYN(N-Glycosylation)K{Fuc:1; Hex:5; HexNAc:4; NeuAc:2}"
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

    def test_n_glycan_sequons(self):
        seq = self.parse_sequence(p2).deglycosylate()
        self.assertTrue(4 in seq.n_glycan_sequon_sites)

    def test_o_glycopeptide(self):
        masses = [
            718.36096708622,
            921.44033960573,
            1067.49824840467,
            1083.49316302423,
            1229.5510718231699
        ]
        for i, stub in enumerate(self.parse_sequence(p4).stub_fragments()):
            self.assertAlmostEqual(masses[i], stub.mass, 3)

    def test_composition(self):
        for p in [p1, p2, p3, p4, p5]:
            seq = self.parse_sequence(p)
            self.assertAlmostEqual(seq.total_mass, seq.total_composition().mass, 4)

    def test_repr_(self):
        seq = self.parse_sequence(p5)
        self.assertTrue("(Carbamidomethyl)" in (repr(seq).split('-')[0]))

    def test_glycan_composition(self):
        seq = self.parse_sequence(p3)
        self.assertEqual(seq.glycan_composition, GlycanComposition(Hex=5, HexNAc=4, NeuAc=2))

    def test_delete(self):
        seq = self.parse_sequence(p1)
        ref = self.parse_sequence(p1)
        seq.delete(0)
        proline = R("P")
        self.assertAlmostEqual(
            seq.total_mass, ref.total_mass - proline.mass, 4)

    def test_fucosylated_stubs(self):
        seq = self.parse_sequence(p6)

        had_fucose = False
        for stub in seq.stub_fragments(extended=1):
            if "Fuc" in stub.name:
                had_fucose = True
        self.assertTrue(had_fucose)


class TestPeptideSequence(PeptideSequenceSuiteBase, unittest.TestCase):
    def parse_sequence(self, seqstr):
        return sequence.PeptideSequence(seqstr)


class TestNamedSequence(PeptideSequenceSuiteBase, unittest.TestCase):
    def parse_sequence(self, seqstr):
        return sequence.NamedSequence("spam", seqstr)


if __name__ == '__main__':
    unittest.main()
