import unittest
import pickle

from glycopeptidepy.structure import sequence, modification, residue, composition, parser
from glypy import GlycanComposition, MonosaccharideResidue


HYDROGEN = composition.Composition("H").mass
R = residue.Residue
Modification = modification.Modification

p1 = "PEPTIDE"
p2 = "YPVLN(N-Glycosylation)VTMPN(Deamidation)NGKFDK{Hex:9; HexNAc:2}"
p3 = "NEEYN(N-Glycosylation)K{Hex:5; HexNAc:4; NeuAc:2}"
p4 = "TVDGT(O-Glycosylation)AR{Fuc:1; Hex:1; HexNAc:1; Neu5Ac:1}"
p5 = "(Carbamidomethyl)-FFYFTPNK"
p6 = "NEEYN(N-Glycosylation)K{Fuc:1; Hex:5; HexNAc:4; NeuAc:2}"
p7 = 'ISASGVEDIS(GAG-Linker)R{Xyl:1; a,enHex:1; aHex:1; Hex:1; HexS:1; HexNAc(S):1}'
p8 = ('QQQHLFGSN(#:iupac,glycosylation_type=n_linked:b-D-Galp-(1-4)-b-D-Glcp2NAc-(1-2)'
      '-a-D-Manp-(1-6)-[b-D-Galp-(1-4)-b-D-Glcp2NAc-(1-4)-a-D-Manp-(1-3)]b-D-Manp-(1-4)'
      '-b-D-Glcp2NAc-(1-4)-?-D-Glcp2NAc)VTDC(Carbamidomethyl)SGNFCLFR')
p9 = "N(N-Glycosylation)ITEIVYLTN(N-Glycosylation)TTIEK{Hex:10; HexNAc:8}"
p10 = ('N(#:iupac,glycosylation_type=n_linked:b-D-Galp-(1-4)-b-D-Glcp2NAc-(1-2)-a-D-Manp'
       '-(1-6)-[b-D-Galp-(1-4)-b-D-Glcp2NAc-(1-4)-a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp2NAc'
       '-(1-4)-?-D-Glcp2NAc)ITEIVYLTN(#:iupac,glycosylation_type=n_linked:b-D-Galp-(1-4)'
       '-b-D-Glcp2NAc-(1-2)-a-D-Manp-(1-6)-[b-D-Galp-(1-4)-b-D-Glcp2NAc-(1-4)-a-D-Manp-(1-3)]'
       'b-D-Manp-(1-4)-b-D-Glcp2NAc-(1-4)-?-D-Glcp2NAc)TTIEK')
p11 = ("N(#:glycosylation_type=n_linked:{Hex:5; HexNAc:4})ITEI"
       "VYLTN(#:glycosylation_type=n_linked:{Hex:5; HexNAc:4})TTIEK")
hexnac_mass = MonosaccharideResidue.from_iupac_lite("HexNAc").mass()
hexose_mass = MonosaccharideResidue.from_iupac_lite("Hex").mass()


def make_fragment(self, key):
    for group in self.get_fragments(key[0], compute_compositions=True):
        for frag in group:
            if frag.name == key:
                return frag


class TestSequenceParser(unittest.TestCase):
    def test_parser(self):
        chunks, mods, glycan, n_term, c_term = parser.sequence_tokenizer(p1)
        # self.assertEqual(len(mods), 0)
        self.assertEqual(len(chunks), len(p1))
        self.assertEqual(glycan, None)

        chunks, mods, glycan, n_term, c_term = parser.sequence_tokenizer(p2)
        self.assertEqual(GlycanComposition.parse("{Hex:9; HexNAc:2}"), glycan)
        # self.assertEqual(len(mods), 2)
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
            "z3": 360.1527 - HYDROGEN,
        }

        for fragment, mass in mapping.items():
            fmass = make_fragment(seq, fragment).mass
            self.assertAlmostEqual(mass, fmass, 2, (mass, fmass, fragment))
            self.assertAlmostEqual(
                mass, make_fragment(seq, fragment).total_composition().mass, 2, (mass, fmass, fragment))

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

    def test_stub_ions_extended(self):
        peptide = self.parse_sequence(p3)
        stubs = sorted({f for f in peptide.stub_fragments(extended=True)
                       if f.glycosylation}, key=lambda x: x.glycosylation.mass())
        for stub in stubs:
            if sum(stub.glycosylation.values()) > 5:
                assert stub.is_extended
            else:
                assert not stub.is_extended

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
        peptide = self.parse_sequence(p8)
        self.assertEqual(peptide, peptide.clone())
        self.assertEqual(peptide.total_mass, peptide.clone().total_mass)
        self.assertEqual(peptide.total_composition(), peptide.clone().total_composition())
        self.assertTrue(peptide.full_structure_equality(peptide.clone()))

    def test_add_remove_modification(self):
        p = self.parse_sequence("ENGTISR")
        ref_mass = p.mass
        assert not p.has_modification(1, "Deamidated")
        p.add_modification(1, "Deamidated")
        assert p.has_modification(1, "Deamidated")
        self.assertAlmostEqual(self.parse_sequence(str(p)).mass, p.mass, 5)

        p.drop_modification(1, "Deamidated")
        assert not p.has_modification(1, "Deamidated")
        self.assertAlmostEqual(p.mass, ref_mass, 5)

    def test_picklability(self):
        for s in [p3, p8]:
            original = self.parse_sequence(s)
            duplicate = pickle.loads(pickle.dumps(original))
            self.assertEqual(
                original.total_composition(),
                duplicate.total_composition())
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

    def test_sequon_finding(self):
        seq = self.parse_sequence("PEPTIDES")
        self.assertIn(3, seq.o_glycan_sequon_sites)
        self.assertNotIn(3, seq.glycosaminoglycan_sequon_sites)

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

    def test_multiply_glycosylated_stubs(self):
        seq = self.parse_sequence("N(N-Glycosylation)ITEIVYLTN(N-Glycosylation)TTIEK{Hex:14; HexNAc:4}")
        stubs = list(seq.stub_fragments(extended=True))
        name_set = {f.name: f for f in stubs}
        self.assertTrue("peptide+Hex4HexNAc2" in name_set)
        f = name_set['peptide+Hex4HexNAc2']
        assert f.is_glycosylated
        assert f.glycosylation_size == 6
        self.assertAlmostEqual(f.glycosylation.mass(), 1072.3806, 2)

    def test_gag_linker_oxonium_ions(self):
        seq = self.parse_sequence(p7)
        self.assertAlmostEqual(seq.total_mass, 2285.766, 2)
        ox_map = {f.name: f for f in seq.glycan_fragments()}
        self.assertAlmostEqual(ox_map["HexNAca,en-Hex"].mass, 361.10089, 3)

    def test_gag_linker_stubs(self):
        seq = self.parse_sequence(p7)
        self.assertAlmostEqual(seq.total_mass, 2285.766, 2)
        peptide_mass = seq.peptide_composition().mass
        stub_map = {f.name: f for f in seq.stub_fragments()}
        bare_peptide = stub_map['peptide']
        self.assertAlmostEqual(bare_peptide.mass, peptide_mass, 2)

    def test_cad_fragmentation(self):
        t1 = self.parse_sequence(p10)
        low_energy_frags = list(t1.glycan_fragments(oxonium=False, all_series=True))
        assert len(low_energy_frags) == 431

    def test_glycan_representations(self):
        t1 = self.parse_sequence(p9)
        t2 = self.parse_sequence(p10)
        t3 = self.parse_sequence(p11)
        self.assertAlmostEqual(t1.total_mass, t2.total_mass, 4)
        self.assertAlmostEqual(t1.total_mass, t3.total_mass, 4)
        self.assertEqual(t1.clone().deglycosylate(), t2.clone().deglycosylate())
        self.assertEqual(t1.clone().deglycosylate(), t3.clone().deglycosylate())

    def test_glycan_tracking(self):
        t1 = self.parse_sequence("PEPT(HexNAc)IDE")
        t2 = self.parse_sequence("PEPT(N-Glycosylation)IDE{HexNAc:1}")
        self.assertAlmostEqual(t1.total_mass, t2.total_mass, 3)
        self.assertNotAlmostEqual(t1.peptide_composition().mass, t2.peptide_composition().mass, 3)

    def test_equalities(self):
        t1 = self.parse_sequence("PEPT(HexNAc)IDE")
        t2 = self.parse_sequence("PEPT(N-Glycosylation)IDE{HexNAc:1}")
        self.assertTrue(t1.base_sequence_equality(t2))
        self.assertFalse(t1.modified_sequence_equality(t2))

    def test_terminal_modification(self):
        base = self.parse_sequence("QDQC(Carbamidomethyl)IYNTTYLNVQR")
        self.assertAlmostEqual(base.mass, 1914.8894, 3)
        modified = base.clone()
        mt = modification.ModificationTable()
        pyro_glu_from_q = mt['Pyro-glu from Q']
        sites = pyro_glu_from_q.find_valid_sites(modified)
        modified.add_modification(sites[0], pyro_glu_from_q())
        self.assertAlmostEqual(modified.mass, 1897.86282353323, 3)

        base = self.parse_sequence("QDQC(Carbamidomethyl)IYNTTYLNVQR")
        p1 = base.clone()
        p1.add_modification(0, modification.Modification("TMT6plex"))
        p2 = base.clone()
        p2.n_term = base.n_term.modify(modification.Modification("TMT6plex"))
        self.assertAlmostEqual(p1.total_mass - p2.total_mass, 0)


class TestPeptideSequence(PeptideSequenceSuiteBase, unittest.TestCase):
    def parse_sequence(self, seqstr):
        return sequence.PeptideSequence(seqstr)


class TestNamedSequence(PeptideSequenceSuiteBase, unittest.TestCase):
    def parse_sequence(self, seqstr):
        return sequence.NamedSequence("spam", seqstr)


class TestAnnotatedSequence(PeptideSequenceSuiteBase, unittest.TestCase):
    def parse_sequence(self, seqstr):
        return sequence.AnnotatedSequence("spam", seqstr, annotations={"color": "green"})


if __name__ == '__main__':
    unittest.main()
