import unittest
import pickle

import glypy

from glycopeptidepy.structure import sequence, fragment


class TestPeptideFragment(unittest.TestCase):
    def test_fragment_hashing(self):
        seq = sequence.PeptideSequence("LVPVPITN(N-Glycosylation)ATLDQITGK{Hex:5; HexNAc:4; Neu5Ac:2}")
        frag_map = dict()
        for position in seq.get_fragments('b'):
            for frag in position:
                frag_map[frag] = 1
        seq = sequence.PeptideSequence("LVPVPITN(N-Glycosylation)ATLDQITGK{Hex:5; HexNAc:4; Neu5Ac:2}")
        for position in seq.get_fragments('b'):
            for frag in position:
                self.assertEqual(frag_map[frag], 1)
                if "HexNAc" in frag.name:
                    assert frag.glycosylation_size > 0
                    assert frag.is_glycosylated
                else:
                    assert not frag.is_glycosylated
                    assert frag.name == frag.base_name()

    def test_pickle(self):
        seq = sequence.PeptideSequence("LVPVPITN(N-Glycosylation)ATLDQITGK{Hex:5; HexNAc:4; Neu5Ac:2}")
        for position in seq.get_fragments('b'):
            for frag in position:
                assert frag == frag.clone()
                assert frag == pickle.loads(pickle.dumps(frag))


class TestStubFragment(unittest.TestCase):
    def test_pickle(self):
        seq = sequence.PeptideSequence("LVPVPITN(N-Glycosylation)ATLDQITGK{Hex:5; HexNAc:4; Neu5Ac:2}")
        for frag in seq.stub_fragments():
            assert frag == frag.clone()
            assert frag == pickle.loads(pickle.dumps(frag))

    def test_naming(self):
        gc = glypy.GlycanComposition(Hex=5, HexNAc=4, NeuAc=2)
        gc.composition_offset = glypy.Composition()
        peptide = sequence.PeptideSequence("LVPVPITNATLDQITGK")
        name = fragment.StubFragment.build_name_from_composition(gc)
        frag = fragment.StubFragment(
            name, peptide.mass + gc.mass(),
            fragment.IonSeries.stub_glycopeptide,
            peptide.total_composition() + gc.total_composition(),
            glycosylation=gc)
        assert frag.name == "peptide+Hex5HexNAc4Neu5Ac2"
        assert frag.base_name() == "peptide+Hex5HexNAc4Neu5Ac2"
        assert abs(frag.mass - (peptide.mass + gc.mass())) < 1e-3
        shift = fragment.ChemicalShift("Na1H-1", glypy.Composition("Na1H-1"))
        frag.chemical_shift = shift
        assert frag.name == "peptide+Hex5HexNAc4Neu5Ac2Na1H-1"
        assert frag.base_name() == "peptide+Hex5HexNAc4Neu5Ac2"
        assert abs(frag.mass - (peptide.mass + gc.mass() + shift.mass)) < 1e-3


if __name__ == '__main__':
    unittest.main()
