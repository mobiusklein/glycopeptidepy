import unittest
import pickle

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


if __name__ == '__main__':
    unittest.main()
