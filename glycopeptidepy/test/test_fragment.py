import unittest

from glycopeptidepy.structure import sequence, modification, fragment, composition


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


if __name__ == '__main__':
    unittest.main()
