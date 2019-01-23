import unittest

from glycopeptidepy.structure import sequence, modification, residue, sequence_composition


class TestSequenceComposition(unittest.TestCase):
    def test_parse(self):
        composition_string = '{N:2; K:1; E:2; Y:1}'
        composition = sequence_composition.SequenceComposition.parse(composition_string)
        seq = sequence.PeptideSequence("NEEYNK")
        self.assertAlmostEqual(seq.mass, composition.mass, 3)

    def test_serialize(self):
        composition = sequence_composition.SequenceComposition()
        composition['N'] = 2
        composition['K'] = 1
        composition['E'] = 2
        composition['Y'] = 1
        assert composition.serialize() == '{N:2; K:1; E:2; Y:1}'

    def test_clone(self):
        composition_string = '{N:2; K:1; E:2; Y:1}'
        composition = sequence_composition.SequenceComposition.parse(composition_string)
        assert composition == composition.clone()

    def test_modified(self):
        composition_string = '{N:1; N(N-Glycosylation):1; K:1; E:2; Y:1}'
        composition = sequence_composition.SequenceComposition.parse(composition_string)
        self.assertAlmostEqual(composition.mass, 998.4192, 3)


if __name__ == '__main__':
    unittest.main()
