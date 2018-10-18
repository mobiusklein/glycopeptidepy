import unittest
import pickle

from glycopeptidepy.structure import sequence, modification, residue, composition, fragmentation_strategy
from glypy import GlycanComposition, Glycan, MonosaccharideResidue


class StubGlycopeptideStrategyTest(unittest.TestCase):
    def test_equivalence(self):
        gp2 = sequence.PeptideSequence('YLGN(#:glycosylation_type=N-Linked:{Hex:5; HexNAc:4; Neu5Ac:1})ATAIFFLPDEGK')
        gp3 = sequence.PeptideSequence((
            'YLGN(#:iupac,glycosylation_type=N-Linked:?-?-Hexp-(?-?)-?-?-Hexp2NAc-(?-?)-a-D-Manp-(1-6)-'
            '[a-D-Neup5Ac-(?-?)-?-?-Hexp-(?-?)-?-?-Hexp2NAc-(?-?)-a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp'
            '2NAc-(1-4)-b-D-Glcp2NAc)ATAIFFLPDEGK'))
        ref = list(fragmentation_strategy.StubGlycopeptideStrategy(gp2, extended=True))
        case = list(fragmentation_strategy.StubGlycopeptideStrategy(gp3, extended=True))
        self.assertEqual(ref, case)


if __name__ == '__main__':
    unittest.main()
