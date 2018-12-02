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

    def test_fucosylation(self):
        gp = sequence.PeptideSequence(
            'YLGN(#:glycosylation_type=N-Linked:{Fuc:2; Hex:5; HexNAc:4; Neu5Ac:1})ATAIFFLPDEGK')
        fucosylated = []
        for frag in gp.stub_fragments(extended=True, extended_fucosylation=False):
            if frag.glycosylation['Fuc'] > 0:
                fucosylated.append(frag)
        assert fucosylated

        fucosylated = []
        for frag in gp.stub_fragments(extended=True, extended_fucosylation=True):
            if frag.glycosylation['Fuc'] > 1:
                fucosylated.append(frag)
        assert fucosylated

    def test_fucosylation_mechanism_limited(self):
        gp = sequence.PeptideSequence("VLGFKPKPPKN(N-Glycosylation)ESLETYPLMMK{Fuc:2; Hex:8; HexNAc:6; Neu5Ac:4}")
        strategy = fragmentation_strategy.StubGlycopeptideStrategy(gp, extended=True)
        assert len(strategy.n_glycan_composition_fragments(strategy.glycan_composition(), 1, 0)) == 69
        assert len(strategy.n_glycan_composition_fragments(strategy.glycan_composition(), 1, 1)) == 69
        assert len(strategy.n_glycan_composition_fragments(strategy.glycan_composition(), 1, 2)) == 35


if __name__ == '__main__':
    unittest.main()
