import unittest
import pickle

from glycopeptidepy.structure import sequence, modification, residue, composition, fragmentation_strategy
from glycopeptidepy.utils.collectiontools import groupby
from glypy import GlycanComposition, Glycan, MonosaccharideResidue


class OxoniumIonStrategyTest(unittest.TestCase):
    def test_no_multi_neuac(self):
        gp2 = sequence.PeptideSequence(
            'YLGN(#:glycosylation_type=N-Linked:{Hex:5; HexNAc:4; Neu5Ac:2})ATAIFFLPDEGK')
        for oxonium_ion in gp2.glycan_fragments():
            assert "Neu5AcNeu5Ac" not in oxonium_ion.name


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

    def test_gag_linker_stubs(self):
        gp = sequence.PeptideSequence("ISAS(GAG-Linker)GVEDISR{Xyl:1; a,en-Hex:1; Hex:2; a-Hex:1; HexNAc(S):1}")
        frags = list(gp.stub_fragments())
        assert len(frags) == 5

    def test_fucosylation_mechanism_limited(self):
        gp = sequence.PeptideSequence("VLGFKPKPPKN(N-Glycosylation)ESLETYPLMMK{Fuc:2; Hex:8; HexNAc:6; Neu5Ac:4}")
        strategy = fragmentation_strategy.StubGlycopeptideStrategy(gp, extended=True)
        assert len(strategy.n_glycan_composition_fragments(strategy.glycan_composition(), 1, 0)) == 69
        assert len(strategy.n_glycan_composition_fragments(strategy.glycan_composition(), 1, 1)) == 69
        assert len(strategy.n_glycan_composition_fragments(strategy.glycan_composition(), 1, 2)) == 35


class EXDFragmentationStrategyTest(unittest.TestCase):
    def test_approximated(self):
        gp2 = sequence.PeptideSequence('YLGN(#:glycosylation_type=N-Linked:{Hex:5; HexNAc:4; Neu5Ac:1})ATAIFFLPDEGK')
        gp3 = sequence.PeptideSequence((
            'YLGN(#:iupac,glycosylation_type=N-Linked:?-?-Hexp-(?-?)-?-?-Hexp2NAc-(?-?)-a-D-Manp-(1-6)-'
            '[a-D-Neup5Ac-(?-?)-?-?-Hexp-(?-?)-?-?-Hexp2NAc-(?-?)-a-D-Manp-(1-3)]b-D-Manp-(1-4)-b-D-Glcp'
            '2NAc-(1-4)-b-D-Glcp2NAc)ATAIFFLPDEGK'))
        seq2 = list(fragmentation_strategy.EXDFragmentationStrategy(gp2, 'c'))
        seq3 = list(fragmentation_strategy.EXDFragmentationStrategy(gp3, 'c'))
        assert len(seq2) == len(seq3)

        for exact_series, approximate_series in zip(seq3, seq2):
            exact_groups = groupby(exact_series, lambda x: round(x.mass))
            missed = 0
            for approx in approximate_series:
                key = round(approx.mass)
                missed += key not in exact_groups
            assert missed < 4


if __name__ == '__main__':
    unittest.main()
