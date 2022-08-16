import unittest

from glycopeptidepy.structure import sequence

from glycopeptidepy.test.sequence_test_suite import PeptideSequenceSuiteBase, TestSequenceParser


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
