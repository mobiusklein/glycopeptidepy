import os
import unittest

from glycopeptidepy.io import fasta

from .common import datafile


n_proteins_in_fasta = 244


class TestFastaIO(unittest.TestCase):
    def open_stream(self):
        return open(datafile("proteins.fa"))

    def test_fasta_parser(self):
        with self.open_stream() as f:
            parser = fasta.ProteinFastaFileParser(f)
            i = 0
            for protein in parser:
                i += 1
                self.assertEqual(protein.name.db, 'sp')
            self.assertEqual(i, n_proteins_in_fasta)
