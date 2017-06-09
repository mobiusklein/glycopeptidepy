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


class TestHeaderParsing(unittest.TestCase):
    def test_proper_uniprot(self):
        header = ">sp|P16066|ANPRA_HUMAN Atrial natriuretic peptide receptor 1"
        result = fasta.default_parser(header)
        self.assertEqual(result.accession, "P16066")

    def test_without_prefix(self):
        header = "sp|P16066|ANPRA_HUMAN Atrial natriuretic peptide receptor 1"
        result = fasta.default_parser(header)
        self.assertEqual(result.accession, "P16066")

    def test_partial_uniprot(self):
        header = "P16066|ANPRA_HUMAN"
        result = fasta.default_parser(header)
        self.assertEqual(result.accession, "P16066")

    def test_accession_alone(self):
        header = "P16066"
        result = fasta.default_parser(header)
        self.assertRaises(AttributeError, lambda: result.accession)
        self.assertEqual(result[0], "P16066")
