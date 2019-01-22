import unittest
import warnings
from io import BytesIO
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from glycopeptidepy.io import fasta

from glycopeptidepy.test.common import datafile


class TestFastaIO(unittest.TestCase):
    def open_stream(self):
        return open(datafile("proteins.fa"), 'rb')

    def test_fasta_parser(self):
        n_proteins_in_fasta = 243
        with self.open_stream() as f:
            parser = fasta.ProteinFastaFileReader(f, index=True)
            assert len(parser.index) == n_proteins_in_fasta
            i = 0
            for protein in parser:
                i += 1
                self.assertEqual(protein.name.db, 'sp')
            assert i == n_proteins_in_fasta
            protein = parser['sp|Q9NYU2|UGGG1_HUMAN UDP-glucose:glycoprotein glucosyltransferase 1']
            assert str(protein).endswith("KREEL")

    def test_fasta_writer(self):
        outbuffer = BytesIO()
        writer = fasta.FastaFileWriter(outbuffer)
        writer.write("description", "PEPTIDE")
        writer.write("description2", "PROTEIN")
        content = outbuffer.getvalue()
        reference = b'''>description
PEPTIDE

>description2
PROTEIN

'''
        assert content == reference
        with self.open_stream() as f:
            parser = fasta.ProteinFastaFileReader(f)
            proteins = list(parser)
        outbuffer = BytesIO()
        writer = fasta.ProteinFastaFileWriter(outbuffer)
        writer.writelines(proteins)
        outbuffer.seek(0)
        new_proteins = list(fasta.ProteinFastaFileReader(outbuffer))
        assert proteins == new_proteins

    def test_replaces_unknown(self):
        reader = fasta.ProteinFastaFileReader(
            datafile("unknown_amino_acid.fa"), replace_unknown=True)
        contents = list(reader)
        assert len(contents) == 1
        reader = fasta.ProteinFastaFileReader(
            datafile("unknown_amino_acid.fa"), replace_unknown=False)
        with warnings.catch_warnings(record=True) as warning_log:
            warnings.simplefilter("always")
            should_fail = list(reader)
            assert len(should_fail) == 0
            assert len(warning_log) == 1
            assert "B" in str(warning_log[0].message)


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


class TestPEFF(unittest.TestCase):
    def open_stream(self, mode='rb'):
        return open(datafile("SmallTestDB-PEFF1.0.peff"), mode)

    def test_peff_parser(self):
        with self.open_stream() as f:
            parser = fasta.PEFFReader(f)
            assert parser.number_of_entries == 29
            proteins = list(parser)
            n = len(proteins)
            # this file was manually corrupted by its creator
            # to include a space in the first sequence
            assert n + 1 == parser.number_of_entries
            annotations = proteins[0].annotations
            assert len(annotations) == 7
            assert annotations['Length'] == 265


class DictWrapperSuiteBase(object):
    def _make_instance(self):
        raise NotImplementedError()

    def _make_reference(self):
        raise NotImplementedError()

    def test_getitem(self):
        inst = self._make_instance()
        ref = self._make_reference()
        for key in ref:
            assert inst[key] == ref[key]

    def test_getattr(self):
        inst = self._make_instance()
        ref = self._make_reference()
        for key in ref:
            assert getattr(inst, key) == ref[key]

    def test_contains(self):
        inst = self._make_instance()
        ref = self._make_reference()
        for key in ref:
            assert key in inst

    def test_keys(self):
        inst = self._make_instance()
        ref = self._make_reference()
        assert set(inst.keys()) == set(ref.keys())

    def test_values(self):
        inst = self._make_instance()
        ref = self._make_reference()
        assert sorted(inst.values(), key=str) == sorted(ref.values(), key=str)

    def test_items(self):
        inst = self._make_instance()
        ref = self._make_reference()
        assert sorted(inst.items()) == sorted(ref.items())

    def test_len(self):
        inst = self._make_instance()
        ref = self._make_reference()
        assert len(inst) == len(ref)


class TestFastaHeader(DictWrapperSuiteBase, unittest.TestCase):
    defline = ">sp|Q9NYU2|UGGG1_HUMAN UDP-glucose:glycoprotein glucosyltransferase 1"

    def _make_instance(self):
        return fasta.default_parser(self.defline)

    def _make_reference(self):
        ref = {'accession': 'Q9NYU2',
               'db': 'sp',
               'description': 'UDP-glucose:glycoprotein glucosyltransferase 1',
               'name': 'UGGG1_HUMAN'}
        return ref


class TestPEFFHeaderBlock(DictWrapperSuiteBase, unittest.TestCase):
    block = b'''# PEFF 1.0
# GeneralComment=this is a hand crafted peff file
# GeneralComment=This is a second comment line.
# //
# DbName=UniProtKB/Swiss-Prot-extract
# Prefix=sp
# DbDescription=selected entries from UniProtKB
# Decoy=false
# DbSource=ftp://ftp.uniprot.org//ftp/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
# DbVersion=55.5
# DbDate=20080610
# NumberOfEntries=28
# Conversion=manual extraction of entries, conversion with prepareDb.pl and then manually edited
# SequenceType=AA
# SpecificKey=3D-Status:status of a 3-D structure
# SpecificValue=3D-Status:(available|unsure|not available)
# SpecificKey=isoform:description of a specific isoform
# SpecificValue=isoform:(xsd:type=string)
# //'''

    def _make_instance(self):
        buff = BytesIO(self.block)
        reader = fasta.PEFFReader(buff)
        return reader.blocks[0]

    def _make_reference(self):
        return {u'Conversion': u'manual extraction of entries, conversion with prepareDb.pl and then manually edited',
                u'DbDate': u'20080610',
                u'DbDescription': u'selected entries from UniProtKB',
                u'DbName': u'UniProtKB/Swiss-Prot-extract',
                u'DbSource': (u'ftp://ftp.uniprot.org//ftp/pub/databases/uniprot/'
                              u'current_release/knowledgebase/complete/uniprot_sprot.dat.gz'),
                u'DbVersion': u'55.5',
                u'Decoy': u'false',
                u'NumberOfEntries': u'28',
                u'Prefix': u'sp',
                u'SequenceType': u'AA',
                u'SpecificKey': [u'3D-Status:status of a 3-D structure',
                                 u'isoform:description of a specific isoform'],
                u'SpecificValue': [u'3D-Status:(available|unsure|not available)',
                                   u'isoform:(xsd:type=string)']}


if __name__ == '__main__':
    unittest.main()
