import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from glycopeptidepy import enzyme
from glycopeptidepy.io.fasta import ProteinFastaFileParser

_protein_fasta = '''
>sp|Q9Y251|HPSE_HUMAN Heparanase OS=Homo sapiens GN=HPSE PE=1 SV=2
MLLRSKPALPPPLMLLLLGPLGPLSPGALPRPAQAQDVVDLDFFT
QEPLHLVSPSFLSVTIDANLATDPRFLILLGSPKLRTLARGLSPA
YLRFGGTKTDFLIFDPKKESTFEERSYWQSQVNQDICKYGSIPPD
VEEKLRLEWPYQEQLLLREHYQKKFKNSTYSRSSVDVLYTFANCS
GLDLIFGLNALLRTADLQWNSSNAQLLLDYCSSKGYNISWELGNE
PNSFLKKADIFINGSQLGEDFIQLHKLLRKSTFKNAKLYGPDVGQ
PRRKTAKMLKSFLKAGGEVIDSVTWHHYYLNGRTATKEDFLNPDV
LDIFISSVQKVFQVVESTRPGKKVWLGETSSAYGGGAPLLSDTFA
AGFMWLDKLGLSARMGIEVVMRQVFFGAGNYHLVDENFDPLPDYW
LSLLFKKLVGTKVLMASVQGSKRRKLRVYLHCTNTDNPRYKEGDL
TLYAINLHNVTKYLRLPYPFSNKQVDKYLLRPLGPHGLLSKSVQL
NGLTLKMVDDQTLPPLMEKPLRPGSSLGLPAFSYSFFVIRNAKVA
ACI
'''

heparanase = next(iter(ProteinFastaFileParser(StringIO(_protein_fasta))))


class TestProtease(unittest.TestCase):
    def test_digest(self):
        trypsin = enzyme.Protease("trypsin")
        for peptide, start, stop, missed in trypsin.cleave(heparanase, 2):
            assert missed < 3
            if peptide == "KFKNSTYSR":
                break
        else:
            raise AssertionError("Did not produce overlapped digest")

    def test_combine(self):
        trypsin = enzyme.Protease("trypsin")
        gluc = enzyme.Protease("glutamyl endopeptidase")
        seq = "PEPTIDER"
        self.assertEqual(trypsin.missed_cleavages(seq), 0)
        self.assertEqual(gluc.missed_cleavages(seq), 2)
        merged = enzyme.Protease.combine(trypsin, gluc)
        self.assertEqual(merged.missed_cleavages(seq), 2)


if __name__ == '__main__':
    unittest.main()
