import unittest

try:
    from StringIO import StringIO
except:
    from io import StringIO

from glycopeptidepy.structure import sequence, modification, residue, composition
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
        for peptide, start, stop in trypsin.cleave(heparanase, 2):
            if peptide == "KFKNSTYSR":
                break
        else:
            raise AssertionError("Did not produce overlapped digest")


if __name__ == '__main__':
    unittest.main()
