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

>sp|P02763|A1AG1_HUMAN Alpha-1-acid glycoprotein 1 OS=Homo sapiens GN=ORM1 PE=1 SV=1
MALSWVLTVLSLLPLLEAQIPLCANLVPVPITNATLDQITGKWFYIASAFRNEEYNKSVQ
EIQATFFYFTPNKTEDTIFLREYQTRQDQCIYNTTYLNVQRENGTISRYVGGQEHFAHLL
ILRDTKTYMLAFDVNDEKNWGLSVYADKPETTKEQLGEFYEALDCLRIPKSDVVYTDWKK
DKCEPLEKQHEKERKQEEGES

'''

parser = ProteinFastaFileParser(StringIO(_protein_fasta))

heparanase = next(parser)
agp1 = next(parser)


class TestProtease(unittest.TestCase):
    def test_digest(self):
        trypsin = enzyme.Protease("trypsin")
        example_found = False
        for peptide, start, stop, missed in trypsin.cleave(heparanase, 2):
            assert missed < 3
            if peptide == "KFKNSTYSR":
                example_found = True
        if not example_found:
            raise AssertionError("Did not produce overlapped digest")

    def test_semispecific(self):
        trypsin = enzyme.Protease("trypsin")
        c_term_example = False
        n_term_example = False
        for peptide, start, end, missed in trypsin.cleave(agp1, semispecific=True):
            assert missed < 3
            if peptide == "LVPVPITNATLDQITGK":
                c_term_example = True
            if peptide == 'ENGT':
                n_term_example = True

        if not (n_term_example and c_term_example):
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
