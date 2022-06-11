from glypy import Composition

from glycopeptidepy.structure.modification import ModificationRule, ModificationTarget, SequenceLocation


rules = []

# [[[cog
# import cog
# from glycopeptidepy.io.cv import uniprot_ptm
# from glycopeptidepy.structure.modification.data.psimod import *
# terms = uniprot_ptm.UniProtPTMListParser().build_table().values()
# cog.out("rules = [\n")
# for src in filter(bool, map(uniprot_ptm.rule_to_source, terms)):
#     cog.out(f"    {src},\n")
# cog.out("]")
# ]]]
rules = [
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0450', monoisotopic_mass=-17.026549, composition=Composition(
        {'H': -3, 'N': -1}), alt_names={'RESID:AA0441', 'PSI-MOD:01624', 'PTM-0450', '(2-aminosuccinimidyl)acetic acid (Asn-Gly)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0312', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'RESID:AA0441', 'PSI-MOD:00952', '(2-aminosuccinimidyl)acetic acid (Asp-Gly)', 'PTM-0312'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0476', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'RESID:AA0601', '(3R)-3-hydroxyarginine', 'PTM-0476', 'PSI-MOD:01956'}),
    ModificationRule([ModificationTarget({'N'}, SequenceLocation['anywhere'])], 'PTM-0369', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PSI-MOD:00035', 'RESID:AA0026', 'PTM-0369', '(3R)-3-hydroxyasparagine'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0371', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PSI-MOD:00036', 'PTM-0371', '(3R)-3-hydroxyaspartate', 'RESID:AA0027'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0477', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'(3S)-3-hydroxyhistidine', 'PTM-0477', 'PSI-MOD:01920'}),
    ModificationRule([ModificationTarget({'P'}, SequenceLocation['anywhere'])], 'PTM-0368', monoisotopic_mass=31.989829,
                     composition=Composition({'O': 2}), alt_names={'PSI-MOD:01402', '(3R,4R)-3,4-dihydroxyproline', 'RESID:AA0479', 'PTM-0368'}),
    ModificationRule([ModificationTarget({'I'}, SequenceLocation['anywhere'])], 'PTM-0336', monoisotopic_mass=31.989829, composition=Composition(
        {'O': 2}), alt_names={'PTM-0336', '(3R,4R)-4,5-dihydroxyisoleucine', 'PSI-MOD:01378', 'RESID:AA0449'}),
    ModificationRule([ModificationTarget({'P'}, SequenceLocation['anywhere'])], 'PTM-0306', monoisotopic_mass=31.989829,
                     composition=Composition({'O': 2}), alt_names={'PTM-0306', 'RESID:AA0282', '(3R,4S)-3,4-dihydroxyproline', 'PSI-MOD:00287'}),
    ModificationRule([ModificationTarget({'I'}, SequenceLocation['anywhere'])], 'PTM-0337', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'RESID:AA0448', '(3R,4S)-4-hydroxyisoleucine', 'PTM-0337', 'PSI-MOD:01377'}),
    ModificationRule([ModificationTarget({'N'}, SequenceLocation['anywhere'])], 'PTM-0370', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PTM-0370', 'PSI-MOD:01401', '(3S)-3-hydroxyasparagine', 'RESID:AA0478'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0473', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'(3S)-3-hydroxyaspartate', 'PTM-0473', 'PSI-MOD:01919'}),
    ModificationRule([ModificationTarget({'I'}, SequenceLocation['anywhere'])], 'PTM-0343', monoisotopic_mass=31.989829, composition=Composition(
        {'O': 2}), alt_names={'PTM-0343', '(3S,4R)-3,4-dihydroxyisoleucine', 'RESID:AA0447', 'PSI-MOD:01376'}),
    ModificationRule([ModificationTarget({'L'}, SequenceLocation['anywhere'])], 'PTM-0491', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PSI-MOD:01373', 'PTM-0491', 'RESID:AA0443', '(4R)-5-hydroxyleucine'}),
    ModificationRule([ModificationTarget({'L'}, SequenceLocation['anywhere'])], 'PTM-0340', monoisotopic_mass=31.989829,
                     composition=Composition({'O': 2}), alt_names={'PSI-MOD:01375', '(4R)-4,5-dihydroxyleucine', 'RESID:AA0445', 'PTM-0340'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0001', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'RESID:AA0365', 'PTM-0001', '(E)-2,3-didehydrotyrosine', 'PSI-MOD:00370'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0002', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'PTM-0002', 'PSI-MOD:00191', 'RESID:AA0183', '(Z)-2,3-didehydrotyrosine'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['c_term'])], 'PTM-0379', monoisotopic_mass=-46.005479, composition=Composition(
        {'C': -1, 'H': -2, 'O': -2}), alt_names={'PTM-0379', 'PSI-MOD:01433', 'RESID:AA0486', '1-amino-2-propanone'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0003', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'PTM-0003', 'PSI-MOD:00275', 'RESID:AA0270', "1'-histidyl-3'-tyrosine (His-Tyr)"}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['anywhere'])], 'PTM-0004', monoisotopic_mass=15.977156,
                     composition=Composition({'O': -1, 'S': 1}), alt_names={'RESID:AA0265', 'PSI-MOD:00270', '1-thioglycine', 'PTM-0004'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0445', monoisotopic_mass=-2.01565, composition=Composition(
        {'H': -2}), alt_names={'1-(tryptophan-3-yl)-tryptophan (Trp-Trp) (interchain with W-...)', 'RESID:AA0544', 'PTM-0445', 'PSI-MOD:01818'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0005', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'PTM-0005', 'RESID:AA0109', "2'-(S-cysteinyl)-histidine (Cys-His)", 'PSI-MOD:00118'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0468', monoisotopic_mass=-33.987721, composition=Composition(
        {'H': -2, 'S': -1}), alt_names={'2,3-didehydroalanine (Cys)', 'PSI-MOD:00793', 'RESID:AA0181', 'PTM-0468'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0006', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'RESID:AA0181', 'PTM-0006', '2,3-didehydroalanine (Ser)', 'PSI-MOD:00189'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0441', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'(E)-2,3-didehydrobutyrine', 'RESID:AA0547', 'PSI-MOD:01841', 'PTM-0441'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0440', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'PTM-0440', '(Z)-2,3-didehydrobutyrine', 'RESID:AA0182', 'PSI-MOD:00190'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0007', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'PTM-0007', 'RESID:AA0547', '2,3-didehydrobutyrine', 'RESID:AA0182', 'PSI-MOD:00190', 'PSI-MOD:01841'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0008', monoisotopic_mass=-
                     2.01565, composition=Composition({'H': -2}), alt_names={'PSI-MOD:00706', '2,3-didehydrotyrosine', 'PTM-0008'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0009', monoisotopic_mass=29.974179,
                     composition=Composition({'H': -2, 'O': 2}), alt_names={'PSI-MOD:00156', 'PTM-0009', 'RESID:AA0147', "2',4',5'-topaquinone"}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0338', monoisotopic_mass=29.974179, composition=Composition(
        {'H': -2, 'O': 2}), alt_names={"2'-cysteinyl-6'-hydroxytryptophan sulfoxide (Trp-Cys)", 'RESID:AA0451', 'PTM-0338', 'PSI-MOD:01380'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0448', monoisotopic_mass=-5.062935, composition=Composition(
        {'H': -7, 'N': -1, 'O': 1}), alt_names={'RESID:AA0550', '2-(3-methylbutanoyl)-5-hydroxyoxazole-4-carbothionic acid (Leu-Cys)', 'PSI-MOD:01844', 'PTM-0448'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0457', monoisotopic_mass=-6.04695, composition=Composition(
        {'H': -6}), alt_names={'2-(4-guanidinobutanoyl)-5-hydroxyimidazole-4-carbothionic acid (Arg-Cys)', 'PSI-MOD:01877', 'PTM-0457', 'RESID:AA0553'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0424', monoisotopic_mass=167.982375, composition=Composition(
        {'C': 3, 'H': 5, 'O': 6, 'P': 1}), alt_names={'PSI-MOD:00797', 'RESID:AA0391', '2-(S-cysteinyl)pyruvic acid O-phosphothioketal', 'PTM-0424'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0010', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'PSI-MOD:00347', '2-cysteinyl-D-allo-threonine (Cys-Thr)', 'RESID:AA0342', 'PTM-0010'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0011', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'PTM-0011', 'RESID:AA0341', '2-cysteinyl-D-phenylalanine (Cys-Phe)', 'PSI-MOD:00346'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0012', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'2-cysteinyl-L-phenylalanine (Cys-Phe)', 'PSI-MOD:00345', 'PTM-0012', 'RESID:AA0340'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0495', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'RESID:AA0559', '2-(S-cysteinyl)-methionine (Cys-Met)', 'PSI-MOD:01857', 'PTM-0495'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0013', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'2-iminomethyl-5-imidazolinone (Gln-Gly)', 'PTM-0013', 'PSI-MOD:00197', 'RESID:AA0189'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0014', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'PTM-0014', '2-iminomethyl-5-imidazolinone (Glu-Gly)', 'RESID:AA0378', 'PSI-MOD:00383'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0015', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'RESID:AA0379', 'PTM-0015', 'PSI-MOD:00384', '2-iminomethyl-5-imidazolinone (Met-Gly)'}),
    ModificationRule([ModificationTarget({'Q'}, SequenceLocation['anywhere'])], 'PTM-0016', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'RESID:AA0273', 'PSI-MOD:00278', '2-methylglutamine', 'PTM-0016'}),
    ModificationRule([ModificationTarget({'W'}, SequenceLocation['anywhere'])], 'PTM-0304', monoisotopic_mass=77.97755, composition=Composition(
        {'C': 1, 'H': 2, 'O': 2, 'S': 1}), alt_names={"2'-methylsulfonyltryptophan", 'RESID:AA0450', 'PSI-MOD:01379', 'PTM-0304'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['n_term'])], 'PTM-0017', monoisotopic_mass=-17.026549,
                     composition=Composition({'H': -3, 'N': -1}), alt_names={'PSI-MOD:00138', 'RESID:AA0129', '2-oxobutanoic acid', 'PTM-0017'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0018', monoisotopic_mass=-37.052764, composition=Composition(
        {'H': -7, 'N': -1, 'O': -1}), alt_names={'PTM-0018', 'RESID:AA0382', '2-tetrahydro-2-pyridyl-5-imidazolinone (Lys-Gly)', 'PSI-MOD:00387'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0390', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'PTM-0390', 'RESID:AA0490', "3-(O4'-tyrosyl)-valine (Val-Tyr)", 'PSI-MOD:01442'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0019', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={"3'-(S-cysteinyl)-tyrosine (Cys-Tyr)", 'RESID:AA0113', 'PTM-0019', 'PSI-MOD:00122'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0020', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'RESID:AA0396', 'PTM-0020', 'PSI-MOD:00803', '3-(S-cysteinyl)-tyrosine (Cys-Tyr)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0021', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'RESID:AA0113', 'PSI-MOD:00122', 'PTM-0021', "3'-(S-cysteinyl)-tyrosine (Tyr-Cys)"}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0022', monoisotopic_mass=31.989829,
                     composition=Composition({'O': 2}), alt_names={'PSI-MOD:00374', '3,4-dihydroxyarginine', 'PTM-0022', 'RESID:AA0369'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0023', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'RESID:AA0146', 'PSI-MOD:00155', "3',4'-dihydroxyphenylalanine", 'PTM-0023'}),
    ModificationRule([ModificationTarget({'P'}, SequenceLocation['anywhere'])], 'PTM-0024', monoisotopic_mass=31.989829,
                     composition=Composition({'O': 2}), alt_names={'RESID:AA0282', 'PTM-0024', 'PSI-MOD:00287', '3,4-dihydroxyproline'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0025', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'RESID:AA0314', 'PSI-MOD:00319', 'PTM-0025', '3-cysteinyl-aspartic acid (Cys-Asp)'}),
    ModificationRule([ModificationTarget({'W'}, SequenceLocation['anywhere'])], 'PTM-0026', monoisotopic_mass=136.125201, composition=Composition(
        {'C': 10, 'H': 16}), alt_names={"3'-geranyl-2',N2-cyclotryptophan", 'PTM-0026', 'PSI-MOD:00817', 'RESID:AA0408'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0027', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={"3'-histidyl-3-tyrosine (His-Tyr)", 'PSI-MOD:00255', 'PTM-0027', 'RESID:AA0250'}),
    ModificationRule([ModificationTarget({'N'}, SequenceLocation['anywhere'])], 'PTM-0028', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PSI-MOD:00035', 'RESID:AA0026', '3-hydroxyasparagine', 'PTM-0028'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0029', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'RESID:AA0027', 'PTM-0029', 'PSI-MOD:00036', '3-hydroxyaspartate'}),
    ModificationRule([ModificationTarget({'F'}, SequenceLocation['anywhere'])], 'PTM-0346', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'3-hydroxyphenylalanine', 'RESID:AA0462', 'PSI-MOD:01385', 'PTM-0346'}),
    ModificationRule([ModificationTarget({'P'}, SequenceLocation['anywhere'])], 'PTM-0030', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PSI-MOD:00038', 'PTM-0030', '3-hydroxyproline', 'RESID:AA0029'}),
    ModificationRule([ModificationTarget({'W'}, SequenceLocation['anywhere'])], 'PTM-0031', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PTM-0031', 'RESID:AA0322', '3-hydroxytryptophan', 'PSI-MOD:00327'}),
    ModificationRule([ModificationTarget({'V'}, SequenceLocation['anywhere'])], 'PTM-0347', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'3-hydroxyvaline', 'RESID:AA0463', 'PTM-0347', 'PSI-MOD:01386'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0032', monoisotopic_mass=45.987721,
                     composition=Composition({'C': 1, 'H': 2, 'S': 1}), alt_names={'RESID:AA0232', 'PSI-MOD:00237', 'PTM-0032', '3-methylthioaspartic acid'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0434', monoisotopic_mass=44.985078, composition=Composition(
        {'H': -1, 'N': 1, 'O': 2}), alt_names={'RESID:AA0537', 'PTM-0434', "3'-nitrotyrosine", 'PSI-MOD:01786'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0033', monoisotopic_mass=-17.992806, composition=Composition(
        {'H': -2, 'O': 1, 'S': -1}), alt_names={'PSI-MOD:00193', '3-oxoalanine (Cys)', 'PTM-0033', 'RESID:AA0185'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0034', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'3-oxoalanine (Ser)', 'PSI-MOD:00835', 'PTM-0034', 'RESID:AA0185'}),
    ModificationRule([ModificationTarget({'F'}, SequenceLocation['n_term'])], 'PTM-0035', monoisotopic_mass=0.984016, composition=Composition(
        {'H': -1, 'N': -1, 'O': 1}), alt_names={'3-phenyllactic acid', 'PTM-0035', 'RESID:AA0128', 'PSI-MOD:00137'}),
    ModificationRule([ModificationTarget({'L'}, SequenceLocation['anywhere'])], 'PTM-0372', monoisotopic_mass=47.984744,
                     composition=Composition({'O': 3}), alt_names={'PTM-0372', 'RESID:AA0480', 'PSI-MOD:01403', "4,5,5'-trihydroxyleucine"}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0036', monoisotopic_mass=31.989829,
                     composition=Composition({'O': 2}), alt_names={'PTM-0036', 'RESID:AA0370', '4,5-dihydroxylysine', 'PSI-MOD:00375'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0038', monoisotopic_mass=79.966331,
                     composition=Composition({'H': 1, 'O': 3, 'P': 1}), alt_names={'RESID:AA0033', '4-aspartylphosphate', 'PSI-MOD:00042', 'PTM-0038'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0039', monoisotopic_mass=43.989829,
                     composition=Composition({'C': 1, 'O': 2}), alt_names={'PTM-0039', 'PSI-MOD:00041', 'RESID:AA0032', '4-carboxyglutamate'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0040', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'4-cysteinyl-glutamic acid (Cys-Glu)', 'RESID:AA0315', 'PSI-MOD:00320', 'PTM-0040'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0041', monoisotopic_mass=27.958529, composition=Composition(
        {'H': -4, 'O': 2}), alt_names={"4'-cysteinyl-tryptophylquinone (Cys-Trp)", 'RESID:AA0313', 'PSI-MOD:00318', 'PTM-0041'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0042', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PTM-0042', 'PSI-MOD:00220', '4-hydroxyarginine', 'RESID:AA0215'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0453', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'4-hydroxyglutamate', 'PTM-0453', 'RESID:AA0487', 'PSI-MOD:01434'}),
    ModificationRule([ModificationTarget({'P'}, SequenceLocation['anywhere'])], 'PTM-0043', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PTM-0043', '4-hydroxyproline', 'PSI-MOD:00039', 'RESID:AA0030'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0373', monoisotopic_mass=24.0, composition=Composition(
        {'C': 2}), alt_names={'PSI-MOD:01391', 'RESID:AA0468', '5-(methoxymethyl)thiazole-4-carboxylic acid (Val-Cys)', 'PTM-0373'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['c_term'])], 'PTM-0406', monoisotopic_mass=143.058243, composition=Composition(
        {'C': 6, 'H': 9, 'N': 1, 'O': 3}), alt_names={'PTM-0406', '5-glutamyl 2-aminoadipic acid', 'RESID:AA0502', 'PSI-MOD:01605'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0403', monoisotopic_mass=197.045309, composition=Composition(
        {'C': 5, 'H': 12, 'N': 1, 'O': 5, 'P': 1}), alt_names={'PSI-MOD:00179', 'PTM-0403', '5-glutamyl glycerylphosphorylethanolamine', 'RESID:AA0170'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['c_term'])], 'PTM-0479', monoisotopic_mass=129.042593, composition=Composition(
        {'C': 5, 'H': 7, 'N': 1, 'O': 3}), alt_names={'5-glutamyl glutamate', 'PSI-MOD:01970', 'PTM-0479', 'RESID:AA0612'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['c_term'])], 'PTM-0407', monoisotopic_mass=128.094963, composition=Composition(
        {'C': 6, 'H': 12, 'N': 2, 'O': 1}), alt_names={'PTM-0407', '5-glutamyl N2-lysine', 'PSI-MOD:01608', 'RESID:AA0505'}),
    ModificationRule([ModificationTarget({'I'}, SequenceLocation['anywhere'])], 'PTM-0466', monoisotopic_mass=13.979265, composition=Composition(
        {'H': -2, 'O': 1}), alt_names={'PSI-MOD:01897', 'RESID:AA0473', '5-hydroxy-3-methylproline (Ile)', 'PTM-0466'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0044', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PTM-0044', 'RESID:AA0028', 'PSI-MOD:00037', '5-hydroxylysine'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0471', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'RESID:AA0028', 'PTM-0471', '(5R)-5-hydroxylysine', 'PSI-MOD:01925'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0472', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PSI-MOD:01918', '(5S)-5-hydroxylysine', 'PTM-0472'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0045', monoisotopic_mass=-18.010565,
                     composition=Composition({'H': -2, 'O': -1}), alt_names={'5-imidazolinone (Ala-Gly)', 'RESID:AA0187', 'PTM-0045', 'PSI-MOD:00195'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0046', monoisotopic_mass=-18.010565,
                     composition=Composition({'H': -2, 'O': -1}), alt_names={'PSI-MOD:00385', 'RESID:AA0380', 'PTM-0046', '5-imidazolinone (Asn-Gly)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0047', monoisotopic_mass=-18.010565,
                     composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0047', '5-imidazolinone (Cys-Gly)', 'RESID:AA0188', 'PSI-MOD:00196'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0048', monoisotopic_mass=-18.010565,
                     composition=Composition({'H': -2, 'O': -1}), alt_names={'5-imidazolinone (Lys-Gly)', 'PTM-0048', 'RESID:AA0381', 'PSI-MOD:00386'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0049', monoisotopic_mass=-18.010565,
                     composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0049', 'PSI-MOD:00192', 'RESID:AA0184', '5-imidazolinone (Ser-Gly)'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0050', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'PTM-0050', 'PSI-MOD:00277', 'RESID:AA0272', '5-methylarginine'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0461', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'PSI-MOD:01900', 'RESID:AA0571', 'PTM-0461', '5-methyloxazole-4-carboxylic acid (Cys-Thr)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0386', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'PSI-MOD:01397', 'RESID:AA0474', '5-methyloxazole-4-carboxylic acid (Ser-Thr)', 'PTM-0386'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0462', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'PSI-MOD:01901', '5-methyloxazole-4-carboxylic acid (Thr-Thr)', 'RESID:AA0572', 'PTM-0462'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0352', monoisotopic_mass=-6.010565, composition=Composition(
        {'C': 1, 'H': -2, 'O': -1}), alt_names={'RESID:AA0469', 'PTM-0352', '5-methylthiazole-4-carboxylic acid (Asn-Cys)', 'PSI-MOD:01392'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0465', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'PSI-MOD:01904', 'RESID:AA0575', '5-methyloxazoline-4-carboxylic acid (Ser-Thr)', 'PTM-0465'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0320', monoisotopic_mass=12.995249, composition=Composition(
        {}), alt_names={'PSI-MOD:01787', "5'-tyrosyl-5'-aminotyrosine (Tyr-Tyr) (interchain with Y-...)", 'RESID:AA0459', 'PTM-0320'}),
    ModificationRule([ModificationTarget({'W'}, SequenceLocation['anywhere'])], 'PTM-0051', monoisotopic_mass=77.910512,
                     composition=Composition({'Br': 1, 'H': -1}), alt_names={'PSI-MOD:00188', "6'-bromotryptophan", 'RESID:AA0179', 'PTM-0051'}),
    ModificationRule([ModificationTarget({'W'}, SequenceLocation['anywhere'])], 'PTM-0444', monoisotopic_mass=33.961028,
                     composition=Composition({'Cl': 1, 'H': -1}), alt_names={'RESID:AA0180', 'PSI-MOD:00886', "5'-chlorotryptophan", 'PTM-0444'}),
    ModificationRule([ModificationTarget({'W'}, SequenceLocation['anywhere'])], 'PTM-0052', monoisotopic_mass=33.961028,
                     composition=Composition({'Cl': 1, 'H': -1}), alt_names={'PTM-0052', "6'-chlorotryptophan", 'PSI-MOD:00886', 'RESID:AA0180'}),
    ModificationRule([ModificationTarget({'W'}, SequenceLocation['anywhere'])], 'PTM-0427', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={"7'-hydroxytryptophan", 'RESID:AA0520', 'PTM-0427', 'PSI-MOD:01664'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0053', monoisotopic_mass=541.061109, composition=Composition(
        {'C': 15, 'H': 21, 'N': 5, 'O': 13, 'P': 2}), alt_names={'PSI-MOD:00177', 'PTM-0053', 'ADP-ribosylarginine', 'RESID:AA0168'}),
    ModificationRule([ModificationTarget({'N'}, SequenceLocation['anywhere'])], 'PTM-0054', monoisotopic_mass=541.061109, composition=Composition(
        {'C': 15, 'H': 21, 'N': 5, 'O': 13, 'P': 2}), alt_names={'PTM-0054', 'PSI-MOD:00236', 'ADP-ribosylasparagine', 'RESID:AA0231'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0497', monoisotopic_mass=541.061109,
                     composition=Composition({'C': 15, 'H': 21, 'N': 5, 'O': 13, 'P': 2}), alt_names={'PTM-0497', 'ADP-ribosyl aspartic acid'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0055', monoisotopic_mass=541.061109, composition=Composition(
        {'C': 15, 'H': 21, 'N': 5, 'O': 13, 'P': 2}), alt_names={'PTM-0055', 'PSI-MOD:00178', 'ADP-ribosylcysteine', 'RESID:AA0169'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0056', monoisotopic_mass=541.061109, composition=Composition(
        {'C': 15, 'H': 21, 'N': 5, 'O': 13, 'P': 2}), alt_names={'PSI-MOD:00242', 'ADP-ribosylserine', 'PTM-0056', 'RESID:AA0237'}),
    ModificationRule([ModificationTarget({'A'}, SequenceLocation['c_term'])], 'PTM-0057', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'PSI-MOD:00090', 'PTM-0057', 'RESID:AA0081', 'Alanine amide'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0418', monoisotopic_mass=-17.026549, composition=Composition(
        {'H': -3, 'N': -1}), alt_names={'RESID:AA0515', 'PSI-MOD:01618', 'PTM-0418', 'Alanine isoaspartyl cyclopeptide (Ala-Asn)'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0059', monoisotopic_mass=-1.031634,
                     composition=Composition({'H': -3, 'N': -1, 'O': 1}), alt_names={'PSI-MOD:00130', 'RESID:AA0121', 'PTM-0059', 'Allysine'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0321', monoisotopic_mass=13.979265, composition=Composition(
        {'H': -2, 'O': 1}), alt_names={'Aminomalonic acid (Ser)', 'PTM-0321', 'RESID:AA0458', 'PSI-MOD:01384'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['c_term'])], 'PTM-0060', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'Arginine amide', 'RESID:AA0082', 'PSI-MOD:00091', 'PTM-0060'}),
    ModificationRule([ModificationTarget({'N'}, SequenceLocation['c_term'])], 'PTM-0062', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'Asparagine amide', 'PSI-MOD:00092', 'PTM-0062', 'RESID:AA0083'}),
    ModificationRule([ModificationTarget({'N'}, SequenceLocation['c_term'])], 'PTM-0335', monoisotopic_mass=386.110369, composition=Composition(
        {'C': 13, 'H': 19, 'N': 6, 'O': 6, 'P': 1}), alt_names={'PSI-MOD:00333', "Aspartic acid 1-[(3-aminopropyl)(5'-adenosyl)phosphono]amide", 'RESID:AA0328', 'PTM-0335'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['c_term'])], 'PTM-0063', monoisotopic_mass=-0.984016,
                     composition=Composition({'H': 1, 'N': 1, 'O': -1}), alt_names={'Aspartic acid 1-amide', 'PSI-MOD:00093', 'RESID:AA0084', 'PTM-0063'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0064', monoisotopic_mass=-15.994915,
                     composition=Composition({'O': -1}), alt_names={'PSI-MOD:00378', 'RESID:AA0373', 'PTM-0064', 'Aspartyl aldehyde'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0489', monoisotopic_mass=-17.026549, composition=Composition(
        {'H': -3, 'N': -1}), alt_names={'Isoaspartyl glycine isopeptide (Asn-Gly)', 'PSI-MOD:00135', 'RESID:AA0126', 'PTM-0489'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0490', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'RESID:AA0126', 'Isoaspartyl glycine isopeptide (Asp-Gly)', 'PTM-0490', 'PSI-MOD:01805'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0066', monoisotopic_mass=28.0313, composition=Composition(
        {'C': 2, 'H': 4}), alt_names={'Asymmetric dimethylarginine', 'PTM-0066', 'RESID:AA0068', 'PSI-MOD:00077'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0314', monoisotopic_mass=-43.989829,
                     composition=Composition({'C': -1, 'O': -2}), alt_names={'Beta-decarboxylated aspartate', 'RESID:AA0001', 'PTM-0314', 'PSI-MOD:00869'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0067', monoisotopic_mass=-18.010565,
                     composition=Composition({'H': -2, 'O': -1}), alt_names={'RESID:AA0112', 'Beta-methyllanthionine (Cys-Thr)', 'PTM-0067', 'PSI-MOD:00121'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0068', monoisotopic_mass=-18.010565,
                     composition=Composition({'H': -2, 'O': -1}), alt_names={'RESID:AA0112', 'Beta-methyllanthionine (Thr-Cys)', 'PSI-MOD:00121', 'PTM-0068'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0069', monoisotopic_mass=-2.01565, composition=Composition(
        {'H': -2}), alt_names={'Beta-methyllanthionine sulfoxide (Thr-Cys)', 'PTM-0069', 'RESID:AA0330', 'PSI-MOD:00335'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0089', monoisotopic_mass=77.910512,
                     composition=Composition({'Br': 1, 'H': -1}), alt_names={'RESID:AA0173', 'PTM-0089', 'PSI-MOD:00182', 'Bromohistidine'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['c_term'])], 'PTM-0090', monoisotopic_mass=368.344301, composition=Composition(
        {'C': 27, 'H': 44}), alt_names={'Cholesterol glycine ester', 'RESID:AA0309', 'PTM-0090', 'PSI-MOD:00314'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0091', monoisotopic_mass=294.183109, composition=Composition(
        {'C': 17, 'H': 26, 'O': 4}), alt_names={'PSI-MOD:00321', 'Cis-14-hydroxy-10,13-dioxo-7-heptadecenoic acid aspartate ester', 'PTM-0091', 'RESID:AA0316'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0092', monoisotopic_mass=0.984016,
                     composition=Composition({'H': -1, 'N': -1, 'O': 1}), alt_names={'RESID:AA0214', 'PSI-MOD:00219', 'Citrulline', 'PTM-0092'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0093', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0093', 'Cyclopeptide (Ala-Arg)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0315', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'Cyclopeptide (Ala-Ile)', 'PTM-0315'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0307', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0307', 'Cyclopeptide (Ala-Pro)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0095', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'Cyclopeptide (Asn-Gly)', 'PTM-0095'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0096', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0096', 'Cyclopeptide (Asp-Asn)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0496', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'Cyclopeptide (Cys-Ile)', 'PTM-0496'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0098', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'Cyclopeptide (Gly-Asn)', 'PTM-0098'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0099', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'Cyclopeptide (Gly-Asp)', 'PTM-0099'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0436', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'Cyclopeptide (His-Asn)', 'PTM-0436'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0437', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0437', 'Cyclopeptide (His-Asp)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0339', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'Cyclopeptide (Ile-Pro)', 'PTM-0339'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0353', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0353', 'Cyclopeptide (Leu-Leu)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0319', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0319', 'Cyclopeptide (Leu-Trp)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0100', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0100', 'Cyclopeptide (Lys-Asp)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0316', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0316', 'Cyclopeptide (Pro-Met)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0317', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0317', 'Cyclopeptide (Pro-Tyr)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0101', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0101', 'Cyclopeptide (Ser-Asn)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0318', monoisotopic_mass=-
                     18.010565, composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0318', 'Cyclopeptide (Ser-Gly)'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['c_term'])], 'PTM-0433', monoisotopic_mass=103.009185, composition=Composition(
        {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 1}), alt_names={'RESID:AA0529', 'CysO-cysteine adduct', 'PSI-MOD:01778', 'PTM-0433'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['c_term'])], 'PTM-0102', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'PSI-MOD:00094', 'Cysteine amide', 'PTM-0102', 'RESID:AA0085'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0104', monoisotopic_mass=45.987721, composition=Composition(
        {'C': 1, 'H': 2, 'S': 1}), alt_names={'PTM-0104', 'Cysteine methyl disulfide', 'RESID:AA0101', 'PSI-MOD:00110'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['c_term'])], 'PTM-0105', monoisotopic_mass=14.01565, composition=Composition(
        {'C': 1, 'H': 2}), alt_names={'PSI-MOD:00114', 'Cysteine methyl ester', 'PTM-0105', 'RESID:AA0105'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0106', monoisotopic_mass=31.972071,
                     composition=Composition({'S': 1}), alt_names={'Cysteine persulfide', 'RESID:AA0269', 'PSI-MOD:00274', 'PTM-0106'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0107', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'RESID:AA0205', 'PSI-MOD:00210', 'Cysteine sulfenic acid (-SOH)', 'PTM-0107'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0108', monoisotopic_mass=31.989829,
                     composition=Composition({'O': 2}), alt_names={'PTM-0108', 'RESID:AA0262', 'Cysteine sulfinic acid (-SO2H)', 'PSI-MOD:00267'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0109', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'PSI-MOD:00363', 'Cysteinyl-selenocysteine (Cys-Sec)', 'PTM-0109', 'RESID:AA0358'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0110', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'PSI-MOD:00363', 'Cysteinyl-selenocysteine (Sec-Cys)', 'RESID:AA0358', 'PTM-0110'}),
    ModificationRule([ModificationTarget({'V'}, SequenceLocation['anywhere'])], 'PTM-0111', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'RESID:AA0388', 'PTM-0111', 'D-4-hydroxyvaline', 'PSI-MOD:00756'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0113', monoisotopic_mass=-15.994915,
                     composition=Composition({'O': -1}), alt_names={'PSI-MOD:00858', 'D-alanine (Ser)', 'RESID:AA0191', 'PTM-0113'}),
    ModificationRule([ModificationTarget({'N'}, SequenceLocation['anywhere'])], 'PTM-0116', monoisotopic_mass=0.984016, composition=Composition(
        {'H': -1, 'N': -1, 'O': 1}), alt_names={'Deamidated asparagine', 'PTM-0116', 'RESID:AA0004', 'PSI-MOD:00684'}),
    ModificationRule([ModificationTarget({'Q'}, SequenceLocation['anywhere'])], 'PTM-0117', monoisotopic_mass=0.984016, composition=Composition(
        {'H': -1, 'N': -1, 'O': 1}), alt_names={'PSI-MOD:00685', 'RESID:AA0006', 'PTM-0117', 'Deamidated glutamine'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['c_term'])], 'PTM-0354', monoisotopic_mass=-43.989829, composition=Composition(
        {'C': -1, 'O': -2}), alt_names={'RESID:AA0465', 'PSI-MOD:01388', 'PTM-0354', 'Decarboxylated threonine'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0408', monoisotopic_mass=251.793295,
                     composition=Composition({'H': -2, 'I': 2}), alt_names={'Diiodotyrosine', 'PSI-MOD:01613', 'RESID:AA0510', 'PTM-0408'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0341', monoisotopic_mass=28.0313,
                     composition=Composition({'C': 2, 'H': 4}), alt_names={'Dimethylated arginine', 'PSI-MOD:00783', 'PTM-0341'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0118', monoisotopic_mass=143.11789, composition=Composition(
        {'C': 7, 'H': 15, 'N': 2, 'O': 1}), alt_names={'RESID:AA0040', 'Diphthamide', 'PTM-0118', 'PSI-MOD:00049'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0309', monoisotopic_mass=-15.977156,
                     composition=Composition({'O': 1, 'S': -1}), alt_names={'PTM-0309', 'PSI-MOD:00892', 'D-serine (Cys)', 'RESID:AA0195'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0125', monoisotopic_mass=438.09405, composition=Composition(
        {'C': 17, 'H': 19, 'N': 4, 'O': 8, 'P': 1}), alt_names={'PSI-MOD:00355', 'FMN phosphoryl serine', 'RESID:AA0350', 'PTM-0125'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0126', monoisotopic_mass=438.09405, composition=Composition(
        {'C': 17, 'H': 19, 'N': 4, 'O': 8, 'P': 1}), alt_names={'PTM-0126', 'RESID:AA0349', 'FMN phosphoryl threonine', 'PSI-MOD:00354'}),
    ModificationRule([ModificationTarget({'Q'}, SequenceLocation['anywhere'])], 'PTM-0127', monoisotopic_mass=14.999666, composition=Composition(
        {'C': 1, 'H': 1, 'N': -1, 'O': 1}), alt_names={'PTM-0127', 'Glutamate methyl ester (Gln)', 'PSI-MOD:00657', 'RESID:AA0072'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0128', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'Glutamate methyl ester (Glu)', 'PSI-MOD:00081', 'PTM-0128', 'RESID:AA0072'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['c_term'])], 'PTM-0129', monoisotopic_mass=-0.984016,
                     composition=Composition({'H': 1, 'N': 1, 'O': -1}), alt_names={'PTM-0129', 'PSI-MOD:00096', 'RESID:AA0087', 'Glutamic acid 1-amide'}),
    ModificationRule([ModificationTarget({'Q'}, SequenceLocation['c_term'])], 'PTM-0130', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'PSI-MOD:00095', 'RESID:AA0086', 'PTM-0130', 'Glutamine amide'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['c_term'])], 'PTM-0132', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'PSI-MOD:00097', 'RESID:AA0088', 'PTM-0132', 'Glycine amide'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['c_term'])], 'PTM-0409', monoisotopic_mass=329.05252, composition=Composition(
        {'C': 10, 'H': 12, 'N': 5, 'O': 6, 'P': 1}), alt_names={'RESID:AA0511', 'PSI-MOD:01614', 'PTM-0409', 'Glycyl adenylate'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0413', monoisotopic_mass=13.961506, composition=Composition(
        {'H': -2, 'O': -1, 'S': 1}), alt_names={'PTM-0413', 'RESID:AA0512', 'Glycyl cysteine dithioester (Cys-Gly) (interchain with G-...)', 'PSI-MOD:01615'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0410', monoisotopic_mass=13.961506, composition=Composition(
        {'H': -2, 'O': -1, 'S': 1}), alt_names={'Glycyl cysteine dithioester (Gly-Cys) (interchain with C-...)', 'RESID:AA0512', 'PTM-0410', 'PSI-MOD:01615'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0422', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'RESID:AA0495', 'Glycyl serine ester (Ser-Gly) (interchain with G-...)', 'PTM-0422', 'PSI-MOD:01585'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0425', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'PTM-0425', 'RESID:AA0495', 'Glycyl serine ester (Gly-Ser) (interchain with S-...)', 'PSI-MOD:01585'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0423', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'PSI-MOD:01586', 'RESID:AA0496', 'PTM-0423', 'Glycyl threonine ester (Thr-Gly) (interchain with G-...)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0426', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'PSI-MOD:01586', 'PTM-0426', 'RESID:AA0496', 'Glycyl threonine ester (Gly-Thr) (interchain with T-...)'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['c_term'])], 'PTM-0148', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'PTM-0148', 'PSI-MOD:00098', 'Histidine amide', 'RESID:AA0089'}),
    ModificationRule([ModificationTarget({'P'}, SequenceLocation['anywhere'])], 'PTM-0149', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PTM-0149', 'Hydroxyproline', 'PSI-MOD:00678'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0150', monoisotopic_mass=87.068414, composition=Composition(
        {'C': 4, 'H': 9, 'N': 1, 'O': 1}), alt_names={'PTM-0150', 'RESID:AA0116', 'PSI-MOD:00125', 'Hypusine'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0411', monoisotopic_mass=125.896648,
                     composition=Composition({'H': -1, 'I': 1}), alt_names={'PTM-0411', 'PSI-MOD:01612', 'RESID:AA0509', 'Iodotyrosine'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0151', monoisotopic_mass=-17.026549, composition=Composition(
        {'H': -3, 'N': -1}), alt_names={'PSI-MOD:00221', 'PTM-0151', 'RESID:AA0216', 'Isoaspartyl cysteine isopeptide (Cys-Asn)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0152', monoisotopic_mass=-17.026549, composition=Composition(
        {'H': -3, 'N': -1}), alt_names={'PSI-MOD:00135', 'Isoaspartyl glycine isopeptide (Gly-Asn)', 'RESID:AA0126', 'PTM-0152'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0330', monoisotopic_mass=-17.026549, composition=Composition(
        {'H': -3, 'N': -1}), alt_names={'Isoaspartyl lysine isopeptide (Lys-Asn)', 'PSI-MOD:00299', 'RESID:AA0294', 'PTM-0330'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0486', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'PTM-0486', 'Isoaspartyl lysine isopeptide (Lys-Asp)', 'RESID:AA0294', 'PSI-MOD:01917'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0155', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'PSI-MOD:00373', 'RESID:AA0368', 'PTM-0155', 'Isodityrosine (Tyr-Tyr)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0156', monoisotopic_mass=-17.026549, composition=Composition(
        {'H': -3, 'N': -1}), alt_names={'PTM-0156', 'RESID:AA0108', 'PSI-MOD:00117', 'Isoglutamyl cysteine thioester (Cys-Gln)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0157', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'PSI-MOD:00365', 'PTM-0157', 'Isoglutamyl glycine isopeptide (Gly-Glu)', 'RESID:AA0360'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0159', monoisotopic_mass=-17.026549, composition=Composition(
        {'H': -3, 'N': -1}), alt_names={'PSI-MOD:00133', 'PTM-0159', 'Isoglutamyl lysine isopeptide (Lys-Gln)', 'RESID:AA0124'}),
    ModificationRule([ModificationTarget({'I'}, SequenceLocation['c_term'])], 'PTM-0161', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'PSI-MOD:00099', 'RESID:AA0090', 'PTM-0161', 'Isoleucine amide'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['n_term'])], 'PTM-0163', monoisotopic_mass=-15.010899,
                     composition=Composition({'H': -1, 'N': -1}), alt_names={'RESID:AA0186', 'PSI-MOD:00194', 'PTM-0163', 'Lactic acid'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0164', monoisotopic_mass=-18.010565,
                     composition=Composition({'H': -2, 'O': -1}), alt_names={'PSI-MOD:00120', 'PTM-0164', 'RESID:AA0110', 'Lanthionine (Cys-Ser)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0165', monoisotopic_mass=-18.010565,
                     composition=Composition({'H': -2, 'O': -1}), alt_names={'PSI-MOD:00120', 'PTM-0165', 'RESID:AA0110', 'Lanthionine (Ser-Cys)'}),
    ModificationRule([ModificationTarget({'L'}, SequenceLocation['c_term'])], 'PTM-0166', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'Leucine amide', 'RESID:AA0091', 'PTM-0166', 'PSI-MOD:00100'}),
    ModificationRule([ModificationTarget({'L'}, SequenceLocation['c_term'])], 'PTM-0167', monoisotopic_mass=14.01565, composition=Composition(
        {'C': 1, 'H': 2}), alt_names={'PTM-0167', 'PSI-MOD:00304', 'Leucine methyl ester', 'RESID:AA0299'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['c_term'])], 'PTM-0168', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'PTM-0168', 'RESID:AA0092', 'Lysine amide', 'PSI-MOD:00101'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['c_term'])], 'PTM-0170', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'RESID:AA0318', 'PSI-MOD:00323', 'Lysine methyl ester', 'PTM-0170'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0171', monoisotopic_mass=11.963614,
                     composition=Composition({'H': -4, 'O': 1}), alt_names={'Lysine tyrosylquinone (Lys-Tyr)', 'PTM-0171', 'PSI-MOD:00238', 'RESID:AA0233'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0323', monoisotopic_mass=11.963614,
                     composition=Composition({'H': -4, 'O': 1}), alt_names={'Lysine tyrosylquinone (Tyr-Lys)', 'PTM-0323', 'PSI-MOD:00238', 'RESID:AA0233'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0439', monoisotopic_mass=87.032028, composition=Composition(
        {'C': 3, 'H': 5, 'N': 1, 'O': 2}), alt_names={'PSI-MOD:01838', 'PTM-0439', 'RESID:AA0123', 'Lysino-D-alanine (Lys)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0172', monoisotopic_mass=-18.010565,
                     composition=Composition({'H': -2, 'O': -1}), alt_names={'PTM-0172', 'PSI-MOD:00132', 'RESID:AA0123', 'Lysinoalanine (Ser-Lys)'}),
    ModificationRule([ModificationTarget({'M'}, SequenceLocation['c_term'])], 'PTM-0173', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'RESID:AA0093', 'PSI-MOD:00102', 'Methionine amide', 'PTM-0173'}),
    ModificationRule([ModificationTarget({'M'}, SequenceLocation['anywhere'])], 'PTM-0175', monoisotopic_mass=31.989829,
                     composition=Composition({'O': 2}), alt_names={'PSI-MOD:00256', 'RESID:AA0251', 'PTM-0175', 'Methionine sulfone'}),
    ModificationRule([ModificationTarget({'M'}, SequenceLocation['anywhere'])], 'PTM-0480', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PTM-0480', 'RESID:AA0581', 'Methionine (R)-sulfoxide', 'PSI-MOD:00720'}),
    ModificationRule([ModificationTarget({'M'}, SequenceLocation['anywhere'])], 'PTM-0481', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PTM-0481', 'PSI-MOD:00721', 'Methionine (S)-sulfoxide'}),
    ModificationRule([ModificationTarget({'M'}, SequenceLocation['anywhere'])], 'PTM-0469', monoisotopic_mass=15.994915,
                     composition=Composition({'O': 1}), alt_names={'PSI-MOD:00719', 'Methionine sulfoxide', 'PTM-0469'}),
    ModificationRule([ModificationTarget({'L'}, SequenceLocation['anywhere'])], 'PTM-0492', monoisotopic_mass=13.979265,
                     composition=Composition({'H': -2, 'O': 1}), alt_names={'(4R)-5-oxoleucine', 'RESID:AA0444', 'PSI-MOD:01374', 'PTM-0492'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0176', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'PSI-MOD:00661', 'Methylhistidine', 'PTM-0176'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0451', monoisotopic_mass=-2.01565, composition=Composition(
        {'H': -2}), alt_names={'PTM-0451', 'PSI-MOD:01859', 'N,N-(cysteine-1,S-diyl)phenylalanine (Cys-Phe)', 'RESID:AA0562'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0037', monoisotopic_mass=-2.01565,
                     composition=Composition({'H': -2}), alt_names={'N,N-(cysteine-1,S-diyl)serine (Cys-Ser)', 'PTM-0037', 'RESID:AA0344', 'PSI-MOD:00349'}),
    ModificationRule([ModificationTarget({'A'}, SequenceLocation['n_term'])], 'PTM-0177', monoisotopic_mass=43.054227, composition=Composition(
        {'C': 3, 'H': 7}), alt_names={'PSI-MOD:00071', 'PTM-0177', 'N,N,N-trimethylalanine', 'RESID:AA0062'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['n_term'])], 'PTM-0485', monoisotopic_mass=43.054775, composition=Composition(
        {'C': 3, 'H': 7}), alt_names={'PTM-0485', 'RESID:AA0619', 'N,N,N-trimethylglycine', 'PSI-MOD:01982'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['n_term'])], 'PTM-0430', monoisotopic_mass=43.054227, composition=Composition(
        {'C': 3, 'H': 7}), alt_names={'N,N,N-trimethylserine', 'PTM-0430', 'PSI-MOD:01784', 'RESID:AA0535'}),
    ModificationRule([ModificationTarget({'A'}, SequenceLocation['n_term'])], 'PTM-0178', monoisotopic_mass=28.0313, composition=Composition(
        {'C': 2, 'H': 4}), alt_names={'RESID:AA0433', 'PTM-0178', 'PSI-MOD:01179', 'N,N-dimethylalanine'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['n_term'])], 'PTM-0484', monoisotopic_mass=28.0313, composition=Composition(
        {'C': 2, 'H': 4}), alt_names={'N,N-dimethylglycine', 'PSI-MOD:01983', 'PTM-0484', 'RESID:AA0620'}),
    ModificationRule([ModificationTarget({'L'}, SequenceLocation['n_term'])], 'PTM-0435', monoisotopic_mass=28.0313, composition=Composition(
        {'C': 2, 'H': 4}), alt_names={'PSI-MOD:01806', 'N,N-dimethylleucine', 'RESID:AA0538', 'PTM-0435'}),
    ModificationRule([ModificationTarget({'P'}, SequenceLocation['n_term'])], 'PTM-0179', monoisotopic_mass=29.038577,
                     composition=Composition({'C': 2, 'H': 5}), alt_names={'PTM-0179', 'N,N-dimethylproline', 'RESID:AA0066', 'PSI-MOD:00075'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['n_term'])], 'PTM-0431', monoisotopic_mass=28.0313,
                     composition=Composition({'C': 2, 'H': 4}), alt_names={'PTM-0431', 'N,N-dimethylserine', 'RESID:AA0534', 'PSI-MOD:01783'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['n_term'])], 'PTM-0180', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'PSI-MOD:00359', 'RESID:AA0354', 'N2-acetylarginine', 'PTM-0180'}),
    ModificationRule([ModificationTarget({'W'}, SequenceLocation['n_term'])], 'PTM-0181', monoisotopic_mass=100.016044, composition=Composition(
        {'C': 4, 'H': 4, 'O': 3}), alt_names={'N2-succinyltryptophan', 'RESID:AA0130', 'PSI-MOD:00139', 'PTM-0181'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['n_term'])], 'PTM-0459', monoisotopic_mass=28.0313, composition=Composition(
        {'C': 2, 'H': 4}), alt_names={'RESID:AA0569', 'N2,N2-dimethylarginine', 'PSI-MOD:01898', 'PTM-0459'}),
    ModificationRule([ModificationTarget({'N'}, SequenceLocation['anywhere'])], 'PTM-0182', monoisotopic_mass=28.0313, composition=Composition(
        {'C': 2, 'H': 4}), alt_names={'PTM-0182', 'RESID:AA0311', 'PSI-MOD:00316', 'N4,N4-dimethylasparagine'}),
    ModificationRule([ModificationTarget({'N'}, SequenceLocation['anywhere'])], 'PTM-0183', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'N4-methylasparagine', 'RESID:AA0070', 'PTM-0183', 'PSI-MOD:00079'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0184', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'N5-methylarginine', 'PSI-MOD:00310', 'PTM-0184', 'RESID:AA0305'}),
    ModificationRule([ModificationTarget({'Q'}, SequenceLocation['anywhere'])], 'PTM-0185', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'RESID:AA0071', 'PSI-MOD:00080', 'PTM-0185', 'N5-methylglutamine'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0355', monoisotopic_mass=541.061109, composition=Composition(
        {'C': 15, 'H': 21, 'N': 5, 'O': 13, 'P': 2}), alt_names={'PSI-MOD:01399', 'RESID:AA0476', 'N6-(ADP-ribosyl)lysine', 'PTM-0355'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0387', monoisotopic_mass=229.014009, composition=Composition(
        {'C': 8, 'H': 8, 'N': 1, 'O': 5, 'P': 1}), alt_names={'N6-(pyridoxal phosphate)lysine', 'RESID:AA0119', 'PSI-MOD:00128', 'PTM-0387'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0388', monoisotopic_mass=266.203451, composition=Composition(
        {'C': 20, 'H': 26}), alt_names={'N6-(retinylidene)lysine', 'PTM-0388', 'PSI-MOD:00129', 'RESID:AA0120'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0186', monoisotopic_mass=59.049141, composition=Composition(
        {'C': 3, 'H': 7, 'O': 1}), alt_names={'N6,N6,N6-trimethyl-5-hydroxylysine', 'PTM-0186', 'PSI-MOD:00364', 'RESID:AA0359'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0187', monoisotopic_mass=43.054227, composition=Composition(
        {'C': 3, 'H': 7}), alt_names={'N6,N6,N6-trimethyllysine', 'PSI-MOD:00083', 'RESID:AA0074', 'PTM-0187'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0188', monoisotopic_mass=28.0313,
                     composition=Composition({'C': 2, 'H': 4}), alt_names={'PSI-MOD:00084', 'N6,N6-dimethyllysine', 'PTM-0188', 'RESID:AA0075'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0189', monoisotopic_mass=72.021129, composition=Composition(
        {'C': 3, 'H': 4, 'O': 2}), alt_names={'N6-1-carboxyethyl lysine', 'RESID:AA0115', 'PTM-0189', 'PSI-MOD:00124'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0190', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'PTM-0190', 'N6-acetyllysine', 'PSI-MOD:00064', 'RESID:AA0055'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0382', monoisotopic_mass=226.077599, composition=Composition(
        {'C': 10, 'H': 14, 'N': 2, 'O': 2, 'S': 1}), alt_names={'RESID:AA0117', 'PSI-MOD:00126', 'PTM-0382', 'N6-biotinyllysine'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0191', monoisotopic_mass=43.989829,
                     composition=Composition({'C': 1, 'O': 2}), alt_names={'PTM-0191', 'N6-carboxylysine', 'PSI-MOD:00123', 'RESID:AA0114'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0475', monoisotopic_mass=68.026215, composition=Composition(
        {'C': 4, 'H': 4, 'O': 1}), alt_names={'RESID:AA0567', 'N6-crotonyllysine', 'PTM-0475', 'PSI-MOD:01892'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0192', monoisotopic_mass=27.994915,
                     composition=Composition({'C': 1, 'O': 1}), alt_names={'N6-formyllysine', 'RESID:AA0211', 'PTM-0192', 'PSI-MOD:00216'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0383', monoisotopic_mass=188.032957, composition=Composition(
        {'C': 8, 'H': 12, 'O': 1, 'S': 2}), alt_names={'RESID:AA0118', 'N6-lipoyllysine', 'PSI-MOD:00127', 'PTM-0383'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0467', monoisotopic_mass=86.000394, composition=Composition(
        {'C': 3, 'H': 2, 'O': 3}), alt_names={'PTM-0467', 'N6-malonyllysine', 'PSI-MOD:01893', 'RESID:AA0568'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0194', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'PTM-0194', 'RESID:AA0076', 'PSI-MOD:00085', 'N6-methyllysine'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0196', monoisotopic_mass=210.198365, composition=Composition(
        {'C': 14, 'H': 26, 'O': 1}), alt_names={'PTM-0196', 'RESID:AA0078', 'N6-myristoyl lysine', 'PSI-MOD:00087'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0197', monoisotopic_mass=238.229666, composition=Composition(
        {'C': 16, 'H': 30, 'O': 1}), alt_names={'RESID:AA0077', 'N6-palmitoyl lysine', 'PSI-MOD:00086', 'PTM-0197'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0198', monoisotopic_mass=426.440996, composition=Composition(
        {'C': 24, 'H': 54, 'N': 6}), alt_names={'N6-poly(methylaminopropyl)lysine', 'PTM-0198', 'PSI-MOD:00283', 'RESID:AA0278'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0438', monoisotopic_mass=100.016044, composition=Composition(
        {'C': 4, 'H': 4, 'O': 3}), alt_names={'N6-succinyllysine', 'PSI-MOD:01819', 'RESID:AA0545', 'PTM-0438'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0487', monoisotopic_mass=114.031694,
                     composition=Composition({'C': 5, 'H': 6, 'O': 3}), alt_names={'N6-glutaryllysine', 'PTM-0487'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0499', monoisotopic_mass=87.044604,
                     composition=Composition({'C': 4, 'H': 7, 'O': 2}), alt_names={'PTM-0499', 'N6-(beta-hydroxybutyrate)lysine'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['n_term'])], 'PTM-0419', monoisotopic_mass=226.19328, composition=Composition(
        {'C': 14, 'H': 26, 'O': 2}), alt_names={'RESID:AA0516', 'PSI-MOD:01690', 'PTM-0419', 'N-[(12R)-12-hydroxymyristoyl]cysteine'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['n_term'])], 'PTM-0420', monoisotopic_mass=224.17763, composition=Composition(
        {'C': 14, 'H': 24, 'O': 2}), alt_names={'PTM-0420', 'PSI-MOD:01691', 'N-(12-oxomyristoyl)cysteine', 'RESID:AA0517'}),
    ModificationRule([ModificationTarget({'A'}, SequenceLocation['n_term'])], 'PTM-0199', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'RESID:AA0041', 'PSI-MOD:00050', 'N-acetylalanine', 'PTM-0199'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['n_term'])], 'PTM-0200', monoisotopic_mass=42.010565,
                     composition=Composition({'C': 2, 'H': 2, 'O': 1}), alt_names={'PTM-0200', 'PSI-MOD:00051', 'N-acetylaspartate', 'RESID:AA0042'}),
    ModificationRule([ModificationTarget({'K'}, SequenceLocation['anywhere'])], 'PTM-0482', monoisotopic_mass=42.010565,
                     composition=Composition({'C': 2, 'H': 2, 'O': 1}), alt_names={'PTM-0482', 'N-acetylated lysine', 'PSI-MOD:00723'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['n_term'])], 'PTM-0201', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'RESID:AA0043', 'N-acetylcysteine', 'PSI-MOD:00052', 'PTM-0201'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['n_term'])], 'PTM-0202', monoisotopic_mass=42.010565,
                     composition=Composition({'C': 2, 'H': 2, 'O': 1}), alt_names={'PTM-0202', 'N-acetylglutamate', 'PSI-MOD:00053', 'RESID:AA0044'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['n_term'])], 'PTM-0203', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'N-acetylglycine', 'PSI-MOD:00055', 'PTM-0203', 'RESID:AA0046'}),
    ModificationRule([ModificationTarget({'I'}, SequenceLocation['n_term'])], 'PTM-0204', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'PSI-MOD:00056', 'N-acetylisoleucine', 'RESID:AA0047', 'PTM-0204'}),
    ModificationRule([ModificationTarget({'M'}, SequenceLocation['n_term'])], 'PTM-0205', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'RESID:AA0049', 'PSI-MOD:00058', 'N-acetylmethionine', 'PTM-0205'}),
    ModificationRule([ModificationTarget({'P'}, SequenceLocation['n_term'])], 'PTM-0206', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'RESID:AA0050', 'PSI-MOD:00059', 'PTM-0206', 'N-acetylproline'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['n_term'])], 'PTM-0207', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'PSI-MOD:00060', 'PTM-0207', 'N-acetylserine', 'RESID:AA0051'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['n_term'])], 'PTM-0208', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'N-acetylthreonine', 'RESID:AA0052', 'PSI-MOD:00061', 'PTM-0208'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['n_term'])], 'PTM-0209', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'PSI-MOD:00062', 'N-acetyltyrosine', 'RESID:AA0053', 'PTM-0209'}),
    ModificationRule([ModificationTarget({'V'}, SequenceLocation['n_term'])], 'PTM-0210', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'PSI-MOD:00063', 'RESID:AA0054', 'PTM-0210', 'N-acetylvaline'}),
    ModificationRule([ModificationTarget({'A'}, SequenceLocation['n_term'])], 'PTM-0374', monoisotopic_mass=43.005814, composition=Composition(
        {'C': 1, 'H': 1, 'N': 1, 'O': 1}), alt_names={'PTM-0374', 'RESID:AA0343', 'N-carbamoylalanine', 'PSI-MOD:00348'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['n_term'])], 'PTM-0331', monoisotopic_mass=176.032088, composition=Composition(
        {'C': 6, 'H': 8, 'O': 6}), alt_names={'PTM-0331', 'N-D-glucuronoyl glycine', 'PSI-MOD:00067', 'RESID:AA0058'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['n_term'])], 'PTM-0211', monoisotopic_mass=27.994915,
                     composition=Composition({'C': 1, 'O': 1}), alt_names={'PSI-MOD:00066', 'N-formylglycine', 'PTM-0211', 'RESID:AA0057'}),
    ModificationRule([ModificationTarget({'M'}, SequenceLocation['n_term'])], 'PTM-0212', monoisotopic_mass=27.994915,
                     composition=Composition({'C': 1, 'O': 1}), alt_names={'N-formylmethionine', 'PSI-MOD:00482', 'RESID:AA0021', 'PTM-0212'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0213', monoisotopic_mass=44.985078,
                     composition=Composition({'H': -1, 'N': 1, 'O': 2}), alt_names={'PSI-MOD:01352', 'Nitrated tyrosine', 'PTM-0213'}),
    ModificationRule([ModificationTarget({'A'}, SequenceLocation['n_term'])], 'PTM-0214', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'N-methylalanine', 'PSI-MOD:00070', 'RESID:AA0061', 'PTM-0214'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['anywhere'])], 'PTM-0483', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'RESID:AA0063', 'PSI-MOD:00072', 'N-methylglycine', 'PTM-0483'}),
    ModificationRule([ModificationTarget({'I'}, SequenceLocation['n_term'])], 'PTM-0215', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'PTM-0215', 'PSI-MOD:00341', 'N-methylisoleucine', 'RESID:AA0336'}),
    ModificationRule([ModificationTarget({'L'}, SequenceLocation['n_term'])], 'PTM-0216', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'PTM-0216', 'PSI-MOD:00342', 'RESID:AA0337', 'N-methylleucine'}),
    ModificationRule([ModificationTarget({'M'}, SequenceLocation['n_term'])], 'PTM-0217', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'PTM-0217', 'PSI-MOD:00073', 'N-methylmethionine', 'RESID:AA0064'}),
    ModificationRule([ModificationTarget({'F'}, SequenceLocation['n_term'])], 'PTM-0218', monoisotopic_mass=14.01565, composition=Composition(
        {'C': 1, 'H': 2}), alt_names={'RESID:AA0065', 'PSI-MOD:00074', 'PTM-0218', 'N-methylphenylalanine'}),
    ModificationRule([ModificationTarget({'P'}, SequenceLocation['n_term'])], 'PTM-0219', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'RESID:AA0419', 'N-methylproline', 'PTM-0219', 'PSI-MOD:00830'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['n_term'])], 'PTM-0432', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'RESID:AA0533', 'PSI-MOD:01782', 'N-methylserine', 'PTM-0432'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['n_term'])], 'PTM-0220', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'N-methyltyrosine', 'RESID:AA0338', 'PTM-0220', 'PSI-MOD:00343'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['n_term'])], 'PTM-0221', monoisotopic_mass=210.198365, composition=Composition(
        {'C': 14, 'H': 26, 'O': 1}), alt_names={'PSI-MOD:00068', 'PTM-0221', 'N-myristoyl glycine', 'RESID:AA0059'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['n_term'])], 'PTM-0222', monoisotopic_mass=238.229666, composition=Composition(
        {'C': 16, 'H': 30, 'O': 1}), alt_names={'RESID:AA0060', 'PSI-MOD:00069', 'N-palmitoyl cysteine', 'PTM-0222'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['n_term'])], 'PTM-0223', monoisotopic_mass=238.229666, composition=Composition(
        {'C': 16, 'H': 30, 'O': 1}), alt_names={'PSI-MOD:00344', 'PTM-0223', 'N-palmitoyl glycine', 'RESID:AA0339'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['n_term'])], 'PTM-0224', monoisotopic_mass=70.005479, composition=Composition(
        {'C': 3, 'H': 2, 'O': 2}), alt_names={'PSI-MOD:00279', 'PTM-0224', 'RESID:AA0274', 'N-pyruvate 2-iminyl-cysteine'}),
    ModificationRule([ModificationTarget({'V'}, SequenceLocation['n_term'])], 'PTM-0225', monoisotopic_mass=70.005479, composition=Composition(
        {'C': 3, 'H': 2, 'O': 2}), alt_names={'RESID:AA0275', 'PTM-0225', 'PSI-MOD:00280', 'N-pyruvate 2-iminyl-valine'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0399', monoisotopic_mass=123.00853, composition=Composition(
        {'C': 2, 'H': 6, 'N': 1, 'O': 3, 'P': 1}), alt_names={'PTM-0399', 'RESID:AA0497', 'PSI-MOD:01587', 'O-(2-aminoethylphosphoryl)serine'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0400', monoisotopic_mass=166.062756, composition=Composition(
        {'C': 5, 'H': 13, 'N': 1, 'O': 3, 'P': 1}), alt_names={'PSI-MOD:01588', 'O-(2-cholinephosphoryl)serine', 'RESID:AA0498', 'PTM-0400'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0391', monoisotopic_mass=340.085794, composition=Composition(
        {'C': 11, 'H': 21, 'N': 2, 'O': 6, 'P': 1, 'S': 1}), alt_names={'PSI-MOD:00159', "O-(pantetheine 4'-phosphoryl)serine", 'PTM-0391', 'RESID:AA0150'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0389', monoisotopic_mass=881.146903, composition=Composition(
        {'C': 26, 'H': 42, 'N': 7, 'O': 19, 'P': 3, 'S': 1}), alt_names={'PTM-0389', 'PSI-MOD:00176', 'RESID:AA0167', 'O-(phosphoribosyl dephospho-coenzyme A)serine'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0230', monoisotopic_mass=154.00311, composition=Composition(
        {'C': 3, 'H': 7, 'O': 5, 'P': 1}), alt_names={'O-(sn-1-glycerophosphoryl)serine', 'RESID:AA0264', 'PSI-MOD:00269', 'PTM-0230'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0231', monoisotopic_mass=783.141485, composition=Composition(
        {'C': 27, 'H': 31, 'N': 9, 'O': 15, 'P': 2}), alt_names={'RESID:AA0145', 'O-8alpha-FAD tyrosine', 'PSI-MOD:00154', 'PTM-0231'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0232', monoisotopic_mass=42.010565,
                     composition=Composition({'C': 2, 'H': 2, 'O': 1}), alt_names={'PSI-MOD:00369', 'PTM-0232', 'O-acetylserine', 'RESID:AA0364'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0233', monoisotopic_mass=42.010565, composition=Composition(
        {'C': 2, 'H': 2, 'O': 1}), alt_names={'PTM-0233', 'PSI-MOD:01171', 'O-acetylthreonine', 'RESID:AA0423'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0332', monoisotopic_mass=329.05252, composition=Composition(
        {'C': 10, 'H': 12, 'N': 5, 'O': 6, 'P': 1}), alt_names={'RESID:AA0203', 'PTM-0332', 'O-AMP-tyrosine', 'PSI-MOD:00208'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0393', monoisotopic_mass=329.05252, composition=Composition(
        {'C': 10, 'H': 12, 'N': 5, 'O': 6, 'P': 1}), alt_names={'PSI-MOD:00272', 'RESID:AA0267', 'PTM-0393', 'O-AMP-threonine'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0235', monoisotopic_mass=154.135765, composition=Composition(
        {'C': 10, 'H': 18, 'O': 1}), alt_names={'PTM-0235', 'RESID:AA0387', 'O-decanoyl threonine', 'PSI-MOD:00392'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0234', monoisotopic_mass=154.135765, composition=Composition(
        {'C': 10, 'H': 18, 'O': 1}), alt_names={'PSI-MOD:00390', 'O-decanoyl serine', 'RESID:AA0385', 'PTM-0234'}),
    ModificationRule([ModificationTarget({'Q'}, SequenceLocation['anywhere'])], 'PTM-0236', monoisotopic_mass=760.730862, composition=Composition(
        {'C': 50, 'H': 96, 'O': 4}), alt_names={'PSI-MOD:00352', 'RESID:AA0347', 'PTM-0236', 'Omega-hydroxyceramide glutamate ester'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0237', monoisotopic_mass=14.01565, composition=Composition(
        {'C': 1, 'H': 2}), alt_names={'PSI-MOD:00078', 'RESID:AA0069', 'Omega-N-methylarginine', 'PTM-0237'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0356', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'RESID:AA0464', 'O-methylthreonine', 'PSI-MOD:01387', 'PTM-0356'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0239', monoisotopic_mass=126.104465, composition=Composition(
        {'C': 8, 'H': 14, 'O': 1}), alt_names={'O-octanoyl serine', 'RESID:AA0290', 'PSI-MOD:00295', 'PTM-0239'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0240', monoisotopic_mass=126.104465, composition=Composition(
        {'C': 8, 'H': 14, 'O': 1}), alt_names={'O-octanoyl threonine', 'RESID:AA0386', 'PSI-MOD:00391', 'PTM-0240'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0241', monoisotopic_mass=238.229666, composition=Composition(
        {'C': 16, 'H': 30, 'O': 1}), alt_names={'RESID:AA0080', 'PSI-MOD:00089', 'O-palmitoyl serine', 'PTM-0241'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0242', monoisotopic_mass=238.229666, composition=Composition(
        {'C': 16, 'H': 30, 'O': 1}), alt_names={'PSI-MOD:00088', 'RESID:AA0079', 'PTM-0242', 'O-palmitoyl threonine'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0333', monoisotopic_mass=306.025302, composition=Composition(
        {'C': 9, 'H': 11, 'N': 2, 'O': 8, 'P': 1}), alt_names={'PSI-MOD:00261', 'RESID:AA0256', 'O-UMP-tyrosine', 'PTM-0333'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0376', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'RESID:AA0238', 'PTM-0376', 'Oxazole-4-carboxylic acid (Cys-Ser)', 'PSI-MOD:00243'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0377', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'PTM-0377', 'Oxazole-4-carboxylic acid (Gly-Ser)', 'PSI-MOD:00245', 'RESID:AA0240'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0463', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'Oxazole-4-carboxylic acid (Ile-Ser)', 'RESID:AA0573', 'PSI-MOD:01902', 'PTM-0463'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0464', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'PTM-0464', 'RESID:AA0574', 'PSI-MOD:01903', 'Oxazole-4-carboxylic acid (Ser-Ser)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0381', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'Oxazoline-4-carboxylic acid (Cys-Ser)', 'PTM-0381', 'RESID:AA0239', 'PSI-MOD:00244'}),
    ModificationRule([ModificationTarget({'F'}, SequenceLocation['c_term'])], 'PTM-0248', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'PTM-0248', 'RESID:AA0094', 'Phenylalanine amide', 'PSI-MOD:00103'}),
    ModificationRule([ModificationTarget({'G'}, SequenceLocation['c_term'])], 'PTM-0249', monoisotopic_mass=699.52029, composition=Composition(
        {'C': 39, 'H': 74, 'N': 1, 'O': 7, 'P': 1}), alt_names={'RESID:AA0346', 'Phosphatidylethanolamine amidated glycine', 'PSI-MOD:00351', 'PTM-0249'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0250', monoisotopic_mass=79.966331, composition=Composition(
        {'H': 1, 'O': 3, 'P': 1}), alt_names={'Phosphoarginine', 'RESID:AA0222', 'PSI-MOD:00227', 'PTM-0250'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0251', monoisotopic_mass=79.966331, composition=Composition(
        {'H': 1, 'O': 3, 'P': 1}), alt_names={'PTM-0251', 'RESID:AA0034', 'PSI-MOD:00043', 'Phosphocysteine'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0252', monoisotopic_mass=79.966331,
                     composition=Composition({'H': 1, 'O': 3, 'P': 1}), alt_names={'PSI-MOD:00890', 'PTM-0252', 'Phosphohistidine'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0253', monoisotopic_mass=79.966331,
                     composition=Composition({'H': 1, 'O': 3, 'P': 1}), alt_names={'Phosphoserine', 'PTM-0253', 'PSI-MOD:00046', 'RESID:AA0037'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0254', monoisotopic_mass=79.966331, composition=Composition(
        {'H': 1, 'O': 3, 'P': 1}), alt_names={'RESID:AA0038', 'Phosphothreonine', 'PTM-0254', 'PSI-MOD:00047'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0255', monoisotopic_mass=79.966331, composition=Composition(
        {'H': 1, 'O': 3, 'P': 1}), alt_names={'PSI-MOD:00048', 'PTM-0255', 'RESID:AA0039', 'Phosphotyrosine'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0449', monoisotopic_mass=-4.0313, composition=Composition(
        {'H': -4}), alt_names={'PSI-MOD:01845', 'RESID:AA0551', 'PTM-0449', 'Proline 5-hydroxy-oxazole-4-carbothionic acid (Pro-Cys)'}),
    ModificationRule([ModificationTarget({'P'}, SequenceLocation['c_term'])], 'PTM-0257', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'PTM-0257', 'PSI-MOD:00104', 'Proline amide', 'RESID:AA0095'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0258', monoisotopic_mass=783.141485, composition=Composition(
        {'C': 27, 'H': 31, 'N': 9, 'O': 15, 'P': 2}), alt_names={'PSI-MOD:00153', 'PTM-0258', 'Pros-8alpha-FAD histidine', 'RESID:AA0144'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0259', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'Pros-methylhistidine', 'RESID:AA0073', 'PSI-MOD:00082', 'PTM-0259'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0260', monoisotopic_mass=79.966331, composition=Composition(
        {'H': 1, 'O': 3, 'P': 1}), alt_names={'Pros-phosphohistidine', 'PTM-0260', 'PSI-MOD:00045', 'RESID:AA0036'}),
    ModificationRule([ModificationTarget({'Q'}, SequenceLocation['n_term'])], 'PTM-0261', monoisotopic_mass=-17.026549, composition=Composition(
        {'H': -3, 'N': -1}), alt_names={'PTM-0261', 'Pyrrolidone carboxylic acid', 'PSI-MOD:00040', 'RESID:AA0031'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['n_term'])], 'PTM-0262', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'RESID:AA0031', 'PTM-0262', 'PSI-MOD:00420', 'Pyrrolidone carboxylic acid (Glu)'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0263', monoisotopic_mass=37.906494, composition=Composition(
        {'H': -10, 'O': 3}), alt_names={'RESID:AA0283', 'PSI-MOD:00288', 'Pyrroloquinoline quinone (Glu-Tyr)', 'PTM-0263'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['n_term'])], 'PTM-0265', monoisotopic_mass=-33.003705, composition=Composition(
        {'H': -3, 'N': -1, 'O': 1, 'S': -1}), alt_names={'Pyruvic acid (Cys)', 'PTM-0265', 'PSI-MOD:00136', 'RESID:AA0127'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['n_term'])], 'PTM-0266', monoisotopic_mass=-17.026549,
                     composition=Composition({'H': -3, 'N': -1}), alt_names={'Pyruvic acid (Ser)', 'PTM-0266', 'PSI-MOD:00807', 'RESID:AA0127'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0447', monoisotopic_mass=316.203845, composition=Composition(
        {'C': 20, 'H': 28, 'O': 3}), alt_names={'PTM-0447', 'RESID:AA0426', 'PSI-MOD:01174', 'S-(15-deoxy-Delta12,14-prostaglandin J2-9-yl)cysteine'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0267', monoisotopic_mass=-64.016044, composition=Composition(
        {'C': -1, 'H': -4, 'O': -3}), alt_names={'S-(2-aminovinyl)-3-methyl-D-cysteine (Thr-Cys)', 'PTM-0267', 'PSI-MOD:00258', 'RESID:AA0253'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0446', monoisotopic_mass=-79.9932, composition=Composition(
        {'C': -1, 'H': -4, 'O': -2, 'S': -1}), alt_names={'S-(2-aminovinyl)-D-cysteine (Cys-Cys)', 'RESID:AA0204', 'PTM-0446', 'PSI-MOD:01849'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0268', monoisotopic_mass=-64.016044, composition=Composition(
        {'C': -1, 'H': -4, 'O': -3}), alt_names={'S-(2-aminovinyl)-D-cysteine (Ser-Cys)', 'RESID:AA0204', 'PTM-0268', 'PSI-MOD:00209'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0443', monoisotopic_mass=-79.9932, composition=Composition(
        {'C': -1, 'H': -4, 'O': -2, 'S': -1}), alt_names={'RESID:AA0548', 'PTM-0443', 'S-(2-aminovinyl)-L-cysteine (Cys-Cys)', 'PSI-MOD:01842'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0414', monoisotopic_mass=146.036779, composition=Composition(
        {'C': 9, 'H': 6, 'O': 2}), alt_names={'RESID:AA0207', 'PSI-MOD:00212', 'PTM-0414', 'S-(4-hydroxycinnamyl)cysteine'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0428', monoisotopic_mass=421.142641, composition=Composition(
        {'C': 26, 'H': 19, 'N': 3, 'O': 3}), alt_names={'PTM-0428', 'S-(coelenterazin-3a-yl)cysteine', 'PSI-MOD:01694'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0421', monoisotopic_mass=418.137616, composition=Composition(
        {'C': 20, 'H': 22, 'N': 2, 'O': 8}), alt_names={'PTM-0421', 'S-(dipyrrolylmethanemethyl)cysteine', 'RESID:AA0252', 'PSI-MOD:00257'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0269', monoisotopic_mass=220.182715, composition=Composition(
        {'C': 15, 'H': 24, 'O': 1}), alt_names={'PSI-MOD:00112', 'S-12-hydroxyfarnesyl cysteine', 'RESID:AA0103', 'PTM-0269'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0270', monoisotopic_mass=456.104615, composition=Composition(
        {'C': 17, 'H': 21, 'N': 4, 'O': 9, 'P': 1}), alt_names={'PTM-0270', 'RESID:AA0351', 'PSI-MOD:00356', 'S-4a-FMN cysteine'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0271', monoisotopic_mass=454.088965, composition=Composition(
        {'C': 17, 'H': 19, 'N': 4, 'O': 9, 'P': 1}), alt_names={'PTM-0271', 'RESID:AA0220', 'S-6-FMN cysteine', 'PSI-MOD:00225'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0272', monoisotopic_mass=783.141485, composition=Composition(
        {'C': 27, 'H': 31, 'N': 9, 'O': 15, 'P': 2}), alt_names={'RESID:AA0143', 'PSI-MOD:00152', 'S-8alpha-FAD cysteine', 'PTM-0272'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['n_term'])], 'PTM-0273', monoisotopic_mass=634.662782, composition=Composition(
        {'C': 43, 'H': 86, 'O': 2}), alt_names={'PTM-0273', 'RESID:AA0223', 'S-archaeol cysteine', 'PSI-MOD:00228'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0452', monoisotopic_mass=396.083866, composition=Composition(
        {'C': 13, 'H': 20, 'N': 2, 'O': 10, 'S': 1}), alt_names={'PTM-0452', 'RESID:AA0563', 'PSI-MOD:01860', 'S-bacillithiol cysteine disulfide'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0324', monoisotopic_mass=13.979265, composition=Composition(
        {'H': -2, 'O': 1}), alt_names={'S-cysteinyl 3-(oxidosulfanyl)alanine (Cys-Cys)', 'PSI-MOD:01383', 'RESID:AA0457', 'PTM-0324'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0415', monoisotopic_mass=119.004099, composition=Composition(
        {'C': 3, 'H': 5, 'N': 1, 'O': 2, 'S': 1}), alt_names={'S-cysteinyl cysteine', 'PSI-MOD:00765', 'PTM-0415', 'RESID:AA0025'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['n_term'])], 'PTM-0274', monoisotopic_mass=576.511761, composition=Composition(
        {'C': 37, 'H': 68, 'O': 4}), alt_names={'PTM-0274', 'S-diacylglycerol cysteine', 'PSI-MOD:00116', 'RESID:AA0107'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['c_term'])], 'PTM-0275', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'RESID:AA0096', 'PSI-MOD:00105', 'Serine amide', 'PTM-0275'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['c_term'])], 'PTM-0276', monoisotopic_mass=831.197041, composition=Composition(
        {'C': 36, 'H': 37, 'N': 3, 'O': 20}), alt_names={'RESID:AA0374', 'Serine microcin E492 siderophore ester', 'PSI-MOD:00379', 'PTM-0276'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0277', monoisotopic_mass=204.187801,
                     composition=Composition({'C': 15, 'H': 24}), alt_names={'S-farnesyl cysteine', 'PTM-0277', 'PSI-MOD:00111', 'RESID:AA0102'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0278', monoisotopic_mass=272.250401, composition=Composition(
        {'C': 20, 'H': 32}), alt_names={'PSI-MOD:00113', 'RESID:AA0104', 'PTM-0278', 'S-geranylgeranyl cysteine'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0311', monoisotopic_mass=305.068156, composition=Composition(
        {'C': 10, 'H': 15, 'N': 3, 'O': 6, 'S': 1}), alt_names={'PTM-0311', 'S-glutathionyl cysteine', 'PSI-MOD:00234', 'RESID:AA0229'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0279', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'RESID:AA0234', 'S-methylcysteine', 'PTM-0279', 'PSI-MOD:00239'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0280', monoisotopic_mass=28.990164, composition=Composition(
        {'H': -1, 'N': 1, 'O': 1}), alt_names={'S-nitrosocysteine', 'PTM-0280', 'PSI-MOD:00235', 'RESID:AA0230'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0281', monoisotopic_mass=238.229666, composition=Composition(
        {'C': 16, 'H': 30, 'O': 1}), alt_names={'RESID:AA0106', 'PSI-MOD:00115', 'PTM-0281', 'S-palmitoyl cysteine'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0282', monoisotopic_mass=79.916521,
                     composition=Composition({'Se': 1}), alt_names={'PSI-MOD:00282', 'PTM-0282', 'RESID:AA0277', 'S-selanylcysteine'}),
    ModificationRule([ModificationTarget({'C'}, SequenceLocation['anywhere'])], 'PTM-0283', monoisotopic_mass=266.260966, composition=Composition(
        {'C': 18, 'H': 34, 'O': 1}), alt_names={'S-stearoyl cysteine', 'PSI-MOD:00816', 'RESID:AA0407', 'PTM-0283'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0284', monoisotopic_mass=79.956815,
                     composition=Composition({'O': 3, 'S': 1}), alt_names={'PSI-MOD:00366', 'Sulfoserine', 'PTM-0284', 'RESID:AA0361'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['anywhere'])], 'PTM-0285', monoisotopic_mass=79.956815,
                     composition=Composition({'O': 3, 'S': 1}), alt_names={'PTM-0285', 'RESID:AA0362', 'Sulfothreonine', 'PSI-MOD:00367'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0286', monoisotopic_mass=79.956815,
                     composition=Composition({'O': 3, 'S': 1}), alt_names={'Sulfotyrosine', 'RESID:AA0172', 'PTM-0286', 'PSI-MOD:00181'}),
    ModificationRule([ModificationTarget({'R'}, SequenceLocation['anywhere'])], 'PTM-0287', monoisotopic_mass=28.0313, composition=Composition(
        {'C': 2, 'H': 4}), alt_names={'RESID:AA0067', 'Symmetric dimethylarginine', 'PTM-0287', 'PSI-MOD:00076'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0416', monoisotopic_mass=90.031694, composition=Composition(
        {'C': 3, 'H': 6, 'O': 3}), alt_names={'PSI-MOD:01177', 'Tele-(1,2,3-trihydroxypropan-2-yl)histidine', 'PTM-0416', 'RESID:AA0431'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0288', monoisotopic_mass=783.141485, composition=Composition(
        {'C': 27, 'H': 31, 'N': 9, 'O': 15, 'P': 2}), alt_names={'PTM-0288', 'Tele-8alpha-FAD histidine', 'PSI-MOD:00226', 'RESID:AA0221'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0289', monoisotopic_mass=454.088965, composition=Composition(
        {'C': 17, 'H': 19, 'N': 4, 'O': 9, 'P': 1}), alt_names={'PSI-MOD:00357', 'RESID:AA0352', 'Tele-8alpha-FMN histidine', 'PTM-0289'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0290', monoisotopic_mass=14.01565,
                     composition=Composition({'C': 1, 'H': 2}), alt_names={'Tele-methylhistidine', 'RESID:AA0317', 'PTM-0290', 'PSI-MOD:00322'}),
    ModificationRule([ModificationTarget({'H'}, SequenceLocation['anywhere'])], 'PTM-0325', monoisotopic_mass=79.966331, composition=Composition(
        {'H': 1, 'O': 3, 'P': 1}), alt_names={'Tele-phosphohistidine', 'RESID:AA0035', 'PSI-MOD:00044', 'PTM-0325'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0460', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'PSI-MOD:01899', 'Thiazole-4-carboxylic acid (Arg-Cys)', 'RESID:AA0570', 'PTM-0460'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0456', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'RESID:AA0541', 'PTM-0456', 'Thiazole-4-carboxylic acid (Glu-Cys)', 'PSI-MOD:01815'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0378', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'Thiazole-4-carboxylic acid (Gly-Cys)', 'RESID:AA0241', 'PTM-0378', 'PSI-MOD:00246'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0361', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'Thiazole-4-carboxylic acid (Ile-Cys)', 'PTM-0361', 'PSI-MOD:01389', 'RESID:AA0466'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0363', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'Thiazole-4-carboxylic acid (Ser-Cys)', 'PSI-MOD:00247', 'RESID:AA0242', 'PTM-0363'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0364', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'Thiazole-4-carboxylic acid (Thr-Cys)', 'PSI-MOD:01406', 'PTM-0364', 'RESID:AA0483'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0365', monoisotopic_mass=-20.026215, composition=Composition(
        {'H': -4, 'O': -1}), alt_names={'Thiazole-4-carboxylic acid (Val-Cys)', 'PTM-0365', 'RESID:AA0467', 'PSI-MOD:01390'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0366', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'Thiazoline-4-carboxylic acid (Phe-Cys)', 'PSI-MOD:01407', 'RESID:AA0484', 'PTM-0366'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0392', monoisotopic_mass=-18.010565, composition=Composition(
        {'H': -2, 'O': -1}), alt_names={'RESID:AA0485', 'PSI-MOD:01408', '(4S)-thiazoline-4-carboxylic acid (Thr-Cys)', 'PTM-0392'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0458', monoisotopic_mass=-4.0313, composition=Composition(
        {'H': -4}), alt_names={'RESID:AA0554', 'Threonine 5-hydroxy-oxazole-4-carbonthionic acid (Thr-Cys)', 'PTM-0458', 'PSI-MOD:01878'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['c_term'])], 'PTM-0293', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'Threonine amide', 'PSI-MOD:00106', 'RESID:AA0097', 'PTM-0293'}),
    ModificationRule([ModificationTarget({'T'}, SequenceLocation['c_term'])], 'PTM-0412', monoisotopic_mass=14.01565, composition=Composition(
        {'C': 1, 'H': 2}), alt_names={'PSI-MOD:01610', 'Threonine methyl ester', 'PTM-0412', 'RESID:AA0507'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0294', monoisotopic_mass=595.612805,
                     composition=Composition({'C': 6, 'I': 4, 'O': 1}), alt_names={'PSI-MOD:00187', 'RESID:AA0178', 'PTM-0294', 'Thyroxine'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['anywhere'])], 'PTM-0295', monoisotopic_mass=469.716158, composition=Composition(
        {'C': 6, 'H': 1, 'I': 3, 'O': 1}), alt_names={'PSI-MOD:00186', 'RESID:AA0177', 'PTM-0295', 'Triiodothyronine'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0417', monoisotopic_mass=93.900563,
                     composition=Composition({'H': -2, 'S': 3}), alt_names={'Trithiocysteine (Cys-Cys)', 'PTM-0417', 'RESID:AA0513', 'PSI-MOD:01616'}),
    ModificationRule([ModificationTarget({'W'}, SequenceLocation['c_term'])], 'PTM-0296', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'Tryptophan amide', 'RESID:AA0098', 'PTM-0296', 'PSI-MOD:00107'}),
    ModificationRule([ModificationTarget(set(), SequenceLocation['anywhere'])], 'PTM-0298', monoisotopic_mass=27.958529, composition=Composition(
        {'H': -4, 'O': 2}), alt_names={'Tryptophan tryptophylquinone (Trp-Trp)', 'RESID:AA0149', 'PSI-MOD:00158', 'PTM-0298'}),
    ModificationRule([ModificationTarget({'W'}, SequenceLocation['anywhere'])], 'PTM-0299', monoisotopic_mass=29.974179,
                     composition=Composition({'H': -2, 'O': 2}), alt_names={'PTM-0299', 'RESID:AA0148', 'PSI-MOD:00157', 'Tryptophylquinone'}),
    ModificationRule([ModificationTarget({'Y'}, SequenceLocation['c_term'])], 'PTM-0302', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'PTM-0302', 'PSI-MOD:00108', 'Tyrosine amide', 'RESID:AA0099'}),
    ModificationRule([ModificationTarget({'V'}, SequenceLocation['c_term'])], 'PTM-0303', monoisotopic_mass=-0.984016, composition=Composition(
        {'H': 1, 'N': 1, 'O': -1}), alt_names={'PTM-0303', 'PSI-MOD:00109', 'RESID:AA0100', 'Valine amide'}),
    ModificationRule([ModificationTarget({'Q'}, SequenceLocation['anywhere'])], 'PTM-0488', monoisotopic_mass=94.053098, composition=Composition(
        {'C': 5, 'H': 6, 'N': 2}), alt_names={'RESID:AA0596', 'PTM-0488', 'PSI-MOD:01950', 'L-isoglutamyl histamine'}),
    ModificationRule([ModificationTarget({'S'}, SequenceLocation['anywhere'])], 'PTM-0493', monoisotopic_mass=236.214016,
                     composition=Composition({}), alt_names={'PTM-0493', 'RESID:AA0455', 'PSI-MOD:01381', 'O-palmitoleyl serine'}),
]
# [[[end]]]
