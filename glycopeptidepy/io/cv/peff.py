from .term import Term


class PEFFTerm(Term):
    __slots__ = ()


peff_cv_terms = []

# [[[cog
# import cog
# from glycopeptidepy.io.cv.term import render_list
# render_list('PEFF CV Term', list_name="peff_cv_terms", term_cls_name="PEFFTerm", writer=cog.out)
# ]]]
peff_cv_terms = [
    PEFFTerm(u'Individual Sequence Entries Section term', u'PEFF:0000003',
             u'CV term that may appear in a PEFF Individual Entry Section.',       'PEFF CV Term', [u'PEFF CV term']),
    PEFFTerm(u'File Header Section term', u'PEFF:0000002',
             u'CV term that may appear in a PEFF File Header Section.',       'PEFF CV Term', [u'PEFF CV term']),
    PEFFTerm(u'EV', u'PEFF:0001008', u'Entry version.',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'PE', u'PEFF:0001009', u'Protein Evidence; A UniprotKB code.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'DbUniqueId', u'PEFF:0001001', u'Sequence Database unique identifier.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'PName', u'PEFF:0001002', u'Protein Name, description.',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'NcbiTaxId', u'PEFF:0001003', u'NCBI taxonomy identifier.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'TaxName', u'PEFF:0001004', u'Taxonomy name (latin or common name).',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'GName', u'PEFF:0001005', u'Gene name.',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'Length', u'PEFF:0001006', u'Sequence length.',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'SV', u'PEFF:0001007', u'Sequence version.',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'XRef', u'PEFF:0001019', u'Cross-reference to an external resource.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'GO', u'PEFF:0001018', u'Gene Ontology code.',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'ModRes', u'PEFF:0001013', u'Modified residue without PSI-MOD identifier.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'ModResPsi', u'PEFF:0001012', u'Modified residue with PSI-MOD identifier.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'Variant', u'PEFF:0001011', u'DEPRECATED in favor of VariantSimple and VariantComplex. Former definition: Sequence variation (substitution, insertion, deletion).',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'Processed', u'PEFF:0001010', u'Processed Molecule.',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'KW', u'PEFF:0001017', u'Entry associated keyword(s).',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'CC', u'PEFF:0001016', u'Entry associated comment.',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'SeqStatus', u'PEFF:0001015', u'Sequence Status. Complete or Fragment.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'AltAC', u'PEFF:0001014', u'Alternative Accession Code.',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'ID', u'PEFF:0001026', u'UniProtKB specific Protein identifier ID; a UniProtKB term.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'ModResUnimod', u'PEFF:0001027', u'Modified residue with Unimod identifier.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'Crc64', u'PEFF:0001024', u'Sequence checksum in crc64.',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'Domain', u'PEFF:0001025', u'Sequence range of a domain.',       'PEFF CV Term', [
             u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'Transit', u'PEFF:0001022', u'Sequence range of transit peptide.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'Conflict', u'PEFF:0001023', u'Sequence conflict; a UniProtKB term.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'Chain', u'PEFF:0001020', u'Sequence range of active processed polypeptide.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'Signal', u'PEFF:0001021', u'Sequence range of signal peptide.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'VariantSimple', u'PEFF:0001028', u'Simple sequence variation of a single amino acid change. A change to a stop codon is permitted with a * symbol. More complex variations must be encoded with the VariantComplex term.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'VariantComplex', u'PEFF:0001029', u'Simple sequence variation of a single amino acid change. A change to a stop codon is permitted with a * symbol. More complex variations must be encoded with the VariantComplex term.',
             'PEFF CV Term', [u'Individual Sequence Entries Section term', u'PEFF CV term']),
    PEFFTerm(u'SequenceType', u'PEFF:0000017', u'Molecular type of the sequences.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'GeneralComment', u'PEFF:0000021', u'PEFF file general comment.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'DatabaseDescription', u'PEFF:0000020', u'Short Description of the PEFF.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'CustomTag', u'PEFF:0000023', u'A tag (short string) used to categorize a sequence annotation (variant or modification).',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'ProteoformDb', u'PEFF:0000022', u'Proteoform database flag.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'Prefix', u'PEFF:0000009', u'Sequence Database Prefix.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'DbName', u'PEFF:0000008', u'Sequence Database Name.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'NumberOfEntries', u'PEFF:0000015', u'Number of sequence entries in the database.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'DbDescription', u'PEFF:0000010', u'Sequence Database Short description.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'Decoy', u'PEFF:0000011', u'Specifies whether the Sequence Database is a Decoy.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'DbVersion', u'PEFF:0000013', u'Database version (release date) according to database provider.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'DbDate', u'PEFF:0000014', u'Database date (release or file date of the source) according to database provider.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'SpecificValue', u'PEFF:0000019', u'PEFF specific values for a defined key.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'SpecificKey', u'PEFF:0000018', u'Db specific information (not included in the current list of allowed keys).',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'DbSource', u'PEFF:0000012', u'Source of the database file.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
    PEFFTerm(u'Conversion', u'PEFF:0000016', u'Description of the conversion from original format to this current one.',
             'PEFF CV Term', [u'File Header Section term', u'PEFF CV term']),
]
# [[[end]]]

peff_cv_terms_by_name = {c.name: c for c in peff_cv_terms}


def peff_cv_term(name, strict=False):
    try:
        return peff_cv_terms_by_name[name]
    except KeyError:
        if not strict:
            return PEFFTerm(name, name, name, name, [name])
        else:
            raise
