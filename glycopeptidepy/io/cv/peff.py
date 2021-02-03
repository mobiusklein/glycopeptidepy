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
    PEFFTerm(u'PEFF file sequence entry term', u'PEFF:0000003',
             u'CV term that may appear in a description line of a PEFF file individual sequence entry.',       'PEFF CV Term', [u'PEFF CV term']),
    PEFFTerm(u'PEFF file header section term', u'PEFF:0000002',
             u'CV term that may appear in a PEFF file header section.',       'PEFF CV Term', [u'PEFF CV term']),
    PEFFTerm(u'NcbiTaxId', u'PEFF:0001003', u'PEFF keyword for the NCBI taxonomy identifier.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'TaxName', u'PEFF:0001004', u'PEFF keyword for the taxonomy name (latin or common name).',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'Length', u'PEFF:0001006', u'PEFF keyword for the sequence length.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'SV', u'PEFF:0001007', u'PEFF keyword for the sequence version.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'EV', u'PEFF:0001008', u'PEFF keyword for the entry version.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'PE', u'PEFF:0001009', u'PEFF keyword for the Protein Evidence; A UniProtKB code 1-5.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'DbUniqueId', u'PEFF:0001001', u'OBSOLETE Sequence database unique identifier.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'PName', u'PEFF:0001002', u'PEFF keyword for the protein full name.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'XRef', u'PEFF:0001019', u'PEFF keyword for the cross-reference to an external resource.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'GO', u'PEFF:0001018', u'PEFF keyword for the Gene Ontology code.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'ModRes', u'PEFF:0001013', u'PEFF keyword for the modified residue without aPSI-MOD or UniMod identifier.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'ModResPsi', u'PEFF:0001012', u'PEFF keyword for the modified residue with PSI-MOD identifier.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'Variant', u'PEFF:0001011', u'OBSOLETE Sequence variation (substitution, insertion, deletion).',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'Processed', u'PEFF:0001010', u'PEFF keyword for information on how the full length original protein sequence can be processed into shorter components such as signal peptides and chains.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'KW', u'PEFF:0001017', u'PEFF keyword for the entry associated keyword(s).',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'CC', u'PEFF:0001016', u'PEFF keyword for the entry associated comment.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'SeqStatus', u'PEFF:0001015', u'PEFF keyword for the sequence status. Complete or Fragment.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'AltAC', u'PEFF:0001014', u'PEFF keyword for the Alternative Accession Code.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'ID', u'PEFF:0001026', u'PEFF keyword for the UniProtKB specific Protein identifier ID; a UniProtKB term.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'ModResUnimod', u'PEFF:0001027', u'PEFF keyword for the modified residue with UniMod identifier.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'Crc64', u'PEFF:0001024', u'PEFF keyword for the Sequence checksum in crc64.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'Domain', u'PEFF:0001025', u'PEFF keyword for the sequence range of a domain.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'Conflict', u'PEFF:0001023', u'PEFF keyword for the sequence conflict; a UniProtKB term.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'VariantSimple', u'PEFF:0001028', u'PEFF keyword for the simple sequence variation of a single amino acid change. A change to a stop codon is permitted with a * symbol. More complex variations must be encoded with the VariantComplex term.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'VariantComplex', u'PEFF:0001029', u'PEFF keyword for a sequence variation that is more complex than a single amino acid change or change to a stop codon.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'GName', u'PEFF:0001005', u'PEFF keyword for the gene name.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'DisulfideBond', u'PEFF:0001031', u'PEFF keyword for the disulfide bonds in this protein, constructed as a sets of annotation identifiers of two half-cystine modifications.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'Proteoform', u'PEFF:0001030', u'PEFF keyword for the proteoforms of this protein, constructed as a set of annotation identifiers.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'Comment', u'PEFF:0001033', u'PEFF keyword for the individual protein entry comment. It is discouraged to put parsable information here. This is only for free-text commentary.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'PEFF molecule processing keyword', u'PEFF:0001032', u'PEFF keyword describing the type of processing event being described.',
             'PEFF CV Term', [u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'ProteoformDb', u'PEFF:0000022', u"PEFF keyword that when set to 'true' indicates that the database contains complete proteoforms.",
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'SequenceType', u'PEFF:0000017', u'PEFF keyword for the molecular type of the sequences.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'GeneralComment', u'PEFF:0000021', u'PEFF keyword for a general comment.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'DatabaseDescription', u'PEFF:0000020', u'PEFF keyword for the short description of the PEFF file.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'OptionalTagDef', u'PEFF:0000023', u'PEFF keyword for the short tag (abbreviation) and longer definition used to annotate a sequence annotation (such as variant or modification) in the OptionalTag location.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'HasAnnotationIdentifiers', u'PEFF:0000024', u"PEFF keyword that when set to 'true' indicates that entries in the database have identifiers for each annotation.",
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'Prefix', u'PEFF:0000009', u'PEFF keyword for the sequence database prefix.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'DbName', u'PEFF:0000008', u'PEFF keyword for the sequence database name.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'NumberOfEntries', u'PEFF:0000015', u'PEFF keyword for the sumber of sequence entries in the database.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'DbDescription', u'PEFF:0000010', u'PEFF keyword for the sequence database short description.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'Decoy', u'PEFF:0000011', u'PEFF keyword for the specifying whether the sequence database is a decoy database.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'DbVersion', u'PEFF:0000013', u'PEFF keyword for the database version (release date) according to database provider.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'DbDate', u'PEFF:0000014', u'OBSOLETE PEFF keyword for the database date (release or file date of the source) according to database provider.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'SpecificValue', u'PEFF:0000019', u'PEFF keyword for the specific values for a custom key.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'SpecificKey', u'PEFF:0000018', u'PEFF keyword for database specific keywords not included in the current controlled vocabulary.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'DbSource', u'PEFF:0000012', u'PEFF keyword for the source of the database file.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'Conversion', u'PEFF:0000016', u'PEFF keyword for the description of the conversion from original format to this current one.',
             'PEFF CV Term', [u'PEFF file header section term', u'PEFF CV term']),
    PEFFTerm(u'transit peptide', u'PEFF:0001022', u'Short peptide present at the N-terminus of a newly synthesized protein that helps the protein through the membrane of its destination organelle.',
             'PEFF CV Term', [u'PEFF molecule processing keyword', u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'mature protein', u'PEFF:0001020', u'Portion of a newly synthesized protein that contributes to a final structure after other components such as signal peptides are removed.',
             'PEFF CV Term', [u'PEFF molecule processing keyword', u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'signal peptide', u'PEFF:0001021', u'Short peptide present at the N-terminus of a newly synthesized protein that is cleaved off and is not part of the final mature protein.',
             'PEFF CV Term', [u'PEFF molecule processing keyword', u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'initiator methionine', u'PEFF:0001035', u'N-terminal methionine residue of a protein that can be co-translationally cleaved.',
             'PEFF CV Term', [u'PEFF molecule processing keyword', u'PEFF file sequence entry term', u'PEFF CV term']),
    PEFFTerm(u'propeptide', u'PEFF:0001034', u'Short peptide that is cleaved off a newly synthesized protein and generally immediately degraded in the process of protein maturation, and is not a signal peptide or transit peptide.',
             'PEFF CV Term', [u'PEFF molecule processing keyword', u'PEFF file sequence entry term', u'PEFF CV term']),
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
