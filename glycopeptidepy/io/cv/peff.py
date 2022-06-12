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
    PEFFTerm(
        "PEFF file header section term",
        "PEFF:0000002",
        "CV term that may appear in a PEFF file header section.",
        "PEFF CV Term",
        ["PEFF CV term"],
    ),
    PEFFTerm(
        "PEFF file sequence entry term",
        "PEFF:0000003",
        "CV term that may appear in a description line of a PEFF file individual sequence entry.",
        "PEFF CV Term",
        ["PEFF CV term"],
    ),
    PEFFTerm(
        "DbName",
        "PEFF:0000008",
        "PEFF keyword for the sequence database name.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "Prefix",
        "PEFF:0000009",
        "PEFF keyword for the sequence database prefix.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "DbDescription",
        "PEFF:0000010",
        "PEFF keyword for the sequence database short description.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "Decoy",
        "PEFF:0000011",
        "PEFF keyword for the specifying whether the sequence database is a decoy database.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "DbSource",
        "PEFF:0000012",
        "PEFF keyword for the source of the database file.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "DbVersion",
        "PEFF:0000013",
        "PEFF keyword for the database version (release date) according to database provider.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "DbDate",
        "PEFF:0000014",
        "OBSOLETE PEFF keyword for the database date (release or file date of the source) according to database provider.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "NumberOfEntries",
        "PEFF:0000015",
        "PEFF keyword for the sumber of sequence entries in the database.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "Conversion",
        "PEFF:0000016",
        "PEFF keyword for the description of the conversion from original format to this current one.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "SequenceType",
        "PEFF:0000017",
        "PEFF keyword for the molecular type of the sequences.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "SpecificKey",
        "PEFF:0000018",
        "PEFF keyword for database specific keywords not included in the current controlled vocabulary.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "SpecificValue",
        "PEFF:0000019",
        "PEFF keyword for the specific values for a custom key.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "DatabaseDescription",
        "PEFF:0000020",
        "PEFF keyword for the short description of the PEFF file.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "GeneralComment",
        "PEFF:0000021",
        "PEFF keyword for a general comment.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "ProteoformDb",
        "PEFF:0000022",
        "PEFF keyword that when set to 'true' indicates that the database contains complete proteoforms.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "OptionalTagDef",
        "PEFF:0000023",
        "PEFF keyword for the short tag (abbreviation) and longer definition used to annotate a sequence annotation (such as variant or modification) in the OptionalTag location.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "HasAnnotationIdentifiers",
        "PEFF:0000024",
        "PEFF keyword that when set to 'true' indicates that entries in the database have identifiers for each annotation.",
        "PEFF CV Term",
        ["PEFF file header section term", "PEFF CV term"],
    ),
    PEFFTerm(
        "DbUniqueId",
        "PEFF:0001001",
        "OBSOLETE Sequence database unique identifier.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "PName",
        "PEFF:0001002",
        "PEFF keyword for the protein full name.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "NcbiTaxId",
        "PEFF:0001003",
        "PEFF keyword for the NCBI taxonomy identifier.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "TaxName",
        "PEFF:0001004",
        "PEFF keyword for the taxonomy name (latin or common name).",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "GName",
        "PEFF:0001005",
        "PEFF keyword for the gene name.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "Length",
        "PEFF:0001006",
        "PEFF keyword for the sequence length.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "SV",
        "PEFF:0001007",
        "PEFF keyword for the sequence version.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "EV",
        "PEFF:0001008",
        "PEFF keyword for the entry version.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "PE",
        "PEFF:0001009",
        "PEFF keyword for the Protein Evidence; A UniProtKB code 1-5.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "Processed",
        "PEFF:0001010",
        "PEFF keyword for information on how the full length original protein sequence can be processed into shorter components such as signal peptides and chains.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "Variant",
        "PEFF:0001011",
        "OBSOLETE Sequence variation (substitution, insertion, deletion).",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "ModResPsi",
        "PEFF:0001012",
        "PEFF keyword for the modified residue with PSI-MOD identifier.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "ModRes",
        "PEFF:0001013",
        "PEFF keyword for the modified residue without aPSI-MOD or UniMod identifier.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "AltAC",
        "PEFF:0001014",
        "PEFF keyword for the Alternative Accession Code.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "SeqStatus",
        "PEFF:0001015",
        "PEFF keyword for the sequence status. Complete or Fragment.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "CC",
        "PEFF:0001016",
        "PEFF keyword for the entry associated comment.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "KW",
        "PEFF:0001017",
        "PEFF keyword for the entry associated keyword(s).",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "GO",
        "PEFF:0001018",
        "PEFF keyword for the Gene Ontology code.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "XRef",
        "PEFF:0001019",
        "PEFF keyword for the cross-reference to an external resource.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "Conflict",
        "PEFF:0001023",
        "PEFF keyword for the sequence conflict; a UniProtKB term.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "Crc64",
        "PEFF:0001024",
        "PEFF keyword for the Sequence checksum in crc64.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "Domain",
        "PEFF:0001025",
        "PEFF keyword for the sequence range of a domain.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "ID",
        "PEFF:0001026",
        "PEFF keyword for the UniProtKB specific Protein identifier ID; a UniProtKB term.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "ModResUnimod",
        "PEFF:0001027",
        "PEFF keyword for the modified residue with UniMod identifier.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "VariantSimple",
        "PEFF:0001028",
        "PEFF keyword for the simple sequence variation of a single amino acid change. A change to a stop codon is permitted with a * symbol. More complex variations must be encoded with the VariantComplex term.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "VariantComplex",
        "PEFF:0001029",
        "PEFF keyword for a sequence variation that is more complex than a single amino acid change or change to a stop codon.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "Proteoform",
        "PEFF:0001030",
        "PEFF keyword for the proteoforms of this protein, constructed as a set of annotation identifiers.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "DisulfideBond",
        "PEFF:0001031",
        "PEFF keyword for the disulfide bonds in this protein, constructed as a sets of annotation identifiers of two half-cystine modifications.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "PEFF molecule processing keyword",
        "PEFF:0001032",
        "PEFF keyword describing the type of processing event being described.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "Comment",
        "PEFF:0001033",
        "PEFF keyword for the individual protein entry comment. It is discouraged to put parsable information here. This is only for free-text commentary.",
        "PEFF CV Term",
        ["PEFF file sequence entry term", "PEFF CV term"],
    ),
    PEFFTerm(
        "mature protein",
        "PEFF:0001020",
        "Portion of a newly synthesized protein that contributes to a final structure after other components such as signal peptides are removed.",
        "PEFF CV Term",
        [
            "PEFF molecule processing keyword",
            "PEFF file sequence entry term",
            "PEFF CV term",
        ],
    ),
    PEFFTerm(
        "signal peptide",
        "PEFF:0001021",
        "Short peptide present at the N-terminus of a newly synthesized protein that is cleaved off and is not part of the final mature protein.",
        "PEFF CV Term",
        [
            "PEFF molecule processing keyword",
            "PEFF file sequence entry term",
            "PEFF CV term",
        ],
    ),
    PEFFTerm(
        "transit peptide",
        "PEFF:0001022",
        "Short peptide present at the N-terminus of a newly synthesized protein that helps the protein through the membrane of its destination organelle.",
        "PEFF CV Term",
        [
            "PEFF molecule processing keyword",
            "PEFF file sequence entry term",
            "PEFF CV term",
        ],
    ),
    PEFFTerm(
        "propeptide",
        "PEFF:0001034",
        "Short peptide that is cleaved off a newly synthesized protein and generally immediately degraded in the process of protein maturation, and is not a signal peptide or transit peptide.",
        "PEFF CV Term",
        [
            "PEFF molecule processing keyword",
            "PEFF file sequence entry term",
            "PEFF CV term",
        ],
    ),
    PEFFTerm(
        "initiator methionine",
        "PEFF:0001035",
        "N-terminal methionine residue of a protein that can be co-translationally cleaved.",
        "PEFF CV Term",
        [
            "PEFF molecule processing keyword",
            "PEFF file sequence entry term",
            "PEFF CV term",
        ],
    ),
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
