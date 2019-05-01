from glypy.utils import Enum


class SequenceLocation(Enum):
    anywhere = None
    n_term = -1
    c_term = -2
    protein_n_term = -3
    protein_c_term = -4


SequenceLocation.n_term.add_name("N-term")
SequenceLocation.n_term.add_name("N-Term")
SequenceLocation.n_term.add_name("n-term")
SequenceLocation.n_term.add_name("N_term")
SequenceLocation.c_term.add_name("C-term")
SequenceLocation.c_term.add_name("C-Term")
SequenceLocation.c_term.add_name("c-term")
SequenceLocation.c_term.add_name("C_term")
SequenceLocation.anywhere.add_name("Anywhere")
SequenceLocation.protein_n_term.add_name("Protein N-term")
SequenceLocation.protein_c_term.add_name("Protein C-term")
SequenceLocation.protein_n_term.add_name("Protein N-Term")
SequenceLocation.protein_c_term.add_name("Protein C-Term")


class ModificationCategory(Enum):
    unknown = None
    glycosylation = 1
    artefact = 2
    substitution = 3
    chemical_derivative = 4
    non_standard_residue = 5
    isotopic_label = 6
    post_translational = 7
    other = 8
    multiple = 10
    pre_translational = 11
    co_translational = 12
    synthetic_peptide_protect = 13
    cross_link = 14
    cid_cleavable_cross_link = 15
    other_cleavable_cross_link = 16

ModificationCategory.substitution.add_name("AA substitution")
ModificationCategory.glycosylation.add_name("other_glycosylation")
ModificationCategory.glycosylation.add_name("N-linked glycosylation")
ModificationCategory.glycosylation.add_name("O-linked glycosylation")
ModificationCategory.glycosylation.add_name("Other glycosylation")
ModificationCategory.chemical_derivative.add_name("Chemical derivative")
ModificationCategory.post_translational.add_name("Post-translational")
ModificationCategory.multiple.add_name("Multiple")
ModificationCategory.artefact.add_name("Artefact")
ModificationCategory.isotopic_label.add_name("Isotopic label")
ModificationCategory.other.add_name("Other")
ModificationCategory.pre_translational.add_name("Pre-translational")
ModificationCategory.non_standard_residue.add_name("Non-standard residue")
ModificationCategory.co_translational.add_name("Co-translational")
ModificationCategory.synthetic_peptide_protect.add_name("Synth. pep. protect. gp.")
ModificationCategory.cross_link.add_name("Cross-link")
ModificationCategory.cid_cleavable_cross_link.add_name("CID cleavable cross-link")
ModificationCategory.other_cleavable_cross_link.add_name("Other cleavable cross-link")