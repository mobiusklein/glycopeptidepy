from .structure.sequence import PeptideSequence, ProteinSequence
from .structure.modification import (
    Modification, ModificationTable, RestrictedModificationTable,
    ModificationRule, ModificationTarget, AnonymousModificationRule,
    Glycosylation, SequenceLocation, ModificationCategory)
from .structure.fragment import (
    IonSeries, PeptideFragment, ChemicalShift,
    SimpleFragment, StubFragment)
from .structure.residue import (
    AminoAcidResidue, get_all_residues, register_degenerate, register_residue)
from .enzyme import Protease, cleave
from .structure.glycan import (
    HashableGlycanComposition, GlycosylationType, TypedGlycanComposition)
from .structure.composition import (Composition, ChemicalCompositionError)


def parse(string, peptide_class=PeptideSequence):
    return peptide_class(string)


__all__ = [
    "PeptideSequence", "ProteinSequence", "parse",

    "AminoAcidResidue", "get_all_residues", "register_residue", "register_degenerate",

    "Modification", "ModificationTable", "RestrictedModificationTable",
    "ModificationRule", "ModificationTarget", "AnonymousModificationRule",

    "IonSeries", "PeptideFragment", "ChemicalShift", "SimpleFragment", "StubFragment",

    "Protease", "cleave",

    "HashableGlycanComposition", "GlycosylationType", "TypedGlycanComposition",
    "Glycosylation",

    "Composition", "ChemicalCompositionError",
]
