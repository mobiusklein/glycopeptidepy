cimport cython

from glycopeptidepy._c.structure.base cimport (
    PeptideSequenceBase,
    SequencePosition, TerminalGroup, AminoAcidResidueBase,
    ModificationBase)


cdef AminoAcidResidueBase _parse_residue(basestring residue_string)

cdef TerminalGroup _make_terminal_group(basestring base_composition_formula, ModificationBase modification=*)


cpdef _init_from_parsed_string(PeptideSequenceBase self, list seq_list, glycan=*, n_term=*, c_term=*)


cpdef _init_from_components(PeptideSequenceBase self, list seq_list, glycan_composition=*,
                            ModificationBase n_term=*, ModificationBase c_term=*)


cpdef _init_termini(PeptideSequenceBase self, TerminalGroup n_term, TerminalGroup c_term)