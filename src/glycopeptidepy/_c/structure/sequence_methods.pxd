cimport cython

from glypy.composition.ccomposition cimport CComposition
from glycopeptidepy._c.structure.base cimport (
    PeptideSequenceBase,
    SequencePosition, TerminalGroup, AminoAcidResidueBase,
    ModificationBase)


cdef AminoAcidResidueBase _parse_residue(basestring residue_string)

cdef TerminalGroup _make_terminal_group(basestring base_composition_formula, ModificationBase modification=*)


cdef class _PeptideSequenceCore(PeptideSequenceBase):

    cpdef _init_from_parsed_string(self, list seq_list, glycan=*, n_term=*, c_term=*)

    cpdef _init_from_string(self, sequence, parser_function)

    cpdef _init_from_components(self, list seq_list, glycan_composition=*, ModificationBase n_term=*, ModificationBase c_term=*)

    cpdef _init_termini(self, TerminalGroup n_term, TerminalGroup c_term)

    cpdef _initialize_fields(self)

    cdef size_t get_size(self)
    cdef TerminalGroup get_n_term(self)
    cdef void set_n_term(self, _value)
    cdef TerminalGroup get_c_term(self)
    cdef void set_c_term(self, _value)

    cpdef basestring get_sequence(self, bint include_glycan=*, bint include_termini=*, str implicit_n_term=*, str implicit_c_term=*)
    cpdef clone(self)

    cpdef bint base_sequence_equality(self, _PeptideSequenceCore other)
    cpdef bint modified_sequence_equality(self, _PeptideSequenceCore other)
    cpdef bint full_structure_equality(self, _PeptideSequenceCore other)

    cpdef CComposition peptide_composition(self)
    cpdef CComposition total_composition(self)

    cdef double get_peptide_backbone_mass(self)
    cdef double get_total_mass(self)

cpdef add_modification(_PeptideSequenceCore self, position, modification_type)