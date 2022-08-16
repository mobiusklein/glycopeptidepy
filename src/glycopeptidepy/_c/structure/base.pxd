cimport cython
from glypy.composition.ccomposition cimport CComposition


cdef class PeptideSequenceBase(object):
    cdef:
        public list sequence
        public double _mass
        public object _glycosylation_manager
        public TerminalGroup _n_term
        public TerminalGroup _c_term

        public CComposition _total_composition
        public CComposition _peptide_composition

    cdef SequencePosition get(self, ssize_t i)

    cpdef _invalidate(self)


cdef class TerminalGroup(object):
    cdef:
        public CComposition base_composition
        public double mass
        public ModificationBase _modification

    cdef ModificationBase get_modification(self)
    cdef CComposition get_composition(self)
    cpdef TerminalGroup  modify(self, ModificationBase modification)
    cdef void set_modification(self, ModificationBase)
    cdef bint equal_to(self, object other)
    cpdef str serialize(self)


cdef class AminoAcidResidueBase(object):
    cdef:
        public basestring name
        public basestring symbol
        public double mass
        public CComposition composition
        public object neutral_loss
        public Py_hash_t _hash

    @staticmethod
    cdef AminoAcidResidueBase _create(str name, str symbol, double mass, CComposition composition)

    cdef bint equal_to(self, AminoAcidResidueBase other)


cdef class ModificationBase(object):
    cdef:
        public basestring name
        public double mass
        public CComposition composition

    cpdef bint is_a(self, object category)

    cpdef basestring serialize(self)

    cpdef bint is_tracked_for(self, object category)


IF int == long:
    DEF PY_VERSION = 3
ELSE:
    DEF PY_VERSION = 2


IF PY_VERSION == 3:
    ctypedef fused modification_or_string:
        str
        bytes
        ModificationBase

ELSE:
    ctypedef fused modification_or_string:
        str
        unicode
        ModificationBase


cdef class SequencePosition(object):
    cdef:
        public AminoAcidResidueBase amino_acid
        public list modifications

    @staticmethod
    cdef SequencePosition _create(AminoAcidResidueBase amino_acid, list modifications)

    cdef inline double get_mass(self)

    cpdef bint has_modification(self, modification)
    cpdef bint is_modified(self)
    cpdef drop_modification(self, modification)
    cpdef add_modification(self, modification)

    cdef inline size_t get_modification_count(self)
    cdef inline ModificationBase get_modification(self, size_t i)