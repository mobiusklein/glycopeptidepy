cimport cython
from glypy.composition.ccomposition cimport CComposition


@cython.freelist(1000000)
cdef class PeptideSequenceBase(object):
    cdef:
        public list sequence
        public double _mass
        public object _glycosylation_manager
        public TerminalGroup _n_term
        public TerminalGroup _c_term

        public object _fragment_index
        public object _fragments_map
        public CComposition _total_composition
        public CComposition _peptide_composition

    cdef SequencePosition get(self, ssize_t i)

    cpdef _invalidate(self)


@cython.freelist(100)
cdef class TerminalGroup(object):
    cdef:
        public CComposition base_composition
        public double mass
        public ModificationBase _modification

    cdef ModificationBase get_modification(self)
    cdef CComposition get_composition(self)
    cpdef TerminalGroup  modify(self, ModificationBase modification)
    cdef void set_modification(self, ModificationBase)


@cython.freelist(100)
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



@cython.freelist(100000)
cdef class ModificationBase(object):
    cdef:
        public basestring name
        public double mass
        public CComposition composition

    cpdef bint is_a(self, object category)

    cpdef basestring serialize(self)

    cpdef bint is_tracked_for(self, object category)


cdef class SequencePosition(object):
    cdef:
        public AminoAcidResidueBase amino_acid
        public list modifications

    @staticmethod
    cdef SequencePosition _create(AminoAcidResidueBase amino_acid, list modifications)
