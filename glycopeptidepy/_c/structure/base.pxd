cimport cython
from glypy.composition.ccomposition cimport CComposition


@cython.freelist(100)
cdef class AminoAcidResidueBase(object):
    cdef:
        public basestring name
        public basestring symbol
        public double mass
        public CComposition composition
        public object neutral_loss

    @staticmethod
    cdef AminoAcidResidueBase _create(str name, str symbol, double mass, CComposition composition)



@cython.freelist(100000)
cdef class ModificationBase(object):
    cdef:
        public basestring name
        public double mass
        public CComposition composition

    cpdef bint is_a(self, object category)


@cython.freelist(100000)
cdef class SequencePosition(object):
    cdef:
        public AminoAcidResidueBase amino_acid
        public list modifications

    @staticmethod
    cdef SequencePosition _create(AminoAcidResidueBase amino_acid, list modifications)
