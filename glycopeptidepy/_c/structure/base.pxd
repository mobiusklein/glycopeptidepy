from glypy.composition.ccomposition cimport CComposition


cdef class AminoAcidResidueBase(object):
    cdef:
        public str name
        public str symbol
        public double mass
        public CComposition composition
        public object neutral_loss

    @staticmethod
    cdef AminoAcidResidueBase _create(str name, str symbol, double mass, CComposition composition)



cdef class ModificationBase(object):
    cdef:
        public str name
        public double mass
        public CComposition composition





cdef class SequencePosition(object):
    cdef:
        public AminoAcidResidueBase amino_acid
        public list modifications

    @staticmethod
    cdef SequencePosition _create(AminoAcidResidueBase amino_acid, list modifications)
