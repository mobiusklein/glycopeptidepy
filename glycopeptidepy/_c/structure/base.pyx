from glypy.composition.ccomposition cimport CComposition


cdef class MoleculeBase(object):

    def __copy__(self):
        return self.clone()


cdef class ResidueBase(MoleculeBase):
    cdef:
        public str name
        public str symbol
        public double mass
        public CComposition composition
        public object neutral_loss


cdef class ModificationBase(MoleculeBase):
    cdef:
        public str name
        public double mass
        public CComposition composition

