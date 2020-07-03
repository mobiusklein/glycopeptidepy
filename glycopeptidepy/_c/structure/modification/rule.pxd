from glypy.composition.ccomposition cimport CComposition
from glycopeptidepy._c.structure.base cimport ModificationBase


cdef class ModificationRuleBase(ModificationBase):
    cdef:
        public set names
        public list categories

        public basestring preferred_name
        public basestring title
        public basestring common_name
        public basestring unimod_name

        public Py_hash_t _hash


cdef class NeutralLossBase(object):
    cdef:
        public CComposition composition
        public double mass
        public object label


cpdef bint is_tracked_for_glycosylation(self, object category)