from glycopeptidepy._c.structure.base cimport ModificationBase


cdef class ModificationRuleBase(ModificationBase):
    cdef:
        public set names
        public list categories

        public basestring preferred_name
        public basestring title
        public basestring common_name
        public basestring unimod_name
