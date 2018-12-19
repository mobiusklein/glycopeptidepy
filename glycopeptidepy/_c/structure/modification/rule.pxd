from glycopeptidepy._c.structure.base cimport ModificationBase


cdef class ModificationRuleBase(ModificationBase):
    cdef:
        public set names
        public set categories

        public str preferred_name
        public str title
        public str common_name
        public str unimod_name