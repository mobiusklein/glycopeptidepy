from glycopeptidepy._c.structure.base cimport ModificationBase
from glycopeptidepy._c.structure.modification.rule cimport ModificationRuleBase


cdef class ModificationInstanceBase(ModificationBase):
    cdef:
        public ModificationRuleBase rule
        public Py_hash_t _hash

    cpdef _init_from_rule(self, ModificationRuleBase rule)

    cpdef bint is_tracked_for(self, category)
