from glycopeptidepy._c.structure.base cimport ModificationBase
from glycopeptidepy._c.structure.modification.rule cimport ModificationRuleBase


cdef class ModificationInstanceBase(ModificationBase):
    cdef:
        public ModificationRuleBase rule
        public object _hash

    cpdef _init_from_rule(self, ModificationRuleBase rule)
