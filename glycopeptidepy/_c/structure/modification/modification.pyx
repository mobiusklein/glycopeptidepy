from glycopeptidepy._c.structure.base cimport ModificationBase
from glycopeptidepy._c.structure.modification.rule cimport ModificationRuleBase


cdef class ModificationInstanceBase(ModificationBase):

    cpdef bint is_a(self, object category):
        return self.rule.is_a(category)
