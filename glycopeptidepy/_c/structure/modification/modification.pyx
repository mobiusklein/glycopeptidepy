from glycopeptidepy._c.structure.base cimport ModificationBase
from glycopeptidepy._c.structure.modification.rule cimport ModificationRuleBase


cdef class ModificationInstanceBase(ModificationBase):

    cpdef bint is_a(self, object category):
        return self.rule.is_a(category)

    cpdef basestring serialize(self):
        cdef:
            basestring rep
        if self.rule.is_standard:
            rep = self.name
        else:
            rep = self.rule.serialize()
        return rep