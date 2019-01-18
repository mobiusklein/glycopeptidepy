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

    cpdef _init_from_rule(self, ModificationRuleBase rule):
        self.name = rule.name
        self.mass = rule.mass
        self.rule = rule
        self._hash = self.rule._hash

        try:
            self.composition = rule.composition
        except AttributeError:
            self.composition = None

    cpdef bint is_tracked_for(self, category):
        return self.rule.is_tracked_for(category)