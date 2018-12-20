
cdef class ModificationRuleBase(ModificationBase):
    cpdef bint is_a(self, object category):
        return category in self.categories