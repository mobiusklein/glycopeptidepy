
cdef class ModificationRuleBase(ModificationBase):
    cpdef bint is_a(self, object category):
        return category in self.categories
    
    def __hash__(self):
        return self._hash
