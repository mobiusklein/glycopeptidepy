
cdef class ModificationRuleBase(ModificationBase):
    cpdef bint is_a(self, object category):
        return category in self.categories
    
    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        if self is other:
            return True
        if not isinstance(self, ModificationRuleBase):
            temp = self
            self = <ModificationRuleBase>other
            other = temp
        idents = (<ModificationRuleBase>self).names
        if isinstance(other, ModificationRuleBase):
            other_idents = (<ModificationRuleBase>other).names
        else:
            other_idents = {other}
        return bool(idents & other_idents)

    def __ne__(self, other):
        return not self == other

