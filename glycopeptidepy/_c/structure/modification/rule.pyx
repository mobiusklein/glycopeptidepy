
cdef class ModificationRuleBase(ModificationBase):
    cpdef bint is_a(self, object category):
        '''Returns whether or not this :class:`ModificationRule` object belongs to
        the specified :class:`~.ModificationCategory`.

        Returns
        -------
        bool
        '''
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


cdef class NeutralLossBase(object):
    def __eq__(self, other):
        if other is None:
            return False
        return self.label == (<NeutralLossBase>other).label and self.mass == (<NeutralLossBase>other).mass

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.label)
