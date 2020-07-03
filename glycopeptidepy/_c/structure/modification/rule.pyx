cimport cython
from glycopeptidepy.structure.modification.descriptors import ModificationCategory

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


cdef object ModificationCategory_glycosylation = ModificationCategory.glycosylation


@cython.binding(True)
cpdef bint is_tracked_for_glycosylation(self, object category):
    """Determine if this :class:`ModificationBase` is tracked by a particular
    behavioral pattern associated with a :class:`~.ModificationCategory`.

    This relationship is distinct from :meth:`is_a` which merely observes that
    the semantic relationship holds, not that any actual behavior is available.

    Parameters
    ----------
    category : :class:`~.ModificationCategory`
        The category to check

    Returns
    -------
    bool
    """
    return category == ModificationCategory_glycosylation