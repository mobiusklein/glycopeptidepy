cimport cython
from glycopeptidepy.structure.modification.descriptors import ModificationCategory
from glypy.utils.cenum cimport EnumValue

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

    cpdef basestring serialize(self):
        '''A string representation for inclusion in sequences'''
        return self.name

    @property
    def is_standard(self):
        """Indicates whether the :class:`ModificationRule` describes a modification which
        has been specified from a reference source, and may be fully reconstructed from
        its name alone.

        This property is not strictly enforced, and does not cover all aspects of the object's
        state. It does however provide a general guide for whether or not a modification is able
        to be translated from :mod:`glycopeptidepy`'s internal representation into something that
        another program *should* be able to understand.

        Returns
        -------
        bool
        """
        return True


cdef class NeutralLossBase(object):
    def __eq__(self, other):
        if other is None:
            return False
        return self.label == (<NeutralLossBase>other).label and self.mass == (<NeutralLossBase>other).mass

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.label)


cdef EnumValue ModificationCategory_glycosylation = ModificationCategory.glycosylation


cdef class GlycosylationBase(object):
    cdef bint get_is_composition(self):
        return self._is_composition

    cdef void set_is_composition(self, bint value):
        self._is_composition = value

    cdef bint get_is_core(self):
        return self._is_core

    cdef void set_is_core(self, bint value):
        self._is_core = value

    cdef EnumValue get_glycosylation_type(self):
        return self._glycosylation_type

    cdef void set_glycosylation_type(self, EnumValue glycosylation_type):
        self._glycosylation_type = glycosylation_type

    cpdef bint is_tracked_for(self, object category):
        return ModificationCategory_glycosylation == category

    @property
    def is_core(self):
        return self.get_is_core()

    @is_core.setter
    def is_core(self, bint value):
        self.set_is_core(value)

    @property
    def is_composition(self):
        return self.get_is_composition()

    @is_composition.setter
    def is_composition(self, bint value):
        self.set_is_composition(value)

    @property
    def glycosylation_type(self):
        return self.get_glycosylation_type()

    @glycosylation_type.setter
    def glycosylation_type(self, value):
        self.set_glycosylation_type(value)


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