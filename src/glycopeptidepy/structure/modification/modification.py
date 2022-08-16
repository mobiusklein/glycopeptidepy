from six import string_types as basestring

from glypy.utils import classproperty

from ..base import ModificationBase

from .source import ModificationTable
from .rule import ModificationRule


def _Modifcation_reconstructor():
    return Modification.__new__(Modification)


class ModificationInstanceBase(ModificationBase):
    def is_a(self, category):
        '''Returns whether or not this :class:`ModificationRule` object belongs to
        the specified :class:`~.ModificationCategory`.

        Returns
        -------
        bool
        '''
        return self.rule.is_a(category)

    def serialize(self):
        if self.rule.is_standard:
            rep = str(self.name)
        else:
            rep = self.rule.serialize()
        return rep

    def _init_from_rule(self, rule):
        self.name = rule.name
        self.mass = rule.mass
        self.rule = rule
        self._hash = self.rule._hash

        try:
            self.composition = rule.composition
        except AttributeError:
            self.composition = None

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        try:
            return other.name in self.rule.names
        except AttributeError:
            return other in self.rule.names

    def __ne__(self, other):
        return not self == other

    def is_tracked_for(self, category):
        """Determine if :attr:`rule` is tracked by a particular
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
        return self.rule.is_tracked_for(category)


try:
    _has_c = True
    from glycopeptidepy._c.structure.modification.modification import ModificationInstanceBase
except ImportError:
    _has_c = False


class Modification(ModificationInstanceBase):
    """Represents a peptide modification, with a fixed mass and optionally
    a composition.

    :class:`Modification` objects are usually bound to a :class:`~.ModificationRule`
    object which defines their properties. Though the essential properties are copied
    into the :class:`Modification` instance, the rule still holds additional metadata
    and behaviors.

    Modifications are equal to anything that their :attr:`rule` is equal to by name,
    and hash based upon their :attr:`name`.

    Attributes
    ----------
    name: str
        The preferred name of the modification
    mass: float
        The monosisotopic mass of the modification
    composition: :class:`~.Composition`
        The elemental composition of the modification. This
        attribute may be missing when the modification is not
        defined from a molecule but an arbitrary mass shift.
    rule: :class:`~.ModificationRule`
        The rule that defines the modification type this object represents
        an instance of.
    """

    # If the C extension base class is available, do not declare redundant slot
    # descriptors.
    if not _has_c:
        __slots__ = ["name", "mass", "rule", "composition", "_hash"]
    else:
        __slots__ = []

    _table = ModificationTable()

    @classmethod
    def register_new_rule(cls, rule):
        ModificationTable.register_new_rule(rule)

    def _resolve_name(self, rule_string):
        """Look up a :class:`~.ModificationRule` by name in :attr:`modification_table`.

        If the modification rule has not been registered with :attr:`modification_table`, and `rule_string`
        is not a full modification specification string, this may fail.

        Parameters
        ----------
        rule_string : :class:`str`
            The name or specification of a :class:`~.ModificationRule`
        """
        return self._table.resolve(rule_string)

    @classproperty
    def modification_table(self):
        """Access the class-wide :class:`~.ModificationTable` instance which all dynamic name
        lookups pass through.

        Returns
        -------
        :class:`~.ModificationTable`
        """
        return self._table

    @classmethod
    def from_rule(cls, rule):
        '''Create a new :class:`Modification` instance from a :class:`~.ModificationRule`
        instance directly.

        This method is marginally faster than the default initialization behavior which
        checks to see if the rule is a string.

        Parameters
        ----------
        rule: :class:`~.ModificationRule`
            The rule this modification will follow

        Returns
        -------
        :class:`Modification`
        '''
        inst = cls.__new__(cls)
        inst._init_from_rule(rule)
        return inst

    def __init__(self, rule):
        """Instantiate a :class:`Modification` instance from a template.

        Parameters
        ----------
        rule : :class:`str` or :class:`~.ModificationRule`
            Either the name of a modification rule as a string (to be resolved with
            :meth:`_resolve_name`) or an actual :class:`~.ModificationRule`, which will
            be used to populate all other attributes.
        """
        if(isinstance(rule, basestring)):
            rule = self._resolve_name(rule)

        self._init_from_rule(rule)

    @property
    def names(self):
        '''Returns the alternative names for this :class:`Modification`'s :attr:`rule`

        Returns
        -------
        frozenset
        '''
        return frozenset(self.rule.names)

    def __repr__(self):
        return self.serialize()

    def __str__(self):
        return self.serialize()

    def __reduce__(self):
        return _Modifcation_reconstructor, (), self.__getstate__()

    def __getstate__(self):
        return [self.name, self.mass, self.rule, self.composition]

    def __setstate__(self, state):
        self.name, self.mass, self.rule, self.composition = state
        self._hash = self.rule._hash

    def clone(self):
        return self.__class__(self.rule)

    def get_fragments(self, *args, **kwargs):
        return self.rule.get_fragments(*args, **kwargs)


# Late bind the ModificationRule's modification_tp field to avoid
# cyclic dependencies
ModificationRule.modification_tp = Modification
