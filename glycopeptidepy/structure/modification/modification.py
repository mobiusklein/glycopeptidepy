from six import string_types as basestring

from ..base import ModificationBase

from .source import ModificationTable
from .rule import ModificationRule


def _Modifcation_reconstructor():
    return Modification.__new__(Modification)


class ModificationInstanceBase(ModificationBase):
    def is_a(self, category):
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


try:
    _has_c = True
    from glycopeptidepy._c.structure.modification.modification import ModificationInstanceBase
except ImportError:
    _has_c = False
    pass


class Modification(ModificationInstanceBase):

    """Represents a molecular modification, which may be bound at a given position,
    or may be present at multiple locations. This class pulls double duty,
    and when describing the total number of modifications of a type on a molecule, its
    position attributes are unused."""

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
        return self._table.resolve(rule_string)

    def __init__(self, rule):
        if(isinstance(rule, basestring)):
            rule = self._resolve_name(rule)

        self._init_from_rule(rule)

    def __repr__(self):
        return self.serialize()

    def __str__(self):
        return self.serialize()

    def __eq__(self, other):
        try:
            return other.name in self.rule.names
        except AttributeError:
            return other in self.rule.names

    def __ne__(self, other):
        return not self == other

    def __reduce__(self):
        return _Modifcation_reconstructor, (), self.__getstate__()

    def __getstate__(self):
        return [self.name, self.mass, self.rule, self.composition]

    def __setstate__(self, state):
        self.name, self.mass, self.rule, self.composition = state
        self._hash = self.rule._hash

    def clone(self):
        return self.__class__(self.rule)

    def is_tracked_for(self, category):
        return self.rule.is_tracked_for(category)

    def get_fragments(self, *args, **kwargs):
        return self.rule.get_fragments(*args, **kwargs)


# Late bind the ModificationRule's modification_tp field to avoid
# cyclic dependencies
ModificationRule.modification_tp = Modification
