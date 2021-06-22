from glycopeptidepy.utils.memoize import memoize

from .base import MoleculeBase
from .modification import Modification
from .composition import Composition, formula

from six import string_types as basestring


class TerminalGroup(MoleculeBase):
    __slots__ = ("base_composition", "mass", "_modification")

    def __init__(self, base_composition, modification=None):
        if not isinstance(base_composition, Composition):
            base_composition = Composition(base_composition)
        self.mass = 0
        self.base_composition = base_composition
        self._modification = None
        if modification is not None:
            self.modification = modification
        self.mass = self._calculate_mass()

    def _calculate_mass(self):
        base_mass = self.base_composition.mass
        mod = self.modification
        if mod is not None:
            base_mass += mod.mass
        return base_mass

    def clone(self):
        return self.__class__(self.base_composition, self.modification)

    def __reduce__(self):
        return self.__class__, (self.base_composition, self.modification)

    @property
    def modification(self):
        return self._modification

    @modification.setter
    def modification(self, value):
        if value is not None:
            new_mass = value.mass
        else:
            new_mass = 0
        if self._modification is not None:
            old_mass = self._modification.mass
        else:
            old_mass = 0
        self.mass += new_mass - old_mass
        self._modification = value

    def modify(self, modification):
        return self.__class__(self.base_composition, modification)

    @property
    def composition(self):
        modification = self.modification
        if modification is None:
            return self.base_composition
        mod_comp = modification.composition
        return self.base_composition + mod_comp

    def __repr__(self):
        template = "{self.__class__.__name__}({self.base_composition}, {self.modification})"
        return template.format(self=self)

    def __str__(self):
        if self.modification is not None:
            return str(self.modification)
        return formula(self.base_composition)

    def __eq__(self, other):
        if other is None:
            return False
        try:
            return (self.base_composition == other.base_composition) and (self.modification == other.modification)
        except AttributeError:
            if isinstance(other, basestring):
                return self.composition == Modification(other).composition
            else:
                try:
                    return self.composition == other.composition
                except AttributeError:
                    return NotImplemented

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(formula(self.composition))

    def serialize(self):
        return str(self)


try:
    from glycopeptidepy._c.structure.base import TerminalGroup
except ImportError:
    pass


@memoize(100)
def _make_terminal_group(base_composition_formula, modification=None):
    return TerminalGroup(base_composition_formula, modification)
