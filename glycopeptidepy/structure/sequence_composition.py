# pragma: no cover

import itertools

try:
    basestring
except Exception:
    from six import string_types as basestring

from glypy import MonosaccharideResidue

from .sequence import PeptideSequenceBase
from .residue import Residue as AminoAcidResidue, memoize, get_all_residues
from .composition import Composition
from .modification import Modification
from .base import SequencePosition

class AminoAcidSequenceBuildingBlock(object):
    @classmethod
    def get_all_common_residues(cls):
        return map(cls, get_all_residues())

    @classmethod
    def get_all_sequencing_residues(cls):
        residues = get_all_residues()
        for residue in list(residues):
            degenerate = residue.is_degenerate
            if degenerate:
                for degen in degenerate:
                    residues.remove(degen)
        return map(cls, residues)

    def __init__(self, residue_, modifications=None, mass=None):
        if modifications is None:
            modifications = ()
        self.residue = residue_
        self.modifications = tuple(modifications)
        if mass is None:
            mass = residue_.mass + sum(m.mass for m in modifications)
        self.mass = mass
        self._string = None

    def __iter__(self):
        yield self.residue
        yield self.modifications

    def _reset(self):
        self._string = None

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        try:
            return self.residue is other.residue and self.modifications == other.modifications
        except AttributeError:
            return other == str(self)

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return "({}, {}):{:.2f}".format(self.residue.symbol, self.modifications, self.mass)

    def __str__(self):
        if self._string is None:
            self._string = "{}{}".format(
                self.residue.symbol,
                "({})".format(self.modifications[0]) if len(self.modifications) > 0 else "")
        return self._string

    def __contains__(self, item):
        return self == item

    @classmethod
    @memoize()
    def from_str(cls, string):
        # A hack to get ModificationBuildingBlock into the same flow of execution
        # Move this logic to a pure function later.
        if string.startswith("@"):
            return ModificationBuildingBlock(string)
        parts = string.split("(")
        aa = AminoAcidResidue(parts[0])
        mod = tuple()
        if len(parts) > 1:
            mod = (Modification(parts[1][:-1]),)
        return cls(aa, mod)

    def orderings(self):
        yield (self)

    def clone(self):
        return self.__class__(self.residue, tuple(self.modifications))

    def as_sequence_position(self):
        return SequencePosition([self.residue, list(self.modifications)])

    @classmethod
    def from_sequence_position(cls, position):
        return cls(position.amino_acid, tuple(position.modifications))


class ModificationBuildingBlock(object):
    # Not all SequenceComposition methods are compatible with these?
    sigil = '@'

    def __init__(self, name):
        try:
            if name.startswith(self.sigil):
                self.name = name
                self.modification = Modification(name[1:])
            else:
                self.name = self.sigil + name
                self.modification = Modification(name)
        except AttributeError:
            self.name = self.sigil + str(name)
            self.modification = name
        self.mass = self.modification.mass
        self._hash = hash(self.modification)

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        try:
            return self.modification == other.modification
        except AttributeError:
            return str(self) == str(other)

    def __ne__(self, other):
        try:
            return self.modification != other.modification
        except AttributeError:
            return str(self) != str(other)

    def __hash__(self):
        return self._hash

    def clone(self):
        return self.__class__(self.name)


n_hexnac = AminoAcidSequenceBuildingBlock.from_str("N(HexNAc)")


class MonosaccharideResidueAdapter(object):
    def __init__(self, residue):
        if isinstance(residue, basestring):
            residue = MonosaccharideResidue.from_iupac_lite(residue)
        self.residue = residue
        self.mass = residue.mass()
        self.symbol = "(%s)" % str(residue)

    def __hash__(self):
        return hash(self.residue)

    def __eq__(self, other):
        if isinstance(other, (MonosaccharideResidue, AminoAcidResidue)):
            return self.symbol == other.symbol
        else:
            return self.symbol == other

    def __ne__(self, other):
        return not self == other


class SequenceComposition(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self)
        self._mass = None
        self._composition_offset = Composition("H2O")
        self.extend(*args)

    def clone(self):
        inst = self.__class__()
        inst += self
        return inst

    def __setitem__(self, key, value):
        if key is None:
            return
        if isinstance(key, basestring):
            key = AminoAcidSequenceBuildingBlock.from_str(key)
        dict.__setitem__(self, key, value)
        self._mass = None

    def __getitem__(self, key):
        if isinstance(key, basestring):
            key = AminoAcidSequenceBuildingBlock.from_str(key)
        return dict.__getitem__(self, key)

    def __contains__(self, key):
        if isinstance(key, dict):
            for k, v in key.items():
                if self[k] - v < 0:
                    return False
            return True
        else:
            if isinstance(key, basestring):
                key = AminoAcidSequenceBuildingBlock.from_str(key)
            return dict.__contains__(self, key)

    @property
    def mass(self):
        if self._mass is not None:
            return self._mass
        mass = self._composition_offset.mass
        for residue_type, count in self.items():
            mass += residue_type.mass * count
        self._mass = mass
        return mass

    def extend(self, *args):
        for residue in args:
            if isinstance(residue, PeptideSequenceBase):
                for pos in residue:
                    self[AminoAcidSequenceBuildingBlock.from_sequence_position(pos)] += 1
            else:
                self[residue] += 1

    def __iadd__(self, other):
        for elem, cnt in (other.items()):
            self[elem] += cnt
        return self

    def __add__(self, other):
        result = self.clone()
        for elem, cnt in other.items():
            result[elem] += cnt
        return result

    def __radd__(self, other):
        return self + other

    def __isub__(self, other):
        for elem, cnt in other.items():
            self[elem] -= cnt
        return self

    def __sub__(self, other):
        result = self.clone()
        for elem, cnt in other.items():
            result[elem] -= cnt
        return result

    def __rsub__(self, other):
        return (self - other) * (-1)

    def __mul__(self, other):
        if not isinstance(other, int):
            raise TypeError(
                'Cannot multiply Composition by non-integer',
                other)
        prod = SequenceComposition()
        for k, v in self.items():
            prod[k] = v * other

        return prod

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        if not isinstance(other, dict):
            return False
        self_items = set([i for i in self.items() if i[1]])
        other_items = set([i for i in other.items() if i[1]])
        return self_items == other_items

    def __neg__(self):
        return -1 * self

    def __missing__(self, key):
        return 0

    @property
    def composition_offset(self):
        return self._composition_offset

    @composition_offset.setter
    def composition_offset(self, value):
        self._mass = None
        self._composition_offset = value

    def flatten(self):
        for k in self:
            for i in range(0, self[k]):
                yield k

    def serialize(self):
        return "{%s}" % '; '.join("{}:{}".format(str(k), v) for k, v in sorted(
            self.items(), key=lambda x: x[0].mass))

    __str__ = serialize

    def __hash__(self):
        return hash(str(self))

    @classmethod
    def parse(cls, string):
        inst = cls()
        tokens = string[1:-1].split('; ')
        for token in tokens:
            residue, count = token.split(":")
            inst[AminoAcidSequenceBuildingBlock.from_str(residue)] = int(count)
        return inst

    def maybe_n_glycosylation(self):
        n_N = self["N"] + self[n_hexnac]
        n_S = self["S"]
        n_T = self["T"]
        used = 0
        possible_sequons = 0
        while(n_N > 0):
            if n_S > 0 or n_T > 0:
                if sum(self.values()) - used >= 3:
                    n_N -= 1
                    if n_S > 0:
                        n_S -= 1
                    elif n_T > 0:
                        n_T -= 1
                    else:
                        break
                    used += 3
                    possible_sequons += 1
                    continue
            break
        return possible_sequons

    def orderings(self):
        for ordering in itertools.permutations(self.flatten()):
            yield ordering


def all_compositions(blocks, size=20):
    components = [None] + list(blocks)
    for case in itertools.combinations_with_replacement(components, size):
        sc = SequenceComposition(*case)
        if len(sc) == 0:
            continue
        yield sc
