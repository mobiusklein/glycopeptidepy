from collections import defaultdict

try:
    from collections.abc import Mapping
except ImportError:
    from collections import Mapping

from glycopeptidepy.utils.collectiontools import decoratordict

from glypy.utils import Enum

from glypy.structure.glycan import NamedGlycan
from glypy.structure.glycan_composition import (
    FrozenMonosaccharideResidue,
    FrozenGlycanComposition, GlycanComposition, HashableGlycanComposition)
from glypy import Composition


class allset(object):
    def __contains__(self, x):
        return True


class GlycosylationType(Enum):
    n_linked = 1
    o_linked = 2
    glycosaminoglycan = 3
    other = 4


GlycosylationType.n_linked.add_name("N-Linked")
GlycosylationType.n_linked.add_name("N-linked")
GlycosylationType.n_linked.add_name("n-linked")
GlycosylationType.n_linked.add_name("N-Glycan")
GlycosylationType.n_linked.add_name("n-glycan")
GlycosylationType.n_linked.add_name("n_glycan")
GlycosylationType.n_linked.name = "N-Linked"

GlycosylationType.o_linked.add_name("O-Linked")
GlycosylationType.o_linked.add_name("O-linked")
GlycosylationType.o_linked.add_name("o-linked")
GlycosylationType.o_linked.add_name("O-Glycan")
GlycosylationType.o_linked.add_name("o-glycan")
GlycosylationType.o_linked.add_name("o_glycan")
GlycosylationType.o_linked.name = "O-Linked"

GlycosylationType.glycosaminoglycan.add_name("Glycosaminoglycan")
GlycosylationType.glycosaminoglycan.add_name("GAG-Linker")
GlycosylationType.glycosaminoglycan.add_name("gag-linker")
GlycosylationType.glycosaminoglycan.add_name("gag_linker")
GlycosylationType.glycosaminoglycan.name = "Glycosaminoglycan"


glycosylation_site_detectors = decoratordict({
    GlycosylationType.n_linked: lambda x: [],
    GlycosylationType.o_linked: lambda x: [],
    GlycosylationType.glycosaminoglycan: lambda x: []
})


class TypedGlycanComposition(HashableGlycanComposition):

    def __init__(self, glycosylation_type=None, *args, **kwargs):
        if isinstance(glycosylation_type, GlycanComposition):
            temp = args[0]
            args = (glycosylation_type,) + args[1:]
            glycosylation_type = temp
        self.glycosylation_type = GlycosylationType[glycosylation_type]
        super(TypedGlycanComposition, self).__init__(*args, **kwargs)

    @classmethod
    def _empty(cls):
        inst = super(cls, cls)._empty()
        inst.glycosylation_type = None
        return inst

    def clone(self, propogate_composition_offset=True):
        dup = self.__class__(self.glycosylation_type, self)
        if not propogate_composition_offset:
            dup.composition_offset = Composition('H2O')
        return dup

    def is_type(self, glycosylation_type):
        return self.glycosylation_type is GlycosylationType[glycosylation_type]

    def as_modification(self):
        from .modification import Glycosylation
        return Glycosylation(self, metadata={
            "glycosylation_type": self.glycosylation_type.name
        })


class GlycanCombination(Mapping):  # pragma: no cover
    def __init__(self, components):
        self.components = tuple(components)
        self._glycan_composition = self._sum_components(self.components)

    def _sum_components(self, components):
        acc = HashableGlycanComposition()
        for x in components:
            acc += x
        acc.composition_offset -= Composition({"H": 2, "O": 1})
        return acc

    def mass(self, *args, **kwargs):
        return self._glycan_composition.mass(*args, **kwargs)

    def total_composition(self, *args, **kwargs):
        return self._glycan_composition.total_composition(*args, **kwargs)

    def __iter__(self):
        return iter(self._glycan_composition)

    def __getitem__(self, key):
        return self._glycan_composition[key]

    def keys(self):
        return self._glycan_composition.keys()

    def values(self):
        return self._glycan_composition.values()

    def items(self):
        return self._glycan_composition.items()

    def clone(self, *args, **kwargs):
        return self._glycan_composition.clone(*args, **kwargs)

    def serialize(self):
        parts = [str(x)[1:-1] for x in self.components]
        return '{%s}' % '|'.join(parts)

    def __repr__(self):
        return self.serialize()

    def __str__(self):
        return self.serialize()

    @classmethod
    def parse(cls, string):
        string = str(string)[1:-1]
        parts = string.split("|")
        components = []
        for part in parts:
            component = HashableGlycanComposition()
            tokens = part.split("; ")
            for token in tokens:
                try:
                    residue, count = token.split(":")
                except ValueError:
                    if string == "{}":
                        return cls([])
                    else:
                        raise ValueError("Malformed Token, %s" % (token,))
                component[FrozenMonosaccharideResidue.from_iupac_lite(residue)] = int(count)
            components.append(component)
        return cls(components)

    def __hash__(self):
        return hash(self._glycan_composition)

    def __eq__(self, other):
        try:
            other = other._glycan_composition
        except AttributeError:
            pass
        return self._glycan_composition == other

    def __ne__(self, other):
        try:
            other = other._glycan_composition
        except AttributeError:
            pass
        return (self._glycan_composition != other)

    def __add__(self, other):
        return self._glycan_composition + other

    def __sub__(self, other):
        return self._glycan_composition - other

    def __mul__(self, other):
        return self._glycan_composition * other

    def __len__(self):
        return len(self._glycan_composition)


class TypedGlycan(NamedGlycan):

    def __init__(self, glycosylation_type=None, *args, **kwargs):
        self.glycosylation_type = GlycosylationType[glycosylation_type]
        super(TypedGlycan, self).__init__(*args, **kwargs)

    def clone(self, index_method='dfs', cls=None):
        if cls is None:
            cls = TypedGlycan
        inst = super(TypedGlycan, self).clone(index_method=index_method, cls=cls)
        inst.glycosylation_type = self.glycosylation_type
        return inst

    def __repr__(self):
        rep = super(TypedGlycan, self).__repr__()
        return "%s\n%s" % (self.glycosylation_type, rep)

    def is_type(self, glycosylation_type):
        return self.glycosylation_type == GlycosylationType[glycosylation_type]

    def as_modification(self, encoding_format=None):
        from .modification import Glycosylation
        return Glycosylation(self, encoding_format, metadata={
            "glycosylation_type": self.glycosylation_type.name
        })


class GlycanCompositionProxy(Mapping):
    '''A mapping-like object that imitates the GlycanComposition interface in
    a read-only fashion.
    '''

    def __init__(self, glycan_composition):
        self.obj = glycan_composition

    def mass(self, *args, **kwargs):
        return self.obj.mass(*args, **kwargs)

    def total_composition(self, *args, **kwargs):
        return self.obj.total_composition(*args, **kwargs)

    def __iter__(self):
        return iter(self.obj)

    def __getitem__(self, key):
        return self.obj[key]

    def _getitem_fast(self, key):
        return self.obj._getitem_fast(key)

    def __contains__(self, key):
        return key in self.obj

    def keys(self):
        return self.obj.keys()

    def values(self):
        return self.obj.values()

    def items(self):
        return self.obj.items()

    def clone(self, *args, **kwargs):
        return self.obj.clone(*args, **kwargs)

    def query(self, query, exact=True, **kwargs):
        return self.obj.query(query, exact=exact, **kwargs)

    def __repr__(self):
        return repr(self.obj)

    def __str__(self):
        return str(self.obj)

    def __hash__(self):
        return hash(self.obj)

    def __eq__(self, other):
        try:
            other = other.obj
        except AttributeError:
            pass
        return self.obj == other

    def __ne__(self, other):
        try:
            other = other.obj
        except AttributeError:
            pass
        return (self.obj != other)

    def __add__(self, other):
        return self.obj + other

    def __sub__(self, other):
        return self.obj - other

    def __mul__(self, other):
        return self.obj * other

    def __len__(self):
        return len(self.obj)


class _GlycanCompositionWithOffsetProxyBase(object):
    __slots__ = ('composition_offset', )

    def __init__(self, obj, offset=None):
        if offset is None:
            offset = Composition()
        self.obj = obj
        self.composition_offset = offset


try:
    from glycopeptidepy._c.structure.glycan import GlycanCompositionWithOffsetProxyBase as _GlycanCompositionWithOffsetProxyBase
except ImportError:
    pass


class GlycanCompositionWithOffsetProxy(_GlycanCompositionWithOffsetProxyBase, GlycanCompositionProxy):
    """A :class:`GlycanCompositionProxy` that pretends to allow you to
    mutate the :attr:`composition_offset` attribute of the underlying
    composition, but instead stores that within the proxy, leaving the
    underlying compostion unchanged.

    This type should be useful for re-using the same :class:`FrozenGlycanComposition`
    object without repeatedly copying and invalidating its cached properties. When
    this object is copied, the underlying GlycanComposition is not copied, only
    the proxy.
    """

    def clone(self, *args, **kwargs):
        return self.__class__(self.obj, self.composition_offset.copy())

    def mass(self, *args, **kwargs):
        mass = self.obj.mass(*args, **kwargs)
        return mass + self.composition_offset.mass

    def total_composition(self):
        comp = self.obj.total_composition()
        return comp + self.composition_offset


WATER_OFFSET = Composition({"H": 2, "O": 1})

class GlycosylationManager(object):
    def __init__(self, parent, aggregate=None):
        self.mapping = {}
        self.parent = parent
        self._aggregate = None
        self._proxy = None
        self._type_track = None
        if aggregate is not None:
            self.set_aggregate(aggregate)
        self._total_glycosylation_size = None

    def clear(self):
        self.invalidate()
        self.aggregate = None
        self.mapping.clear()

    def __setitem__(self, key, value):
        self.invalidate()
        self.mapping[key] = value

    def __delitem__(self, key):
        self.invalidate()
        del self.mapping[key]

    def pop(self, key, default=None):
        self.invalidate()
        self.mapping.pop(key, default)

    def update(self, base):
        self.invalidate()
        self.mapping.update(base)

    def keys(self):
        return self.mapping.keys()

    def values(self):
        return self.mapping.values()

    def items(self):
        return self.mapping.items()

    def __getitem__(self, key):
        return self.mapping[key]

    def __len__(self):
        return len(self.mapping)

    def invalidate(self):
        self._proxy = None
        self._type_track = None
        self._total_glycosylation_size = None

    @property
    def glycosylation_types(self):
        if self._type_track is None:
            self._type_track = self._update_track()
        return self._type_track

    def count_glycosylation_type(self, glycosylation_type):
        return len(self.glycosylation_types[glycosylation_type])

    def _update_track(self):
        track = defaultdict(list)
        for position, value in self.items():
            glycosylation = value.rule
            track[glycosylation.glycosylation_type].append((position, value))
        return track

    setdefault = None
    popitem = None

    def copy(self):
        inst = self.__class__(
            self.parent,
            self.aggregate.clone() if self.aggregate is not None else None)
        inst.update(self)
        return inst

    @property
    def aggregate(self):
        return self._aggregate

    @aggregate.setter
    def aggregate(self, value):
        self._aggregate = value
        if self._aggregate is not None:
            self._patch_aggregate()
        self.invalidate()

    def _patch_aggregate(self):
        self.aggregate.composition_offset -= WATER_OFFSET

    @property
    def glycan_composition(self):
        if self._proxy is None:
            self._proxy = self._make_glycan_composition_proxy()
        return self._proxy

    def total_glycosylation_size(self):
        if self._total_glycosylation_size is None:
            self._total_glycosylation_size = sum(self.aggregate.values())
        return self._total_glycosylation_size

    def clone(self):
        return self.copy()

    def total_composition(self):
        total = Composition()
        has_aggregate = self.aggregate is not None
        for key, value in self.items():
            if has_aggregate and value.rule.is_core:
                continue
            total += value.composition
        if has_aggregate:
            total += self.aggregate.total_composition()
        return total

    def mass(self):
        total = 0.
        has_aggregate = self.aggregate is not None
        for key, value in self.items():
            if has_aggregate and value.rule.is_core:
                continue
            total += value.mass
        if has_aggregate:
            total += self.aggregate.mass()
        return total

    def is_fully_specified_topologies(self):
        is_fully_specified = len(self) > 0
        for key, value in self.items():
            if value.rule.is_core:
                is_fully_specified = False
                break
            elif value.rule.is_composition:
                is_fully_specified = False
                break
        return is_fully_specified

    def _make_glycan_composition_proxy(self):
        if self.aggregate is not None:
            base = self.aggregate.clone()
        else:
            base = HashableGlycanComposition()
            # Represent the initial amide bond between the peptide
            # and the first glycan. Subsequent glycans do not need
            # further chemical losses because of the dehyration built
            # directly into the Residue abstraction.
            base.composition_offset -= WATER_OFFSET
        for key, value in self.items():
            if value.rule.is_core:
                continue
            elif value.rule.is_composition:
                base += value.rule.glycan
            else:
                # Convert Glycan object into a composition, using the original
                # detatched topology to omit the "aglycone" group which represents
                # the connection between the glycan and the peptide, which penalizes
                # the composition by H2O. This H2O is lost when that bond is formed,
                # but doesn't need to be explicitly included as the loss is tracked
                # when initializing the base above.
                gc = HashableGlycanComposition.from_glycan(
                    value.rule._original)
                base += gc
        return GlycanCompositionProxy(base)


try:
    from glycopeptidepy._c.structure.glycan import GlycosylationManager, set_glycan_composition_proxy_type
    set_glycan_composition_proxy_type(GlycanCompositionProxy)
    del set_glycan_composition_proxy_type
except ImportError:
    pass
