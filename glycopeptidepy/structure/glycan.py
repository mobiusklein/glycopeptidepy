from collections import Mapping

from glycopeptidepy.utils.collectiontools import decoratordict

from glypy.utils import Enum

from glypy.structure.glycan import NamedGlycan
from glypy import Composition
from glypy.composition.glycan_composition import HashableGlycanComposition


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
GlycosylationType.n_linked.name = "N-Linked"

GlycosylationType.o_linked.add_name("O-Linked")
GlycosylationType.o_linked.add_name("O-linked")
GlycosylationType.o_linked.add_name("o-linked")
GlycosylationType.o_linked.name = "O-Linked"

GlycosylationType.glycosaminoglycan.add_name("Glycosaminoglycan")
GlycosylationType.glycosaminoglycan.name = "Glycosaminoglycan"


glycosylation_site_detectors = decoratordict({
    GlycosylationType.n_linked: lambda x: [],
    GlycosylationType.o_linked: lambda x: [],
    GlycosylationType.glycosaminoglycan: lambda x: []
})


class TypedGlycanComposition(HashableGlycanComposition):

    def __init__(self, glycosylation_type=None, *args, **kwargs):
        self.glycosylation_type = GlycosylationType[glycosylation_type]
        super(TypedGlycanComposition, self).__init__(*args, **kwargs)

    def is_type(self, glycosylation_type):
        return self.glycosylation_type is GlycosylationType[glycosylation_type]


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


class GlycanCompositionProxy(Mapping):
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

    def keys(self):
        return self.obj.keys()

    def values(self):
        return self.obj.values()

    def items(self):
        return self.obj.items()

    def clone(self, *args, **kwargs):
        return self.obj.clone(*args, **kwargs)

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


class GlycosylationManager(dict):
    def __init__(self, parent, aggregate=None):
        self.parent = parent
        self._aggregate = None
        self.aggregate = aggregate
        self._proxy = None

    def invalidate(self):
        self._proxy = None

    def clear(self):
        self.invalidate()
        self.aggregate = None
        dict.clear(self)

    def __setitem__(self, key, value):
        self.invalidate()
        dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        self.invalidate()
        dict.__delitem__(self, key)

    def pop(self, key, default=None):
        self.invalidate()
        dict.pop(self, key, default)

    def update(self, base):
        self.invalidate()
        dict.update(self, base)

    setdefault = None
    popitem = None

    def copy(self):
        inst = self.__class__(
            self.parent,
            self.aggregate.clone() if self.aggregate is not None else None)
        inst.update(self)
        return inst

    def clone(self):
        return self.copy()

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
        offset = Composition({"H": 2, "O": 1})
        self.aggregate.composition_offset -= offset

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

    def _make_glycan_composition_proxy(self):
        if self.aggregate is not None:
            base = self.aggregate.clone()
        else:
            base = HashableGlycanComposition()
            # Represent the initial amide bond between the peptide
            # and the first glycan. Subsequent glycans do not need
            # further chemical losses because of the dehyration built
            # directly into the Residue abstraction.
            base.composition_offset -= Composition({"H": 2, "O": 1})
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
                gc = HashableGlycanComposition.from_glycan(value.rule._original)
                base += gc
        return GlycanCompositionProxy(base)

    @property
    def glycan_composition(self):
        if self._proxy is None:
            self._proxy = self._make_glycan_composition_proxy()
        return self._proxy
