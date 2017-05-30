from glycopeptidepy.utils.collectiontools import decoratordict

from glypy.utils import Enum
from .modification import AnonymousModificationRule

import glypy
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

    def __init__(self, glycosylation_type, *args, **kwargs):
        self.glycosylation_type = GlycosylationType[glycosylation_type]
        super(TypedGlycanComposition).__init__(self, *args, **kwargs)

    def as_modification(self):
        return AnonymousModificationRule("Glycan[{0}]".format(';'.join(
            self.serialize(),
            self.mass,
            self.glycosylation_type.name)))()

    def is_type(self, glycosylation_type):
        return self.glycosylation_type is GlycosylationType[glycosylation_type]


class GlycanCompositionProxy(object):
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

    def __repr__(self):
        return repr(self.obj)

    def __str__(self):
        return str(self.obj)

    def __hash__(self):
        return hash(self.obj)

    def __eq__(self, other):
        return self.obj == other

    def __ne__(self, other):
        return not (self.obj == other)


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

    def as_modification(self):
        return AnonymousModificationRule("Glycan[{0}]".format(';'.join(self.name, self.mass)))()


class GlycosylationSite(object):

    attachment_loss = Composition("H2O")

    def __init__(self, parent, position, glycosylation=None, site_types=None):
        if site_types is None:
            site_types = [GlycosylationType.n_linked]
        self.parent = parent
        self.position = position
        self.glycosylation = glycosylation
        self.site_types = [GlycosylationType[site_type] for site_type in site_types]

    @property
    def mass(self):
        if self.glycosylation is not None:
            return self.glycosylation.mass() - self.attachment_loss.mass
        else:
            return 0

    def _set_glycosylation(self, glycan):
        self._glycosylation = glycan

    def _get_glycosylation(self):
        return self._glycosylation

    glycosylation = property(_get_glycosylation, _set_glycosylation)

    def add_site_type(self, site_type):
        self.site_types.append(site_type)

    def total_composition(self):
        if self._glycosylation is not None:
            return self._glycosylation.total_composition() - self.attachment_loss
        return Composition()

    def __repr__(self):
        rep = 'GlycosylationSite(%d, %s, %s)' % (self.position, self.glycosylation, ', '.join(
            site_type.name for site_type in self.site_types))
        return rep

    def __eq__(self, other):
        try:
            return (self.position == other.position) and (self.glycosylation == other.glycosylation)
        except:
            return False

    def __ne__(self, other):
        try:
            return (self.position != other.position) and (self.glycosylation != other.glycosylation)
        except:
            return True

    def __hash__(self):
        return hash(self.position)


class MassableValueDict(dict):
    def total_composition(self):
        total = Composition()
        for key, value in self.items():
            total += value.total_composition()
        return total

    def mass(self):
        total = 0.
        for key, value in self.items():
            total += value.mass
        return total


class GlycosylationManager(MassableValueDict):
    def __init__(self, parent, *glycosites):
        self.parent = parent
        for glycosite in glycosites:
            self[glycosite.position] = glycosite

    def glycosylation_sites(self, glycosylation_type=None):
        if glycosylation_type is None:
            return list(self.keys())
        else:
            glycosylation_type = GlycosylationType[glycosylation_type]
            results = []
            for position, site in self.items():
                if glycosylation_type in site.site_types:
                    results.append(position)
            return results

    def add_site(self, position, glycosylation=None, site_types=(GlycosylationType.n_linked,)):
        site = GlycosylationSite(self.parent, position, glycosylation, site_types)
        self[position] = site

    @classmethod
    def extract_sites(cls, peptide):
        site_map = {}
        for glycosylation_type, site_fn in glycosylation_site_detectors.items():
            for site in site_fn(peptide):
                if site in site_map:
                    obj = site_map[site]
                    obj.add_site_type(glycosylation_type)
                else:
                    site_map[site] = GlycosylationSite(peptide, site, site_types=[glycosylation_type])
        return cls(peptide, *site_map.values())
