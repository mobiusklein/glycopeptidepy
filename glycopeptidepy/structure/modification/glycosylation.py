import re

from functools import partial

from six import string_types as basestring

from glypy import Substituent, Glycan
from glypy.io import glycoct, iupac, linear_code, wurcs
from glypy.structure.glycan_composition import (
    FrozenMonosaccharideResidue, FrozenGlycanComposition)


from ..composition import Composition
from ..glycan import TypedGlycan, TypedGlycanComposition, GlycosylationType

from .target import ModificationTarget
from .rule import ModificationRule
from .descriptors import ModificationCategory


_hexnac = FrozenMonosaccharideResidue.from_iupac_lite("HexNAc")
_hexose = FrozenMonosaccharideResidue.from_iupac_lite("Hex")
_xylose = FrozenMonosaccharideResidue.from_iupac_lite("Xyl")


hexnac_modification = ModificationRule.from_unimod({
    "title": "HexNAc",
    "composition": {
        "H": 13,
        "C": 8,
        "O": 5,
        "N": 1
    },
    "record_id": 43,
    "mono_mass": 203.079373,
    "full_name": "N-Acetylhexosamine",
    "specificity": [{
        "position": "Anywhere",
        "hidden": True,
        "site": "T",
        "classification": "Other glycosylation",
        "spec_group": 3
    }, {
        "position": "Anywhere",
        "hidden": True,
        "site": "S",
        "classification": "Other glycosylation",
        "spec_group": 2
    }, {
        "position": "Anywhere",
        "hidden": True,
        "site": "N",
        "classification": "N-linked glycosylation",
        "spec_group": 1
    }],
})


xylose_modification = ModificationRule.from_unimod({
    "title": "Xyl",
    "composition": dict(_xylose.total_composition()),
    "mono_mass": _xylose.mass(),
    "full_name": "Xylose",
    "specificity": [{
        "position": "Anywhere",
        "hidden": True,
        "site": "S",
        "classification": "O-linked glycosylation",
        "spec_group": 2
    }],
})


glycan_resolvers = dict()
glycan_resolvers['glycoct'] = partial(glycoct.loads, structure_class=TypedGlycan)
glycan_resolvers['iupac'] = partial(iupac.loads, structure_class=TypedGlycan)
glycan_resolvers['iupac_simple'] = partial(iupac.loads, structure_class=TypedGlycan, dialect='simple')
glycan_resolvers['linear_code'] = partial(linear_code.loads, structure_class=TypedGlycan)
glycan_resolvers['wurcs2'] = partial(wurcs.loads, structure_class=TypedGlycan)
glycan_resolvers['iupaclite'] = TypedGlycanComposition.parse


def parse_glycan(glycan_format, encoded_string):
    try:
        parser = glycan_resolvers[glycan_format]
        return parser(encoded_string)
    except KeyError:
        raise KeyError("Could not resolve glycan parser for %r (%r)" % (glycan_format, encoded_string))


def _Glycosylation_reconstructor(string):
    return Glycosylation.try_parse(string)


@ModificationRule.resolve.register
class Glycosylation(ModificationRule):
    """
    Incubator Idea - Represent occupied glycosylation sites
    using the Modification interface.

    Attributes
    ----------
    categories : list
        Description
    mass : TYPE
        Description
    name : TYPE
        Description
    names : TYPE
        Description
    options : dict
        Description
    parser : TYPE
        Description
    title : TYPE
        Description
    """
    @classmethod
    def _parse(cls, rule_string):
        if rule_string.startswith("#"):
            rule_string = rule_string[1:]
        format_type = "iupaclite"
        glycan_definition = rule_string
        metadata = {}
        if rule_string.startswith(":"):
            match = re.search(r":([^:]*?):(.+)", rule_string, re.DOTALL)
            if match:
                metadata = match.group(1).split(",")
                if '=' not in metadata[0]:
                    format_type = metadata[0]
                    if not format_type:
                        format_type = 'iupaclite'
                    metadata = metadata[1:]
                metadata = dict([token.split("=") for token in metadata])
                glycan_definition = match.group(2)
            else:
                raise ValueError("Cannot recognize glycan format %r" % (rule_string,))
        return parse_glycan(format_type, glycan_definition), format_type, metadata

    @classmethod
    def try_parse(cls, rule_string):
        try:
            glycan, encoding_format, metadata = cls._parse(rule_string)
            return cls(glycan, encoding_format, metadata)
        except Exception:
            return None

    def __init__(self, glycan, encoding_format=None, metadata=None):
        if metadata is None:
            metadata = {}
        if isinstance(glycan, basestring):
            if encoding_format is None:
                glycan, encoding_format, _metadata = self._parse(glycan)
                metadata.update(_metadata)
            else:
                glycan = parse_glycan(encoding_format, glycan)
        else:
            glycan = glycan.clone()

        if isinstance(glycan, Glycan):
            self._is_composition = False
            if encoding_format is None:
                encoding_format = "glycoct"
        else:
            self._is_composition = True
            encoding_format = "iupaclite"

        # Prepare information to encode the rule
        self.encoding_format = encoding_format
        self.metadata = metadata
        self.glycan = glycan
        try:
            self.glycan.glycosylation_type = GlycosylationType[metadata['glycosylation_type']]
        except KeyError:
            pass
        self._original = glycan.clone()
        self.common_name = self._make_string()

        # Only after the common name is set, prepare the patch
        self._patch_dehydration()
        self.mass = self.glycan.mass()
        self.composition = self.glycan.total_composition()

        self._simple_configuration()

    def __reduce__(self):
        return _Glycosylation_reconstructor, (self.name, )

    def _simple_configuration(self):
        self.title = self.common_name
        self.name = self.common_name
        self.unimod_name = self.common_name
        self.categories = [ModificationCategory.glycosylation]
        self.names = {self.common_name, self.title, self.name}
        self.aliases = set()
        self.options = {}
        self._hash = hash(self.name)
        self.neutral_losses = []
        self._n_term_targets = None
        self._c_term_targets = None
        self.targets = set()

    def _make_string(self):
        template = "#:{}:{}"
        metadata = ",".join("%s=%s" % (k, v) for k, v in self.metadata.items())
        if metadata:
            if self.encoding_format != 'iupaclite':
                metadata = "{},{}".format(self.encoding_format, metadata)
        else:
            if self.encoding_format != 'iupaclite':
                metadata = self.encoding_format
        if self.is_composition:
            return template.format(metadata, self.glycan.serialize())
        else:
            return template.format(
                metadata, self.glycan.serialize(self.encoding_format))

    @property
    def glycosylation_type(self):
        return self.glycan.glycosylation_type

    @property
    def is_composition(self):
        return self._is_composition

    @property
    def is_core(self):
        return False

    def _patch_dehydration(self):
        if self.is_composition:
            self.glycan.composition_offset = Composition()
        else:
            self.glycan.root.add_substituent(
                Substituent("aglycone", composition=Composition()),
                position=1, parent_loss=Composition("HO"))

    def clone(self):
        return self.__class__(
            self._original, self.encoding_format, self.metadata)

    def is_tracked_for(self, category):
        return category == ModificationCategory.glycosylation

    def get_fragments(self, *args, **kwargs):
        if self.is_composition:
            raise TypeError("Cannot generate fragments from composition")
        for frag in self.glycan.fragments(*args, **kwargs):
            yield frag

    def total_composition(self):
        return self.glycan.total_composition()


class CoreGlycosylation(Glycosylation):
    @property
    def is_composition(self):
        return True

    @property
    def is_core(self):
        return True

    def get_fragments(self, *args, **kwargs):
        return []

    def _common_init(self):
        self.names = {self.unimod_name, self.title, self.name, self.common_name}
        self.options = {}
        self.aliases = set()
        self.categories = [ModificationCategory.glycosylation]
        self._hash = hash(self.name)
        self.neutral_losses = []
        self._n_term_targets = None
        self._c_term_targets = None
        self.targets = set()

    def clone(self):
        return self.__class__(self.mass, )

    def __reduce__(self):
        return self.__class__, (self.mass, )


class NGlycanCoreGlycosylation(CoreGlycosylation):
    mass_ladder = {k: FrozenGlycanComposition.parse(k).total_composition() - Composition("H2O") for k in {
        "{HexNAc:1}",
        "{HexNAc:2}",
        "{HexNAc:2; Hex:1}",
        "{HexNAc:2; Hex:2}",
        "{HexNAc:2; Hex:3}",
    }}

    def __init__(self, base_mass=_hexnac.mass()):
        self.common_name = "NGlycanCoreGlycosylation"
        self.mass = base_mass
        self.title = "N-Glycan Core Glycosylation"
        self.unimod_name = "N-Glycosylation"
        self.name = self.unimod_name
        self.targets = [(ModificationTarget("N"))]
        self.composition = _hexnac.total_composition().clone()
        self._common_init()

    @property
    def glycosylation_type(self):
        return GlycosylationType.n_linked

    def clone(self):
        return self.__class__(self.mass)


class MucinOGlycanCoreGlycosylation(CoreGlycosylation):
    mass_ladder = {k: FrozenGlycanComposition.parse(k).total_composition() - Composition("H2O") for k in {
        "{HexNAc:1}",
        "{HexNAc:1; Hex:1}",
    }}

    def __init__(self, base_mass=_hexnac.mass()):
        self.common_name = "MucinOGlycanCoreGlycosylation"
        self.mass = base_mass
        self.title = "Mucin O-Glycan Core Glycosylation"
        self.unimod_name = "O-Glycosylation"
        self.name = self.unimod_name
        self.targets = [ModificationTarget("S"), ModificationTarget("T")]
        self.composition = _hexnac.total_composition().clone()
        self._common_init()

    @property
    def glycosylation_type(self):
        return GlycosylationType.o_linked

    def clone(self):
        return self.__class__(self.mass)


OGlycanCoreGlycosylation = MucinOGlycanCoreGlycosylation


class GlycosaminoglycanLinkerGlycosylation(CoreGlycosylation):
    mass_ladder = {k: FrozenGlycanComposition.parse(k).total_composition() - Composition("H2O") for k in {
        "{Xyl:1}",
        "{Xyl:1; Hex:1}",
        "{Xyl:1; Hex:2}",
        "{Xyl:1; Hex:2; HexA:1}",
    }}

    def __init__(self, base_mass=_xylose.mass()):
        self.common_name = "GlycosaminoglycanLinkerGlycosylation"
        self.mass = base_mass
        self.title = "Glycosaminoglycan Linker Glycosylation"
        self.unimod_name = "GAG-Linker"
        self.name = self.unimod_name
        self.targets = [ModificationTarget("S")]
        self.composition = _hexnac.total_composition().clone()
        self._common_init()

    @property
    def glycosylation_type(self):
        return GlycosylationType.glycosaminoglycan

    def clone(self):
        return self.__class__(self.mass)


class OGlcNAcylation(CoreGlycosylation):
    mass_ladder = {
        "{GlcNAc:1}": _hexnac.total_composition()
    }

    def __init__(self, base_mass=_hexnac.mass()):
        self.common_name = "O-GlcNAc"
        self.name = self.common_name
        self.mass = base_mass
        self.title = self.common_name
        self.unimod_name = self.common_name
        self.targets = [ModificationTarget("S"), ModificationTarget("T")]
        self.composition = _hexnac.total_composition().clone()
        self._common_init()

    def get_fragments(self, *args, **kwargs):
        for label_loss in self.mass_ladder.items():
            yield label_loss

    def clone(self):
        return self.__class__(self.mass)


def _GlycanFragment_reconstructor():
    return GlycanFragment.__new__(GlycanFragment)


class GlycanFragment(Glycosylation):
    @property
    def is_composition(self):
        return True

    @property
    def is_core(self):
        return False

    def get_fragments(self, *args, **kwargs):
        return []

    def __init__(self, fragment):
        self.mass = fragment.mass
        self.name = self.unimod_name = self.title = self.common_name = fragment.name
        self.names = {self.name, }
        self.composition = fragment.composition

        self.categories = [ModificationCategory.glycosylation]
        self.aliases = set()
        self.targets = set()
        self.options = {}
        self._hash = hash(self.name)

    def __reduce__(self):
        return _GlycanFragment_reconstructor, (), self.__getstate__()

    def __getstate__(self):
        state = {
            'name': self.name,
            'composition': self.composition,
            'mass': self.mass,
        }
        return state

    def __setstate__(self, state):
        self.mass = state['mass']
        self.name = self.unimod_name = self.title = self.common_name = state['name']
        self.names = {self.name, }
        self.composition = state['composition']

        self.categories = [ModificationCategory.glycosylation]
        self.aliases = set()
        self.targets = set()
        self.options = {}
        self._hash = hash(self.name)
