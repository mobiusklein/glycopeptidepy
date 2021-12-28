'''Represents fragments of peptides and glycopeptides.
'''
import re
from collections import defaultdict
from six import add_metaclass
from .modification import (
    Modification, NGlycanCoreGlycosylation, OGlycanCoreGlycosylation,
    GlycosaminoglycanLinkerGlycosylation, ModificationCategory)
from .glycan import HashableGlycanComposition
from .composition import Composition, formula

_n_glycosylation = NGlycanCoreGlycosylation()
_o_glycosylation = OGlycanCoreGlycosylation()
_gag_linker_glycosylation = GlycosaminoglycanLinkerGlycosylation()
_modification_hexnac = Modification("HexNAc").rule
_modification_xylose = Modification("Xyl").rule

fragment_pairing = {
    "a": "x",
    "b": "y",
    "c": "z",
    "x": "a",
    "y": "b",
    "z": "c",
}

fragment_shift_composition = {
    'a': -Composition("CHO"),
    'b': -Composition({"H": 1}),
    'c': Composition("NH2"),
    'x': Composition("CO2") - Composition("H"),
    'y': Composition('H'),
    'z': (-Composition("NH")),
    'zp': (-Composition("NH2")),
}

fragment_shift = {
    k: v.mass for k, v in fragment_shift_composition.items()
}


fragment_direction = {
    "a": 1,
    "b": 1,
    "c": 1,
    "x": -1,
    "y": -1,
    "z": -1,
    "zp": -1,
}

generic_chemical_shifts_composition = {
    "-NH3": -Composition("NH3"),
    "-H2O": -Composition("H2O"),
    # "-NH3-NH3": -Composition("(NH3)2"),
    "-H2O-H2O": -Composition("(H2O)2"),
    "-NH3-H2O": -Composition("NH3H2O")
}


def format_negative_composition(composition):
    '''Format a chemical composition as a formula with negative
    element counts.
    '''
    return "-%s" % formula({k: -v for k, v in composition.items()})

# pylint: disable=access-member-before-definition,assigning-non-slot,invalid-unary-operand-type


try:
    from glycopeptidepy._c.structure.fragment import ChemicalShiftBase as _ChemicalShiftBase
except ImportError:
    _ChemicalShiftBase = object


class ChemicalShift(_ChemicalShiftBase):
    '''Represents a chemical shift, a gain or loss of chemical composition.

    Attributes
    ----------
    name: str
        The name of the modification
    composition: Composition
        The chemical composition difference
    mass: float
        The mass of the chemical composition shift
    '''
    def __init__(self, name, composition=None):
        if composition is None:
            composition = generic_chemical_shifts_composition[name]
        self.name = name
        self.composition = composition
        self.mass = composition.mass

    def clone(self):
        '''Create a copy of this object

        Returns
        -------
        :class:`ChemicalShift`
        '''
        return self.__class__(self.name, self.composition.clone())

    def __str__(self):
        return self.name

    def __repr__(self):
        return "%s(name=%r)" % (self.__class__.__name__, self.name)

    def is_loss(self):
        '''Returns :const:`True` if the chemical shift results in a net loss in mass

        Returns
        -------
        :class:`bool`
        '''
        return self.mass < 0

    AllLosses = object()


NeutralLoss = ChemicalShift


class FragmentBase(object):
    """Base class for all Fragment types. Defines basic
    name generation and neutral loss handling functions.

    Attributes
    ----------
    chemical_shift : ChemicalShift
        The ChemicalShift associated with this fragment, or None.
        If a ChemicalShift, its composition and mass are subtracted
        from this object's composition and mass attributes.
    name: str
        The human readable description of this fragment
    series: IonSeries
        The ion ladder this fragment is derived from

    """

    __slots__ = ("_chemical_shift", "_name", "_hash")

    def get_series(self):
        '''Return the :class:`IonSeries` this fragment belongs to.

        Returns
        -------
        :class:`IonSeries`
        '''
        raise NotImplementedError()

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(self.name)
        return self._hash

    def __eq__(self, other):
        try:
            return self.name == other.name and abs(self.mass - other.mass) < 1e-5
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self == other

    @property
    def series(self):
        '''The :class:`IonSeries` this fragment is drawn from.

        Returns
        -------
        :class:`IonSeries`
        '''
        return self.get_series()

    def clone(self):
        '''Create a copy of this object

        Returns
        -------
        :class:`FragmentBase`
        '''
        raise NotImplementedError()

    def get_chemical_shift(self):
        '''Returns the chemical shift associated with the fragment.

        Returns
        -------
        :class:`ChemicalShift`
        '''
        return self._chemical_shift

    def set_chemical_shift(self, chemical_shift):
        '''Sets the chemical shift associated with the fragment, updating
        the :attr:`mass` and :attr:`chemical_shift`.
        '''
        if self._chemical_shift is not None:
            self.mass -= self._chemical_shift.mass
        self._chemical_shift = chemical_shift
        if chemical_shift is not None:
            self.mass += chemical_shift.mass

    def get_fragment_name(self):
        parts = [self._name]
        chemical_shift = self.chemical_shift
        if chemical_shift is not None:
            parts.append(str(chemical_shift))
        return ''.join(parts)

    def base_name(self):
        return self._name

    @property
    def name(self):
        if self._name is None:
            self._name = self.get_fragment_name()
        return self._name

    @name.setter
    def set_name(self, name):
        self._name = name

    chemical_shift = property(get_chemical_shift, set_chemical_shift)

    def generate_chemical_shifts(self, losses=ChemicalShift.AllLosses):
        if losses is ChemicalShift.AllLosses:
            losses = generic_chemical_shifts_composition.keys()

        for loss in losses:
            frag = self.clone()
            if loss is not None:
                frag.chemical_shift = ChemicalShift(loss)
            yield frag


class PeptideFragment(FragmentBase):
    '''Represents a peptide backbone fragment, such as a peptide b or y ion.

    Attributes
    ----------
    series: :class:`IonSeries`
        The ion series this fragment is from
    bare_mass: float
        The unmodified mass of this fragment
    mass: float
        The modified, correct mass of this fragment
    position: int
        The peptide backbone position index this fragment came from, relative to the
        starting terminal
    modification_dict: :class:`~.ModificationIndex`
        A count of the modifications on this fragment
    flanking_amino_acids: list
        The adjacent amino acids of the amide bond that broke to form this fragment
    glycosylation: dict
        A position to glycan mapping for glycan structures
    chemical_shift: :class:`ChemicalShift`
        An associated chemical shift
    composition: :class:`~.Composition`
        The elemental composition of this fragment
    '''
    concerned_modifications = set(
        [_n_glycosylation,
         _modification_hexnac,
         _o_glycosylation,
         _gag_linker_glycosylation,
         _modification_xylose])

    __slots__ = ("kind", "position", "modification_dict", "bare_mass",
                 "flanking_amino_acids", "glycosylation",
                 "_chemical_shift", "composition", 'mass')

    def __init__(self, kind, position, modification_dict, mass,
                 flanking_amino_acids=None, glycosylation=None, chemical_shift=None,
                 composition=None):
        self.kind = kind

        # The mass value is the bare backbone's mass
        self.bare_mass = mass
        self.modification_dict = modification_dict
        self.mass = mass
        self.composition = composition
        self._chemical_shift = None
        self._name = None
        self._hash = None

        self.flanking_amino_acids = flanking_amino_acids
        self.position = position

        self.glycosylation = glycosylation
        self.chemical_shift = chemical_shift

        self._update_mass_with_modifications()

    def _update_mass_with_modifications(self):
        for key, value in self.modification_dict.items():
            self.mass += (key).mass * value
        chemical_shift = self.get_chemical_shift()
        if chemical_shift is not None:
            self.mass += chemical_shift.mass

    def get_series(self):
        return self.kind

    def clone(self):
        return self.__class__(
            self.series, self.position, dict(self.modification_dict),
            self.bare_mass, tuple(
                self.flanking_amino_acids),
            self.glycosylation.clone() if self.glycosylation is not None else None,
            self._chemical_shift.clone() if self._chemical_shift is not None else None,
            self.composition.clone())

    def total_composition(self):
        composition = self.composition.clone()
        chemical_shift = self.chemical_shift
        if chemical_shift is not None:
            composition += chemical_shift.composition
        return composition

    def __reduce__(self):
        return self.__class__, (
            self.kind, self.position, self.modification_dict, self.bare_mass,
            self.flanking_amino_acids, self.glycosylation,
            self.chemical_shift, self.composition)

    def base_name(self):
        """Simply return string like b2, y3 with no modification information."""
        fragment_name = []
        fragment_name.append(str(self.series))
        fragment_name.append(str(self.position))
        return ''.join(fragment_name)

    def get_fragment_name(self):
        """Build a complete fragment name listing any glycosylation modifications
        """
        fragment_name = []
        fragment_name.append(str(self.series))
        fragment_name.append(str(self.position))

        # Only concerned modifications are reported.
        for mod_rule, count in self.modification_dict.items():
            if mod_rule in self.concerned_modifications or mod_rule.is_a(ModificationCategory.glycosylation):
                if count > 1:
                    fragment_name.extend(
                        ['+', str(count), str(mod_rule.name)])
                elif count == 1:
                    fragment_name.extend(['+', str(mod_rule.name)])
                else:
                    pass

        if self.chemical_shift is not None:
            fragment_name.append(str(self.chemical_shift))
        name = ''.join(fragment_name)
        return name

    @property
    def is_glycosylated(self):
        if self.glycosylation:
            return True
        else:
            for mod in self.modification_dict:
                if mod in self.concerned_modifications or mod.is_a(ModificationCategory.glycosylation):
                    return True
        return False

    @property
    def glycosylation_size(self):
        if self.glycosylation is not None:
            raise NotImplementedError()
        size = 0
        size += self.modification_dict.get(_modification_hexnac, 0)
        size += self.modification_dict.get(_modification_xylose, 0)
        return size

    def __repr__(self):
        template = ("PeptideFragment({self.series}, {self.position}, "
                    "{self.mass}, {modification_dict}, {self.flanking_amino_acids},"
                    " {self.chemical_shift})")
        return template.format(self=self, modification_dict=dict(self.modification_dict))



try:
    _PeptideFragment = PeptideFragment
    from glycopeptidepy._c.structure.fragment import PeptideFragment as _CPeptideFragment
    PeptideFragment = _CPeptideFragment
except ImportError:
    pass


class SimpleFragment(FragmentBase):
    __slots__ = ["name", "mass", "kind", "composition",
                 "_chemical_shift", "is_glycosylated"]

    def __init__(self, name, mass, kind, composition, chemical_shift=None, is_glycosylated=False):
        self._name = None
        self._hash = None
        self._chemical_shift = None
        self.name = name
        self.mass = mass
        self.kind = kind
        self.chemical_shift = chemical_shift
        self.composition = composition
        self.is_glycosylated = is_glycosylated

    def clone(self):
        corrected_mass = self.mass
        if self._chemical_shift is not None:
            corrected_mass -= self.chemical_shift.mass
        return self.__class__(self.name, corrected_mass, self.kind,
                              self.composition,
                              self._chemical_shift.clone() if self._chemical_shift is not None else None,
                              self.is_glycosylated)

    def __reduce__(self):
        corrected_mass = self.mass
        if self._chemical_shift is not None:
            corrected_mass -= self.chemical_shift.mass
        return self.__class__, (self.name, corrected_mass, self.kind, self.composition,
                                self.chemical_shift, self.is_glycosylated)

    def __repr__(self):
        return ("{self.__class__.__name__}(name={self.name}, "
                "mass={self.mass:.04f}, series={self.kind})").format(self=self)

    def get_series(self):
        return self.kind


try:
    _SimpleFragment = SimpleFragment
    from glycopeptidepy._c.structure.fragment import SimpleFragment as _CSimpleFragment
    SimpleFragment = _CSimpleFragment
except ImportError:
    pass


try:
    from glycopeptidepy._c.structure.fragment import _NameTree
except ImportError:

    class _NameTree(object):
        def __init__(self, name=None, children=None):
            if children is None:
                children = defaultdict(_NameTree)
            self.name = name
            self.children = children

        def __getitem__(self, key):
            return self.children[key]

        def traverse(self, parts):
            node = self
            for kv in parts:
                node = node[kv]
            return node


class StubFragment(SimpleFragment):
    __slots__ = ['glycosylation', 'is_extended']

    _name_cache = _NameTree()

    def __init__(self, name, mass, kind, composition, chemical_shift=None, is_glycosylated=False,
                 glycosylation=None, is_extended=False):
        if glycosylation is None:
            glycosylation = HashableGlycanComposition()
        super(StubFragment, self).__init__(name, mass, kind, composition, chemical_shift, is_glycosylated)
        self.glycosylation = glycosylation
        self.is_extended = is_extended

    def clone(self):
        dup = super(StubFragment, self).clone()
        dup.glycosylation = self.glycosylation
        dup.is_extended = self.is_extended
        return dup

    @property
    def glycosylation_size(self):
        return sum(self.glycosylation.values())

    @classmethod
    def build_name_from_composition(cls, glycan_composition):
        name = 'peptide'
        root = cls._name_cache
        parts = list(glycan_composition.items())
        # traverse the tree to find the node referring to
        # this combination of monosaccharides and counts
        node = root.traverse(parts)
        extended_key = node.name
        if extended_key is None:
            parts.sort(key=lambda x: x[0].mass())
            node.name = extended_key = ''.join("%s%d" % kv for kv in parts)
        if extended_key:
            name = "%s+%s" % (name, extended_key)
        return name

    def __reduce__(self):
        proto = list(super(StubFragment, self).__reduce__())
        proto[1] = proto[1] + (self.glycosylation, self.is_extended)
        return tuple(proto)

    def base_name(self):
        if self.chemical_shift is None:
            return self._name
        else:
            return self.build_name_from_composition(self.glycosylation)


try:
    from glycopeptidepy._c.structure.fragment import StubFragment
except ImportError:
    pass


class MemoizedIonSeriesMetaclass(type):

    def __call__(cls, name=None, *args, **kwargs):
        if not hasattr(cls, "_cache"):
            cls._cache = dict()
        try:
            if name is not None:
                return cls._cache[name]
            else:
                raise Exception("Must provide a name parameter")
        except KeyError:
            if name is not None:
                inst = type.__call__(cls, name=name, *args, **kwargs)
                cls._cache[inst.name] = inst
                return inst
            else:
                raise KeyError("Cannot find an IonSeries for %r" % (name))


try:
    from glycopeptidepy._c.structure.fragment import IonSeriesBase as _IonSeriesBase
except ImportError:
    class _IonSeriesBase(object):
        __slots__ = ('name', 'direction', 'includes_peptide', 'mass_shift', 'regex', 'composition_shift',
                     '_hash', 'int_code')

        def __eq__(self, other):
            try:
                return self.name == other.name
            except AttributeError:
                return self.name == other

        def __ne__(self, other):
            try:
                return self.name != other.name
            except AttributeError:
                return self.name != other

        def __hash__(self):
            return self._hash


@add_metaclass(MemoizedIonSeriesMetaclass)
class IonSeries(_IonSeriesBase):
    '''An ion series type, reflecting a terminal type and a bond cleavage type.

    Attributes
    ----------
    name: str
        The name of the ion series, like 'b' or 'y'
    direction: int
        Whether the fragments from this ion series arise from the N-terminal (1),
        C-terminal (-1), or neither (0).
    includes_peptide: bool
        Whether the fragment includes the peptide backbone or not.
    mass_shift: float
        The mass delta applied to the raw components of the fragment.
    composition_shift: :class:`~.Composition`
        The elemental composition delta applied to the raw components of the fragment.
    '''

    @classmethod
    def get(cls, name):
        '''Get a particular cached :class:`IonSeries` by name

        Returns
        -------
        :class:`IonSeries`
        '''
        return cls(name)

    def __init__(self, name, direction=None, includes_peptide=True, mass_shift=None, regex=None,
                 composition_shift=None, int_code=-1):
        if direction is None:
            if name in fragment_direction:
                direction = fragment_direction[name]
            else:
                direction = 0
        if mass_shift is None:
            if name in fragment_shift:
                mass_shift = fragment_shift[name]
            else:
                mass_shift = 0.
        self.name = name
        self.direction = direction
        self.includes_peptide = includes_peptide
        self.mass_shift = mass_shift
        self.regex = re.compile(regex) if regex is not None else regex
        if composition_shift is None:
            self.composition_shift = fragment_shift_composition.get(self.name, Composition())
        else:
            self.composition_shift = composition_shift
        self._hash = hash(self.name)
        self.int_code = int_code

    def __reduce__(self):
        return self.__class__, (self.name, self.direction, self.includes_peptide, self.mass_shift,
                                self.regex, self.composition_shift, self.int_code)

    def __repr__(self):
        template = ("{self.__class__.__name__}({self.name}, "
                    "direction={self.direction}, mass_shift={self.mass_shift})")
        return template.format(self=self)

    def __str__(self):
        return str(self.name)

    def is_member(self, key):
        '''Test if ``key`` matches this series name.

        Parameters
        ----------
        key: str
            The name of a fragment to test for membership

        Returns
        -------
        bool
        '''
        if self.regex is None:
            return key.startswith(self.name)
        else:
            return self.regex.search(key)

    __call__ = is_member

    def __radd__(self, other):
        return other + self.mass_shift

    def __rsub__(self, other):
        return other - self.mass_shift

    def __neg__(self):
        return -self.mass_shift

    def __pos__(self):
        return self.mass_shift

    def __add__(self, other):
        return self.mass_shift + other

    def __sub__(self, other):
        return self.mass_shift - other


IonSeries.b = IonSeries("b")
IonSeries.y = IonSeries("y")
IonSeries.c = IonSeries("c")
IonSeries.z = IonSeries("z")
IonSeries.zp = IonSeries("zp")
IonSeries.z1 = IonSeries.zp

IonSeries.precursor = IonSeries("precursor")
IonSeries.oxonium_ion = IonSeries("oxonium_ion", includes_peptide=False)
IonSeries.stub_glycopeptide = IonSeries("stub_glycopeptide")
IonSeries.other = IonSeries('other')

for i, series in enumerate((IonSeries.b, IonSeries.y, IonSeries.c, IonSeries.z, IonSeries.zp, IonSeries.precursor,
                            IonSeries.oxonium_ion, IonSeries.stub_glycopeptide), 1):
    series.int_code = i
