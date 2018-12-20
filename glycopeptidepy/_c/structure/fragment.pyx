cimport cython

from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_Append, PyList_GetItem, PyList_SetItem, PyList_New
from cpython.dict cimport PyDict_SetItem, PyDict_Keys, PyDict_Values, PyDict_Items
from cpython.int cimport PyInt_AsLong
from cpython.float cimport PyFloat_AsDouble
from cpython.tuple cimport PyTuple_GetItem


from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.structure.base cimport ModificationBase

from glycopeptidepy.structure.modification import (
    Modification, NGlycanCoreGlycosylation, OGlycanCoreGlycosylation,
    GlycosaminoglycanLinkerGlycosylation, ModificationCategory)

_n_glycosylation = NGlycanCoreGlycosylation()
_o_glycosylation = OGlycanCoreGlycosylation()
_gag_linker_glycosylation = GlycosaminoglycanLinkerGlycosylation()
_modification_hexnac = Modification("HexNAc").rule
_modification_xylose = Modification("Xyl").rule


@cython.freelist(100)
cdef class ChemicalShiftBase(object):

    def clone(self):
        return self.__class__(self.name, self.composition.clone())

    def __str__(self):
        return self.name

    def __repr__(self):
        return "%s(name=%r)" % (self.__class__.__name__, self.name)

    def is_loss(self):
        return self.mass < 0



cdef class IonSeriesBase(object):

    def __eq__(self, other):
        try:
            return (self is other) or (self.name == other.name)
        except AttributeError:
            return self.name == other

    def __ne__(self, other):
        return not self == other


@cython.freelist(1000000)
cdef class FragmentBase(object):
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

    cpdef IonSeriesBase get_series(self):
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
        return self.get_series()

    def clone(self):
        raise NotImplementedError()

    cpdef ChemicalShiftBase get_chemical_shift(self):
        return self._chemical_shift

    cpdef set_chemical_shift(self, ChemicalShiftBase chemical_shift):
        if self._chemical_shift is not None:
            self.mass += self._chemical_shift.mass
        self._chemical_shift = chemical_shift
        if chemical_shift is not None:
            self.mass -= chemical_shift.mass

    chemical_shift = property(get_chemical_shift, set_chemical_shift)

    cpdef str get_fragment_name(self):
        parts = [self._name]
        chemical_shift = self.chemical_shift
        if chemical_shift is not None:
            parts.append(str(chemical_shift))
        return ''.join(parts)

    @property
    def name(self):
        if self._name is None:
            self._name = self.get_fragment_name()
        return self._name

    @name.setter
    def name(self, name):
        self._name = name


cdef object ModificationCategory_glycosylation = ModificationCategory.glycosylation


cdef class PeptideFragment(FragmentBase):
    concerned_modifications = set(
        [_n_glycosylation,
         _modification_hexnac,
         _o_glycosylation,
         _gag_linker_glycosylation,
         _modification_xylose])

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

        self.flanking_amino_acids = flanking_amino_acids
        self.position = position

        self.glycosylation = glycosylation
        self.set_chemical_shift(chemical_shift)
        self._update_mass_with_modifications()

        self._name = self.get_fragment_name()

        self._hash = hash(self._name)


    cdef void _update_mass_with_modifications(self):
        cdef:
            list modifications
            ChemicalShiftBase chemical_shift
            size_t i, n
            ModificationBase mod
            object v
            tuple mod_count

        modifications = PyDict_Items(self.modification_dict)
        n = PyList_GET_SIZE(modifications)
        for i in range(n):
            mod_count = <tuple>PyList_GetItem(modifications, i)
            mod = <ModificationBase>PyTuple_GetItem(mod_count, 0)
            v = <object>PyTuple_GetItem(mod_count, 1)
            self.mass += mod.mass * PyInt_AsLong(v)

        chemical_shift = self.get_chemical_shift()
        if chemical_shift is not None:
            self.mass += chemical_shift.mass

    cpdef IonSeriesBase get_series(self):
        return self.kind

    cpdef clone(self):
        return self.__class__(
            self.series, self.position, dict(self.modification_dict),
            self.bare_mass, list(
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

    cpdef str get_fragment_name(self):
        cdef:
            list fragment_name, mod_rule_counts
            tuple mod_rule_count
            ModificationBase mod_rule
            int count
            str name
            ChemicalShiftBase chemical_shift

        fragment_name = []
        fragment_name.append(self.get_series().name)
        fragment_name.append(str(self.position))

        mod_rule_counts = PyDict_Items(self.modification_dict)
        n = PyList_GET_SIZE(mod_rule_counts)

        for i in range(n):
            mod_rule_count = <tuple>PyList_GetItem(mod_rule_counts, i)
            mod_rule = <ModificationBase>PyTuple_GetItem(mod_rule_count, 0)
            if mod_rule.is_a(ModificationCategory_glycosylation):
                count = PyInt_AsLong(<object>PyTuple_GetItem(mod_rule_count, 1))
                if count > 1:
                    fragment_name.extend(
                        ['+', str(count), mod_rule.name])
                elif count == 1:
                    fragment_name.extend(['+', mod_rule.name])
                else:
                    pass

        chemical_shift = self.get_chemical_shift()

        if chemical_shift is not None:
            fragment_name.append(str(chemical_shift))
        name = ''.join(fragment_name)
        return name

    @property
    def is_glycosylated(self):
        if self.glycosylation:
            return True
        else:
            for mod in self.modification_dict:
                if mod.is_a(ModificationCategory_glycosylation):
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
        return ("PeptideFragment(%(type)s %(position)s %(mass)s "
                "%(modification_dict)s %(flanking_amino_acids)s %(chemical_shift)r)") % {
            "type": self.series, "position": self.position, "mass": self.mass,
            "modification_dict": self.modification_dict, "flanking_amino_acids": self.flanking_amino_acids,
            "chemical_shift": self.chemical_shift
        }
