# cython: embedsignature=True

from libc.stdlib cimport malloc, free, realloc
from libc.string cimport strcpy, memcpy

cimport cython
from cpython cimport PyObject
from cpython.ref cimport Py_INCREF
from cpython.object cimport PyObject_Str
from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_Append, PyList_GetItem, PyList_SetItem, PyList_New
from cpython.dict cimport PyDict_SetItem, PyDict_Keys, PyDict_Values, PyDict_Items, PyDict_Next
from cpython.int cimport PyInt_AsLong
from cpython.float cimport PyFloat_AsDouble
from cpython.tuple cimport PyTuple_GetItem

from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.structure.base cimport ModificationBase

from glycopeptidepy._c.compat cimport PyStr_AsUTF8AndSize, PyStr_FromStringAndSize

from glycopeptidepy._c.count_table cimport CountTable, CountTableIterator

from glycopeptidepy.structure.modification import (
    Modification, NGlycanCoreGlycosylation, OGlycanCoreGlycosylation,
    GlycosaminoglycanLinkerGlycosylation, ModificationCategory)

cdef ModificationBase _n_glycosylation = NGlycanCoreGlycosylation()
cdef ModificationBase _o_glycosylation = OGlycanCoreGlycosylation()
cdef ModificationBase _gag_linker_glycosylation = GlycosaminoglycanLinkerGlycosylation()
cdef ModificationBase _modification_hexnac = Modification("HexNAc").rule
cdef ModificationBase _modification_xylose = Modification("Xyl").rule


cdef struct string_cell:
    char* string
    size_t size

DEF ARRAY_SIZE = 2 ** 12
DEF DEFAULT_FRAGMENT_NAME_BUFFER_SIZE = 128

cdef string_cell[ARRAY_SIZE] str_ints
cdef int i
cdef Py_ssize_t z
cdef str pa
cdef char* a
cdef char* atemp

for i in range(ARRAY_SIZE):
    pa = str(i)
    a = <char*>malloc(sizeof(char) * z)
    atemp = PyStr_AsUTF8AndSize(pa, &z)
    strcpy(a, atemp)
    a[z] = "\0"
    str_ints[i].string = a
    str_ints[i].size = z


@cython.freelist(100)
cdef class ChemicalShiftBase(object):

    cpdef ChemicalShiftBase clone(self):
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
        chemical_shift = self.get_chemical_shift()
        if chemical_shift is not None:
            parts.append(chemical_shift.name)
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

    @staticmethod
    cdef PeptideFragment _create(IonSeriesBase kind, int position, CountTable modification_dict, double mass,
                                 list flanking_amino_acids=None, object glycosylation=None,
                                 ChemicalShiftBase chemical_shift=None, CComposition composition=None):
        cdef PeptideFragment self = PeptideFragment.__new__(PeptideFragment)
        self.kind = kind
        self.position = position

        self.bare_mass = mass
        self.mass = mass
        self.modification_dict = modification_dict
        self.composition = composition
        self._chemical_shift = None
        
        self.flanking_amino_acids = flanking_amino_acids
        self.glycosylation = glycosylation
        self.set_chemical_shift(chemical_shift)

        self._update_mass_with_modifications()

        self._name = self.get_fragment_name()
        self._hash = hash(self._name)
        return self

    def __init__(self, kind, position, modification_dict, mass, flanking_amino_acids=None,
                 glycosylation=None, chemical_shift=None, composition=None):
        self.kind = kind

        # The mass value is the bare backbone's mass
        self.bare_mass = PyFloat_AsDouble(mass)
        self.modification_dict = modification_dict
        self.mass = self.bare_mass
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
            ChemicalShiftBase chemical_shift
            ModificationBase mod
            CountTable modification_dict
            CountTableIterator modification_iterator
            PyObject *pkey
            long count
            Py_ssize_t ppos = 0

        modification_dict = self.modification_dict
        modification_iterator = CountTableIterator._create(modification_dict)

        while modification_iterator.has_more():
            ppos = modification_iterator.get_next_value(&pkey, &count)
            if ppos != 0:
                break
            mod = <ModificationBase>pkey
            self.mass += mod.mass * count

        chemical_shift = self.get_chemical_shift()
        if chemical_shift is not None:
            self.mass += chemical_shift.mass

    cpdef IonSeriesBase get_series(self):
        return self.kind

    cpdef clone(self):
        return PeptideFragment._create(
            self.series, self.position, self.modification_dict.copy(),
            self.bare_mass, list(self.flanking_amino_acids),
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
        fragment_name.append(self.get_series().name)
        fragment_name.append(str(self.position))
        return ''.join(fragment_name)

    cpdef str get_fragment_name(self):
        cdef:
            ModificationBase mod_rule
            long count
            ChemicalShiftBase chemical_shift
            CountTable modification_dict
            CountTableIterator modification_iterator
            str name
            PyObject *pkey
            PyObject *pvalue
            Py_ssize_t index, size_ref, buffer_size
            Py_ssize_t ppos = 0

            string_cell int_conv

            char* tmp_buffer
            char[DEFAULT_FRAGMENT_NAME_BUFFER_SIZE] default_name_buffer
            char* name_buffer
            char* oversized_buffer
            bint needs_free

        name_buffer = <char*>default_name_buffer
        buffer_size = DEFAULT_FRAGMENT_NAME_BUFFER_SIZE
        needs_free = False
        index = 0
        tmp_buffer = PyStr_AsUTF8AndSize(self.get_series().name, &size_ref)
        memcpy(&name_buffer[index], tmp_buffer, size_ref)
        index += size_ref

        int_conv = str_ints[self.position]
        memcpy(&name_buffer[index], int_conv.string, int_conv.size)
        index += int_conv.size

        modification_dict = self.modification_dict
        modification_iterator = CountTableIterator._create(modification_dict)
        while modification_iterator.has_more():
            ppos = modification_iterator.get_next_value(&pkey, &count)
            if ppos != 0:
                break
            mod_rule = <ModificationBase>pkey
            if buffer_size - index < 30:
                if needs_free:
                    name_buffer = <char*>realloc(name_buffer, sizeof(char) * buffer_size * 2)
                    buffer_size *= 2
                else:
                    oversized_buffer = <char*>malloc(sizeof(char) * buffer_size * 2)
                    memcpy(oversized_buffer, name_buffer, index)
                    name_buffer = oversized_buffer
                    needs_free = True
                    buffer_size *= 2
            if mod_rule.is_a(ModificationCategory_glycosylation):
                if count > 1:
                    name_buffer[index] = "+"
                    index += 1
                    int_conv = str_ints[count]
                    memcpy(&name_buffer[index], int_conv.string, int_conv.size)
                    index += int_conv.size
                    tmp_buffer = PyStr_AsUTF8AndSize(mod_rule.name, &size_ref)   
                    if size_ref + index > buffer_size:
                        if needs_free:
                            name_buffer = <char*>realloc(name_buffer, sizeof(char) * (buffer_size + size_ref) * 2)
                            buffer_size = (buffer_size + size_ref) * 2
                        else:
                            oversized_buffer = <char*>malloc(sizeof(char) * (buffer_size + size_ref) * 2)
                            memcpy(oversized_buffer, name_buffer, index)
                            name_buffer = oversized_buffer
                            needs_free = True
                            buffer_size = (buffer_size + size_ref) * 2
                    memcpy(&name_buffer[index], tmp_buffer, size_ref)
                    index += size_ref
                elif count == 1:
                    name_buffer[index] = "+"
                    index += 1
                    tmp_buffer = PyStr_AsUTF8AndSize(mod_rule.name, &size_ref)   
                    if size_ref + index > buffer_size:
                        if needs_free:
                            name_buffer = <char*>realloc(name_buffer, sizeof(char) * (buffer_size + size_ref) * 2)
                            buffer_size = (buffer_size + size_ref) * 2
                        else:
                            oversized_buffer = <char*>malloc(sizeof(char) * (buffer_size + size_ref) * 2)
                            memcpy(oversized_buffer, name_buffer, index)
                            name_buffer = oversized_buffer
                            needs_free = True
                            buffer_size = (buffer_size + size_ref) * 2
                    memcpy(&name_buffer[index], tmp_buffer, size_ref)
                    index += size_ref

        chemical_shift = self.get_chemical_shift()
        if chemical_shift is not None:
            tmp_buffer = PyStr_AsUTF8AndSize(chemical_shift.name, &size_ref)
            if size_ref + index > buffer_size:
                if needs_free:
                    name_buffer = <char*>realloc(name_buffer, sizeof(char) * (buffer_size + size_ref))
                    buffer_size = (buffer_size + size_ref)
                else:
                    oversized_buffer = <char*>malloc(sizeof(char) * (buffer_size + size_ref))
                    memcpy(oversized_buffer, name_buffer, index)
                    name_buffer = oversized_buffer
                    needs_free = True
                    buffer_size = (buffer_size + size_ref) 
            memcpy(&name_buffer[index], tmp_buffer, size_ref)
            index += size_ref        
        name = PyStr_FromStringAndSize(&name_buffer[0], index)
        if needs_free:
            free(name_buffer)
        return name

    @property
    def is_glycosylated(self):
        cdef:
            PyObject *pkey
            long value
            Py_ssize_t ppos = 0
            ModificationBase mod
            CountTable modification_dict
            CountTableIterator modification_iterator
        if self.glycosylation:
            return True
        else:
            modification_dict = self.modification_dict
            modification_iterator = CountTableIterator._create(modification_dict)
            while modification_iterator.has_more():
                ppos = modification_iterator.get_next_value(&pkey, &value)
                if ppos != 0:
                    break
                mod = <ModificationBase>pkey
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
