# cython: embedsignature=True

from libc.stdlib cimport malloc, free, realloc
from libc.string cimport strcpy, memcpy

cimport cython
from cpython cimport PyObject
from cpython.ref cimport Py_INCREF
from cpython.object cimport PyObject_Str
from cpython cimport PyErr_SetString
from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_Append, PyList_GetItem, PyList_SetItem, PyList_New
from cpython.sequence cimport PySequence_Fast, PySequence_Fast_GET_ITEM, PySequence_Fast_GET_SIZE
from cpython.dict cimport PyDict_SetItem, PyDict_Keys, PyDict_Values, PyDict_Items, PyDict_Next, PyDict_GetItem
from cpython.int cimport PyInt_AsLong, PyInt_FromLong
from cpython.float cimport PyFloat_AsDouble
from cpython.tuple cimport PyTuple_GetItem

from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.structure.base cimport ModificationBase

from glycopeptidepy._c.compat cimport PyStr_AsUTF8AndSize, PyStr_FromStringAndSize

from glycopeptidepy._c.count_table cimport CountTable, CountTableIterator, count_table_items

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
cdef dict str_int_overflow_intern_cache = {}

cdef void init_str_ints():
    cdef int i
    cdef Py_ssize_t z
    cdef str pa
    cdef char* a
    cdef char* atemp

    for i in range(ARRAY_SIZE):
        pa = str(i)
        atemp = PyStr_AsUTF8AndSize(pa, &z)
        a = <char*>malloc(sizeof(char) * z)
        strcpy(a, atemp)
        a[z] = "\0"
        str_ints[i].string = a
        str_ints[i].size = z


init_str_ints()


cdef int get_str_int(int i, string_cell* out):
    '''Get a `string_cell` for a given integer.

    For integers < 2 ** 12, go through the `str_ints` array
    fast path, otherwise use `str_int_overflow_intern_cache`.

    `str_int_overflow_intern_cache` keeps a mapping from PyInt to
    an immortalized PyStr whose internal buffer is used to serve
    the `char*` for `string_cell` instances larger than  2 ** 12 - 1.
    '''
    cdef:
        PyObject* ptemp
        object py_i
        str pa
        Py_ssize_t z
    if i < ARRAY_SIZE:
        out[0] = str_ints[i]
        return 0

    py_i = PyInt_FromLong(i)
    ptemp = PyDict_GetItem(str_int_overflow_intern_cache, py_i)
    if ptemp == NULL:
        pa = str(py_i)
        Py_INCREF(pa)
        PyDict_SetItem(str_int_overflow_intern_cache, py_i, pa)
    else:
        pa = <str>ptemp
    atemp = PyStr_AsUTF8AndSize(pa, &z)
    out.string = atemp
    out.size = z
    return 0


def _debug_get_str_int(int i):
    cdef:
        string_cell out
        int result
    result = get_str_int(i, &out)
    print(result)
    print(out.string)
    print(out.size)

def _get_string_cell(int i):
    cdef string_cell cell = str_ints[i]
    print(cell.size)
    print(cell.string)


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
        cdef:
            IonSeriesBase other_tp
            object temp
        if not isinstance(self, IonSeriesBase):
            temp = self
            self = <IonSeriesBase>other
            other = temp
        if isinstance(other, IonSeriesBase):
            other_tp = <IonSeriesBase>other
            return self is other_tp or self.name == other_tp.name
        else:
            return self.name == other

    def __ne__(self, other):
        cdef:
            IonSeriesBase other_tp
            object temp
        if not isinstance(self, IonSeriesBase):
            temp = self
            self = <IonSeriesBase>other
            other = temp
        if isinstance(other, IonSeriesBase):
            other_tp = <IonSeriesBase>other
            return self.name != other_tp.name
        else:
            return self.name != other

    def __hash__(self):
        return self._hash


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
        '''Return the :class:`IonSeries` this fragment belongs to.

        Returns
        -------
        :class:`IonSeries`
        '''
        raise NotImplementedError()

    def __hash__(self):
        if self._hash == -1:
            if self._name is None:
                self._update_hash_name()
            else:
                self._hash = hash(self._name)
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

    cpdef ChemicalShiftBase get_chemical_shift(self):
        '''Returns the chemical shift associated with the fragment.

        Returns
        -------
        :class:`ChemicalShift`
        '''
        return self._chemical_shift

    cpdef set_chemical_shift(self, ChemicalShiftBase chemical_shift):
        '''Sets the chemical shift associated with the fragment, updating
        the :attr:`mass` and :attr:`chemical_shift`.
        '''
        if self._chemical_shift is None and chemical_shift is None:
            return
        if self._chemical_shift is not None:
            self.mass -= self._chemical_shift.mass
        self._chemical_shift = chemical_shift
        if chemical_shift is not None:
            self.mass += chemical_shift.mass
        self._name = self.get_fragment_name()
        self._hash = hash(self._name)

    cdef void _update_hash_name(self):
        self._name = self.get_fragment_name()
        self._hash = hash(self._name)

    chemical_shift = property(get_chemical_shift, set_chemical_shift)

    cpdef str get_fragment_name(self):
        parts = [self._name]
        chemical_shift = self.get_chemical_shift()
        if chemical_shift is not None:
            parts.append(chemical_shift.name)
        return ''.join(parts)

    cpdef str base_name(self):
        """Simply return string like b2, y3 with no modification information."""
        return self._name

    cdef str get_name(self):
        if self._name is None:
            self._name = self.get_fragment_name()
        return self._name

    @property
    def name(self):
        return self.get_name()

    @name.setter
    def name(self, name):
        self._name = name


cdef object ModificationCategory_glycosylation = ModificationCategory.glycosylation


cdef class PeptideFragment(FragmentBase):
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
    @staticmethod
    cdef PeptideFragment _create(IonSeriesBase kind, int position, CountTable modification_dict, double mass,
                                 list flanking_amino_acids=None, dict glycosylation=None,
                                 ChemicalShiftBase chemical_shift=None, CComposition composition=None,
                                 double* delta_mass=NULL):
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
        if chemical_shift is not None:
            self.set_chemical_shift(chemical_shift)

        if delta_mass == NULL:
            self._update_mass_with_modifications()
        else:
            self.mass += delta_mass[0]
        # if delta_mass != NULL:
        #     assert abs((self.mass - self.bare_mass) - delta_mass[0]) < 1e-3, ((self.mass - self.bare_mass), delta_mass[0])
        self._name = None
        self._hash = -1
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

        self._name = None
        self._hash = -1

    cdef void _update_mass_with_modifications(self):
        cdef:
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

    cpdef IonSeriesBase get_series(self):
        return self.kind

    cpdef clone(self):
        return PeptideFragment._create(
            self.series, self.position, self.modification_dict.copy(),
            self.bare_mass, list(self.flanking_amino_acids),
            self.glycosylation.copy() if self.glycosylation is not None else None,
            self._chemical_shift.clone() if self._chemical_shift is not None else None,
            self.composition.clone() if self.composition is not None else None)

    def total_composition(self):
        if self.composition is None:
            raise ValueError("Composition data is missing")
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

    cpdef str base_name(self):
        """Simply return string like b2, y3 with no modification information."""
        fragment_name = []
        fragment_name.append(self.get_series().name)
        fragment_name.append(str(self.position))
        return ''.join(fragment_name)

    cpdef str get_fragment_name(self):
        """Build a complete fragment name listing any glycosylation modifications
        """
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

        get_str_int(self.position, &int_conv)
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

    cdef bint _is_glycosylated(self):
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
    def is_glycosylated(self):
        return self._is_glycosylated()

    @property
    def glycosylation_size(self):
        if self.glycosylation is not None:
            raise NotImplementedError()
        size = 0
        size += self.modification_dict.get(_modification_hexnac, 0)
        size += self.modification_dict.get(_modification_xylose, 0)
        return size

    cdef long get_glycosylation_size(self) except -1:
        cdef:
            long size
        if self.glycosylation is not None:
            PyErr_SetString(NotImplementedError, "Method not implemented for glycan collections yet")
            return -1
        size = 0
        size += self.modification_dict.getitem(_modification_hexnac)
        size += self.modification_dict.getitem(_modification_xylose)
        return size

    def __repr__(self):
        template = ("PeptideFragment({self.series}, {self.position}, "
                    "{self.mass}, {modification_dict}, {self.flanking_amino_acids},"
                    " {self.chemical_shift})")
        return template.format(self=self, modification_dict=dict(self.modification_dict))


cdef class SimpleFragment(FragmentBase):

    def __init__(self, name, mass, kind, composition, chemical_shift=None, is_glycosylated=False):
        self._name = name
        self._hash = hash(self._name)
        self._chemical_shift = None

        self.mass = mass
        self.kind = kind
        self.composition = composition

        self.set_chemical_shift(chemical_shift)
        self.is_glycosylated = is_glycosylated

    cpdef clone(self):
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

    cpdef IonSeriesBase get_series(self):
        return self.kind


@cython.final
cdef class _NameTree(object):
    @staticmethod
    cdef _NameTree _create():
        cdef _NameTree self = _NameTree.__new__(_NameTree)
        self.name = None
        self.value = None
        self.children = dict()
        return self

    def __init__(self, name=None, children=None):
        if children is None:
            children = dict()
        self.name = name
        self.value = None
        self.children = children

    def __getitem__(self, key):
        return self.get(key)

    cpdef _NameTree get(self, key):
        cdef:
            PyObject* presult
            _NameTree result
        presult = PyDict_GetItem(self.children, key)
        if presult == NULL:
            result = _NameTree._create()
            PyDict_SetItem(self.children, key, result)
            return result
        else:
            result = <_NameTree>presult
            return result

    cpdef traverse(self, list_or_iterable parts):
        cdef:
            _NameTree node
            size_t i, n
            object parts_fast

        node = self
        if list_or_iterable is list:
            n = PyList_GET_SIZE(parts)
            for i in range(n):
                kv = <object>PyList_GET_ITEM(parts, i)
                node = node.get(kv)
        else:
            parts_fast = PySequence_Fast(parts, "Expected a sequence")
            n = PySequence_Fast_GET_SIZE(parts_fast)
            for i in range(n):
                kv = <object>PySequence_Fast_GET_ITEM(parts_fast, i)
                node = node.get(kv)
        return node


cdef _NameTree stub_fragment_name_cache = _NameTree()


cdef object pair_mass(tuple x):
    tmp = <object>PyTuple_GetItem(x, 0)
    return tmp.mass()


cdef basestring build_name_from_composition(mapping_types glycan_composition):
        cdef:
            basestring name, extended_key
            _NameTree root, node
            list parts
        name = 'peptide'
        root = stub_fragment_name_cache
        if mapping_types is dict:
            parts = PyDict_Items(glycan_composition)
        elif mapping_types is CountTable:
            parts = count_table_items(glycan_composition.table)
        elif mapping_types is object:
            parts = list(glycan_composition.items())
        # traverse the tree to find the node referring to
        # this combination of monosaccharides and counts
        node = root.traverse(parts)
        extended_key = node.name
        if extended_key is None:
            parts.sort(key=pair_mass)
            node.name = extended_key = ''.join("%s%d" % kv for kv in parts)
        if extended_key:
            name = "%s+%s" % (name, extended_key)
        return name


cdef class StubFragment(FragmentBase):

    _name_cache = stub_fragment_name_cache

    @staticmethod
    cdef StubFragment _create(str name, double mass, IonSeriesBase kind, CComposition composition,
                              ChemicalShiftBase chemical_shift, bint is_glycosylated, object glycosylation,
                              bint is_extended, int glycosylation_size=-1):
        cdef:
            StubFragment self

        self = StubFragment.__new__(StubFragment)
        self._name = name
        self._hash = hash(self._name)
        self._chemical_shift = None
        self._glycosylation_size = glycosylation_size

        self.mass = mass
        self.kind = kind
        self.composition = composition
        if chemical_shift is not None:
            StubFragment.set_chemical_shift(self, chemical_shift)

        self.is_glycosylated = is_glycosylated
        self.glycosylation = glycosylation
        self.is_extended = is_extended

        return self

    def __init__(self, name, mass, kind, composition, chemical_shift=None, is_glycosylated=False,
                 glycosylation=None, is_extended=False):
        self._name = name
        self._hash = hash(self._name)
        self._chemical_shift = None

        self.mass = mass
        self.kind = kind
        self.composition = composition

        self.set_chemical_shift(chemical_shift)
        self.is_glycosylated = is_glycosylated

        self.glycosylation = glycosylation
        self.is_extended = is_extended
        self._glycosylation_size = -1

    cpdef clone(self):
        cdef StubFragment dup = StubFragment._create(
            self.name, self.mass, self.kind, self.composition, self._chemical_shift,
            self.is_glycosylated, self.glycosylation, self.is_extended,
            self._glycosylation_size)
        return dup

    cdef int get_glycosylation_size(self):
        if self._glycosylation_size == -1:
            self._glycosylation_size = sum(self.glycosylation.values())
        return self._glycosylation_size

    @property
    def glycosylation_size(self):
        return self.get_glycosylation_size()

    cpdef str base_name(self):
        """Simply return the base fragment's name, omitting any chemical modification,
        e.g. peptide+HexNAc1"""

        # If the chemical shift is :const:`None`, there was no change to :attr:`_name`,
        # so just return it, otherwise we need to recompute the name from scratch using
        # :attr:`glycosylation` and :meth:`build_name_from_composition`
        if self._chemical_shift is None:
            return self._name
        else:
            return self.build_name_from_composition(self.glycosylation)

    @classmethod
    def build_name_from_composition(cls, glycan_composition):
        return build_name_from_composition(glycan_composition)

    def __reduce__(self):
        corrected_mass = self.mass
        if self._chemical_shift is not None:
            corrected_mass -= self.chemical_shift.mass
        return self.__class__, (self.name, corrected_mass, self.kind, self.composition,
                                self.chemical_shift, self.is_glycosylated, self.glycosylation,
                                self.is_extended)

    def __repr__(self):
        return ("{self.__class__.__name__}(name={self.name}, "
                "mass={self.mass:.04f}, series={self.kind})").format(self=self)

    cpdef IonSeriesBase get_series(self):
        return self.kind
