cimport cython
from glypy.composition.ccomposition cimport CComposition
from glycopeptidepy._c.count_table cimport CountTable


cdef class ChemicalShiftBase(object):
    cdef:
        public str name
        public CComposition composition
        public double mass

    cpdef ChemicalShiftBase clone(self)


cdef class IonSeriesBase(object):
    cdef:
        public str name
        public int direction
        public bint includes_peptide
        public double mass_shift
        public CComposition composition_shift
        public Py_hash_t _hash
        public int int_code


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

    cdef:
        public str _name
        public Py_hash_t _hash
        public double mass
        public ChemicalShiftBase _chemical_shift

    cpdef IonSeriesBase get_series(self)

    cpdef ChemicalShiftBase get_chemical_shift(self)

    cpdef set_chemical_shift(self, ChemicalShiftBase chemical_shift)

    cpdef str base_name(self)
    cpdef str get_fragment_name(self)
    cdef void _update_hash_name(self)

    cdef str get_name(self)



cdef class SimpleFragment(FragmentBase):
    cdef:
        public IonSeriesBase kind
        public bint is_glycosylated
        public CComposition composition

    cpdef clone(self)


cdef class StubFragment(FragmentBase):
    cdef:
        public IonSeriesBase kind
        public bint is_glycosylated
        public CComposition composition
        public object glycosylation
        public bint is_extended
        int _glycosylation_size

    cpdef clone(self)

    @staticmethod
    cdef StubFragment _create(str name, double mass, IonSeriesBase kind, CComposition composition,
                              ChemicalShiftBase chemical_shift, bint is_glycosylated, object glycosylation,
                              bint is_extended, int glycosylation_size=*)

    cdef int get_glycosylation_size(self)


ctypedef fused mapping_types:
        CountTable
        dict
        object


cdef class PeptideFragment(FragmentBase):
    cdef:
        public IonSeriesBase kind
        public int position
        public CountTable modification_dict
        public double bare_mass
        public list flanking_amino_acids
        public dict glycosylation
        public CComposition composition

    @staticmethod
    cdef PeptideFragment _create(IonSeriesBase kind, int position, CountTable modification_dict, double mass,
                                 list flanking_amino_acids=*, dict glycosylation=*,
                                 ChemicalShiftBase chemical_shift=*, CComposition composition=?, double* delta_mass=?)

    cpdef clone(self)

    cdef void _update_mass_with_modifications(self)
    cdef bint _is_glycosylated(self)
    cdef long get_glycosylation_size(self) except -1


ctypedef fused list_or_iterable:
    list
    object


@cython.final
cdef class _NameTree(object):
    cdef:
        public object name
        public object value
        public dict children

    @staticmethod
    cdef _NameTree _create()

    cpdef _NameTree get(self, key)
    cpdef traverse(self, list_or_iterable parts)

cdef basestring build_name_from_composition(mapping_types glycan_composition)