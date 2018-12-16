from glypy.composition.ccomposition cimport CComposition


cdef class ChemicalShiftBase(object):
    cdef:
        public str name
        public CComposition composition
        public double mass


cdef class IonSeriesBase(object):
    cdef:
        public str name
        public int direction
        public bint includes_peptide
        public double mass_shift
        public CComposition composition_shift
        public int _hash


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
        public int _hash
        public double mass
        public ChemicalShiftBase _chemical_shift

    cpdef IonSeriesBase get_series(self)

    cpdef ChemicalShiftBase get_chemical_shift(self)

    cpdef set_chemical_shift(self, ChemicalShiftBase chemical_shift)

    cpdef str get_fragment_name(self)


cdef class PeptideFragment(FragmentBase):
    cdef:
        public IonSeriesBase kind
        public int position
        public dict modification_dict
        public double bare_mass
        public list flanking_amino_acids
        public object glycosylation
        public CComposition composition

    cpdef clone(self)

    cdef _update_mass_with_modifications(self)