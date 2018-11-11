from glypy.composition.ccomposition cimport CComposition


cdef class MoleculeBase(object):

    def __copy__(self):
        return self.clone()


cdef class ResidueBase(MoleculeBase):
    cdef:
        public str name
        public str symbol
        public double mass
        public CComposition composition
        public object neutral_loss


cdef class ModificationBase(MoleculeBase):
    cdef:
        public str name
        public double mass
        public CComposition composition


cdef class ChemicalShiftBase(object):
    cdef:
        public str name
        public CComposition composition
        public double mass

    def clone(self):
        return self.__class__(self.name, self.composition.clone())

    def __str__(self):
        return self.name

    def __repr__(self):
        return "%s(name=%r)" % (self.__class__.__name__, self.name)

    def is_loss(self):
        return self.mass < 0


cdef class IonSeriesBase(object):
    cdef:
        public str name
        public int direction
        public bint includes_peptide
        public double mass_shift
        public CComposition composition_shift
        public int _hash

    def __eq__(self, other):
        try:
            return (self is other) or (self.name == other.name)
        except AttributeError:
            return self.name == other

    def __ne__(self, other):
        return not self == other


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

    def get_chemical_shift(self):
        return self._chemical_shift

    def set_chemical_shift(self, ChemicalShiftBase chemical_shift):
        if self._chemical_shift is not None:
            self.mass += self._chemical_shift.mass
        self._chemical_shift = chemical_shift
        if chemical_shift is not None:
            self.mass -= chemical_shift.mass

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
