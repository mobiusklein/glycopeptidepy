from glypy.composition.ccomposition cimport CComposition
from glypy.composition import formula as _formula

from cpython.list cimport PyList_GetItem
from cpython.sequence cimport PySequence_GetItem


cdef object formula = _formula


cdef class PeptideSequenceBase(object):

    cdef SequencePosition get(self, ssize_t i):
        return <SequencePosition>PyList_GetItem(self.sequence, i)

    cpdef _invalidate(self):
        pass


cdef class TerminalGroup(object):
    def __init__(self, base_composition, modification=None):
        if not isinstance(base_composition, CComposition):
            base_composition = CComposition(base_composition)
        self.base_composition = base_composition
        self._modification = None
        if modification is not None:
            self.modification = modification
        self.mass = self._calculate_mass()

    def _calculate_mass(self):
        base_mass = self.base_composition.mass
        mod = self.modification
        if mod is not None:
            base_mass += mod.mass
        return base_mass

    def clone(self):
        return self.__class__(self.base_composition, self.modification)

    def __reduce__(self):
        return self.__class__, (self.base_composition, self.modification)


    cdef ModificationBase get_modification(self):
        return self._modification

    cdef void set_modification(self, ModificationBase value):
        cdef:
            double new_mass, old_mass
        if value is not None:
            new_mass = value.mass
        else:
            new_mass = 0
        if self._modification is not None:
            old_mass = self._modification.mass
        else:
            old_mass = 0
        self.mass += new_mass - old_mass
        self._modification = value

    @property
    def modification(self):
        return self.get_modification()

    @modification.setter
    def modification(self, ModificationBase value):
        self.set_modification(value)

    cpdef TerminalGroup  modify(self, ModificationBase modification):
        return self.__class__(self.base_composition, modification)

    cdef CComposition get_composition(self):
        cdef ModificationBase modification = self.get_modification()
        
        if modification is None:
            return self.base_composition
        mod_comp = modification.composition
        return self.base_composition + mod_comp

    @property
    def composition(self):
        return self.get_composition()

    def __repr__(self):
        template = "{self.__class__.__name__}({self.base_composition}, {self.modification})"
        return template.format(self=self)

    def __str__(self):
        if self.modification is not None:
            return str(self.modification)
        return formula(self.base_composition)

    def __eq__(self, other):
        if other is None:
            return False
        try:
            return (self.base_composition == other.base_composition) and (self.modification == other.modification)
        except AttributeError:
            if isinstance(other, basestring):
                from glycopeptidepy.structure.modification import Modification
                return self.composition == Modification(other).composition
            else:
                try:
                    return self.composition == other.composition
                except AttributeError:
                    return NotImplemented

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(formula(self.get_composition()))

    def serialize(self):
        return str(self)


cdef class AminoAcidResidueBase(object):
    '''
    A base type for classes describing amino acid residues
    '''

    def __copy__(self):
        return self.clone()

    @staticmethod
    cdef AminoAcidResidueBase _create(str name, str symbol, double mass, CComposition composition):
        cdef AminoAcidResidueBase inst = AminoAcidResidueBase.__new__(AminoAcidResidueBase)
        inst.name = name
        inst.symbol = symbol
        inst.mass = mass
        inst.composition = composition
        return inst

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        if self is other:
            return True
        try:
            return self.name == other.name or self.symbol == other.symbol
        except AttributeError:
            return self.name == other or self.symbol == other

    cdef bint equal_to(self, AminoAcidResidueBase other):
        return self.name == other.name or self.symbol == other.symbol

    def __ne__(self, other):
        if self is other:
            return False
        try:
            return self.name != other.name and self.symbol != other.symbol
        except AttributeError:
            return self.name != other and self.symbol != other


cdef class ModificationBase(object):
    
    def __copy__(self):
        return self.clone()

    cpdef bint is_a(self, object category):
        return False

    cpdef basestring serialize(self):
        '''A string representation for inclusion in sequences'''
        return self.name

    cpdef bint is_tracked_for(self, object category):
        return False


cdef class SequencePosition(object):

    def __init__(self, parts):
        self.amino_acid = PySequence_GetItem(parts, 0)
        self.modifications = PySequence_GetItem(parts, 1)

    def __iter__(self):
        yield self.amino_acid
        yield self.modifications

    def __getitem__(self, i):
        if i == 0:
            return self.amino_acid
        elif i == 1:
            return self.modifications
        else:
            raise IndexError(i)

    def __setitem__(self, i, value):
        if i == 0:
            self.amino_acid = value
        elif i == 1:
            self.modifications = value
        else:
            raise IndexError(i)

    def __eq__(self, other):
        cdef:
            SequencePosition other_t
        if isinstance(other, SequencePosition):
            other_t = other
            return self.amino_acid == other_t.amino_acid and self.modifications == other_t.modifications
        else:
            return self.amino_acid == other.amino_acid and self.modifications == other.modifications            

    def __ne__(self, other):
        return not (self == other)

    def __reduce__(self):
        return self.__class__, ((self.amino_acid, self.modifications),)

    def __len__(self):
        return 2

    def __repr__(self):
        return "[%r, %r]" % (self.amino_acid, self.modifications)

    @staticmethod
    cdef SequencePosition _create(AminoAcidResidueBase amino_acid, list modifications):
        cdef SequencePosition inst = SequencePosition.__new__(SequencePosition)
        inst.amino_acid = amino_acid
        inst.modifications = modifications
        return inst

