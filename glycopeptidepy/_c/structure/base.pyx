from glypy.composition.ccomposition cimport CComposition
from cpython.sequence cimport PySequence_GetItem


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


cdef class ModificationBase(object):
    
    def __copy__(self):
        return self.clone()

    cpdef bint is_a(self, object category):
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

    def __eq__(self, other):
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

