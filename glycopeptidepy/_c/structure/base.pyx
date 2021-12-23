cimport cython
from glypy.composition.ccomposition cimport CComposition
from glypy.composition import formula as _formula

from cpython.list cimport PyList_GetItem, PyList_GET_ITEM, PyList_Size, PyList_GET_SIZE
from cpython.sequence cimport PySequence_GetItem


cdef object formula = _formula

DEF SAFE_GET = 1

@cython.freelist(10000)
cdef class PeptideSequenceBase(object):

    cdef SequencePosition get(self, ssize_t i):
        IF SAFE_GET:
            return <SequencePosition>PyList_GetItem(self.sequence, i)
        ELSE:
            return <SequencePosition>PyList_GET_ITEM(self.sequence, i)

    cpdef _invalidate(self):
        pass


@cython.freelist(100)
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

    cpdef TerminalGroup modify(self, ModificationBase modification):
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
        return self.serialize()

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

    cdef bint equal_to(self, object other):
        cdef:
            TerminalGroup other_terminal
        if isinstance(other, TerminalGroup):
            other_terminal = <TerminalGroup>other
            return (self.base_composition == other_terminal.base_composition) and (
                self.modification == other_terminal.modification)
        else:
            if isinstance(other, basestring):
                from glycopeptidepy.structure.modification import Modification
                return self.composition == Modification(other).composition
            elif isinstance(other, ModificationBase):
                return self.composition == (<ModificationBase>other).composition
            else:
                try:
                    return self.composition == other.composition
                except AttributeError:
                    return False

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(formula(self.get_composition()))

    cpdef str serialize(self):
        if self.modification is not None:
            return str(self.modification)
        return formula(self.base_composition)


@cython.freelist(100)
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


@cython.freelist(100000)
cdef class ModificationBase(object):

    def __copy__(self):
        return self.clone()

    cpdef bint is_a(self, object category):
        '''Returns whether or not this :class:`ModificationBase` object belongs to
        the specified :class:`~.ModificationCategory`.

        Returns
        -------
        bool
        '''
        return False

    cpdef basestring serialize(self):
        '''A string representation for inclusion in sequences'''
        return self.name

    cpdef bint is_tracked_for(self, object category):
        """Determine if this :class:`ModificationBase` is tracked by a particular
        behavioral pattern associated with a :class:`~.ModificationCategory`.

        This relationship is distinct from :meth:`is_a` which merely observes that
        the semantic relationship holds, not that any actual behavior is available.

        Parameters
        ----------
        category : :class:`~.ModificationCategory`
            The category to check

        Returns
        -------
        bool
        """
        return False


@cython.final
@cython.freelist(10000000)
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
            size_t n, i
        if isinstance(other, SequencePosition):
            other_t = other
            if not self.amino_acid.equal_to(other_t.amino_acid):
                return False
            n = self.get_modification_count()
            if n != other_t.get_modification_count():
                return False
            if self.modifications != other_t.modifications:
                return False
            return True
        else:
            return self.amino_acid == other.amino_acid and (self.modifications == other.modifications)

    def __ne__(self, other):
        return not (self == other)

    def __reduce__(self):
        return self.__class__, ((self.amino_acid, self.modifications),)

    def __len__(self):
        return 2

    def __repr__(self):
        if self.modifications is None:
            modifications = '-'
        else:
            modifications = "%r" % self.modifications
        return "[%r, %s]" % (self.amino_acid, modifications)

    @staticmethod
    cdef SequencePosition _create(AminoAcidResidueBase amino_acid, list modifications):
        cdef SequencePosition inst = SequencePosition.__new__(SequencePosition)
        inst.amino_acid = amino_acid
        inst.modifications = modifications
        return inst

    cdef inline double get_mass(self):
        cdef:
            double mass
            size_t i, n
            ModificationBase mod

        mass = 0
        mass += self.amino_acid.mass
        n = self.get_modification_count()
        for i in range(n):
            mod = self.get_modification(i)
            mass += mod.mass
        return mass

    cpdef bint has_modification(self, modification):
        cdef:
            size_t i, n
            ModificationBase mod

        n = self.get_modification_count()
        for i in range(n):
            mod = self.get_modification(i)
            if modification == mod:
                return True
        return False

    cpdef bint is_modified(self):
        return self.get_modification_count() > 0

    cpdef add_modification(self, modification):
        if self.modifications is None:
            self.modifications = [modification]
        else:
            self.modifications.append(modification)

    cpdef drop_modification(self, modification):
        cdef:
            size_t i, n,
            int index
            ModificationBase mod

        n = self.get_modification_count()
        index = -1
        for i in range(n):
            mod = <ModificationBase>PyList_GET_ITEM(self.modifications, i)
            if modification == mod:
                index = i
                break
        if index == -1:
            raise ValueError("Modification %s not found" % modification)
        return self.modifications.pop(index)

    @property
    def mass(self):
        return self.get_mass()

    cdef size_t get_modification_count(self):
        if self.modifications is None:
            return 0
        return PyList_GET_SIZE(self.modifications)

    cdef ModificationBase get_modification(self, size_t i):
        return <ModificationBase>PyList_GET_ITEM(self.modifications, i)