cimport cython

from cpython.ref cimport Py_INCREF
from cpython cimport PyObject
from cpython.list cimport PyList_GetItem, PyList_SetItem, PyList_Size, PyList_New, PyList_GET_ITEM, PyList_GET_SIZE
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.int cimport PyInt_AsLong
from cpython.float cimport PyFloat_AsDouble

from glycopeptidepy._c.structure.base cimport (
    PeptideSequenceBase,
    SequencePosition, TerminalGroup, AminoAcidResidueBase,
    ModificationBase)

from glycopeptidepy._c.structure.modification.source cimport ModificationTableBase
from glycopeptidepy._c.structure.modification.modification cimport ModificationInstanceBase

from glycopeptidepy._c.parser import sequence_tokenizer

from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.structure.constants cimport Configuration
from glycopeptidepy._c.structure.glycan cimport GlycosylationManager, GlycanCompositionWithOffsetProxy

from glycopeptidepy.structure.residue import AminoAcidResidue
from glycopeptidepy.structure.modification import Modification, ModificationCategory, SequenceLocation
from glycopeptidepy.structure import constants as _structure_constants

cdef CComposition WATER = CComposition("H2O")

cdef Configuration structure_constants = _structure_constants

cdef object ModificationImpl = Modification
cdef object AminoAcidResidueImpl = AminoAcidResidue

cdef object SequenceLocation_n_term = SequenceLocation.n_term
cdef object SequenceLocation_c_term = SequenceLocation.c_term

cdef object ModificationCategory_glycosylation = ModificationCategory.glycosylation

cdef ModificationTableBase modification_table = ModificationImpl.modification_table


cpdef ModificationBase get_modification_by_name(basestring name):
    cdef:
        ModificationInstanceBase inst

    rule = modification_table.resolve(name)
    inst = <ModificationInstanceBase>ModificationImpl.__new__(ModificationImpl)
    inst._init_from_rule(rule)
    return inst


cdef dict terminal_group_cache = dict()

# Helpers and Caches -----------------

cdef TerminalGroup _make_terminal_group(basestring base_composition_formula, ModificationBase modification=None):
    """A cache for creating :class:`~.TerminalGroup` instances from base formula and
    modification.

    Returns
    -------
    TerminalGroup
    """
    key = (base_composition_formula, modification)
    result = PyDict_GetItem(terminal_group_cache, key)
    if result == NULL:
        inst = TerminalGroup(base_composition_formula, modification)
        PyDict_SetItem(terminal_group_cache, key, inst)
        return inst
    return <TerminalGroup>result


cdef dict amino_acid_cache = dict()


cdef AminoAcidResidueBase _parse_residue(basestring residue_string):
    """A substitute cache for :class:`~.AminoAcidResidue` instantiation.

    Normally, creating an :class:`AminoAcidResidue` instance passes through
    a cache in the class's constructor. This however incurs a lot of Python
    overhead because it must go through the class's ``__call__`` slot *every
    time*, so there are many pointer jumps along the way.

    This function implements a cache in front of that cache because creating
    AminoAcidResidue instances is done so frequently by the :class:`_PeptideSequenceCore`
    initialization code. It just maintains a static mapping between the amino acid symbol
    to the instance of the AminoAcidResidue that represents it, e.g. "N" -> Asparagine.

    Returns
    -------
    AminoAcidResidueBase
    """
    result = PyDict_GetItem(amino_acid_cache, residue_string)
    if result == NULL:
        inst = AminoAcidResidueImpl(residue_string)
        PyDict_SetItem(amino_acid_cache, residue_string, inst)
        return inst
    return <AminoAcidResidueBase>result


cdef dict implicit_terminal_cache = dict()


cdef ModificationBase _get_implicit_terminal(basestring terminal_name):
    """A cache for storing the Modification mapping for implicit N- and C-terminal
    groups

    Returns
    -------
    ModificationBase
    """
    result = PyDict_GetItem(implicit_terminal_cache, terminal_name)
    if result == NULL:
        inst = ModificationImpl(terminal_name)
        PyDict_SetItem(implicit_terminal_cache, terminal_name, inst)
        return inst
    return <ModificationBase>result

# -------------------


cdef CComposition _total_composition(_PeptideSequenceCore sequence):
    cdef:
        CComposition total
        size_t i, j, m, n
        SequencePosition position
        ModificationBase mod

    n = sequence.get_size()
    total = CComposition()
    for i in range(n):
        position = sequence.get(i)
        total += position.amino_acid.composition
        m = SequencePosition.get_modification_count(position)
        for j in range(m):
            mod = SequencePosition.get_modification(position, j)
            if mod.is_tracked_for(ModificationCategory_glycosylation):
                continue
            total += mod.composition

    total += sequence._n_term.get_composition()
    total += sequence._c_term.get_composition()
    gc = sequence.glycan_composition
    if gc:
        total += gc.total_composition()
    return total


cdef class _PeptideSequenceCore(PeptideSequenceBase):
    def __init__(self, sequence=None, parser_function=None, **kwargs):
        if parser_function is None:
            parser_function = sequence_tokenizer
        self._initialize_fields()
        if sequence is not None:
            self._init_from_string(sequence, parser_function)

    cpdef _init_from_parsed_string(self, list seq_list, glycan=None, n_term=None, c_term=None):
        cdef:
            size_t i, j, n, m
            list sequence, item, mods, mod_strs
            double mass
            ModificationBase mod
            AminoAcidResidueBase res
            SequencePosition position
            bint has_glycan
            GlycosylationManager glycosylation_manager

        i = 0
        n = PyList_GET_SIZE(seq_list)
        self.sequence = sequence = PyList_New(n)
        glycosylation_manager = <GlycosylationManager>self._glycosylation_manager
        mass = 0
        for i in range(n):
            item = <list>PyList_GET_ITEM(seq_list, i)
            res = _parse_residue(<basestring>PyList_GET_ITEM(item, 0))
            mass += res.mass
            mod_strs = <list>PyList_GET_ITEM(item, 1)
            mods = None
            if mod_strs is not None:
                mods = []
                m = PyList_GET_SIZE(mod_strs)
                for j in range(m):
                    mod_str = <basestring>PyList_GET_ITEM(mod_strs, j)
                    if mod_str != '':
                        mod = get_modification_by_name(mod_str)
                        if mod.is_tracked_for(ModificationCategory_glycosylation):
                            glycosylation_manager[i] = mod
                        mods.append(mod)
                        mass += mod.mass

            position = SequencePosition._create(res, mods)
            Py_INCREF(position)
            PyList_SetItem(sequence, i, position)

        self._mass = mass
        has_glycan = glycan is not None
        if has_glycan:
            glycosylation_manager.set_aggregate(glycan)

        if isinstance(n_term, basestring):
            if n_term != structure_constants.N_TERM_DEFAULT:
                n_term = get_modification_by_name(n_term)
            else:
                n_term = None
        if isinstance(c_term, basestring):
            if c_term != structure_constants.C_TERM_DEFAULT:
                c_term = get_modification_by_name(c_term)
            else:
                c_term = None
        n_term_group = _make_terminal_group(structure_constants.N_TERM_DEFAULT, n_term)
        c_term_group = _make_terminal_group(structure_constants.C_TERM_DEFAULT, c_term)
        self._init_termini(
            n_term_group,
            c_term_group)

    cpdef _init_from_string(self, sequence, parser_function):
        """Initialize a :class:`PeptideSequence` from a parse-able string.

        Parameters
        ----------
        sequence : :class:`str`
            The string to parse
        parser_function : :class:`Callable`
            The parsing function to use
        """
        if sequence is None or sequence == "":
            pass
        else:
            seq_list, _, glycan, n_term, c_term = <tuple>parser_function(
                sequence)
            self._init_from_parsed_string(seq_list, glycan, n_term, c_term)

    cpdef _init_from_components(self, list seq_list, glycan_composition=None, ModificationBase n_term=None, ModificationBase c_term=None):
        """Initialize a :class:`PeptideSequence`'s state from a sequence of (
        :class:`~.AminoAcidResidue`, [:class:`~.Modification`, ]) pairs, with optional extra information
        defining the glycan composition aggregate and the terminal groups

        Parameters
        ----------
        seq_list : :class:`Sequence`
            The sequence of (:class:`~.AminoAcidResidue`, [:class:`~.Modification`, ])
        glycan_composition : :class:`~.GlycanComposition`, optional
            The aggregate glycan composition associated with the sequence not fully specified
            by the modifications in ``seq_list``.
        n_term : :class:`~.Modification`, optional
            The N-terminal modification, if any
        c_term : :class:`~.Modification`, optional
            The C-terminal modification, if any
        """
        cdef:
            size_t i, j, n, m
            list sequence, mods, mod_strs
            SequencePosition item, position
            double mass
            ModificationBase mod
            AminoAcidResidueBase res
            GlycosylationManager glycosylation_manager

        i = 0
        n = PyList_GET_SIZE(seq_list)
        self.sequence = sequence = PyList_New(n)
        mass = 0
        glycosylation_manager = <GlycosylationManager>self._glycosylation_manager
        for i in range(n):
            item = <SequencePosition>PyList_GET_ITEM(seq_list, i)
            mass += item.amino_acid.mass
            if item.modifications is not None:
                m = SequencePosition.get_modification_count(item)
                mods = PyList_New(m)
                for j in range(m):
                    mod = SequencePosition.get_modification(item, j)
                    Py_INCREF(mod)
                    PyList_SetItem(mods, j, mod)
                    if mod.is_tracked_for(ModificationCategory_glycosylation):
                        glycosylation_manager[i] = mod
                    mass += mod.mass
            else:
                mods = None

            position = SequencePosition._create(item.amino_acid, mods)
            Py_INCREF(position)
            PyList_SetItem(sequence, i, position)

        self._mass = mass
        if glycan_composition is not None:
            glycosylation_manager.set_aggregate(glycan_composition.clone())

        self._init_termini(
            _make_terminal_group(structure_constants.N_TERM_DEFAULT, n_term),
            _make_terminal_group(structure_constants.C_TERM_DEFAULT, c_term))

    cpdef _init_termini(self, TerminalGroup n_term, TerminalGroup c_term):
        '''This method assumes that the terminal groups are not yet initialized,
        and that their mass contributions are not reflected in :attr:`_mass`
        '''
        self._n_term = n_term
        self._c_term = c_term
        self._mass += n_term.mass + c_term.mass

    cpdef _initialize_fields(self):
        """Initialize all mutable fields to their default, empty values.
        """
        self._mass = 0.0
        self.sequence = []

        self._glycosylation_manager = GlycosylationManager._create(self, None)

        self._n_term = None
        self._c_term = None

        self._total_composition = None
        self._peptide_composition = None

    cpdef _invalidate(self):
        self._total_composition = None
        self._peptide_composition = None

    cdef TerminalGroup get_n_term(self):
        return self._n_term

    cdef void set_n_term(self, _value):
        cdef:
            TerminalGroup value
            double reset_mass, new_mass
        if isinstance(_value, ModificationBase):
            value = self._n_term.modify(_value)
        else:
            value = _value

        self._invalidate()
        reset_mass = 0
        try:
            reset_mass = self._n_term.mass
        except AttributeError:
            pass
        self._n_term = value
        new_mass = 0
        try:
            new_mass = value.mass
        except AttributeError:
            pass
        self._mass += new_mass - reset_mass

    cdef TerminalGroup get_c_term(self):
        return self._c_term

    cdef void set_c_term(self, _value):
        cdef:
            TerminalGroup value
            double reset_mass, new_mass
        if isinstance(_value, ModificationBase):
            value = self._c_term.modify(_value)
        else:
            value = _value
        self._invalidate()
        reset_mass = 0
        try:
            reset_mass = self._c_term.mass
        except AttributeError:
            pass
        self._c_term = value
        new_mass = 0
        try:
            new_mass = value.mass
        except AttributeError:
            pass
        self._mass += new_mass - reset_mass

    def __repr__(self):
        n_term = ""
        if self.n_term is not None:
            n_term = "({0})-".format(self.n_term)
        c_term = ""
        if self.c_term is not None:
            c_term = "-({0})".format(self.c_term)
        aggregate = self._glycosylation_manager.aggregate if self._glycosylation_manager.aggregate is not None else ""
        rep = "{n_term}{sequence}{c_term}{glycan_aggregate}[{_mass}]".format(
            n_term=n_term, c_term=c_term,
            glycan_aggregate=aggregate,
            sequence=self.sequence,
            _mass=self._mass)
        return rep

    @property
    def n_term(self):
        return self.get_n_term()

    @n_term.setter
    def n_term(self, _value):
        self.set_n_term(_value)

    @property
    def c_term(self):
        return self.get_c_term()

    @c_term.setter
    def c_term(self, _value):
        self.set_c_term(_value)

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value

    cdef double get_peptide_backbone_mass(self):
        cdef:
            double total
            size_t i, j, m, n
            SequencePosition position
            ModificationBase mod

        n = self.get_size()
        total = 0
        for i in range(n):
            position = self.get(i)
            total += position.amino_acid.mass
            m = SequencePosition.get_modification_count(position)
            for j in range(m):
                mod = SequencePosition.get_modification(position, j)
                if mod.is_tracked_for(ModificationCategory_glycosylation):
                    continue
                total += mod.mass

        total += self._n_term.mass
        total += self._c_term.mass
        return total

    @property
    def peptide_backbone_mass(self):
        return self.get_peptide_backbone_mass()

    cdef double get_total_mass(self):
        cdef:
            double total
        total = self.get_peptide_backbone_mass()
        gc = self.glycan_composition
        if gc:
            total += PyFloat_AsDouble(gc.mass())
        return total

    @property
    def total_mass(self):
        return self.get_total_mass()

    cpdef CComposition peptide_composition(self):
        if self._peptide_composition is None:
            if self.glycan is None:
                glycan_composition = CComposition()
            else:
                glycan_composition = self.glycan_composition.total_composition()
            self._peptide_base_composition = self.total_composition() - glycan_composition
        return self._peptide_base_composition

    cpdef CComposition total_composition(self):
        if self._total_composition is None:
            self._total_composition = _total_composition(self)
        return self._total_composition

    def __iter__(self):
        return iter(self.sequence)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.sequence[index]
        cdef c_index = PyInt_AsLong(index)
        if c_index < 0:
            return self.sequence[index]
        sub = self.get(c_index)
        return sub

    def __setitem__(self, index, value):
        self._invalidate()
        self.sequence[index] = value

    def __len__(self):
        return PyList_Size(self.sequence)

    cdef size_t get_size(self):
        return PyList_Size(self.sequence)

    # equality methods

    def __eq__(self, other):
        if isinstance(other, _PeptideSequenceCore):
            return self.full_structure_equality(<_PeptideSequenceCore>other)
        else:
            return str(self) == str(other)

    cpdef bint base_sequence_equality(self, _PeptideSequenceCore other):
        cdef:
            size_t n, m, i
        n = self.get_size()
        m = other.get_size()
        if n != m:
            return False
        for i in range(n):
            if self.get(i).amino_acid != other.get(i).amino_acid:
                return False
        return True

    cpdef bint modified_sequence_equality(self, _PeptideSequenceCore other):
        cdef:
            size_t n, m, i
            SequencePosition a, b
        n = self.get_size()
        m = other.get_size()
        if n != m:
            return False
        for i in range(n):
            a = self.get(i)
            b = other.get(i)
            if a != b:
                return False
        if self.get_n_term() != other.get_n_term():
            return  False
        if self.get_c_term() != other.get_c_term():
            return False
        return True

    cpdef bint full_structure_equality(self, _PeptideSequenceCore other):
        sequences_equal = self.modified_sequence_equality(other)
        if sequences_equal:
            return self.glycan == other.glycan
        else:
            return False

    def __ne__(self, other):
        return str(self) != str(other)

    def __hash__(self):
        return hash(self.get_sequence())

    cpdef basestring get_sequence(self, bint include_glycan=True, bint include_termini=True, str implicit_n_term=None, str implicit_c_term=None):
        """
        Generate human readable sequence string. Called by :meth:`__str__`

        Parameters
        ----------
        include_glycan: bool
            Whether to include the glycan in the resulting string. Defaults to `True`
        include_termini: bool
            Whether to include the N- and C-termini. Make sure this is `True` if you want non-standard
            termini to be properly propagated.


        Returns
        -------
        str
        """
        if implicit_n_term is None:
            implicit_n_term = structure_constants.N_TERM_DEFAULT
        if implicit_c_term is None:
            implicit_c_term = structure_constants.C_TERM_DEFAULT

        cdef:
            list seq_list,  mod_strs
            size_t i, j, n, m
            SequencePosition position
            ModificationBase mod
            basestring mod_str, n_term, c_term, rep
            bint needs_terminals
            TerminalGroup terminal
            ModificationBase implicit_terminal

        seq_list = []
        n = self.get_size()
        for i in range(n):
            position = self.get(i)
            m = SequencePosition.get_modification_count(position)
            if m > 0:
                if m == 1:
                    mod = SequencePosition.get_modification(position, 0)
                    mod_str = "(%s)" % mod.serialize()
                else:
                    mod_strs = []
                    for j in range(m):
                        mod = SequencePosition.get_modification(position, j)
                        mod_strs.append(mod.serialize())
                    mod_str = '|'.join(mod_strs)
                    mod_str = "(%s)" % mod_str
                seq_list.append(position.amino_acid.symbol + mod_str)
            else:
                seq_list.append(position.amino_acid.symbol)

        rep = ''.join(seq_list)
        if include_termini:
            needs_terminals = False
            # Handle N-terminal
            n_term = ""
            terminal = self.get_n_term()
            implicit_terminal = _get_implicit_terminal(implicit_n_term)
            if terminal is not None and not terminal.equal_to(implicit_terminal):
                n_term = "(%s)-" % terminal.serialize()
                needs_terminals = True

            # Handle C-terminal
            c_term = ""
            terminal = self.get_c_term()
            implicit_terminal = _get_implicit_terminal(implicit_c_term)
            if terminal is not None and not terminal.equal_to(implicit_terminal):
                c_term = "-(%s)" % terminal.serialize()
                needs_terminals = True
            if needs_terminals:
                rep = "".join((n_term, rep, c_term))
        if include_glycan:
            if self._glycosylation_manager.aggregate is not None:
                rep += str(self._glycosylation_manager.aggregate)
        return rep

    def __str__(self):
        return self.get_sequence()

    cpdef clone(self):
        cdef _PeptideSequenceCore inst = self.__class__()
        if self._glycosylation_manager.aggregate is not None:
            glycan = self._glycosylation_manager.aggregate.clone()
            if isinstance(glycan, GlycanCompositionWithOffsetProxy):
                glycan.composition_offset = CComposition._create(None)
            else:
                glycan.composition_offset = CComposition._create(WATER)
        else:
            glycan = None
        inst._init_from_components(
            self.sequence, glycan,
            self.get_n_term().get_modification(),
            self.get_c_term().get_modification())
        return inst


@cython.binding(True)
cpdef add_modification(_PeptideSequenceCore self, position, modification_type):
    '''
    Add a modification by name to a specific residue. If the
    position is the N-term or the C-term, the terminal modification will
    be replaced.

    Parameters
    ----------
    position: int
        The position of the modification to add
    modification_type: str or Modification
        The modification to add
    '''
    cdef:
        long position_
        ModificationBase mod
        SequencePosition pos
    self._invalidate()

    if isinstance(modification_type, Modification):
        mod = modification_type
    else:
        mod = Modification(rule=modification_type)

    if position is SequenceLocation_n_term:
        self.n_term = self.n_term.modify(mod)
    elif position is SequenceLocation_c_term:
        self.c_term = self.c_term.modify(mod)
    else:
        position_ = PyInt_AsLong(position)
        if (position_ == -1) or (position_ >= self.get_size()):
            raise IndexError(
                "Invalid modification position. %s, %s, %s" %
                (position, str(self.sequence), modification_type))
        pos = self.get(position_)
        pos.add_modification(mod)
        self._mass += mod.mass
        if mod.is_tracked_for(ModificationCategory_glycosylation):
            self._glycosylation_manager[position] = mod