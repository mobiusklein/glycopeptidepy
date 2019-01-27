cimport cython

from cpython.ref cimport Py_INCREF
from cpython cimport PyObject
from cpython.list cimport PyList_GetItem, PyList_SetItem, PyList_Size, PyList_New
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem

from glycopeptidepy._c.structure.base cimport (
    PeptideSequenceBase,
    SequencePosition, TerminalGroup, AminoAcidResidueBase,
    ModificationBase)

from glycopeptidepy._c.parser import sequence_tokenizer

from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.structure.constants cimport Configuration

from glycopeptidepy.structure.residue import AminoAcidResidue
from glycopeptidepy.structure.modification import Modification, ModificationCategory
from glycopeptidepy.structure.glycan import GlycosylationManager
from glycopeptidepy.structure import constants as _structure_constants



cdef Configuration structure_constants = _structure_constants
cdef object ModificationImpl = Modification
cdef object AminoAcidResidueImpl = AminoAcidResidue
cdef object GlycosylationManagerImpl = GlycosylationManager

cdef dict terminal_group_cache = dict()


cdef TerminalGroup _make_terminal_group(basestring base_composition_formula, ModificationBase modification=None):
    key = (base_composition_formula, modification)
    result = PyDict_GetItem(terminal_group_cache, key)
    if result == NULL:
        inst = TerminalGroup(base_composition_formula, modification)
        PyDict_SetItem(terminal_group_cache, key, inst)
        return inst
    return <TerminalGroup>result


cdef dict amino_acid_cache = dict()


cdef AminoAcidResidueBase _parse_residue(basestring residue_string):
    result = PyDict_GetItem(amino_acid_cache, residue_string)
    if result == NULL:
        inst = AminoAcidResidueImpl(residue_string)
        PyDict_SetItem(amino_acid_cache, residue_string, inst)
        return inst
    return <AminoAcidResidueBase>result


cdef object ModificationCategory_glycosylation = ModificationCategory.glycosylation


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
        m = PyList_Size(position.modifications)
        for j in range(m):
            mod = <ModificationBase>PyList_GetItem(position.modifications, j)
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

        i = 0
        n = PyList_Size(seq_list)
        self.sequence = sequence = PyList_New(n)
        mass = 0
        for i in range(n):
            item = <list>PyList_GetItem(seq_list, i)
            res = _parse_residue(<basestring>PyList_GetItem(item, 0))
            mass += res.mass
            mods = []
            mod_strs = <list>PyList_GetItem(item, 1)
            m = PyList_Size(mod_strs)
            for j in range(m):
                mod_str = <basestring>PyList_GetItem(mod_strs, j)
                if mod_str != '':
                    mod = ModificationImpl(mod_str)
                    if mod.is_tracked_for(ModificationCategory_glycosylation):
                        self._glycosylation_manager[i] = mod
                    mods.append(mod)
                    mass += mod.mass

            position = SequencePosition._create(res, mods)
            Py_INCREF(position)
            PyList_SetItem(sequence, i, position)

        self._mass = mass
        has_glycan = glycan is not None
        if has_glycan:
            self._glycosylation_manager.aggregate = glycan

        if isinstance(n_term, basestring):
            if n_term != structure_constants.N_TERM_DEFAULT:
                n_term = ModificationImpl(n_term)
            else:
                n_term = None
        if isinstance(c_term, basestring):
            if c_term != structure_constants.C_TERM_DEFAULT:
                c_term = ModificationImpl(c_term)
            else:
                c_term = None
        n_term_group = _make_terminal_group(structure_constants.N_TERM_DEFAULT, n_term)
        c_term_group = _make_terminal_group(structure_constants.C_TERM_DEFAULT, c_term)
        self._init_termini(
            _make_terminal_group(structure_constants.N_TERM_DEFAULT, n_term),
            _make_terminal_group(structure_constants.C_TERM_DEFAULT, c_term))

    cpdef _init_from_string(self, sequence, parser_function):
        """Initialize a :class:`PeptideSequence` from a parse-able string.

        Parameters
        ----------
        sequence : :class:`str`
            The string to parse
        parser_function : :class:`Callable`
            The parsing function to use
        """
        if sequence == "" or sequence is None:
            pass
        else:
            seq_list, modifications, glycan, n_term, c_term = parser_function(
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

        i = 0
        n = PyList_Size(seq_list)
        self.sequence = sequence = PyList_New(n)
        mass = 0

        for i in range(n):
            item = <SequencePosition>PyList_GetItem(seq_list, i)
            mass += item.amino_acid.mass
            m = PyList_Size(item.modifications)
            mods = PyList_New(m)
            for j in range(m):
                mod = <ModificationBase>PyList_GetItem(item.modifications, j)
                Py_INCREF(mod)
                PyList_SetItem(mods, j, mod)
                if mod.is_tracked_for(ModificationCategory_glycosylation):
                    self._glycosylation_manager[i] = mod
                mass += mod.mass
            position = SequencePosition._create(item.amino_acid, mods)
            Py_INCREF(position)
            PyList_SetItem(sequence, i, position)

        self._mass = mass
        if glycan_composition is not None:
            self._glycosylation_manager.aggregate = glycan_composition.clone()
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
        self._fragment_index = None

        self._glycosylation_manager = GlycosylationManagerImpl(self)

        self._n_term = None
        self._c_term = None

        self._fragments_map = {}
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

    @property
    def peptide_backbone_mass(self):
        return self.peptide_composition().mass

    @property
    def total_mass(self):
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
            m = PyList_Size(position.modifications)
            for j in range(m):
                mod = <ModificationBase>PyList_GetItem(position.modifications, j)
                if mod.is_tracked_for(ModificationCategory_glycosylation):
                    continue
                total += mod.mass

        total += self._n_term.mass
        total += self._c_term.mass
        gc = self.glycan_composition
        if gc:
            total += gc.mass()
        return total

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
        sub = self.get(index)
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
        return str(self) == str(other)

    def base_sequence_equality(self, other):
        if len(self) != len(other):
            return False
        for a, b in zip(self, other):
            if a[0] != b[0]:
                return False
        return True

    def modified_sequence_equality(self, other):
        if len(self) != len(other):
            return False
        for a, b in zip(self, other):
            if a[0] != b[0] or set(a[1]) != set(b[1]):
                return False
        return True

    def full_structure_equality(self, other):
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

        seq_list = []
        n = self.get_size()
        for i in range(n):
            position = self.get(i)
            mod_str = ''
            m = PyList_Size(position.modifications)
            if m > 0:
                mod_strs = []
                for j in range(m):
                    mod = <ModificationBase>PyList_GetItem(position.modifications, j)
                    mod_strs.append(mod.serialize())
                mod_str = '|'.join(mod_strs)
                mod_str = "(" + mod_str + ")"
            seq_list.append(position.amino_acid.symbol + mod_str)
        rep = ''.join(seq_list)
        if include_termini:
            needs_terminals = False
            n_term = ""
            if self.n_term is not None and self.get_n_term() != implicit_n_term:
                n_term = "({0})-".format(self.get_n_term().serialize())
                needs_terminals = True
            c_term = ""
            if self.c_term is not None and self.get_c_term() != implicit_c_term:
                c_term = "-({0})".format(self.get_c_term().serialize())
                needs_terminals = True
            if needs_terminals:
                rep = "{}{}{}".format(n_term, rep, c_term)
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
            glycan.composition_offset = CComposition("H2O")
        else:
            glycan = None
        inst._init_from_components(
            self.sequence, glycan,
            self.n_term.modification,
            self.c_term.modification)
        return inst