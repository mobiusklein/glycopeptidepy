cimport cython

from cpython cimport PyObject
from cpython.list cimport PyList_GetItem, PyList_SetItem, PyList_Size
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem

from glycopeptidepy._c.structure.base cimport (
    PeptideSequenceBase,
    SequencePosition, TerminalGroup, AminoAcidResidueBase,
    ModificationBase)

from glycopeptidepy._c.structure.constants cimport Configuration

from glycopeptidepy.structure.residue import AminoAcidResidue
from glycopeptidepy.structure.modification import Modification, ModificationCategory
from glycopeptidepy.structure import constants as _structure_constants


cdef Configuration structure_constants = _structure_constants


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
        inst = AminoAcidResidue(residue_string)
        PyDict_SetItem(amino_acid_cache, residue_string, inst)
        return inst
    return <AminoAcidResidueBase>result


cdef object ModificationCategory_glycosylation = ModificationCategory.glycosylation


@cython.binding
def _init_from_parsed_string(PeptideSequenceBase self, list seq_list, glycan=None, n_term=None, c_term=None, **kwargs):
    cdef:
        size_t i, j, n, m
        list sequence, item, mods, mod_strs
        double mass
        ModificationBase mod
        AminoAcidResidueBase res

    i = 0
    n = PyList_Size(seq_list)
    self.sequence = sequence = [None for _ in seq_list]
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
                mod = Modification(mod_str)
                if mod.is_tracked_for(ModificationCategory_glycosylation):
                    self._glycosylation_manager[i] = mod
                mods.append(mod)
                mass += mod.mass
        self.sequence[i] = SequencePosition._create(res, mods)
    self._mass = mass
    has_glycan = glycan != "" and glycan is not None
    if has_glycan:
        self._glycosylation_manager.aggregate = glycan.clone()

    if isinstance(n_term, basestring):
        if n_term != structure_constants.N_TERM_DEFAULT:
            n_term = Modification(n_term)
        else:
            n_term = None
    if isinstance(c_term, basestring):
        if c_term != structure_constants.C_TERM_DEFAULT:
            c_term = Modification(c_term)
        else:
            c_term = None
    n_term_group = _make_terminal_group(structure_constants.N_TERM_DEFAULT, n_term)
    c_term_group = _make_terminal_group(structure_constants.C_TERM_DEFAULT, c_term)
    self._init_termini(n_term_group, c_term_group)


@cython.binding
def _init_from_components(PeptideSequenceBase self, list seq_list, glycan_composition=None,
                          ModificationBase n_term=None, ModificationBase c_term=None, **kwargs):
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
        SequencePosition item
        double mass
        ModificationBase mod
        AminoAcidResidueBase res

    i = 0
    n = PyList_Size(seq_list)
    self.sequence = sequence = [None for _ in seq_list]
    mass = 0

    for i in range(n):
        item = <SequencePosition>PyList_GetItem(seq_list, i)
        mass += item.amino_acid.mass
        m = PyList_Size(item.modifications)
        mods = []
        for j in range(m):
            mod = <ModificationBase>PyList_GetItem(item.modifications, j)
            mods.append(mod)
            if mod.is_tracked_for(ModificationCategory_glycosylation):
                self._glycosylation_manager[i] = mod
            mass += mod.mass
        self.sequence[i] = SequencePosition._create(item.amino_acid, mods)

    self._mass = mass
    if glycan_composition is not None:
        self._glycosylation_manager.aggregate = glycan_composition.clone()
    self._init_termini(
        _make_terminal_group(structure_constants.N_TERM_DEFAULT, n_term),
        _make_terminal_group(structure_constants.C_TERM_DEFAULT, c_term))
