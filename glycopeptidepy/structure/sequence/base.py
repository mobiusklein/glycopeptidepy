from six import string_types as basestring

from .. import constants as structure_constants
from ..base import PeptideSequenceBase, SequencePosition
from ..composition import Composition
from ..glycan import GlycosylationManager
from ..modification import Modification, ModificationCategory
from ..residue import Residue
from ..terminal_group import _make_terminal_group
from ..parser import sequence_tokenizer


def _total_composition(sequence):
    total = Composition()
    for position in sequence:
        total += position[0].composition
        for mod in position[1]:
            if mod.is_tracked_for(ModificationCategory.glycosylation):
                continue
            total += mod.composition
    total += sequence.n_term.composition
    total += sequence.c_term.composition
    gc = sequence.glycan_composition
    if gc:
        total += gc.total_composition()
    return total


class _PeptideSequenceCore(PeptideSequenceBase):
    def __init__(self, sequence=None, parser_function=None, **kwargs):
        if parser_function is None:
            parser_function = sequence_tokenizer
        self._initialize_fields()
        self._init_from_string(sequence, parser_function, **kwargs)

    def _init_from_parsed_string(self, seq_list, glycan=None, n_term=None, c_term=None, **kwargs):
        i = 0
        self.sequence = [None for _ in seq_list]
        mass = 0
        for item in seq_list:
            res = Residue.parse(item[0])
            mass += res.mass
            mods = []
            for mod in item[1]:
                if mod != '':
                    mod = Modification(mod)
                    if mod.is_tracked_for(ModificationCategory.glycosylation):
                        self._glycosylation_manager[i] = mod
                    mods.append(mod)
                    mass += mod.mass
            self.sequence[i] = SequencePosition([res, mods])
            i += 1
        self.mass = mass
        has_glycan = glycan is not None
        if has_glycan:
            self._glycosylation_manager.aggregate = glycan

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

    def _init_from_string(self, sequence, parser_function, **kwargs):
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
            self._init_from_parsed_string(seq_list, glycan, n_term, c_term, **kwargs)

    def _init_from_components(self, seq_list, glycan_composition=None, n_term=None, c_term=None, **kwargs):
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
        mass = 0
        self.sequence = [None for _ in seq_list]
        i = 0
        for res, mods in seq_list:
            mass += res.mass
            for mod in mods:
                if mod.is_tracked_for(ModificationCategory.glycosylation):
                    self._glycosylation_manager[i] = mod
                mass += mod.mass
            self.sequence[i] = SequencePosition([res, list(mods)])
            i += 1
        self.mass = mass
        if glycan_composition is not None:
            self._glycosylation_manager.aggregate = glycan_composition.clone()
        self._init_termini(
            _make_terminal_group(structure_constants.N_TERM_DEFAULT, n_term),
            _make_terminal_group(structure_constants.C_TERM_DEFAULT, c_term))

    def _init_termini(self, n_term, c_term):
        self._n_term = n_term
        self._c_term = c_term
        self._mass += n_term.mass + c_term.mass

    def _initialize_fields(self):
        """Initialize all mutable fields to their default, empty values.
        """
        self._mass = 0.0
        self.sequence = []
        self._fragment_index = None

        self._glycosylation_manager = GlycosylationManager(self)

        self._n_term = None
        self._c_term = None

        self._total_composition = None
        self._peptide_composition = None

    def _invalidate(self):
        self._total_composition = None
        self._peptide_composition = None

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
        return self._n_term

    @n_term.setter
    def n_term(self, value):
        if isinstance(value, Modification):
            value = self._n_term.modify(value)
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
        self.mass += new_mass - reset_mass

    @property
    def c_term(self):
        return self._c_term

    @c_term.setter
    def c_term(self, value):
        if isinstance(value, Modification):
            value = self._c_term.modify(value)
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
        self.mass += new_mass - reset_mass

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
        total = 0
        for position in self:
            total += position[0].mass
            for mod in position[1]:
                if mod.is_tracked_for(ModificationCategory.glycosylation):
                    continue
                total += mod.mass
        total += self.n_term.mass
        total += self.c_term.mass
        gc = self.glycan_composition
        if gc:
            total += gc.mass()
        return total

    def peptide_composition(self):
        if self._peptide_composition is None:
            if self.glycan is None:
                glycan_composition = Composition()
            else:
                glycan_composition = self.glycan_composition.total_composition()
            self._peptide_base_composition = self.total_composition() - glycan_composition
        return self._peptide_base_composition

    def total_composition(self):
        if self._total_composition is None:
            self._total_composition = _total_composition(self)
        return self._total_composition

    def __iter__(self):
        return iter(self.sequence)

    def __getitem__(self, index):
        sub = self.sequence[index]
        return sub

    def __setitem__(self, index, value):
        self._invalidate()
        self.sequence[index] = value

    def __len__(self):
        return len(self.sequence)

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
        return hash(str(self))

    def get_sequence(self, include_glycan=True, include_termini=True, implicit_n_term=None, implicit_c_term=None):
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

        seq_list = []
        for x, y in self.sequence:
            mod_str = ''
            if y:
                mod_str = '|'.join(mod.serialize() for mod in y)
                mod_str = ''.join(['(', mod_str, ')'])
            seq_list.append(''.join([x.symbol, mod_str]))
        rep = ''.join(seq_list)
        if include_termini:
            n_term = ""
            if self.n_term is not None and self.n_term != implicit_n_term:
                n_term = "({0})-".format(self.n_term.serialize())
            c_term = ""
            if self.c_term is not None and self.c_term != implicit_c_term:
                c_term = "-({0})".format(self.c_term.serialize())
            rep = "{0}{1}{2}".format(n_term, rep, c_term)
        if include_glycan:
            if self._glycosylation_manager.aggregate is not None:
                rep += str(self._glycosylation_manager.aggregate)
        return rep

    def __str__(self):
        return self.get_sequence()

    def clone(self):
        inst = self.__class__()
        if self._glycosylation_manager.aggregate is not None:
            glycan = self._glycosylation_manager.aggregate.clone()
            glycan.composition_offset = Composition("H2O")
        else:
            glycan = None
        inst._init_from_components(
            self.sequence, glycan,
            self.n_term.modification,
            self.c_term.modification)
        return inst


try:
    from glycopeptidepy._c.structure import sequence_methods
    _PeptideSequenceCore = sequence_methods._PeptideSequenceCore
except ImportError as err:
    print(err)
