from six import string_types as basestring

from . import PeptideSequenceBase, MoleculeBase, SequencePosition
from . import constants as structure_constants

from .composition import Composition, formula
from .fragment import (
    IonSeries, _n_glycosylation, _o_glycosylation,
    _gag_linker_glycosylation)
from .modification import (
    Modification, SequenceLocation, ModificationCategory)
from .residue import Residue

from .fragmentation_strategy import (
    HCDFragmentationStrategy,
    CADFragmentationStrategy,
    StubGlycopeptideStrategy,
    OxoniumIonStrategy)

from glypy import GlycanComposition

from .parser import sequence_tokenizer

from .glycan import (
    GlycosylationType, GlycosylationManager,
    glycosylation_site_detectors)

from ..utils.memoize import memoize


b_series = IonSeries.b
y_series = IonSeries.y
oxonium_ion_series = IonSeries.oxonium_ion
stub_glycopeptide_series = IonSeries.stub_glycopeptide


def list_to_sequence(seq_list, wrap=True):
    flat_chunks = []
    for chunk in seq_list:
        if(isinstance(chunk[0], list)):
            flat_chunks.extend(chunk)
        else:
            flat_chunks.append(chunk)
    seq = Sequence.from_iterable(flat_chunks) if wrap else flat_chunks
    return seq


@glycosylation_site_detectors(GlycosylationType.o_linked)
def find_o_glycosylation_sequons(sequence, allow_modified=frozenset()):
    try:
        iter(allow_modified)
        allow_modified = set(allow_modified) | {Modification("HexNAc")}
    except TypeError:
        allow_modified = (Modification("HexNAc"),)

    positions = []
    ser = Residue("S")
    thr = Residue("T")
    asn = Residue("N")

    site_set = (ser, thr)

    if isinstance(sequence, basestring):
        sequence = parse(sequence)

    for i, position in enumerate(sequence):
        if position[0] in site_set:
            if ((len(position[1]) == 0) or position[1][0] in allow_modified) and not (
                    i - 2 >= 0 and sequence[i - 2][0] == asn):
                positions.append(i)
    return positions


@glycosylation_site_detectors(GlycosylationType.n_linked)
def find_n_glycosylation_sequons(sequence, allow_modified=frozenset()):
    try:
        iter(allow_modified)
        allow_modified = set(allow_modified) | {Modification(
            "HexNAc"), Modification(_o_glycosylation)}
    except TypeError:
        allow_modified = (Modification("HexNAc"),
                          Modification(_o_glycosylation))
    state = "seek"  # [seek, n, ^p, st]
    if isinstance(sequence, basestring):
        sequence = PeptideSequence(sequence)

    asn = Residue("Asn")
    pro = Residue("Pro")
    ser = Residue("Ser")
    thr = Residue("Thr")

    i = 0
    positions = []
    n_pos = None
    while(i < len(sequence)):
        next_pos = sequence[i]
        if state == "seek":
            # A sequon starts with an Asn residue without modifications, or for counting
            # purposes one that has already been glycosylated
            if next_pos[0] == asn:
                if ((len(next_pos[1]) == 0) or next_pos[1][0] in allow_modified):
                    n_pos = i
                    state = "n"
        elif state == "n":
            if next_pos[0] != pro:
                state = "^p"
            else:
                state = "seek"
                i = n_pos
                n_pos = None
        elif state == "^p":
            if next_pos[0] in {ser, thr}:
                positions.append(n_pos)
            i = n_pos
            n_pos = None
            state = "seek"
        i += 1
    return positions


@glycosylation_site_detectors(GlycosylationType.glycosaminoglycan)
def find_glycosaminoglycan_sequons(sequence, allow_modified=frozenset()):
    try:
        iter(allow_modified)
        allow_modified = set(allow_modified) | {Modification(
            "Xyl"), Modification(_gag_linker_glycosylation)}
    except TypeError:
        allow_modified = (Modification("Xyl"),
                          Modification(_gag_linker_glycosylation))
    state = "seek"  # [seek, s, g1, x, g2]
    ser = Residue("Ser")
    # gly = Residue("Gly")

    i = 0
    positions = []

    # s_position = None
    if isinstance(sequence, basestring):
        sequence = PeptideSequence(sequence)
    while(i < len(sequence)):
        next_pos = sequence[i]
        if state == "seek":
            if next_pos[0] == ser:
                if ((len(next_pos[1]) == 0) or next_pos[1][0] in allow_modified):
                    positions.append(i)
        i += 1
    return positions


def total_composition(sequence):
    if isinstance(sequence, basestring):
        sequence = parse(sequence)
    return _total_composition(sequence)


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


def _calculate_mass(sequence):
    total = 0
    for position in sequence:
        total += position[0].mass
        for mod in position[1]:
            if mod.is_tracked_for(ModificationCategory.glycosylation):
                continue
            total += mod.mass
    total += sequence.n_term.mass
    total += sequence.c_term.mass
    gc = sequence.glycan_composition
    if gc:
        total += gc.mass()
    return total


class TerminalGroup(MoleculeBase):
    __slots__ = ("base_composition", "mass", "_modification")

    def __init__(self, base_composition, modification=None):
        if not isinstance(base_composition, Composition):
            base_composition = Composition(base_composition)
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

    @property
    def modification(self):
        return self._modification

    @modification.setter
    def modification(self, value):
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

    def modify(self, modification):
        return self.__class__(self.base_composition, modification)

    @property
    def composition(self):
        modification = self.modification
        if modification is None:
            return self.base_composition
        mod_comp = modification.composition
        return self.base_composition + mod_comp

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
                return self.composition == Modification(other).composition
            else:
                try:
                    return self.composition == other.composition
                except AttributeError:
                    return NotImplemented

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(formula(self.composition))

    def serialize(self):
        return str(self)


try:
    from glycopeptidepy._c.structure.base import TerminalGroup
except ImportError:
    pass


@memoize(100)
def _make_terminal_group(base_composition_formula, modification=None):
    return TerminalGroup(base_composition_formula, modification)


class PeptideSequence(PeptideSequenceBase):
    '''
    Represents a peptide that may have post-translational modifications
    including glycosylation.

    Attributes
    ----------
    sequence: list
        The underlying container for positions in the amino acid sequence
    n_term: :class:`TerminalGroup`
    c_term: :class:`TerminalGroup`
        Terminal modifications (N-terminus and C-terminus respectively) which
        default to H and OH respectively.
    glycan: Glycan or GlycanComposition
        The total glycan moiety attached to the molecule. The current semantics
        do not cleanly support more than one glycosylation per sequence for generating
        glycan sequence fragments.
    mass: float
        The pre-calculated monoisotopic mass of the molecule. This quantity is
        assumes that the glycan's glycosidic bonds have been broken, leaving only
        the amide-bound HexNAc as a modification attached to the amino acid backbone
    '''
    position_class = SequencePosition

    @classmethod
    def from_iterable(cls, iterable, glycan_composition=None, n_term=None, c_term=None, text=None, **kwargs):
        seq = cls()
        if text is None:
            iterable = list(iterable)
            try:
                if isinstance(iterable[0][0], basestring):
                    text = True
            except IndexError:
                text = True
        if text:
            seq._init_from_parsed_string(iterable, glycan_composition, n_term, c_term)
        else:
            seq._init_from_components(iterable, glycan_composition, n_term, c_term, **kwargs)
        return seq

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

        self._fragments_map = {}
        self._total_composition = None
        self._peptide_composition = None

    def __init__(self, sequence=None, parser_function=None, **kwargs):
        if parser_function is None:
            parser_function = sequence_tokenizer
        self._initialize_fields()
        self._init_from_string(sequence, parser_function, **kwargs)

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

    def __len__(self):
        return len(self.sequence)

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
        return _calculate_mass(self)

    @property
    def glycan(self):
        return self.glycan_composition

    @glycan.setter
    def glycan(self, value):
        if isinstance(value, GlycanComposition):
            self._glycosylation_manager.aggregate = value
        elif isinstance(value, GlycosylationManager):
            self._glycosylation_manager = value
        self._invalidate()

    @property
    def glycosylation_manager(self):
        return self._glycosylation_manager

    @property
    def glycan_composition(self):
        return self._glycosylation_manager.glycan_composition

    @property
    def total_glycosylation_size(self):
        return self._glycosylation_manager.total_glycosylation_size()

    @property
    def n_term(self):
        return self._n_term

    @n_term.setter
    def n_term(self, value):
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

    def __iter__(self):
        return iter(self.sequence)

    def __getitem__(self, index):
        sub = self.sequence[index]
        return sub

    def __setitem__(self, index, value):
        self._invalidate()
        self.sequence[index] = value

    def subsequence(self, slice_obj):
        sub = self[slice_obj]
        subseq = self.from_iterable(sub)
        if slice_obj.start == 0:
            subseq.n_term = self.n_term
        if slice_obj.stop == len(self):
            subseq.c_term = self.c_term
        return subseq

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

    def _retrack_sequence(self):
        self._glycosylation_manager.clear()
        for i, position in enumerate(self):
            for mod in position[1]:
                if mod.is_tracked_for(ModificationCategory.glycosylation):
                    self._glycosylation_manager[i] = mod

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

    def deglycosylate(self):
        self._invalidate()
        _glycosylation_enum = ModificationCategory.glycosylation
        for i, pos in enumerate(self):
            mods = [mod.name for mod in pos[1] if mod.is_a(
                _glycosylation_enum)]
            for mod in mods:
                self.drop_modification(i, mod)
        self._glycosylation_manager.clear()
        self.glycan = None
        return self

    def break_at(self, idx):
        if self._fragment_index is None:
            self._build_fragment_index()
        return self._fragment_index[idx]

    def get_fragments(self, kind, chemical_shifts=None, strategy=None, **kwargs):
        """Generate fragments from this structure according to the strategy specified
        by ``strategy``, returning an iterator over the sequence of theoretical fragments.

        There may be multiple fragments for a single position, resulting in the iterator
        yielding lists of fragment objects per step.

        Parameters
        ----------
        kind : :class:`~.IonSeries` or :class:`str`
            The name of the ion series to produce fragments from, either as a string or the
            :class:`~.IonSeries` object to use
        chemical_shifts : dict, optional
            A :class:`~.Mapping` between :class:`~.AminoAcidResidue` and a listof acceptable chemical shifts to apply
            to the produced fragments containing that :class:`~.AminoAcidResidue`, to be applied combinatorially.
        strategy : :class:`~.FragmentationStrategyBase` type, optional
            The strategy type to employ when producing fragments. Defaults to :class:`~.HCDFragmentationStrategy`.
        **kwargs
            Passed to ``strategy``

        Returns
        -------
        :class:`~.FragmentationStrategyBase`
            The fragmentation iterator
        """
        if strategy is None:
            strategy = HCDFragmentationStrategy
        losses = kwargs.pop('neutral_losses', [])
        if losses and not chemical_shifts:
            chemical_shifts = losses
        if kind == stub_glycopeptide_series:
            return ([frag] for frag in self.stub_fragments(True))
        else:
            return strategy(self, kind, chemical_shifts=chemical_shifts, **kwargs)

    def drop_modification(self, position, modification_type):
        '''
        Drop a modification by name from a specific residue. If the
        position is the N-term or the C-term, the terminal modification will
        be reset to the default.

        Parameters
        ----------
        position: int
            The position of the modification to drop
        modification_type: str or Modification
            The modification to drop
        '''
        dropped_index = None
        self._invalidate()
        if position is SequenceLocation.n_term:
            self.n_term = _make_terminal_group(structure_constants.N_TERM_DEFAULT)
            return
        elif position is SequenceLocation.c_term:
            self.c_term = _make_terminal_group(structure_constants.C_TERM_DEFAULT)
            return

        for i, mod in enumerate(self.sequence[position][1]):
            if modification_type == mod.rule:
                dropped_index = i
                break
        try:
            drop_mod = self.sequence[position][1].pop(dropped_index)
            self.mass -= drop_mod.mass
            if drop_mod.is_tracked_for(ModificationCategory.glycosylation):
                self._glycosylation_manager.pop(position)
        except (IndexError, ValueError):
            raise ValueError("Modification not found! %s @ %s" %
                             (modification_type, position))

    def add_modification(self, position, modification_type):
        self._invalidate()
        if isinstance(modification_type, Modification):
            mod = modification_type
        else:
            mod = Modification(rule=modification_type)

        if position is SequenceLocation.n_term:
            self.n_term = self.n_term.modify(mod)
        elif position is SequenceLocation.c_term:
            self.c_term = self.c_term.modify(mod)
        else:
            if (position == -1) or (position >= len(self.sequence)):
                raise IndexError(
                    "Invalid modification position. %s, %s, %s" %
                    (position, str(self.sequence), modification_type))

            self.sequence[position][1].append(mod)
            self.mass += mod.mass
            if mod.is_tracked_for(ModificationCategory.glycosylation):
                self._glycosylation_manager[position] = mod

    def _build_fragment_index(self, types=tuple('bycz')):
        self._fragment_index = [[] for i in range(len(self) + 1)]
        for series in types:
            series = IonSeries(series)
            if series.direction > 0:
                g = self.get_fragments(
                    series)
                for frags in g:
                    position = self._fragment_index[frags[0].position]
                    position.append(frags)
            else:
                g = self.get_fragments(
                    series)
                for frags in g:
                    position = self._fragment_index[
                        len(self) - frags[0].position]
                    position.append(frags)

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

    def insert(self, position, residue, modifications=None):
        if modifications is None:
            modifications = []
        self._invalidate()
        self.sequence.insert(position, [residue, modifications])
        self.mass += residue.mass
        for mod in modifications:
            self.mass += mod.mass
        self._retrack_sequence()

    def delete(self, position):
        self._invalidate()
        residue, mods = self.sequence.pop(position)
        self.mass -= residue.mass
        for mod in mods:
            self.mass -= mod.mass
        self._retrack_sequence()

    def substitute(self, position, residue):
        old_residue = self.sequence[position][0]
        self.mass -= old_residue.mass
        self.mass += residue.mass
        self.sequence[position][0] = residue
        self._invalidate()
        self._retrack_sequence()

    def append(self, residue, modification=None):
        self._invalidate()
        self.mass += residue.mass
        next_pos = [residue]
        if modification is None:
            next_pos.append([])
        else:
            next_pos.append([modification])
            self.mass += modification.mass
        self.sequence.append(SequencePosition(next_pos))
        self._retrack_sequence()

    def extend(self, sequence):
        self._invalidate()
        if not isinstance(sequence, PeptideSequenceBase):
            sequence = PeptideSequence(sequence)
        self.sequence.extend(sequence.sequence)
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        self._retrack_sequence()

    def leading_extend(self, sequence):
        self._invalidate()
        if not isinstance(sequence, PeptideSequenceBase):
            sequence = PeptideSequence(sequence)
        self.sequence = sequence.sequence + self.sequence
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        self._retrack_sequence()

    @property
    def n_glycan_sequon_sites(self):
        return find_n_glycosylation_sequons(self, structure_constants.ALLOW_MODIFIED_ASPARAGINE)

    @property
    def o_glycan_sequon_sites(self):
        return find_o_glycosylation_sequons(self)

    @property
    def glycosaminoglycan_sequon_sites(self):
        return find_glycosaminoglycan_sequons(self)

    @property
    def glycosylation_sites(self):
        return self.n_glycan_sequon_sites

    def stub_fragments(self, extended=False, extended_fucosylation=False, strategy=None, **kwargs):
        if strategy is None:
            strategy = StubGlycopeptideStrategy
        return strategy(
            self, extended, extended_fucosylation=extended_fucosylation, **kwargs)

    def glycan_fragments(self, oxonium=True, all_series=False, allow_ambiguous=False,
                         include_large_glycan_fragments=True, maximum_fragment_size=5,
                         strategy=None):
        r'''
        Generate all oxonium ions for the attached glycan, and
        if `all_series` is `True`, then include the B/Y glycan
        ion ladder, with the peptide attached to the Y ion ladder.

        Parameters
        ----------
        all_series: bool
            Generate the B/Y+peptide ion ladder, otherwise just 2-residue
            pairs for all monosaccharides in :attr:`self.glycan`

        Yields
        ------
        SimpleFragment
        '''
        if strategy is None:
            strategy = OxoniumIonStrategy
        return strategy(
            self, oxonium=oxonium, all_series=all_series, allow_ambiguous=allow_ambiguous,
            include_large_glycan_fragments=include_large_glycan_fragments,
            maximum_fragment_size=maximum_fragment_size)

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


try:
    from glycopeptidepy._c.structure import sequence_methods
    PeptideSequence._init_from_parsed_string = sequence_methods._init_from_parsed_string
    PeptideSequence._init_from_components = sequence_methods._init_from_components
except ImportError:
    pass


Sequence = PeptideSequence
parse = Sequence


class NamedSequence(PeptideSequence):

    def __init__(self, name=None, sequence=None, parser_function=None, **kwargs):
        super(NamedSequence, self).__init__(
            sequence, parser_function, **kwargs)
        self.name = name

    def clone(self):
        dup = super(NamedSequence, self).clone()
        dup.name = self.name
        return dup

    def __repr__(self):
        string = super(NamedSequence, self).__str__()
        name = str(self.name)
        if name.startswith(">"):
            name = name[1:]
        return ">%s\n%s" % (name, string)


class AnnotatedSequence(NamedSequence):

    def __init__(self, name=None, sequence=None, parser_function=None, annotations=None,
                 **kwargs):
        super(AnnotatedSequence, self).__init__(
            name, sequence, parser_function=parser_function, **kwargs)
        self.annotations = self._prepare_annotations(annotations)

    @staticmethod
    def _prepare_annotations(annotation_data):
        if annotation_data:
            return dict(annotation_data)
        else:
            return dict()

    def clone(self):
        dup = super(AnnotatedSequence, self).clone()
        dup.annotations = self._prepare_annotations(self.annotations)
        return dup


class ProteinSequence(AnnotatedSequence):
    pass
