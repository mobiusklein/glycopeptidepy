import itertools
from collections import Counter

from six import string_types as basestring

from . import PeptideSequenceBase, MoleculeBase
from . import constants as structure_constants

from .composition import Composition
from .fragment import (
    PeptideFragment,
    SimpleFragment, IonSeries, _n_glycosylation, _o_glycosylation,
    _gag_linker_glycosylation)
from .modification import (
    Modification, SequenceLocation, ModificationCategory,
    ModificationIndex)
from .residue import Residue

from .fragmentation_strategy import HCDFragmentationStrategy

from glypy import GlycanComposition, Glycan, ReducedEnd
from glypy.composition.glycan_composition import (
    FrozenGlycanComposition, FrozenMonosaccharideResidue,
    MonosaccharideResidue)

from .parser import sequence_tokenizer

from .glycan import (
    GlycosylationType, GlycosylationManager,
    glycosylation_site_detectors, GlycanCompositionProxy)

from ..utils.iterators import peekable
from ..utils.collectiontools import (
    descending_combination_counter, _AccumulatorBag)


b_series = IonSeries.b
y_series = IonSeries.y
oxonium_ion_series = IonSeries.oxonium_ion
stub_glycopeptide_series = IonSeries.stub_glycopeptide


def remove_labile_modifications(residue):
    has_copied = False
    try:
        for position, substituent in residue.substituents():
            if substituent.name in ("sulfate", "phosphate"):
                if not has_copied:
                    residue = residue.clone(
                        monosaccharide_type=MonosaccharideResidue)
                residue.drop_substituent(position, substituent)
    except AttributeError:
        pass
    return residue


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


class PeptideSequence(PeptideSequenceBase):
    '''
    Represents a peptide that may have post-translational modifications
    including glycosylation.

    Attributes
    ----------
    sequence: list
        The underlying container for positions in the amino acid sequence
    modification_index: ModificationIndex
        A count of different modifications attached to the amino acid sequence
    n_term: Modification or Composition
    c_term: Modification or Composition
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
    position_class = list

    @classmethod
    def from_iterable(cls, iterable):
        seq = cls(None)
        n_term = structure_constants.N_TERM_DEFAULT
        c_term = structure_constants.C_TERM_DEFAULT
        i = 0
        for pos, next_pos in peekable(iterable):
            i += 1
            try:
                resid, mods = pos
            except ValueError:
                if i == 0:
                    n_term = pos
                elif next_pos is peekable.sentinel:
                    c_term = pos
                else:
                    raise
            if not isinstance(resid, Residue):
                resid = Residue(symbol=resid)
            seq.mass += resid.mass
            mod_list = []
            for mod in mods:
                if mod == "":
                    continue
                if not isinstance(mod, Modification):
                    mod = Modification(mod)
                if mod.is_tracked_for(ModificationCategory.glycosylation):
                    seq._glycosylation_manager[i] = mod
                mod_list.append(mod)
                seq.mass += mod.mass
                seq.modification_index[mod.name] += 1
            seq.sequence.append(cls.position_class([resid, mod_list]))
        if not isinstance(n_term, MoleculeBase):
            n_term = Modification(n_term)
        if not isinstance(c_term, MoleculeBase):
            c_term = Modification(c_term)

        seq.n_term = n_term
        seq.c_term = c_term
        return seq

    def __init__(self, sequence=None, parser_function=None, **kwargs):
        if parser_function is None:
            parser_function = sequence_tokenizer
        self._mass = 0.0
        self.sequence = []
        self.modification_index = ModificationIndex()
        self._fragment_index = None
        self._glycan = None
        self._glycan_composition = None

        self._glycosylation_manager = GlycosylationManager(self)

        self._n_term = None
        self._c_term = None

        self._fragments_map = {}
        self._total_composition = None
        self._peptide_composition = None

        if sequence == "" or sequence is None:
            pass
        else:
            seq_list, modifications, glycan, n_term, c_term = parser_function(
                sequence)
            i = 0
            for item in seq_list:
                res = Residue(item[0])
                self.mass += res.mass
                mods = []
                for mod in item[1]:
                    if mod != '':
                        mod = Modification(mod)
                        if mod.is_tracked_for(ModificationCategory.glycosylation):
                            self._glycosylation_manager[i] = mod
                        mods.append(mod)
                        self.modification_index[mod] += 1
                        self.mass += mod.mass
                i += 1
                self.sequence.append(self.position_class([res, mods]))

            if glycan != "":
                self._glycosylation_manager.aggregate = glycan.clone()

            self.glycan = glycan if glycan != "" else None
            self.n_term = Modification(n_term) if isinstance(
                n_term, basestring) else n_term
            self.c_term = Modification(c_term) if isinstance(
                c_term, basestring) else c_term

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
        rep = "{n_term}{sequence}{c_term}{glycan}[{_mass}]".format(
            n_term=n_term, c_term=c_term,
            glycan=aggregate,
            **self.__dict__)
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
        # if isinstance(value, GlycanComposition):
        #     self._patch_glycan_composition()
        # elif isinstance(value, Glycan):
        #     value.reducing_end = ReducedEnd(-Composition("H2O"), valence=0)
        # self._glycan_composition = None

    @property
    def glycan_composition(self):
        return self._glycosylation_manager.glycan_composition

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
        subseq = Sequence.from_iterable(sub)
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

    def get_fragments(self, kind, neutral_losses=None, **kwargs):
        """Return a list of mass values for each fragment of `kind`"""
        if kind == stub_glycopeptide_series:
            for frag in self.stub_fragments(True):
                yield [frag]
            raise StopIteration()
        else:
            for frag in HCDFragmentationStrategy(self, kind):
                yield frag

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
            self.n_term = Modification("H")
            return
        elif position is SequenceLocation.c_term:
            self.c_term = Modification("OH")
            return

        for i, mod in enumerate(self.sequence[position][1]):
            if modification_type == mod.rule:
                dropped_index = i
                break
        try:
            drop_mod = self.sequence[position][1].pop(dropped_index)
            self.mass -= drop_mod.mass
            self.modification_index[drop_mod.name] -= 1
            if drop_mod.is_tracked_for(ModificationCategory.glycosylation):
                self._glycosylation_manager.pop(position)
        except (IndexError, ValueError):
            raise ValueError("Modification not found! %s @ %s" %
                             (modification_type, position))

    def add_modification(self, position=None, modification_type=None):
        self._invalidate()
        if position is None and isinstance(modification_type, Modification):
            position = modification_type.position
        if isinstance(modification_type, Modification):
            mod = modification_type
        else:
            mod = Modification(rule=modification_type, mod_pos=position)

        if position is SequenceLocation.n_term:
            self.n_term = mod
        elif position is SequenceLocation.c_term:
            self.c_term = mod
        else:
            if (position == -1) or (position >= len(self.sequence)):
                raise IndexError(
                    "Invalid modification position. %s, %s, %s" %
                    (position, str(self.sequence), modification_type))

            self.sequence[position][1].append(mod)
            self.mass += mod.mass
            self.modification_index[mod.name] += 1
            if mod.is_tracked_for(ModificationCategory.glycosylation):
                self._glycosylation_manager[position] = mod

    def fragment(self, key):
        try:
            return self._fragments_map[key]
        except KeyError:
            for group in self.get_fragments(key[0]):
                for frag in group:
                    self._fragments_map[frag.name] = frag
            try:
                return self._fragments_map[key]
            except KeyError:
                raise KeyError("Unknown Fragment %r" % (key,))

    def _build_fragment_index(self, types=tuple('by'), neutral_losses=None):
        self._fragment_index = [[] for i in range(len(self) + 1)]
        for series in types:
            series = IonSeries(series)
            if series.direction > 0:
                g = self.get_fragments(
                    series, neutral_losses=neutral_losses)
                for frags in g:
                    position = self._fragment_index[frags[0].position]
                    position.append(frags)
            else:
                g = self.get_fragments(
                    series, neutral_losses=neutral_losses)
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
        inst = self.__class__(sequence=str(self))
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
            self.modification_index[modification.name] += 1
        self.sequence.append(self.position_class(next_pos))
        self._retrack_sequence()

    def extend(self, sequence):
        self._invalidate()
        if not isinstance(sequence, PeptideSequenceBase):
            sequence = PeptideSequence(sequence)
        self.sequence.extend(sequence.sequence)
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        for mod, count in sequence.modification_index.items():
            self.modification_index[mod] += count
        self._retrack_sequence()

    def leading_extend(self, sequence):
        self._invalidate()
        if not isinstance(sequence, PeptideSequenceBase):
            sequence = PeptideSequence(sequence)
        self.sequence = sequence.sequence + self.sequence
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        for mod, count in sequence.modification_index.items():
            self.modification_index[mod] += count
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

    def stub_fragments(self, extended=False):
        n_glycan = self.modification_index[_n_glycosylation] > 0
        o_glycan = self.modification_index[_o_glycosylation] > 0
        gag_linker = self.modification_index[_gag_linker_glycosylation] > 0
        if sum((n_glycan, o_glycan, gag_linker)) > 2:
            raise ValueError(
                "Does not support mixed-type glycan fragmentation (yet)")
        if n_glycan:
            return self.n_glycan_stub_fragments(extended=extended)
        elif o_glycan:
            return self.o_glycan_stub_fragments(extended=extended)
        elif gag_linker:
            return self.gag_linker_stub_fragments(extended=extended)
        else:
            raise ValueError("No Glycan Class Detected")

    def n_glycan_stub_fragments(self, extended=False):
        glycan = self.glycan_composition
        fucose_count = glycan['Fuc'] or glycan['dHex']
        hexnac_in_aggregate = glycan['HexNAc']
        hexose_in_aggregate = glycan["Hex"]
        core_count = self.modification_index[_n_glycosylation]

        per_site_shifts = []
        hexose = FrozenMonosaccharideResidue.from_iupac_lite("Hex")
        hexnac = FrozenMonosaccharideResidue.from_iupac_lite("HexNAc")
        fucose = FrozenMonosaccharideResidue.from_iupac_lite("Fuc")
        base_composition = self.peptide_composition()
        base_mass = base_composition.mass

        def fucosylate_increment(shift):
            fucosylated = shift.copy()
            fucosylated['key'] = fucosylated['key'].copy()
            fucosylated['mass'] += fucose.mass()
            fucosylated['composition'] = fucosylated[
                'composition'] + fucose.total_composition()
            fucosylated['key']["Fuc"] = 1
            return fucosylated

        for i in range(core_count):
            core_shifts = []
            for hexnac_count in range(min(hexnac_in_aggregate + 1, 3)):
                if hexnac_count == 0:
                    shift = {
                        "mass": 0,
                        "composition": Composition(),
                        "key": {}
                    }
                    core_shifts.append(shift)
                elif hexnac_count == 1:
                    shift = {
                        "mass": (hexnac_count * hexnac.mass()),
                        "composition": hexnac_count * hexnac.total_composition(),
                        "key": {"HexNAc": hexnac_count}
                    }
                    core_shifts.append(shift)
                    if i < fucose_count:
                        fucosylated = fucosylate_increment(shift)
                        core_shifts.append(fucosylated)
                elif hexnac_count == 2:
                    shift = {
                        "mass": (hexnac_count * hexnac.mass()),
                        "composition": hexnac_count * hexnac.total_composition(),
                        "key": {"HexNAc": hexnac_count}
                    }
                    core_shifts.append(shift)

                    if i < fucose_count:
                        fucosylated = fucosylate_increment(shift)
                        core_shifts.append(fucosylated)

                    for hexose_count in range(1, min(hexose_in_aggregate + 1, 4)):
                        shift = {
                            "mass": (hexnac_count * hexnac.mass()) + (hexose_count * hexose.mass()),
                            "composition": (hexnac_count * hexnac.total_composition()) + (
                                hexose_count * hexose.total_composition()),
                            "key": {"HexNAc": hexnac_count, "Hex": hexose_count}
                        }
                        core_shifts.append(shift)
                        if i < fucose_count:
                            fucosylated = fucosylate_increment(shift)
                            core_shifts.append(fucosylated)

                        # After the core motif has been exhausted, speculatively add
                        # on the remaining core monosaccharides sequentially until
                        # exhausted.
                        if hexose_count == 3 and hexnac_in_aggregate >= 2 * core_count and extended:
                            for extra_hexnac_count in range(0, 3):
                                if extra_hexnac_count + hexnac_count > hexnac_in_aggregate:
                                    continue
                                shift = {
                                    "mass": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                        hexose_count * hexose.mass()),
                                    "composition": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                        hexose_count * hexose.total_composition()),
                                    "key": {"HexNAc": hexnac_count + extra_hexnac_count, "Hex": hexose_count}
                                }
                                core_shifts.append(shift)
                                if i < fucose_count:
                                    fucosylated = fucosylate_increment(shift)
                                    core_shifts.append(fucosylated)
                                for extra_hexose_count in range(1, 3):
                                    if extra_hexose_count + hexose_count > hexose_in_aggregate:
                                        continue
                                    shift = {
                                        "mass": (
                                            (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                            (hexose_count + extra_hexose_count) * hexose.mass()),
                                        "composition": (
                                            (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                            (hexose_count + extra_hexose_count) * hexose.total_composition()),
                                        "key": {"HexNAc": hexnac_count + extra_hexnac_count, "Hex": (
                                            hexose_count + extra_hexose_count)}
                                    }
                                    core_shifts.append(shift)
                                    if i < fucose_count:
                                        fucosylated = fucosylate_increment(
                                            shift)
                                        core_shifts.append(fucosylated)
            per_site_shifts.append(core_shifts)
        for positions in itertools.product(*per_site_shifts):
            key_base = 'peptide'
            names = _AccumulatorBag()
            mass = base_mass
            composition = base_composition.clone()
            for site in positions:
                mass += site['mass']
                names += (site['key'])
                composition += site['composition']
            invalid = False
            for key, value in names.items():
                if glycan[key] < value:
                    invalid = True
                    break
            if invalid:
                continue
            extended_key = ''.join("%s%d" % kv for kv in sorted(names.items()))
            if len(extended_key) > 0:
                key_base = "%s+%s" % (key_base, extended_key)
            yield SimpleFragment(name=key_base, mass=mass, composition=composition, kind=stub_glycopeptide_series)

    def o_glycan_stub_fragments(self, extended=False):
        glycan = self.glycan_composition
        fucose_count = glycan['Fuc'] or glycan['dHex']
        core_count = self.modification_index[_o_glycosylation]

        per_site_shifts = []
        hexose = FrozenMonosaccharideResidue.from_iupac_lite("Hex")
        hexnac = FrozenMonosaccharideResidue.from_iupac_lite("HexNAc")
        fucose = FrozenMonosaccharideResidue.from_iupac_lite("Fuc")
        base_composition = self.peptide_composition()
        base_mass = base_composition.mass

        hexnac_in_aggregate = glycan[hexnac]
        hexose_in_aggregate = glycan[hexose]

        def fucosylate_increment(shift):
            fucosylated = shift.copy()
            fucosylated['key'] = fucosylated['key'].copy()
            fucosylated['mass'] += fucose.mass()
            fucosylated['composition'] = fucosylated[
                'composition'] + fucose.total_composition()
            fucosylated['key']["Fuc"] = 1
            return fucosylated

        for i in range(core_count):
            core_shifts = []
            for hexnac_count in range(3):
                if hexnac_in_aggregate < hexnac_count:
                    continue
                if hexnac_count == 0:
                    shift = {
                        "mass": 0,
                        "composition": Composition(),
                        "key": {}
                    }
                    core_shifts.append(shift)
                elif hexnac_count >= 1:
                    shift = {
                        "mass": (hexnac_count * hexnac.mass()),
                        "composition": hexnac_count * hexnac.total_composition(),
                        "key": {"HexNAc": hexnac_count}
                    }
                    core_shifts.append(shift)
                    if i < fucose_count:
                        fucosylated = fucosylate_increment(shift)
                        core_shifts.append(fucosylated)
                    for hexose_count in range(0, 2):
                        if hexose_in_aggregate < hexose_count:
                            continue
                        shift = {
                            "mass": (
                                (hexnac_count) * hexnac.mass()) + (
                                (hexose_count) * hexose.mass()),
                            "composition": (
                                (hexnac_count) * hexnac.total_composition()) + (
                                (hexose_count) * hexose.total_composition()),
                            "key": {"HexNAc": hexnac_count, "Hex": (
                                hexose_count)}
                        }
                        if hexose_count > 0:
                            core_shifts.append(shift)
                            if i < fucose_count:
                                fucosylated = fucosylate_increment(shift)
                                core_shifts.append(fucosylated)
                        # After the core motif has been exhausted, speculatively add
                        # on the remaining core monosaccharides sequentially until
                        # exhausted.
                        if extended and hexnac_in_aggregate - hexnac_count >= 0:
                            for extra_hexnac_count in range(hexnac_in_aggregate - hexnac_count):
                                shift = {
                                    "mass": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                        (hexose_count) * hexose.mass()),
                                    "composition": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                        (hexose_count) * hexose.total_composition()),
                                    "key": {"HexNAc": hexnac_count + extra_hexnac_count, "Hex": (
                                        hexose_count)}
                                }
                                if i < fucose_count:
                                    fucosylated = fucosylate_increment(shift)
                                    core_shifts.append(fucosylated)
                                if hexose_in_aggregate > hexose_count:
                                    for extra_hexose_count in range(hexose_in_aggregate - hexose_count):
                                        shift = {
                                            "mass": (
                                                (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                                (hexose_count + extra_hexose_count) * hexose.mass()),
                                            "composition": (
                                                (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                                (hexose_count + extra_hexose_count) * hexose.total_composition()),
                                            "key": {"HexNAc": hexnac_count + extra_hexnac_count, "Hex": (
                                                hexose_count + extra_hexose_count)}
                                        }
                                        if i < fucose_count:
                                            fucosylated = fucosylate_increment(
                                                shift)
                                            core_shifts.append(fucosylated)

            per_site_shifts.append(core_shifts)
        seen = set()
        for positions in itertools.product(*per_site_shifts):
            key_base = 'peptide'
            names = _AccumulatorBag()
            mass = base_mass
            composition = base_composition.clone()
            for site in positions:
                mass += site['mass']
                names += (site['key'])
                composition += site['composition']
            invalid = False
            for key, value in names.items():
                if glycan[key] < value:
                    invalid = True
                    break
            if invalid:
                continue

            extended_key = ''.join(
                "%s%d" % kv for kv in sorted(names.items()) if kv[1] > 0)
            if len(extended_key) > 0:
                key_base = "%s+%s" % (key_base, extended_key)
            if key_base in seen:
                continue
            seen.add(key_base)
            yield SimpleFragment(name=key_base, mass=mass, composition=composition, kind=stub_glycopeptide_series)

    def gag_linker_stub_fragments(self, extended=False):
        glycan = self.glycan_composition
        # fucose_count = glycan['Fuc'] or glycan['dHex']
        core_count = self.modification_index[_gag_linker_glycosylation]

        per_site_shifts = []
        hexose = FrozenMonosaccharideResidue.from_iupac_lite("Hex")
        # hexnac = FrozenMonosaccharideResidue.from_iupac_lite("HexNAc")
        # fucose = FrozenMonosaccharideResidue.from_iupac_lite("Fuc")
        xyl = FrozenMonosaccharideResidue.from_iupac_lite("Xyl")
        hexa = FrozenMonosaccharideResidue.from_iupac_lite("HexA")

        base_composition = self.peptide_composition()
        base_mass = base_composition.mass
        for i in range(core_count):
            core_shifts = []
            for xyl_count in range(2):
                if xyl_count == 0:
                    shift = {
                        "mass": 0,
                        "composition": Composition(),
                        "key": {}
                    }
                    core_shifts.append(shift)
                else:
                    shift = {
                        "mass": xyl.mass() * xyl_count,
                        "composition": xyl.total_composition() * xyl_count,
                        "key": {
                            xyl: xyl_count
                        }
                    }
                    core_shifts.append(shift)
                if xyl_count > 0:
                    # TODO: Handle modified Hexose residues here too.
                    for hexose_count in range(1, 3):
                        shift = {
                            "mass": ((xyl.mass() * xyl_count) + (
                                hexose.mass() * hexose_count)),
                            "composition": (
                                (xyl.total_composition() * xyl_count) + (
                                    hexose.total_composition() * hexose_count)),
                            "key": {
                                xyl: xyl_count,
                                hexose: hexose_count
                            }
                        }
                        core_shifts.append(shift)
                        if hexose_count == 2:
                            shift = {
                                "mass": ((xyl.mass() * xyl_count) + (
                                    hexose.mass() * hexose_count) + hexa.mass()),
                                "composition": (
                                    (xyl.total_composition() * xyl_count) + (
                                        hexose.total_composition() * hexose_count) + hexa.total_composition()),
                                "key": {
                                    xyl: xyl_count,
                                    hexose: hexose_count,
                                    hexa: 1
                                }
                            }
                            core_shifts.append(shift)
            per_site_shifts.append(core_shifts)
        seen = set()
        for positions in itertools.product(*per_site_shifts):
            key_base = 'peptide'
            names = _AccumulatorBag()
            mass = base_mass
            composition = base_composition.clone()
            for site in positions:
                mass += site['mass']
                names += (site['key'])
                composition += site['composition']
            invalid = False
            for key, value in names.items():
                if glycan[key] < value:
                    invalid = True
                    break
            if invalid:
                continue
            extended_key = ''.join("%s%d" % kv for kv in sorted(names.items()))
            if len(extended_key) > 0:
                key_base = "%s+%s" % (key_base, extended_key)
            if key_base in seen:
                continue
            seen.add(key_base)
            yield SimpleFragment(name=key_base, mass=mass, composition=composition, kind=stub_glycopeptide_series)

    def glycan_fragments(self, oxonium=True, all_series=False, allow_ambiguous=False,
                         include_large_glycan_fragments=True, maximum_fragment_size=5):
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
        water = Composition("H2O")
        water2 = water * 2
        side_chain_plus_carbon = Composition("CH2O")
        water2_plus_sidechain_plus_carbon = water2 + side_chain_plus_carbon
        _hexnac = FrozenMonosaccharideResidue.from_iupac_lite("HexNAc")
        _hexose = FrozenMonosaccharideResidue.from_iupac_lite("Hex")
        _neuac = FrozenMonosaccharideResidue.from_iupac_lite("NeuAc")

        if oxonium:
            glycan = (self.glycan_composition).clone()
            for k in glycan:
                k = remove_labile_modifications(k)
                key = str(k)
                mass = k.mass()
                composition = k.total_composition()
                yield SimpleFragment(
                    name=key, mass=mass,
                    composition=composition,
                    kind=oxonium_ion_series)
                yield SimpleFragment(
                    name=key + "-H2O", mass=mass - water.mass,
                    composition=composition - water,
                    kind=oxonium_ion_series)
                yield SimpleFragment(
                    name=key + "-H4O2", mass=mass - water2.mass,
                    composition=composition - (
                        water2), kind=oxonium_ion_series)
                yield SimpleFragment(
                    name=key + "-CH6O3",
                    mass=mass - water2_plus_sidechain_plus_carbon.mass,
                    composition=composition - water2_plus_sidechain_plus_carbon,
                    kind=oxonium_ion_series)
            for i in range(2, 4):
                for kk in itertools.combinations_with_replacement(sorted(glycan, key=str), i):
                    invalid = False
                    for k, v in Counter(kk).items():
                        if glycan[k] < v:
                            invalid = True
                            break
                    if invalid:
                        continue
                    kk = list(map(remove_labile_modifications, kk))
                    key = ''.join(map(str, kk))
                    mass = sum(k.mass() for k in kk)
                    composition = sum((k.total_composition()
                                       for k in kk), Composition())
                    yield SimpleFragment(
                        name=key, mass=mass, kind=oxonium_ion_series, composition=composition)
                    yield SimpleFragment(
                        name=key + "-H2O", mass=mass - water.mass, kind=oxonium_ion_series,
                        composition=composition - water)
                    yield SimpleFragment(
                        name=key + "-H4O2", mass=mass - water2.mass, kind=oxonium_ion_series,
                        composition=composition - (water2))
                    yield SimpleFragment(
                        name=key + "-CH6O3", mass=mass - water2_plus_sidechain_plus_carbon.mass,
                        kind=oxonium_ion_series,
                        composition=composition - water2_plus_sidechain_plus_carbon)

        if isinstance(self.glycan, Glycan) and all_series:
            glycan = self.glycan
            base_composition = self.total_composition() - self.glycan.total_composition()
            base_mass = base_composition.mass

            for fragment in glycan.fragments("BY"):
                if fragment.is_reducing():
                    # TODO:
                    # When self.glycan is adjusted for the attachment cost with the anchoring
                    # amino acid, this water.mass penalty can be removed
                    yield SimpleFragment(
                        name="peptide+" + fragment.name, mass=base_mass + fragment.mass - water.mass,
                        kind=stub_glycopeptide_series, composition=base_composition + fragment.composition)
                else:
                    yield SimpleFragment(
                        name=fragment.name, mass=fragment.mass, kind=oxonium_ion_series,
                        composition=fragment.composition)
        elif allow_ambiguous and all_series:
            if self.modification_index[_n_glycosylation] > 0:
                _offset = Composition()
                total = (self.glycan_composition).clone()
                total_count = sum(total.values())

                base = FrozenGlycanComposition(Hex=3, HexNAc=2)
                remainder = total - base

                peptide_base_composition = self.total_composition() - self.glycan_composition.total_composition()
                stub_composition = peptide_base_composition + base.total_composition() - water
                stub_mass = stub_composition.mass

                # GlycanComposition's clone semantics do not propagate the
                # composition_offset attribute yet. Should it?
                remainder.composition_offset = _offset
                remainder_elemental_composition = remainder.total_composition()
                remainder_mass = remainder.mass()

                for composition in descending_combination_counter(remainder):
                    frag_size = sum(composition.values())

                    # Don't waste time on compositions that are impossible under
                    # common enzymatic pathways in this already questionable
                    # stop-gap
                    if composition.get(_hexnac, 0) + composition.get(
                            _hexose, 0) < composition.get(_neuac, 0):
                        continue

                    composition = FrozenGlycanComposition(composition)
                    composition.composition_offset = _offset

                    elemental_composition = composition.total_composition()
                    composition_mass = elemental_composition.mass

                    if frag_size > 2 and include_large_glycan_fragments and frag_size < maximum_fragment_size:
                        string_form = composition.serialize()
                        yield SimpleFragment(
                            name=string_form, mass=composition_mass,
                            composition=elemental_composition, kind=oxonium_ion_series)
                        yield SimpleFragment(
                            name=string_form + "-H2O", mass=composition_mass - water.mass,
                            composition=elemental_composition - water, kind=oxonium_ion_series)

                    if (total_count - frag_size) < (maximum_fragment_size + 4):
                        f = SimpleFragment(
                            name="peptide+" + str(total - composition),
                            mass=stub_mass + remainder_mass - composition_mass,
                            composition=stub_composition +
                            remainder_elemental_composition - elemental_composition,
                            kind=stub_glycopeptide_series)
                        yield f
        elif all_series:
            raise TypeError(
                "Cannot generate B/Y fragments from non-Glycan {}".format(self.glycan))

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
        return ">%s\n%s" % (self.name, string)


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
