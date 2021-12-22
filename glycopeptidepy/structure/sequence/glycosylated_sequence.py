from glypy import GlycanComposition

from .. import constants as structure_constants
from ..glycan import (GlycosylationManager, glycosylation_site_detectors, GlycosylationType, GlycanCompositionWithOffsetProxy)
from ..modification import ModificationCategory, Modification
from ..residue import Residue
from ..fragmentation_strategy import OxoniumIonStrategy, StubGlycopeptideStrategy
from ..fragment import (
    _n_glycosylation, _o_glycosylation,
    _gag_linker_glycosylation)


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

    for i, position in enumerate(sequence):
        if position[0] in site_set:
            if (not position.modifications or position.modifications[0] in allow_modified) and not (
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
            if next_pos.amino_acid == asn:
                if (not next_pos.modifications or next_pos.modifications[0] in allow_modified):
                    n_pos = i
                    state = "n"
        elif state == "n":
            if next_pos.amino_acid != pro:
                state = "^p"
            else:
                state = "seek"
                i = n_pos
                n_pos = None
        elif state == "^p":
            if next_pos.amino_acid in {ser, thr}:
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
    while(i < len(sequence)):
        next_pos = sequence[i]
        if state == "seek":
            if next_pos.amino_acid == ser:
                if (not next_pos.modifications or next_pos.modifications[0] in allow_modified):
                    positions.append(i)
        i += 1
    return positions


class GlycosylatedSequenceMixin(object):
    __slots__ = ()

    @property
    def glycan(self):
        return self.glycan_composition

    @glycan.setter
    def glycan(self, value):
        if isinstance(value, (GlycanComposition, GlycanCompositionWithOffsetProxy)):
            self._glycosylation_manager.aggregate = value
        elif isinstance(value, GlycosylationManager):
            self._glycosylation_manager = value # pylint: disable=assigning-non-slot
        elif value is None:
            pass
        else:
            raise TypeError("Cannot convert %s %r" % (type(value), value))
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

    def deglycosylate(self):
        self._invalidate()
        _glycosylation_enum = ModificationCategory.glycosylation
        for i, pos in enumerate(self):
            if pos.modifications:
                mods = [mod.name for mod in pos.modifications if mod.is_a(
                    _glycosylation_enum)]
                for mod in mods:
                    self.drop_modification(i, mod)
        self._glycosylation_manager.clear()
        # This line is a no-op.
        # self.glycan = None
        return self

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

    def stub_fragments(self, extended=False, extended_fucosylation=False, strategy=None, **kwargs):
        if strategy is None:
            strategy = StubGlycopeptideStrategy
        return strategy(
            self, extended, extended_fucosylation=extended_fucosylation, **kwargs)

    def is_multiply_glycosylated(self):
        return len(self.glycosylation_manager) > 1

    def is_o_glycosylated(self):
        return self.glycosylation_manager.count_glycosylation_type(GlycosylationType.o_linked)

    def is_n_glycosylated(self):
        return self.glycosylation_manager.count_glycosylation_type(GlycosylationType.n_linked)

    def is_gag_linker_glycosylated(self):
        return self.glycosylation_manager.count_glycosylation_type(GlycosylationType.glycosaminoglycan)
