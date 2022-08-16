from collections import defaultdict, Counter
from itertools import product, combinations, combinations_with_replacement

from glypy.structure.glycan_composition import (
    FrozenMonosaccharideResidue,
    FrozenGlycanComposition,
    MonosaccharideResidue,
    SubstituentResidue,
    HashableGlycanComposition)

from glycopeptidepy.utils.collectiontools import descending_combination_counter, _AccumulatorBag
from glycopeptidepy.structure.composition import (
    Composition)

from glycopeptidepy.structure.glycan import (
    GlycosylationType,)

from glycopeptidepy.structure.fragment import (
    IonSeries,
    SimpleFragment,
    StubFragment,
    _NameTree)


from .base import FragmentationStrategyBase


from six.moves import range


oxonium_ion_series = IonSeries.oxonium_ion

_substituent_residue_cache = {}
labile_substituents = ('sulfate', 'phosphate')


def remove_labile_modifications(residue):
    has_copied = False
    labile_substituents_removed = None
    try:
        for position, substituent in residue.substituents():
            if substituent.name in labile_substituents:
                if labile_substituents_removed is None:
                    labile_substituents_removed = []
                try:
                    substituent_residue = _substituent_residue_cache[substituent.name]
                except KeyError:
                    substituent_residue = _substituent_residue_cache[substituent.name] = SubstituentResidue(substituent.name)
                labile_substituents_removed.append(substituent_residue)
                if not has_copied:
                    residue = residue.clone(
                        monosaccharide_type=MonosaccharideResidue)
                residue.drop_substituent(position, substituent)
    except AttributeError:
        pass
    return residue, labile_substituents_removed


def separate_labile_modifications(glycan_composition):
    labile_modifications = []
    new_counts = HashableGlycanComposition()
    for residue, count in glycan_composition.items():
        new_residue, separated_modifications = remove_labile_modifications(residue)
        new_counts[new_residue] += count
        if separated_modifications:
            for mod in separated_modifications:
                new_counts[mod] += 1
            labile_modifications.extend(labile_modifications)
    return new_counts, labile_modifications


class _GlycanFragmentingStrategyBase(object):

    def __init__(self, *args, **kwargs):
        super(_GlycanFragmentingStrategyBase, self).__init__(*args, **kwargs)

    def _guess_max_glycan_cleavages(self):
        glycans = self.peptide.glycosylation_manager.items()
        counter = 0
        for pos, glycan in glycans:
            rule = glycan.rule
            if not rule.is_composition and not rule.is_core:
                new_count = rule.glycan.count_branches()
                if new_count > counter:
                    counter = new_count
        return counter


class CADFragmentationStrategy(FragmentationStrategyBase, _GlycanFragmentingStrategyBase):
    """Generate glycopeptide fragments derived from glycosidic
    bond cleavages as in lower energy CAD fragmentation

    Attributes
    ----------
    max_cleavages : int
        The maximum number of glycan cleavages to allow
    """

    def __init__(self, peptide, max_cleavages=None, **kwargs):
        super(CADFragmentationStrategy, self).__init__(peptide, **kwargs)
        if max_cleavages is None:
            max_cleavages = self._guess_max_glycan_cleavages()
        self.max_cleavages = max_cleavages
        self._generator = self._build_fragments(self.max_cleavages)

    def _build_fragments(self, max_cleavages=2):
        glycans = self.peptide.glycosylation_manager.items()
        n = len(glycans)
        base_composition = self.peptide.peptide_composition()
        for i in range(1, n + 1):
            for dissociated in combinations(glycans, i):
                remaining = dict(glycans)
                for position, glycan in dissociated:
                    remaining.pop(position)
                reference_composition = base_composition.copy()
                for position, glycan in remaining.items():
                    reference_composition += glycan.rule.total_composition()

                glycan_B_ions = defaultdict(list)
                glycan_Y_ions = defaultdict(list)
                for position, glycan in dissociated:
                    for f in glycan.rule.get_fragments("BY", max_cleavages=max_cleavages):
                        if 'B' in f.series:
                            glycan_B_ions[position].append(f)
                        else:
                            glycan_Y_ions[position].append(f)
                for b_ion_set in glycan_B_ions.values():
                    for f in b_ion_set:
                        yield f
                key_order = list(glycan_Y_ions.keys())
                for y_ion_set in product(*glycan_Y_ions.values()):
                    name = ",".join("%d:%s" % (k, f.name) for k, f in zip(key_order, y_ion_set))
                    fragment_composition = reference_composition.copy()
                    is_glycosylated = False
                    for f in y_ion_set:
                        fragment_composition += f.composition
                        is_glycosylated = True
                    full_name = "peptide"
                    if name:
                        full_name = "%s+%s" % (full_name, name)
                    f = SimpleFragment(
                        name=full_name,
                        mass=fragment_composition.mass,
                        composition=fragment_composition,
                        is_glycosylated=is_glycosylated,
                        kind=IonSeries.stub_glycopeptide)
                    yield f

    def __next__(self):
        return next(self._generator)

    def next(self):
        return self.__next__()

    def __iter__(self):
        return self


class _MonosaccharideDefinitionCacher(object):
    def __init__(self, *args, **kwargs):
        self._hexose = None
        self._hexnac = None
        self._neuac = None
        self._neugc = None
        self._fucose = None
        self._dhex = None
        self._xylose = None
        self._hexa = None
        super(_MonosaccharideDefinitionCacher, self).__init__(*args, **kwargs)

    def _prepare_monosaccharide(self, name):
        return FrozenMonosaccharideResidue.from_iupac_lite(name)

    @property
    def hexose(self):
        if self._hexose is None:
            self._hexose = self._prepare_monosaccharide("Hex")
        return self._hexose

    @property
    def hexnac(self):
        if self._hexnac is None:
            self._hexnac = self._prepare_monosaccharide("HexNAc")
        return self._hexnac

    @property
    def fucose(self):
        if self._fucose is None:
            self._fucose = self._prepare_monosaccharide("Fuc")
        return self._fucose

    @property
    def dhex(self):
        if self._dhex is None:
            self._dhex = self._prepare_monosaccharide("dHex")
        return self._dhex

    @property
    def neuac(self):
        if self._neuac is None:
            self._neuac = self._prepare_monosaccharide("Neu5Ac")
        return self._neuac

    @property
    def neugc(self):
        if self._neugc is None:
            self._neugc = self._prepare_monosaccharide("Neu5Gc")
        return self._neugc

    @property
    def xylose(self):
        if self._xylose is None:
            self._xylose = self._prepare_monosaccharide("Xyl")
        return self._xylose

    @property
    def hexa(self):
        if self._hexa is None:
            self._hexa = self._prepare_monosaccharide("HexA")
        return self._hexa


_Hex2NAc = FrozenMonosaccharideResidue.from_iupac_lite("Hex2NAc")
_Glc2NAc = FrozenMonosaccharideResidue.from_iupac_lite("Glc2NAc")
_Gal2NAc = FrozenMonosaccharideResidue.from_iupac_lite("Gal2NAc")


class GlycanCompositionFragmentStrategyBase(FragmentationStrategyBase):
    def __init__(self, peptide, use_query=False, *args, **kwargs):
        super(GlycanCompositionFragmentStrategyBase, self).__init__(peptide, *args, **kwargs)
        # to be initialized closer to when the glycan composition will
        # be used. Determine whether to use the fast-path __getitem__ or
        # go through the slower but more general GlycanComposition.query
        # method to look up monosaccharides.
        self._use_query = use_query
        self._generator = None

    def glycan_composition(self):
        return self.peptide.glycan_composition

    def _guess_query_mode(self, glycan_composition):
        # these guesses will work for N-glycans and common types of mucin-type O-glycans
        # and GAG linkers
        flag = glycan_composition._getitem_fast(_Hex2NAc) +\
            glycan_composition._getitem_fast(_Glc2NAc) +\
            glycan_composition._getitem_fast(_Gal2NAc)
        return flag or self._use_query

    def count_glycosylation_type(self, glycotype):
        return self.peptide.glycosylation_manager.count_glycosylation_type(glycotype)


_HEX = FrozenMonosaccharideResidue.from_iupac_lite("Hex")
_HEXNAC = FrozenMonosaccharideResidue.from_iupac_lite("HexNAc")
_FUC = FrozenMonosaccharideResidue.from_iupac_lite("Fuc")
_DHEX = FrozenMonosaccharideResidue.from_iupac_lite("dHex")
_XYL = FrozenMonosaccharideResidue.from_iupac_lite("Xyl")
_AHEX = FrozenMonosaccharideResidue.from_iupac_lite("aHex")


class GlycanCompositionFragment(object):
    __slots__ = ("mass", "composition", "key", "is_extended", "_hash_key")

    def __init__(self, mass, composition, key, is_extended=False):
        self.mass = mass
        self.composition = composition
        self.key = key
        self.is_extended = is_extended
        self._hash_key = -1

    def __getitem__(self, key):
        if key == "mass":
            return self.mass
        elif key == "composition":
            return self.composition
        elif key == "key":
            return self.key
        elif key == 'is_extended':
            return self.is_extended
        else:
            raise KeyError(key)

    def __hash__(self):
        if self._hash_key == -1:
            self._hash_key = hash(self._get_glycan_composition())
        return self._hash_key

    def _get_glycan_composition(self):
        return _prepare_glycan_composition_from_mapping(self.key)

    def __eq__(self, other):
        return self.key == other.key and self.is_extended == other.is_extended


class _CompositionTree(object):
    def __init__(self, root=None):
        if root is None:
            root = _NameTree()
        self.root = root

    def __getitem__(self, key):
        return self.root[key]

    def build(self, pair_sequence):
        node = self.root.traverse(pair_sequence)
        if node.name is None:
            # modify the hashable glycan composition before calculating
            # its hash value
            node.name = HashableGlycanComposition()
            for k, v in pair_sequence:
                node.name._setitem_fast(k, v)
            # Force the population of the _mass and _str caches
            node.name.mass()
            str(node.name)
        return node.name


_composition_tree_root = _CompositionTree()
_composition_name_cache = dict()


def _prepare_glycan_composition_from_mapping(mapping):
    pair_sequence = mapping.items()
    return _composition_tree_root.build(pair_sequence)


class StubGlycopeptideStrategy(GlycanCompositionFragmentStrategyBase, _MonosaccharideDefinitionCacher):

    """A fragmentation strategy that generates intact peptide + glycan Y fragments from
    glycopeptides where the glycan structure is not fully known. When the structure is
    known, :class:`CADFragmentationStrategy` should be used instead

    Attributes
    ----------
    extended : bool
        Whether to fragment beyond the conserved core
    extended_fucosylation : bool
        Whether to consider multiple ``Fuc`` residues per fragment rather than one per core
    """

    def __init__(self, peptide, extended=True, use_query=False, extended_fucosylation=False, **kwargs):
        self.extended = extended
        self.extended_fucosylation = extended_fucosylation
        super(StubGlycopeptideStrategy, self).__init__(peptide, use_query)

    def __next__(self):
        if self._generator is None:
            self._generator = self.stub_fragments()
        return next(self._generator)

    def fucosylate_increment(self, shift):
        fucosylated = shift.copy()
        fucosylated['key'] = fucosylated['key'].copy()
        fucosylated['mass'] += self.fucose.mass()
        if self.compute_compositions:
            fucosylated['composition'] = fucosylated[
                'composition'] + self.fucose.total_composition()
        fucosylated['key'][_FUC] = 1
        return fucosylated

    def xylosylate_increment(self, shift):
        xylosylated = shift.copy()
        xylosylated['key'] = xylosylated['key'].copy()
        xylosylated['mass'] += self.xylose.mass()
        if self.compute_compositions:
            xylosylated['composition'] = xylosylated[
                'composition'] + self.xylose.total_composition()
        xylosylated['key'][_XYL] = 1
        return xylosylated

    def fucosylate_extended(self, shift, fucose_count):
        result = [None for i in range(fucose_count)]
        for i in range(1, fucose_count + 1):
            fucosylated = shift.copy()
            fucosylated['key'] = fucosylated['key'].copy()
            fucosylated['mass'] += self.fucose.mass()
            if self.compute_compositions:
                fucosylated['composition'] = fucosylated[
                    'composition'] + self.fucose.total_composition() * i
            fucosylated['key'][_FUC] = i
            result[i - 1] = fucosylated
        return result

    def modified_increment(self, modified, shift):
        modified = modified.copy()
        modified['key'] = modified['key'].copy()
        modified['mass'] += shift.mass()
        if self.compute_compositions:
            modified['composition'] = modified['composition'] + \
                shift.total_composition()
        modified['key'][shift] = 1
        return modified

    def _validate_glycan_composition(self, aggregate_glycosylation, glycan):
        invalid = False
        if self._use_query:
            for key, value in aggregate_glycosylation.items():
                if glycan.query(key) < value:
                    invalid = True
                    break
        else:
            for key, value in aggregate_glycosylation.items():
                if glycan[key] < value:
                    invalid = True
                    break
        return invalid

    def n_glycan_stub_fragments(self):
        glycan = self.glycan_composition()
        self._use_query = self._guess_query_mode(glycan)
        core_count = self.count_glycosylation_type(GlycosylationType.n_linked)
        per_site_shifts = []
        base_composition = self.peptide_composition()
        base_mass = base_composition.mass
        seen = set()
        for i in range(core_count):
            core_shifts = self.n_glycan_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)
        for positions in product(*per_site_shifts):
            aggregate_glycosylation = _AccumulatorBag()
            mass = base_mass
            composition = base_composition.clone()
            is_extended = False
            n_positions = len(positions)
            for site in positions:
                mass += site['mass']
                if n_positions > 1:
                    aggregate_glycosylation += (site['key'])
                else:
                    aggregate_glycosylation = site['key']
                composition += site['composition']
                is_extended |= site['is_extended']
            is_glycosylated = (mass != base_mass)
            if len(positions) > 1:
                invalid = self._validate_glycan_composition(aggregate_glycosylation, glycan)
                if invalid:
                    continue
            glycosylation = _prepare_glycan_composition_from_mapping(
                aggregate_glycosylation)
            name_key = str(glycosylation)
            try:
                name = _composition_name_cache[name_key]
            except KeyError:
                _composition_name_cache[name_key] = name = StubFragment.build_name_from_composition(
                    aggregate_glycosylation)
            if name in seen:
                continue
            seen.add(name)
            yield StubFragment(
                name=name,
                mass=mass,
                composition=composition,
                is_glycosylated=is_glycosylated,
                kind=IonSeries.stub_glycopeptide,
                glycosylation=glycosylation,
                is_extended=is_extended)

    def n_glycan_composition_fragments(self, glycan, core_count=1, iteration_count=0):
        """Generate theoretical N-glycan Y fragment compositions containing the core motif
        plus a portion of the extended branches if :attr:`extended` is used. Unless :attr:`extend

        Parameters
        ----------
        glycan : :class:`~.GlycanComposition`
            Description
        core_count : int, optional
            The total number of glycans attached to the peptide
        iteration_count : int, optional
            The core index. Unless :attr:`extended_fucosylation` is true, this will
            limit the number of Fucose residues to one per core at most, and if
            the iteration count exceeds the number of cores no Fucose will be added

        Returns
        -------
        list of dict
            The compositions corresponding to the theoretical fragments
        """
        hexose = self.hexose
        hexnac = self.hexnac

        if self._use_query:
            fucose_count = glycan.query('Fuc') + glycan.query('dHex')
            xylose_count = glycan.query('Xyl')
            hexnac_in_aggregate = glycan.query('HexNAc')
            hexose_in_aggregate = glycan.query('Hex')
        else:
            fucose_count = glycan[_FUC] + glycan[_DHEX]
            xylose_count = glycan[_XYL]
            hexnac_in_aggregate = glycan[_HEXNAC]
            hexose_in_aggregate = glycan[_HEX]

        core_shifts = []
        base_hexnac = min(hexnac_in_aggregate + 1, 3)
        for hexnac_count in range(base_hexnac):
            if hexnac_count == 0:
                shift = {
                    "mass": 0,
                    "composition": Composition(),
                    "key": {},
                    'is_extended': False,
                }
                core_shifts.append(shift)
            elif hexnac_count == 1:
                shift = {
                    "mass": (hexnac_count * hexnac.mass()),
                    "composition": hexnac_count * hexnac.total_composition(),
                    "key": {_HEXNAC: hexnac_count},
                    'is_extended': False,
                }
                core_shifts.append(shift)
                if iteration_count < fucose_count:
                    fucosylated = self.fucosylate_increment(shift)
                    core_shifts.append(fucosylated)
            elif hexnac_count == 2:
                shift = {
                    "mass": (hexnac_count * hexnac.mass()),
                    "composition": hexnac_count * hexnac.total_composition(),
                    "key": {_HEXNAC: hexnac_count},
                    'is_extended': False,
                }
                core_shifts.append(shift)

                if not self.extended_fucosylation:
                    if iteration_count < fucose_count:
                        fucosylated = self.fucosylate_increment(shift)
                        core_shifts.append(fucosylated)
                        if iteration_count < xylose_count:
                            xylosylated = self.xylosylate_increment(fucosylated)
                            core_shifts.append(xylosylated)
                elif fucose_count > 0:
                    core_shifts.extend(self.fucosylate_extended(shift, fucose_count))
                    if iteration_count < xylose_count:
                        xylosylated = self.xylosylate_increment(fucosylated)
                        core_shifts.append(xylosylated)

                if iteration_count < xylose_count:
                    xylosylated = self.xylosylate_increment(shift)
                    core_shifts.append(xylosylated)

                for hexose_count in range(1, min(hexose_in_aggregate + 1, 4)):
                    shift = {
                        "mass": (hexnac_count * hexnac.mass()) + (hexose_count * hexose.mass()),
                        "composition": (hexnac_count * hexnac.total_composition()) + (
                            hexose_count * hexose.total_composition()),
                        "key": {_HEXNAC: hexnac_count, _HEX: hexose_count},
                        'is_extended': False,
                    }
                    core_shifts.append(shift)

                    if not self.extended_fucosylation:
                        if iteration_count < fucose_count:
                            fucosylated = self.fucosylate_increment(shift)
                            core_shifts.append(fucosylated)
                            if iteration_count < xylose_count:
                                xylosylated = self.xylosylate_increment(fucosylated)
                                core_shifts.append(xylosylated)
                    elif fucose_count > 0:
                        core_shifts.extend(self.fucosylate_extended(shift, fucose_count))
                        if iteration_count < xylose_count:
                            xylosylated = self.xylosylate_increment(fucosylated)
                            core_shifts.append(xylosylated)

                    if iteration_count < xylose_count:
                        xylosylated = self.xylosylate_increment(shift)
                        core_shifts.append(xylosylated)

                    # After the core motif has been exhausted, speculatively add
                    # on the remaining core monosaccharides sequentially until
                    # exhausted.
                    if hexose_count == 3 and hexnac_in_aggregate >= 2 * core_count and self.extended:
                        for extra_hexnac_count in range(0, hexnac_in_aggregate - hexnac_count + 1):
                            if extra_hexnac_count + hexnac_count > hexnac_in_aggregate:
                                continue
                            if extra_hexnac_count > 0:
                                shift = {
                                    "mass": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                        hexose_count * hexose.mass()),
                                    "composition": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                        hexose_count * hexose.total_composition()),
                                    "key": {_HEXNAC: hexnac_count + extra_hexnac_count, _HEX: hexose_count},
                                    'is_extended': True
                                }
                                core_shifts.append(shift)

                                if not self.extended_fucosylation:
                                    if iteration_count < fucose_count:
                                        fucosylated = self.fucosylate_increment(shift)
                                        core_shifts.append(fucosylated)
                                        if iteration_count < xylose_count:
                                            xylosylated = self.xylosylate_increment(fucosylated)
                                            core_shifts.append(xylosylated)
                                elif fucose_count > 0:
                                    core_shifts.extend(self.fucosylate_extended(shift, fucose_count))

                                if iteration_count < xylose_count:
                                    xylosylated = self.xylosylate_increment(shift)
                                    core_shifts.append(xylosylated)

                            for extra_hexose_count in range(1, hexose_in_aggregate - hexose_count + 1):
                                if extra_hexose_count + hexose_count > hexose_in_aggregate:
                                    continue
                                shift = {
                                    "mass": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                        (hexose_count + extra_hexose_count) * hexose.mass()),
                                    "composition": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                        (hexose_count + extra_hexose_count) * hexose.total_composition()),
                                    "key": {_HEXNAC: hexnac_count + extra_hexnac_count, _HEX: (
                                        hexose_count + extra_hexose_count)},
                                    'is_extended': True
                                }
                                core_shifts.append(shift)

                                if not self.extended_fucosylation:
                                    if iteration_count < fucose_count:
                                        fucosylated = self.fucosylate_increment(shift)
                                        core_shifts.append(fucosylated)
                                        if iteration_count < xylose_count:
                                            xylosylated = self.xylosylate_increment(fucosylated)
                                            core_shifts.append(xylosylated)
                                elif fucose_count > 0:
                                    core_shifts.extend(self.fucosylate_extended(shift, fucose_count))
                                    if iteration_count < xylose_count:
                                        xylosylated = self.xylosylate_increment(fucosylated)
                                        core_shifts.append(xylosylated)

                                if iteration_count < xylose_count:
                                    xylosylated = self.xylosylate_increment(shift)
                                    core_shifts.append(xylosylated)
        return core_shifts

    def o_glycan_stub_fragments(self):
        glycan = self.glycan_composition()
        self._use_query = self._guess_query_mode(glycan)
        core_count = self.count_glycosylation_type(GlycosylationType.o_linked)
        per_site_shifts = []

        base_composition = self.peptide_composition()
        base_mass = base_composition.mass

        for i in range(core_count):
            core_shifts = self.o_glycan_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)
        seen = set()
        for positions in product(*per_site_shifts):
            aggregate_glycosylation = _AccumulatorBag()
            mass = base_mass
            composition = base_composition.clone()
            for site in positions:
                mass += site['mass']
                aggregate_glycosylation += (site['key'])
                composition += site['composition']
            is_glycosylated = (composition != base_composition)
            invalid = False
            if self._use_query:
                for key, value in aggregate_glycosylation.items():
                    if glycan.query(key) < value:
                        invalid = True
                        break
            else:
                for key, value in aggregate_glycosylation.items():
                    if glycan[key] < value:
                        invalid = True
                        break
            if invalid:
                continue
            name = StubFragment.build_name_from_composition(aggregate_glycosylation)
            if name in seen:
                continue
            seen.add(name)
            glycosylation = _prepare_glycan_composition_from_mapping(
                aggregate_glycosylation)
            yield StubFragment(
                name=name,
                mass=mass,
                composition=composition,
                is_glycosylated=is_glycosylated,
                kind=IonSeries.stub_glycopeptide,
                glycosylation=glycosylation)

    def o_glycan_composition_fragments(self, glycan, core_count=1, iteration_count=0):
        hexose = self.hexose
        hexnac = self.hexnac

        if self._use_query:
            fucose_count = glycan.query('Fuc') + glycan.query('dHex')
            hexnac_in_aggregate = glycan.query('HexNAc')
            hexose_in_aggregate = glycan.query('Hex')
        else:
            fucose_count = glycan[_FUC] + glycan[_DHEX]
            hexnac_in_aggregate = glycan[_HEXNAC]
            hexose_in_aggregate = glycan[_HEX]
        core_shifts = []
        for hexnac_count in range(3):
            if hexnac_in_aggregate < hexnac_count:
                continue
            if hexnac_count == 0:
                shift = {
                    "mass": 0,
                    "composition": Composition(),
                    "key": {},
                    "is_extended": False
                }
                core_shifts.append(shift)
            elif hexnac_count >= 1:
                shift = {
                    "mass": (hexnac_count * hexnac.mass()),
                    "composition": hexnac_count * hexnac.total_composition(),
                    "key": {_HEXNAC: hexnac_count},
                    "is_extended": False
                }
                core_shifts.append(shift)
                if iteration_count < fucose_count:
                    fucosylated = self.fucosylate_increment(shift)
                    core_shifts.append(fucosylated)
                for hexose_count in range(0, 2):
                    if hexose_in_aggregate < hexose_count:
                        continue
                    if hexose_count > 0:
                        shift = {
                            "mass": (
                                (hexnac_count) * hexnac.mass()) + (
                                (hexose_count) * hexose.mass()),
                            "composition": (
                                (hexnac_count) * hexnac.total_composition()) + (
                                (hexose_count) * hexose.total_composition()),
                            "key": {_HEXNAC: hexnac_count, _HEX: (
                                hexose_count)},
                            "is_extended": False
                        }
                        core_shifts.append(shift)
                        if iteration_count < fucose_count:
                            fucosylated = self.fucosylate_increment(shift)
                            core_shifts.append(fucosylated)
                    # After the core motif has been exhausted, speculatively add
                    # on the remaining core monosaccharides sequentially until
                    # exhausted.
                    if self.extended and hexnac_in_aggregate - hexnac_count >= 0:
                        for extra_hexnac_count in range(hexnac_in_aggregate - hexnac_count + 1):
                            if extra_hexnac_count:
                                shift = {
                                    "mass": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                        (hexose_count) * hexose.mass()),
                                    "composition": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                        (hexose_count) * hexose.total_composition()),
                                    "key": {_HEXNAC: hexnac_count + extra_hexnac_count, _HEX: (
                                        hexose_count)},
                                    'is_extended': True,
                                }
                                core_shifts.append(shift)
                                if iteration_count < fucose_count:
                                    fucosylated = self.fucosylate_increment(shift)
                                    core_shifts.append(fucosylated)
                            if hexose_in_aggregate > hexose_count and hexose_count > 0:
                                for extra_hexose_count in range(hexose_in_aggregate - hexose_count):
                                    if extra_hexose_count > 0 and extra_hexose_count + hexose_count > 0:
                                        shift = {
                                            "mass": (
                                                (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                                (hexose_count + extra_hexose_count) * hexose.mass()),
                                            "composition": (
                                                (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                                (hexose_count + extra_hexose_count) * hexose.total_composition()),
                                            "key": {_HEXNAC: hexnac_count + extra_hexnac_count, _HEX: (
                                                hexose_count + extra_hexose_count)},
                                            "is_extended": True
                                        }
                                        core_shifts.append(shift)
                                        if iteration_count < fucose_count:
                                            fucosylated = self.fucosylate_increment(
                                                shift)
                                            core_shifts.append(fucosylated)

        return core_shifts

    def gag_linker_stub_fragments(self):
        glycan = self.glycan_composition()
        self._use_query = self._guess_query_mode(glycan)
        core_count = self.count_glycosylation_type(GlycosylationType.glycosaminoglycan)

        per_site_shifts = []

        base_composition = self.peptide_composition()
        base_mass = base_composition.mass
        for i in range(core_count):
            core_shifts = self.gag_linker_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)
        seen = set()
        for positions in product(*per_site_shifts):
            aggregate_glycosylation = _AccumulatorBag()
            mass = base_mass
            composition = base_composition.clone()
            for site in positions:
                mass += site['mass']
                aggregate_glycosylation += (site['key'])
                composition += site['composition']
            is_glycosylated = (composition != base_composition)
            invalid = False
            if self._use_query:
                for key, value in aggregate_glycosylation.items():
                    if glycan.query(key) < value:
                        invalid = True
                        break
            else:
                for key, value in aggregate_glycosylation.items():
                    if glycan[key] < value:
                        invalid = True
                        break
            if invalid:
                continue
            name = StubFragment.build_name_from_composition(aggregate_glycosylation)
            if name in seen:
                continue
            seen.add(name)
            glycosylation = _prepare_glycan_composition_from_mapping(
                aggregate_glycosylation)
            yield StubFragment(
                name=name,
                mass=mass,
                composition=composition,
                is_glycosylated=is_glycosylated,
                kind=IonSeries.stub_glycopeptide,
                glycosylation=glycosylation)

    def gag_linker_composition_fragments(self, glycan, core_count=1, iteration_count=0):
        hexose = self.hexose
        xyl = self.xylose
        hexa = self.hexa

        xyl_in_aggregate = glycan[xyl]
        if xyl_in_aggregate == 0:
            xyl_in_aggregate = glycan.query("Pen", exact=False)
        core_shifts = []
        for xyl_count in range(max(0, xyl_in_aggregate) + 1):
            if xyl_count == 0:
                shift = {
                    "mass": 0,
                    "composition": Composition(),
                    "key": {},
                    "is_extended": False
                }
                core_shifts.append(shift)
            else:
                shift = {
                    "mass": xyl.mass() * xyl_count,
                    "composition": xyl.total_composition() * xyl_count,
                    "key": {
                        _XYL: xyl_count
                    },
                    "is_extended": False
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
                            _XYL: xyl_count,
                            _HEX: hexose_count
                        },
                        "is_extended": False
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
                                _XYL: xyl_count,
                                _HEX: hexose_count,
                                _AHEX: 1
                            },
                            "is_extended": False
                        }
                        core_shifts.append(shift)
        return core_shifts

    def stub_fragments(self):
        n_glycan = self.count_glycosylation_type(GlycosylationType.n_linked) > 0
        o_glycan = self.count_glycosylation_type(GlycosylationType.o_linked) > 0
        gag_linker = self.count_glycosylation_type(GlycosylationType.glycosaminoglycan) > 0
        if (n_glycan + o_glycan + gag_linker) > 2:
            raise ValueError(
                "Does not support mixed-type glycan fragmentation (yet)")
        if n_glycan:
            return self.n_glycan_stub_fragments()
        elif o_glycan:
            return self.o_glycan_stub_fragments()
        elif gag_linker:
            return self.gag_linker_stub_fragments()
        else:
            if len(self.peptide.glycosylation_manager) > 0:
                raise ValueError("Unknown Glycan Class Detected")
        return (a for a in [])


try:
    from glycopeptidepy._c.structure.fragmentation_strategy.glycan import (
        GlycanCompositionFragmentStrategyBase, StubGlycopeptideStrategy, GlycanCompositionFragment)
    _has_c = True
except ImportError:
    _has_c = False


_StubGlycopeptideStrategy = StubGlycopeptideStrategy


class LabileAwareStubGlycopeptideStrategy(_StubGlycopeptideStrategy):
    def __init__(self, peptide, extended=True, use_query=False, extended_fucosylation=False, detatch_substituents=False, **kwargs):
        self.detatch_substituents = detatch_substituents
        self.labile_modifications = None
        super(LabileAwareStubGlycopeptideStrategy,
              self).__init__(peptide, use_query=use_query, extended=extended, extended_fucosylation=extended_fucosylation, **kwargs)

    def glycan_composition(self):
        gc = super(LabileAwareStubGlycopeptideStrategy, self).glycan_composition()
        if self.detatch_substituents:
            gc, self.labile_modifications = separate_labile_modifications(gc)
            return gc
        else:
            return gc


labile_monosaccharides_not_to_duplicate = {
    FrozenMonosaccharideResidue.from_iupac_lite("NeuAc"),
    FrozenMonosaccharideResidue.from_iupac_lite("NeuGc"),
    FrozenMonosaccharideResidue.from_iupac_lite("Neu"),
    FrozenMonosaccharideResidue.from_iupac_lite("Fuc"),
    FrozenMonosaccharideResidue.from_iupac_lite("dHex"),
}

class OxoniumIonStrategy(GlycanCompositionFragmentStrategyBase, _MonosaccharideDefinitionCacher):

    def __init__(self, peptide, use_query=False, oxonium=True, all_series=False, allow_ambiguous=False,
                 include_large_glycan_fragments=True, maximum_fragment_size=5, *args, **kwargs):
        self.oxonium = oxonium
        self.all_series = all_series
        self.allow_ambiguous = allow_ambiguous
        self.maximum_fragment_size = maximum_fragment_size
        self.include_large_glycan_fragments = include_large_glycan_fragments
        super(OxoniumIonStrategy, self).__init__(peptide, use_query=use_query, *args, **kwargs)

    def __next__(self):
        if self._generator is None:
            self._generator = self._build_generator()
        return next(self._generator)

    def _oxonium_fragments_get_monosaccharide_list(self, glycan):
        monosaccharides = dict(glycan)
        for mono, count in list(monosaccharides.items()):
            dissociated, lost_mods = remove_labile_modifications(mono)
            if dissociated != mono:
                monosaccharides[dissociated] = count
                for mod in lost_mods:
                    pass
        return monosaccharides

    def _glycan_structural_dissociation(self, max_cleavages=2):
        return CADFragmentationStrategy(self.peptide, max_cleavages)

    def _oxonium_ions(self):
        water = Composition("H2O")
        water2 = water * 2
        side_chain_plus_carbon = Composition("CH2O")
        two_side_chains_plus_carbon = side_chain_plus_carbon * 2
        water2_plus_sidechain_plus_carbon = water2 + side_chain_plus_carbon
        water_plus_two_side_chains_plus_carbon = water + two_side_chains_plus_carbon

        glycan = (self.glycan_composition()).clone()

        monosaccharides = self._oxonium_fragments_get_monosaccharide_list(glycan)
        for k in monosaccharides:
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
                name=key + "-C2H4O2", mass=mass - two_side_chains_plus_carbon.mass,
                composition=composition - (
                    two_side_chains_plus_carbon), kind=oxonium_ion_series)
            yield SimpleFragment(
                name=key + "-CH6O3",
                mass=mass - water2_plus_sidechain_plus_carbon.mass,
                composition=composition - water2_plus_sidechain_plus_carbon,
                kind=oxonium_ion_series)
            yield SimpleFragment(
                name=key + "-C2H6O3",
                mass=mass - water_plus_two_side_chains_plus_carbon.mass,
                composition=composition - water_plus_two_side_chains_plus_carbon,
                kind=oxonium_ion_series)
        for i in range(2, self.maximum_fragment_size + 1):
            for kk in combinations_with_replacement(sorted(monosaccharides, key=str), i):
                invalid = False
                for k, v in Counter(kk).items():
                    if monosaccharides[k] < v:
                        invalid = True
                        break
                    elif k in labile_monosaccharides_not_to_duplicate and v > 1:
                        invalid = True
                        break
                if invalid:
                    continue
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

    def _ambiguous_all_series(self):
        water = Composition("H2O")
        _hexnac = self.hexnac
        _hexose = self.hexose
        _neuac = self.neuac

        if self.peptide.glycosylation_manager.count_glycosylation_type(GlycosylationType.n_linked) > 0:
            _offset = Composition()
            total = (self.glycan_composition()).clone()
            total_count = sum(total.values())

            base = FrozenGlycanComposition(Hex=3, HexNAc=2)
            remainder = total - base

            peptide_base_composition = self.total_composition() - total.total_composition()
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

                if frag_size > 2 and self.include_large_glycan_fragments and frag_size < self.maximum_fragment_size:
                    string_form = composition.serialize()
                    yield SimpleFragment(
                        name=string_form, mass=composition_mass,
                        composition=elemental_composition, kind=oxonium_ion_series)
                    yield SimpleFragment(
                        name=string_form + "-H2O", mass=composition_mass - water.mass,
                        composition=elemental_composition - water, kind=oxonium_ion_series)

                if (total_count - frag_size) < (self.maximum_fragment_size + 4):
                    f = StubFragment(
                        name="peptide+" + str(total - composition),
                        mass=stub_mass + remainder_mass - composition_mass,
                        composition=stub_composition + remainder_elemental_composition - elemental_composition,
                        is_glycosylated=True,
                        kind=IonSeries.stub_glycopeptide,
                        glycosylation=composition,
                        is_extended=None)
                    yield f

    def _build_generator(self):
        if self.oxonium:
            for f in self._oxonium_ions():
                yield f
        if self.peptide.glycosylation_manager.is_fully_specified_topologies() and self.all_series:
            for f in self._glycan_structural_dissociation():
                yield f
        elif self.allow_ambiguous and self.all_series:
            for f in self._ambiguous_all_series():
                yield f
        elif self.all_series:
            raise TypeError(
                "Cannot generate B/Y fragments from non-Glycan {}".format(self.glycan))
