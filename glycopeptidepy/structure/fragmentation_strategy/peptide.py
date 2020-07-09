# -*- coding: utf-8 -*-
from collections import defaultdict
from itertools import product, combinations

from glypy.structure.glycan_composition import FrozenGlycanComposition

from glycopeptidepy.utils.collectiontools import descending_combination_counter, _AccumulatorBag
from glycopeptidepy.structure import constants as structure_constants
from glycopeptidepy.structure.composition import (
    Composition)
from glycopeptidepy.structure.modification import (
    ModificationCategory,
    GlycanFragment,
    ModificationIndex)
from glycopeptidepy.structure.glycan import (
    GlycosylationManager,
    GlycosylationType)
from glycopeptidepy.structure.fragment import (
    IonSeries,
    SimpleFragment,
    PeptideFragment,
    format_negative_composition,
    ChemicalShift)
from glycopeptidepy.structure.residue import (
    AminoAcidResidue as R,
    residue_to_neutral_loss as _residue_to_neutral_loss)


from .base import (
    FragmentationStrategyBase,
    glycosylation_type_to_core,
    _n_glycosylation,
    _o_glycosylation,
    _gag_linker_glycosylation,
    _modification_xylose,
    _modification_hexnac)

from .glycan import (
    _GlycanFragmentingStrategyBase, StubGlycopeptideStrategy)


ammonia_loss = ChemicalShift("-NH3", -Composition({"N": 1, "H": 3}))
water_loss = ChemicalShift("-H2O", -Composition({"O": 1, "H": 2}))


class PeptideFragmentationStrategyBase(FragmentationStrategyBase):

    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1, include_neutral_losses=False, **kwargs):
        if chemical_shifts is None:
            chemical_shifts = dict()
        super(PeptideFragmentationStrategyBase,
              self).__init__(peptide, **kwargs)
        self.series = IonSeries(series)
        self.chemical_shift_rules = defaultdict(list, chemical_shifts)
        self.max_chemical_shifts = max_chemical_shifts
        self.include_neutral_losses = include_neutral_losses
        self._initialize_fields()

    def _initialize_fields(self):
        self.direction = self.series.direction

        # Null Values
        self.running_mass = 0
        if self.compute_compositions:
            self.running_composition = Composition()
        else:
            self.running_composition = None
        self.index = None
        self.size = len(self.peptide)
        self.modification_index = ModificationIndex()
        self.glycosylation_manager = GlycosylationManager(self.peptide)
        self.amino_acids_counter = _AccumulatorBag()

        self.running_mass += self.series.mass_shift
        if self.compute_compositions:
            self.running_composition += self.series.composition_shift
        self._initialize_start_terminal()

    def _get_viable_chemical_shift_combinations(self):
        shifts = []
        if self.include_neutral_losses:
            shifts.append(ammonia_loss)
            shifts.append(water_loss)
        if not self.chemical_shift_rules:
            return shifts
        for residue, count in self.amino_acids_counter.items():
            loss_composition = self.chemical_shift_rules[residue]
            if loss_composition:
                for i in range(count):
                    shifts.append(loss_composition + [Composition()])
        loss_composition_combinations = {}
        for comb in combinations(shifts, self.max_chemical_shifts):
            for prod in product(*comb):
                loss_composition = sum(prod, Composition())
                loss_composition_combinations[format_negative_composition(loss_composition)] = loss_composition
        return [ChemicalShift(k, v) for k, v in loss_composition_combinations.items() if v]

    def _initialize_start_terminal(self):
        if self.direction > 0:
            self.running_mass += self.peptide.n_term.mass
            if self.compute_compositions:
                self.running_composition += self.peptide.n_term.composition
            self.index = -1
        elif self.direction < 0:
            self.running_mass += self.peptide.c_term.mass
            if self.compute_compositions:
                self.running_composition += self.peptide.c_term.composition
            self.index = len(self.peptide)
        else:
            raise ValueError("Unknown direction %r" % (self.series.direction,))

    def composition_of(self, residue, modifications):
        composition = Composition(residue.composition)
        for mod in modifications:
            composition += mod.composition
        return composition

    def has_more(self):
        if self.direction > 0:
            return self.index < self.size - 2
        else:
            return self.index > 1

    def flanking_residues(self):
        residues = [self.peptide[self.index][0],
                    self.peptide[self.index + self.direction][0]]
        if self.direction < 0:
            residues = residues[::-1]
        return residues

    def name_index_of(self):
        if self.direction > 0:
            return self.index + structure_constants.FRAG_OFFSET
        else:
            return (self.size - self.index - 1) + structure_constants.FRAG_OFFSET

    def __next__(self):
        if self.has_more():
            return self.step()
        else:
            raise StopIteration()

    def track_glycosylation(self, index, glycosylation):
        self.glycosylation_manager[self.index] = glycosylation
        self.modification_index[glycosylation] += 1

    def _update_state(self):
        self.index += self.direction
        residue, modifications = self.peptide[self.index]
        for mod in modifications:
            if mod.is_tracked_for(ModificationCategory.glycosylation):
                self.track_glycosylation(self.index, mod)
            else:
                self.modification_index[mod] += 1
        self.amino_acids_counter[residue] += 1
        self.running_mass += residue.mass
        if self.compute_compositions:
            composition = self.composition_of(residue, modifications)
            self.running_composition += composition

    def _build_fragments(self):
        fragments_from_site = []
        frag = PeptideFragment(
            self.series,
            self.name_index_of(),
            self.modification_index.copy(),
            self.running_mass,
            glycosylation=self.glycosylation_manager.copy(),
            flanking_amino_acids=self.flanking_residues(),
            composition=self.running_composition)
        shifts = self._get_viable_chemical_shift_combinations()
        for fragment in self.partial_loss(frag):
            fragments_from_site.append(fragment)
            for shift in shifts:
                f = fragment.clone()
                f.chemical_shift = shift
                fragments_from_site.append(f)
        return fragments_from_site

    def step(self):
        self._update_state()
        fragments_from_site = self._build_fragments()
        return fragments_from_site

    def __repr__(self):
        return "%s(%s, %r, %r, %0.4f, %d)" % (
            self.__class__.__name__,
            self.peptide, self.series, self.running_composition,
            self.running_mass, self.index)


class HCDFragmentationStrategy(PeptideFragmentationStrategyBase):
    modifications_of_interest = {k.name: k for k in ([
        _n_glycosylation,
        _modification_hexnac,
        _o_glycosylation,
        _gag_linker_glycosylation,
        _modification_xylose
    ])}

    modification_compositions = {
        k.name: k.composition for k in modifications_of_interest.values()
    }

    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1, include_neutral_losses=False, **kwargs):
        super(HCDFragmentationStrategy, self).__init__(
            peptide, series, chemical_shifts, max_chemical_shifts, include_neutral_losses, **kwargs)
        self._last_modification_set = None
        self._last_modification_variants = None

    def _get_core_for(self, glycosylation):
        try:
            if not glycosylation.rule.is_core:
                glycosylation = glycosylation_type_to_core[glycosylation.rule.glycosylation_type]()
            return glycosylation
        except KeyError:
            raise ValueError("Cannot determine which core to use for {}".format(
                glycosylation.rule.glycosylation_type))

    def composition_of(self, residue, modifications):
        composition = Composition(residue.composition)
        for mod in modifications:
            if mod.is_tracked_for(ModificationCategory.glycosylation):
                mod = self._get_core_for(mod)
            composition += mod.composition
        return composition

    def track_glycosylation(self, index, glycosylation):
        # HCD strategy does not track intact topologies
        try:
            glycosylation = self._get_core_for(glycosylation)
        except KeyError:
            raise ValueError("Cannot determine which core to use for {}".format(
                glycosylation.rule.glycosylation_type))
        self.glycosylation_manager[self.index] = glycosylation
        self.modification_index[glycosylation] += 1

    def _get_modifications_of_interest(self, fragment):
        modifications = fragment.modification_dict
        if self.compute_compositions:
            delta_composition = Composition()
        else:
            delta_composition = None
        other_modifications = ModificationIndex()
        modifications_of_interest = ModificationIndex()

        for k, v in modifications.items():
            if k.name in self.modifications_of_interest:
                modifications_of_interest[k] = v
                if self.compute_compositions:
                    delta_composition += self.modification_compositions[k.name] * v
            else:
                other_modifications[k] = v

        return modifications_of_interest, delta_composition, other_modifications

    def _replace_cores(self, modifications_of_interest):
        n_cores = modifications_of_interest.pop(_n_glycosylation, 0)
        o_cores = modifications_of_interest.pop(_o_glycosylation, 0)
        gag_cores = modifications_of_interest.pop(_gag_linker_glycosylation, 0)

        # Convert core glycosylation into the residual monosaccharides remaining
        # after dissociation
        hexnac_cores = n_cores * 2 + o_cores
        if hexnac_cores:
            modifications_of_interest[_modification_hexnac] += hexnac_cores
        if gag_cores:
            modifications_of_interest[_modification_xylose] += gag_cores
        return modifications_of_interest

    def _generate_modification_variants(self, interesting_modifications, other_modifications):
        variants = []
        for varied_modifications in descending_combination_counter(interesting_modifications):
            updated_modifications = ModificationIndex(other_modifications)
            for k, v in varied_modifications.items():
                if v != 0:
                    updated_modifications[k] += v
            if self.compute_compositions:
                extra_composition = Composition()
                for mod, mod_count in varied_modifications.items():
                    extra_composition += self.modification_compositions[mod.name] * mod_count
            variants.append((updated_modifications, extra_composition))
        return variants

    def partial_loss(self, fragment):
        (modifications_of_interest,
         delta_composition,
         other_modifications) = self._get_modifications_of_interest(fragment)

        if self.compute_compositions:
            base_composition = fragment.composition - delta_composition
        else:
            base_composition = None

        self._replace_cores(modifications_of_interest)

        if (modifications_of_interest, other_modifications) == self._last_modification_set:
            variants = self._last_modification_variants
        else:
            self._last_modification_set = modifications_of_interest, other_modifications
            variants = self._last_modification_variants = self._generate_modification_variants(
                modifications_of_interest, other_modifications)

        series = self.series
        fragments = []
        for updated_modifications, extra_composition in variants:
            fragments.append(
                PeptideFragment(
                    series, fragment.position, updated_modifications, fragment.bare_mass,
                    flanking_amino_acids=fragment.flanking_amino_acids,
                    composition=(base_composition + extra_composition) if self.compute_compositions else None))
        return fragments



try:
    _has_c = True
    from glycopeptidepy._c.structure.fragmentation_strategy.peptide import (
        PeptideFragmentationStrategyBase,
        HCDFragmentationStrategy)
except ImportError:
    _has_c = False

# Cooper, H. J., Hudgins, R. R., Håkansson, K., & Marshall, A. G. (2002).
# Characterization of amino acid side chain losses in electron capture dissociation.
# Journal of the American Society for Mass Spectrometry, 13(3), 241–249.
# https://doi.org/10.1016/S1044-0305(01)00357-9
exd_sidechain_losses = defaultdict(list, {
    R("His"): [-Composition("H6C4N2")],
    R("Arg"): [-Composition("CH5N3"),
               -Composition("C4H11N3"),
               -Composition("C1N2H4"),
               -Composition("N2H6")],
    R("Asn"): [-Composition("CONH3")],
    R("Gln"): [-Composition("CONH3")],
    R("Met"): [-Composition("C3H6S")]
})

for key, value in _residue_to_neutral_loss.items():
    exd_sidechain_losses[R(key)].extend(value)


class GlycanCompositionFragment(SimpleFragment):
    __slots__ = ('glycan_composition', )

    def __init__(self, name, mass, kind, composition, glycan_composition, chemical_shift=None):
        super(GlycanCompositionFragment, self).__init__(
            name, mass, kind, composition, chemical_shift=chemical_shift, is_glycosylated=True)
        self.glycan_composition = glycan_composition


class _GlycanFragmentingMixin(_GlycanFragmentingStrategyBase):

    def _make_fragmentation_strategy(self):
        strat = StubGlycopeptideStrategy(self.peptide, True)
        return strat

    def _get_glycan_composition_fragments(self, glycosylation, series, position):
        strat = self._make_fragmentation_strategy()
        gc = strat.glycan_composition()
        if glycosylation.rule.glycosylation_type == GlycosylationType.n_linked:
            gen = strat.n_glycan_composition_fragments(gc)
        elif glycosylation.rule.glycosylation_type == GlycosylationType.o_linked:
            gen = strat.o_glycan_composition_fragments(gc)
        elif glycosylation.rule.glycosylation_type == GlycosylationType.glycosaminoglycan:
            gen = strat.gag_linker_composition_fragments(gc)
        else:
            gen = []
        fragments = []
        for frag_spec in gen:
            key = frag_spec['key']
            if not key:
                continue
            name = ''.join("%s%d" % kv for kv in sorted(
                key.items(), key=lambda x: x[0].mass()))
            mass = frag_spec['mass']
            composition = frag_spec['composition']
            fragments.append(GlycanCompositionFragment(
                name, mass, IonSeries.other, composition, frag_spec['key']))
        return fragments

    def _get_glycan_structure_fragments(self, glycosylation, series, position):
        fragments = sorted(glycosylation.rule.get_fragments(
            series, max_cleavages=self.max_glycan_cleavages),
            key=lambda x: x.mass)
        return fragments

    def _get_glycan_fragments(self, glycosylation, series, position):
        try:
            return self._glycan_fragment_cache[glycosylation, series, position]
        except KeyError:
            if glycosylation.rule.is_composition:
                fragments = self._get_glycan_composition_fragments(glycosylation, series, position)
            else:
                fragments = self._get_glycan_structure_fragments(glycosylation, series, position)
            self._glycan_fragment_cache[glycosylation, series, position] = fragments
            return fragments

    def _wrap_glycan_fragment(self, fragment):
        try:
            return self._glycan_fragment_wrapper_cache[fragment]
        except KeyError:
            result = self._glycan_fragment_wrapper_cache[fragment] = GlycanFragment(
                fragment)
            return result


class EXDFragmentationStrategy(PeptideFragmentationStrategyBase, _GlycanFragmentingMixin):
    glycan_fragment_ladder = "Y"

    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1, max_glycan_cleavages=None,
                 include_neutral_losses=False, **kwargs):
        super(EXDFragmentationStrategy, self).__init__(
            peptide, series, chemical_shifts, max_chemical_shifts, include_neutral_losses, **kwargs)
        if max_glycan_cleavages is None:
            max_glycan_cleavages = self._guess_max_glycan_cleavages()
        self.max_glycan_cleavages = max_glycan_cleavages
        self._glycan_fragment_cache = {}
        self._glycan_fragment_wrapper_cache = {}
        self._glycan_fragment_combination_cache = {}

    def _strip_glycosylation_from_modification_dict(self, modification_dict):
        glycosylations = []
        stripped_modifications = {}
        for position, glycosylation in sorted(self.glycosylation_manager.items()):
            glycosylations.append((position, glycosylation))
        for modification, count in modification_dict.items():
            if modification.is_tracked_for(ModificationCategory.glycosylation):
                continue
            stripped_modifications[modification] = count
        return glycosylations, stripped_modifications

    def _get_glycan_loss_combinations(self, glycosylations, glycan_fragment_ladder):
        glycosylations = tuple(glycosylations)
        try:
            wrapped_fragment_combinations = self._glycan_fragment_combination_cache[
                glycosylations, glycan_fragment_ladder]
        except KeyError:
            glycan_fragments = [
                [None] + self._get_glycan_fragments(glycosylation, glycan_fragment_ladder, position)
                for position, glycosylation in glycosylations
            ]
            fragment_combinations = product(*glycan_fragments)

            has_multiple_nonlocalized_cores = len(glycosylations) > 1 and sum(
                g[1].rule.is_core for g in glycosylations) == len(glycosylations)

            wrapped_fragment_combinations = []
            for fragment_set in fragment_combinations:
                if fragment_set == ():
                    continue
                aggregate_glycan_composition = None
                if has_multiple_nonlocalized_cores:
                    aggregate_glycan_composition = _AccumulatorBag()
                subfragments = []
                for glycan_fragment in fragment_set:
                    if glycan_fragment is None:
                        continue
                    subfragment = self._wrap_glycan_fragment(glycan_fragment)
                    if has_multiple_nonlocalized_cores and glycan_fragment is not None:
                        aggregate_glycan_composition += glycan_fragment.glycan_composition
                    subfragments.append(subfragment)
                if aggregate_glycan_composition is not None:
                    invalid = False
                    for k, v in self.peptide.glycan_composition.items():
                        if aggregate_glycan_composition[k] > v:
                            invalid = True
                            break
                    if invalid:
                        continue
                wrapped_fragment_combinations.append(subfragments)
            self._glycan_fragment_combination_cache[glycosylations,
                                                    glycan_fragment_ladder] = wrapped_fragment_combinations
        return wrapped_fragment_combinations

    def partial_loss(self, fragment, glycan_fragment_ladder=None):
        if glycan_fragment_ladder is None:
            glycan_fragment_ladder = self.glycan_fragment_ladder
        # make a clone in case we mutate anything
        base = fragment.clone()

        (glycosylations,
         stripped_modifications) = self._strip_glycosylation_from_modification_dict(base.modification_dict)
        if self.compute_compositions:
            base_composition = base.composition
            for _glycosylation_position, glycosylation in glycosylations:
                base_composition -= glycosylation.composition
        else:
            base_composition = None
        bare_fragment = PeptideFragment(
            fragment.series, fragment.position, ModificationIndex(stripped_modifications),
            fragment.bare_mass,
            flanking_amino_acids=fragment.flanking_amino_acids,
            composition=base_composition.clone() if self.compute_compositions else None)

        wrapped_fragment_combinations = self._get_glycan_loss_combinations(
            glycosylations, glycan_fragment_ladder)

        results = []
        for fragment_set in wrapped_fragment_combinations:
            if self.compute_compositions:
                delta_composition = Composition()
            else:
                delta_composition = None
            new_modifications = ModificationIndex(stripped_modifications)
            for subfragment in fragment_set:
                new_modifications[subfragment] += 1

            extended_fragment = PeptideFragment(
                bare_fragment.series, bare_fragment.position, new_modifications,
                bare_fragment.bare_mass,
                flanking_amino_acids=bare_fragment.flanking_amino_acids,
                composition=bare_fragment.composition + delta_composition if self.compute_compositions else None)
            results.append(extended_fragment)
        results.append(fragment)
        return results
