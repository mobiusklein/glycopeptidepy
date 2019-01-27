# -*- coding: utf-8 -*-
from collections import defaultdict
from itertools import product, combinations

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


class PeptideFragmentationStrategyBase(FragmentationStrategyBase):

    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1):
        if chemical_shifts is None:
            chemical_shifts = dict()
        super(PeptideFragmentationStrategyBase, self).__init__(peptide)
        self.series = IonSeries(series)
        self.chemical_shift_rules = defaultdict(list, chemical_shifts)
        self.max_chemical_shifts = max_chemical_shifts
        self._initialize_fields()

    def _initialize_fields(self):
        self.direction = self.series.direction

        # Null Values
        self.running_mass = 0
        self.running_composition = Composition()
        self.index = None
        self.size = len(self.peptide)
        self.modification_index = ModificationIndex()
        self.glycosylation_manager = GlycosylationManager(self.peptide)
        self.amino_acids_counter = _AccumulatorBag()

        self.running_mass += self.series.mass_shift
        self.running_composition += self.series.composition_shift
        self._initialize_start_terminal()

    def _get_viable_chemical_shift_combinations(self):
        shifts = []
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
            self.running_composition += self.peptide.n_term.composition
            self.index = -1
        elif self.direction < 0:
            self.running_mass += self.peptide.c_term.mass
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
        composition = self.composition_of(residue, modifications)
        self.running_mass += residue.mass
        self.running_composition += composition

    def _build_fragments(self):
        fragments_from_site = []
        frag = PeptideFragment(
            self.series,
            self.name_index_of(),
            dict(self.modification_index),
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

    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1):
        super(HCDFragmentationStrategy, self).__init__(peptide, series, chemical_shifts, max_chemical_shifts)
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
        modifications = dict(fragment.modification_dict)
        delta_composition = Composition()
        other_modifications = dict()
        modifications_of_interest = defaultdict(int)

        for k, v in modifications.items():
            if k.name in self.modifications_of_interest:
                modifications_of_interest[k] = v
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

            extra_composition = Composition()
            for mod, mod_count in varied_modifications.items():
                extra_composition += self.modification_compositions[mod.name] * mod_count
            variants.append((updated_modifications, extra_composition))
        return variants

    def partial_loss(self, fragment):
        (modifications_of_interest,
         delta_composition,
         other_modifications) = self._get_modifications_of_interest(fragment)

        base_composition = fragment.composition - delta_composition

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
                    series, fragment.position, dict(updated_modifications), fragment.bare_mass,
                    flanking_amino_acids=fragment.flanking_amino_acids,
                    composition=base_composition + extra_composition))
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


class EXDFragmentationStrategy(PeptideFragmentationStrategyBase, _GlycanFragmentingStrategyBase):
    glycan_fragment_ladder = "Y"

    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1, max_glycan_cleavages=None):
        super(EXDFragmentationStrategy, self).__init__(peptide, series, chemical_shifts, max_chemical_shifts)
        if max_glycan_cleavages is None:
            max_glycan_cleavages = self._guess_max_glycan_cleavages()
        self.max_glycan_cleavages = max_glycan_cleavages
        self._glycan_fragment_cache = {}

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

    def _get_glycan_fragments(self, glycosylation, series, position):
        try:
            return self._glycan_fragment_cache[glycosylation, series, position]
        except KeyError:
            if glycosylation.rule.is_composition:
                strat = StubGlycopeptideStrategy(self.peptide, True)
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
                    name = ''.join("%s%d" % kv for kv in sorted(key.items()))
                    mass = frag_spec['mass']
                    composition = frag_spec['composition']
                    fragments.append(SimpleFragment(name, mass, 'Y', composition))
                return fragments
            else:
                fragments = sorted(glycosylation.rule.get_fragments(
                    series, max_cleavages=self.max_glycan_cleavages),
                    key=lambda x: x.mass)
                self._glycan_fragment_cache[glycosylation, series, position] = fragments
                return fragments

    def partial_loss(self, fragment, glycan_fragment_ladder=None):
        if glycan_fragment_ladder is None:
            glycan_fragment_ladder = self.glycan_fragment_ladder
        # make a clone in case we mutate anything
        base = fragment.clone()

        (glycosylations,
         stripped_modifications) = self._strip_glycosylation_from_modification_dict(
            base.modification_dict)
        base_composition = base.composition
        for position, glycosylation in glycosylations:
            base_composition -= glycosylation.composition
        bare_fragment = PeptideFragment(
            fragment.series, fragment.position, dict(stripped_modifications),
            fragment.bare_mass,
            flanking_amino_acids=fragment.flanking_amino_acids,
            composition=base_composition.clone())

        fragments = [
            [None] + self._get_glycan_fragments(glycosylation, glycan_fragment_ladder, position)
            for position, glycosylation in glycosylations
        ]
        fragment_combinations = product(*fragments)
        results = []
        for fragment_set in fragment_combinations:
            if fragment_set == ():
                continue
            delta_composition = Composition()
            new_modifications = ModificationIndex(stripped_modifications)
            for glycan_fragment in fragment_set:
                if glycan_fragment is None:
                    continue
                subfragment = GlycanFragment(glycan_fragment)
                new_modifications[subfragment] += 1
                delta_composition += subfragment.composition

            extended_fragment = PeptideFragment(
                bare_fragment.series, bare_fragment.position, dict(new_modifications),
                bare_fragment.bare_mass,
                flanking_amino_acids=bare_fragment.flanking_amino_acids,
                composition=bare_fragment.composition + delta_composition)
            results.append(extended_fragment)
        # yield the intact fragment
        results.append(fragment)
        return results
