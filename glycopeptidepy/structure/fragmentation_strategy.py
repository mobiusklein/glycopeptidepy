# -*- coding: utf-8 -*-
from collections import defaultdict
from itertools import product

from glycopeptidepy.utils.collectiontools import descending_combination_counter
from glycopeptidepy.structure import constants as structure_constants
from glycopeptidepy.structure.composition import (
    Composition)
from glycopeptidepy.structure.modification import (
    ModificationCategory,
    Modification,
    NGlycanCoreGlycosylation,
    OGlycanCoreGlycosylation,
    GlycosaminoglycanLinkerGlycosylation,
    GlycanFragment,
    ModificationIndex)
from glycopeptidepy.structure.glycan import (
    GlycosylationManager,
    GlycosylationType)
from glycopeptidepy.structure.fragment import (
    IonSeries,
    PeptideFragment,
    NeutralLoss)
from glycopeptidepy.structure.residue import AminoAcidResidue as R


_n_glycosylation = NGlycanCoreGlycosylation()
_o_glycosylation = OGlycanCoreGlycosylation()
_gag_linker_glycosylation = GlycosaminoglycanLinkerGlycosylation()
_modification_hexnac = Modification("HexNAc").rule
_modification_xylose = Modification("Xyl").rule


glycosylation_type_to_core = {
    GlycosylationType.n_linked: _n_glycosylation,
    GlycosylationType.o_linked: _o_glycosylation,
    GlycosylationType.glycosaminoglycan: _gag_linker_glycosylation
}


class SimpleCounter(dict):
    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            return 0


class FragmentationStrategyBase(object):

    def __init__(self, peptide, series):
        self.peptide = peptide
        self.series = IonSeries(series)

        # Null Values
        self.running_mass = 0
        self.running_composition = Composition()
        self.index = None
        self.size = len(self.peptide)
        self.modification_index = ModificationIndex()
        self.glycosylation_manager = GlycosylationManager(self.peptide)
        self.amino_acids_counter = SimpleCounter()

        self.running_mass += self.series.mass_shift
        self.running_composition += self.series.composition_shift
        self._initialize_start_terminal()

    @property
    def direction(self):
        return self.series.direction

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

    def next(self):
        return self.__next__()

    def __iter__(self):
        return self

    def track_glycosylation(self, index, glycosylation):
        self.glycosylation_manager[self.index] = glycosylation
        self.modification_index[glycosylation] += 1

    def _update_state(self):
        self.index += self.direction
        residue, modifications = self.peptide[self.index]
        for mod in modifications:
            if mod.is_a(ModificationCategory.glycosylation):
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
        fragments_from_site.extend(self.partial_loss(frag))
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


class HCDFragmentationStrategy(FragmentationStrategyBase):
    modifications_of_interest = [
        _n_glycosylation,
        _modification_hexnac,
        _o_glycosylation,
        _gag_linker_glycosylation,
        _modification_xylose
    ]

    modification_compositions = {
        k: k.composition for k in modifications_of_interest
    }

    def _get_core_for(self, glycosylation):
        try:
            if not glycosylation.rule.is_core:
                glycosylation = glycosylation_type_to_core[glycosylation.rule.glycosylation_type]
            return glycosylation
        except KeyError:
            raise ValueError("Cannot determine which core to use for {}".format(
                glycosylation.rule.glycosylation_type))

    def composition_of(self, residue, modifications):
        composition = Composition(residue.composition)
        for mod in modifications:
            if mod.is_a(ModificationCategory.glycosylation):
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
        modifications_of_interest = defaultdict(
            int, {k: v for k, v in modifications.items()
                  if k in self.modifications_of_interest})
        delta_composition = sum(
            (self.modification_compositions[k] * v
             for k, v in modifications_of_interest.items()),
            Composition())
        return modifications_of_interest, delta_composition

    def _replace_cores(self, modifications_of_interest):
        n_cores = modifications_of_interest.pop(_n_glycosylation, 0)
        o_cores = modifications_of_interest.pop(_o_glycosylation, 0)
        gag_cores = modifications_of_interest.pop(_gag_linker_glycosylation, 0)

        # Convert core glycosylation into the residual monosaccharides remaining
        # after dissociation
        modifications_of_interest[_modification_hexnac] += n_cores * 2 + o_cores
        modifications_of_interest[_modification_xylose] += gag_cores
        return modifications_of_interest

    def _generate_modification_variants(self, interesting_modifications, other_modifications):
        for varied_modifications in descending_combination_counter(interesting_modifications):
            updated_modifications = other_modifications.copy()
            updated_modifications.update(
                {k: v for k, v in varied_modifications.items() if v != 0})

            extra_composition = Composition()
            for mod, mod_count in varied_modifications.items():
                extra_composition += self.modification_compositions[mod] * mod_count
            yield updated_modifications, extra_composition

    def partial_loss(self, fragment):
        (modifications_of_interest,
         delta_composition) = self._get_modifications_of_interest(fragment)
        base_composition = fragment.composition - delta_composition
        self._replace_cores(modifications_of_interest)
        other_modifications = {
            k: v for k, v in fragment.modification_dict.items()
            if k not in self.modifications_of_interest
        }
        for updated_modifications, extra_composition in self._generate_modification_variants(
                modifications_of_interest, other_modifications):
            yield PeptideFragment(
                fragment.series, fragment.position, dict(updated_modifications), fragment.bare_mass,
                golden_pairs=fragment.golden_pairs,
                flanking_amino_acids=fragment.flanking_amino_acids,
                composition=base_composition + extra_composition)


# Cooper, H. J., Hudgins, R. R., Håkansson, K., & Marshall, A. G. (2002).
# Characterization of amino acid side chain losses in electron capture dissociation.
# Journal of the American Society for Mass Spectrometry, 13(3), 241–249.
# https://doi.org/10.1016/S1044-0305(01)00357-9
exd_sidechain_losses = {
    R("His"): [Composition("H6C4N2")],
    R("Arg"): [Composition("CH5N3"),
               Composition("C4H11N3"),
               Composition("C1N2H4"),
               Composition("N2H6")],
    R("Asn"): [Composition("CONH3")],
    R("Gln"): [Composition("CONH3")],
    R("Met"): [Composition("C3H6S")]
}


class EXDFragmentationStrategy(FragmentationStrategyBase):
    glycan_fragment_ladder = "Y"

    def __init__(self, peptide, series, max_glycan_cleavages=2):
        super(EXDFragmentationStrategy, self).__init__(peptide, series)
        self.max_glycan_cleavages = max_glycan_cleavages

    def _strip_glycosylation_from_modification_dict(self, modification_dict):
        glycosylations = []
        stripped_modifications = {}
        for position, glycosylation in sorted(self.glycosylation_manager.items()):
            glycosylations.append(glycosylation)
        for modification, count in modification_dict.items():
            if modification.is_a(ModificationCategory.glycosylation):
                continue
            stripped_modifications[modification] = count
        return glycosylations, stripped_modifications

    def partial_loss(self, fragment, glycan_fragment_ladder=None):
        if glycan_fragment_ladder is None:
            glycan_fragment_ladder = self.glycan_fragment_ladder
        # make a clone in case we mutate anything
        base = fragment.clone()

        (glycosylations,
         stripped_modifications) = self._strip_glycosylation_from_modification_dict(
            base.modification_dict)
        base_composition = base.composition
        for glycosylation in glycosylations:
            base_composition -= glycosylation.composition
        bare_fragment = PeptideFragment(
            fragment.series, fragment.position, dict(stripped_modifications),
            fragment.bare_mass,
            golden_pairs=fragment.golden_pairs,
            flanking_amino_acids=fragment.flanking_amino_acids,
            composition=base_composition.clone())

        fragments = [
            [None] + sorted(glycosylation.rule.get_fragments(
                glycan_fragment_ladder, max_cleavages=self.max_glycan_cleavages),
                key=lambda x: x.mass)
            for glycosylation in glycosylations if not glycosylation.rule.is_composition
        ]
        fragment_combinations = product(*fragments)
        for fragment_set in fragment_combinations:
            delta_composition = Composition()
            new_modifications = ModificationIndex(stripped_modifications)
            for glycan_fragment in fragment_set:
                if glycan_fragment is None:
                    continue
                subfragment = GlycanFragment(glycan_fragment)
                new_modifications[subfragment] += 1
                delta_composition += subfragment.composition

            extended_fragment = PeptideFragment(
                bare_fragment.series, bare_fragment.position, new_modifications,
                bare_fragment.bare_mass,
                golden_pairs=bare_fragment.golden_pairs,
                flanking_amino_acids=bare_fragment.flanking_amino_acids,
                composition=bare_fragment.composition + delta_composition)
            yield extended_fragment

        # yield the intact fragment
        yield fragment
