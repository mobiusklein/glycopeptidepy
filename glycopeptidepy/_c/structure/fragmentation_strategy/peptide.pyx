# cython: embedsignature=True, profile=True

from collections import defaultdict
from itertools import product, combinations

cimport cython

from cpython.ref cimport Py_INCREF
from cpython cimport PyObject
from cpython.tuple cimport PyTuple_GetItem
from cpython.list cimport PyList_GetItem, PyList_SetItem, PyList_Size, PyList_New
from cpython.dict cimport (PyDict_GetItem, PyDict_SetItem, PyDict_Next,
                           PyDict_Keys, PyDict_Update, PyDict_DelItem, PyDict_Size)

from cpython.int cimport PyInt_AsLong, PyInt_Check, PyInt_FromLong

from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.collectiontools cimport descending_combination_counter
from glycopeptidepy._c.structure.constants cimport Configuration
from glycopeptidepy._c.structure.fragment cimport IonSeriesBase, PeptideFragment, ChemicalShiftBase
from glycopeptidepy._c.structure.sequence_methods cimport _PeptideSequenceCore
from glycopeptidepy._c.structure.base cimport ModificationBase, AminoAcidResidueBase, SequencePosition
from glycopeptidepy._c.structure.modification.modification cimport ModificationInstanceBase

from glycopeptidepy.structure.glycan import GlycosylationManager
from glycopeptidepy.structure import constants as _structure_constants
from glycopeptidepy.structure.modification import (
    ModificationCategory,
    ModificationIndex)

from glycopeptidepy._c.structure.fragmentation_strategy.base cimport (
    FragmentationStrategyBase,
    _n_glycosylation,
    _modification_hexnac,
    _o_glycosylation,
    _gag_linker_glycosylation,
    _modification_xylose,
    glycosylation_type_to_core)

from glycopeptidepy.utils.collectiontools import _AccumulatorBag

from glycopeptidepy.structure.modification import (
    Modification,
    NGlycanCoreGlycosylation,
    OGlycanCoreGlycosylation,
    GlycosaminoglycanLinkerGlycosylation)

from glycopeptidepy.structure.fragment import (
    ChemicalShift, IonSeries, format_negative_composition)

from glycopeptidepy.structure.glycan import (GlycosylationType)


cdef Configuration structure_constants = _structure_constants
cdef ModificationCategory_glycosylation = ModificationCategory.glycosylation


cdef class PeptideFragmentationStrategyBase(FragmentationStrategyBase):

    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1):
        if chemical_shifts is None:
            chemical_shifts = dict()
        super(PeptideFragmentationStrategyBase, self).__init__(peptide)
        self.series = IonSeries(series)
        self.chemical_shift_rules = defaultdict(list, chemical_shifts)
        self.max_chemical_shifts = max_chemical_shifts
        self._initialize_fields()

    cpdef _initialize_fields(self):
        self.direction = self.series.direction

        # Null Values
        self.running_mass = 0
        self.running_composition = CComposition()
        self.index = -1
        self.size = self.peptide.get_size()

        self.modification_index = ModificationIndex()
        self.glycosylation_manager = GlycosylationManager(self.peptide)
        self.amino_acids_counter = _AccumulatorBag()

        self.running_mass += self.series.mass_shift
        self.running_composition.add_from(self.series.composition_shift)
        self._initialize_start_terminal()

    cpdef list _get_viable_chemical_shift_combinations(self):
        shifts = []
        if not self.chemical_shift_rules:
            return shifts
        for residue, count in self.amino_acids_counter.items():
            loss_composition = self.chemical_shift_rules[residue]
            if loss_composition:
                for i in range(count):
                    shifts.append(loss_composition + [CComposition()])
        loss_composition_combinations = {}
        for comb in combinations(shifts, self.max_chemical_shifts):
            for prod in product(*comb):
                loss_composition = sum(prod, CComposition())
                loss_composition_combinations[format_negative_composition(loss_composition)] = loss_composition
        return [ChemicalShift(k, v) for k, v in loss_composition_combinations.items() if v]

    cpdef _initialize_start_terminal(self):
        if self.direction > 0:
            self.running_mass += self.peptide.get_n_term().mass
            self.running_composition.add_from(self.peptide.get_n_term().get_composition())
            self.index = -1
        elif self.direction < 0:
            self.running_mass += self.peptide.get_c_term().mass
            self.running_composition.add_from(self.peptide.get_c_term().get_composition())
            self.index = self.peptide.get_size()
        else:
            raise ValueError("Unknown direction %r" % (self.series.direction,))

    cpdef CComposition composition_of(self, SequencePosition position):
        cdef:
            size_t i, n
            CComposition composition
            ModificationBase mod

        n = PyList_Size(position.modifications)
        composition = position.amino_acid.composition.clone()
        for i in range(n):
            mod = <ModificationBase>PyList_GetItem(position.modifications, i)
            composition.add_from(mod.composition)
        return composition

    cpdef bint has_more(self):
        if self.direction > 0:
            return self.index < self.size - 2
        else:
            return self.index > 1

    cpdef list flanking_residues(self):
        residues = [self.peptide[self.index][0],
                    self.peptide[self.index + self.direction][0]]
        if self.direction < 0:
            residues = residues[::-1]
        return residues

    cpdef long name_index_of(self):
        if self.direction > 0:
            return self.index + structure_constants.FRAG_OFFSET
        else:
            return (self.size - self.index - 1) + structure_constants.FRAG_OFFSET

    def __next__(self):
        if self.has_more():
            return self.step()
        else:
            raise StopIteration()

    cpdef track_glycosylation(self, index, glycosylation):
        self.glycosylation_manager[self.index] = glycosylation
        self.modification_index[glycosylation] += 1

    cpdef _update_state(self):
        cdef:
            SequencePosition position
            ModificationBase mod
            size_t i, n
        self.index += self.direction
        position = self.peptide[self.index]
        n = PyList_Size(position.modifications)
        for i in range(n):
            mod = <ModificationBase>PyList_GetItem(position.modifications, i)
            if mod.is_tracked_for(ModificationCategory_glycosylation):
                self.track_glycosylation(self.index, mod)
            else:
                self.modification_index[mod] += 1
        self.amino_acids_counter[position.amino_acid] += 1
        composition = self.composition_of(position)
        self.running_mass += position.amino_acid.mass
        self.running_composition.add_from(composition)

    cpdef list _build_fragments(self):
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

    cpdef list partial_loss(self, PeptideFragment fragment):
        return [fragment]

    cpdef object step(self):
        self._update_state()
        fragments_from_site = self._build_fragments()
        return fragments_from_site

    def __repr__(self):
        return "%s(%s, %r, %r, %0.4f, %d)" % (
            self.__class__.__name__,
            self.peptide, self.series, self.running_composition,
            self.running_mass, self.index)


hcd_modifications_of_interest = {k.name: k for k in [
        _n_glycosylation,
        _modification_hexnac,
        _o_glycosylation,
        _gag_linker_glycosylation,
        _modification_xylose
    ]}

hcd_modification_compositions = {
    k.name: k.composition for k in hcd_modifications_of_interest.values()
}


cdef class ModificationConfiguration(object):

    @staticmethod
    cdef ModificationConfiguration _create(object modifications_of_interest, object other_modifications, CComposition delta_composition):
        cdef ModificationConfiguration inst = ModificationConfiguration.__new__(ModificationConfiguration)
        inst.modifications_of_interest = modifications_of_interest
        inst.other_modifications = other_modifications
        inst.delta_composition = delta_composition
        return inst

    def __init__(self, modifications_of_interest, other_modifications, delta_composition):
        self.modifications_of_interest = modifications_of_interest
        self.other_modifications = other_modifications
        self.delta_composition = delta_composition

    def __iter__(self):
        yield self.modifications_of_interest
        yield self.delta_composition
        yield self.other_modifications

    def __repr__(self):
        return "{self.__class__.__name__}({self.modifications_of_interest}, {self.other_modifications})".format(
            self=self)

    cdef bint equal_to(self, ModificationConfiguration other):
        if other is None:
            return False
        if self.modifications_of_interest != other.modifications_of_interest:
            return False
        if self.other_modifications != other.other_modifications:
            return False
        return True


cdef class HCDFragmentationStrategy(PeptideFragmentationStrategyBase):
    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1):
        super(HCDFragmentationStrategy, self).__init__(peptide, series, chemical_shifts, max_chemical_shifts)
        self._last_modification_set = None
        self._last_modification_variants = None

    cpdef _get_core_for(self, ModificationInstanceBase glycosylation):
        try:
            if not glycosylation.rule.is_core:
                glycosylation = glycosylation_type_to_core[glycosylation.rule.glycosylation_type]()
            return glycosylation
        except KeyError:
            raise ValueError("Cannot determine which core to use for {}".format(
                glycosylation.rule.glycosylation_type))

    cpdef CComposition composition_of(self, SequencePosition position):
        cdef:
            CComposition composition
            ModificationInstanceBase mod
            size_t i, n
        composition = position.amino_acid.composition.clone()
        n = PyList_Size(position.modifications)
        for i in range(n):
            mod = <ModificationInstanceBase>PyList_GetItem(position.modifications, i)
            if mod.is_tracked_for(ModificationCategory_glycosylation):
                mod = self._get_core_for(mod)
            composition.add_from(mod.composition)
        return composition

    cpdef track_glycosylation(self, index, glycosylation):
        # HCD strategy does not track intact topologies
        try:
            glycosylation = self._get_core_for(glycosylation)
        except KeyError:
            raise ValueError("Cannot determine which core to use for {}".format(
                glycosylation.rule.glycosylation_type))
        self.glycosylation_manager[self.index] = glycosylation
        self.modification_index[glycosylation] += 1

    cpdef ModificationConfiguration _get_modifications_of_interest(self, PeptideFragment fragment):
        cdef:
            ModificationConfiguration result
            CComposition delta_composition
            dict modifications, other_modifications
            object modifications_of_interest
            Py_ssize_t pos
            PyObject *key
            PyObject *value
            ModificationBase mod
            long count


        modifications = dict(fragment.modification_dict)
        delta_composition = CComposition()
        other_modifications = dict()
        modifications_of_interest = defaultdict(int)
        pos = 0
        while PyDict_Next(modifications, &pos, &key, &value):
            mod = <ModificationBase>key
            count = PyInt_AsLong(<object>value)
            if mod.name in hcd_modifications_of_interest:
                modifications_of_interest[mod] = count
                delta_composition.add_from(mod.composition * count)
            else:
                other_modifications[mod] = <object>value

        result = ModificationConfiguration._create(
            modifications_of_interest, other_modifications, delta_composition)
        return result

    cpdef _replace_cores(self, modifications_of_interest):
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

    cpdef list _generate_modification_variants(self, interesting_modifications, other_modifications):
        cdef:
            list variants, variant_modification_list
            dict varied_modifications
            object updated_modifications
            size_t i, n
            Py_ssize_t pos
            PyObject *key
            PyObject *value
            ModificationBase mod
            long count

        variants = []
        variant_modification_list = list(descending_combination_counter(interesting_modifications))
        n = PyList_Size(variant_modification_list)
        for i in range(n):
            varied_modifications = <dict>PyList_GetItem(variant_modification_list, i)
            updated_modifications = ModificationIndex(other_modifications)
            extra_composition = CComposition()
            pos = 0
            while PyDict_Next(varied_modifications, &pos, &key, &value):
                mod = <ModificationBase>key
                count = PyInt_AsLong(<object>value)
                if count != 0:
                    updated_modifications[mod] += <object>value
                    extra_composition.add_from(mod.composition * count)

            variants.append((updated_modifications, extra_composition))
        return variants

    cpdef list partial_loss(self, PeptideFragment fragment):
        cdef:
            CComposition base_composition, delta_composition
            ModificationConfiguration mod_config
            list variants, fragments
            IonSeriesBase series
            tuple variant_pair
            size_t i, n

        mod_config = self._get_modifications_of_interest(fragment)

        base_composition = fragment.composition.clone()
        base_composition.subtract_from(mod_config.delta_composition)

        self._replace_cores(mod_config.modifications_of_interest)

        if mod_config.equal_to(self._last_modification_set):
            variants = self._last_modification_variants
        else:
            self._last_modification_set = mod_config
            variants = self._last_modification_variants = self._generate_modification_variants(
                mod_config.modifications_of_interest,
                mod_config.other_modifications)

        series = self.series
        fragments = []
        n = PyList_Size(variants)
        for i in range(n):
            variant_pair = <tuple>PyList_GetItem(variants, i)
            updated_modifications = <object>PyTuple_GetItem(variant_pair, 0)
            extra_composition = <CComposition>PyTuple_GetItem(variant_pair, 1)
            fragments.append(
                PeptideFragment(
                    series, fragment.position, dict(updated_modifications), fragment.bare_mass,
                    flanking_amino_acids=fragment.flanking_amino_acids,
                    composition=base_composition + extra_composition))
        return fragments
