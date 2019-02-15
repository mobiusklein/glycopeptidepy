# cython: embedsignature=True

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
from glycopeptidepy._c.count_table cimport CountTable, CountTableIterator
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

cdef ChemicalShiftBase ammonia_loss = ChemicalShift("-NH3", -CComposition({"N": 1, "H": 3}))
cdef ChemicalShiftBase water_loss = ChemicalShift("-H2O", -CComposition({"O": 1, "H": 2}))


cdef class PeptideFragmentationStrategyBase(FragmentationStrategyBase):

    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1, include_neutral_losses=False):
        if chemical_shifts is None:
            chemical_shifts = dict()
        super(PeptideFragmentationStrategyBase, self).__init__(peptide)
        self.series = IonSeries(series)
        self.chemical_shift_rules = defaultdict(list, chemical_shifts)
        self.max_chemical_shifts = max_chemical_shifts
        self.include_neutral_losses = include_neutral_losses
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
        cdef:
            bint has_neutral_loss_rules
            list shifts

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
        cdef:
            list residues
            SequencePosition pos_a, pos_b

        residues = PyList_New(2)
        pos_a = self.peptide.get(self.index)
        pos_b = self.peptide.get(self.index + self.direction)
        if self.direction < 0:
            Py_INCREF(pos_b.amino_acid)
            PyList_SetItem(residues, 0, pos_b.amino_acid)
            Py_INCREF(pos_a.amino_acid)
            PyList_SetItem(residues, 1, pos_a.amino_acid)
        else:
            Py_INCREF(pos_a.amino_acid)
            PyList_SetItem(residues, 0, pos_a.amino_acid)
            Py_INCREF(pos_b.amino_acid)
            PyList_SetItem(residues, 1, pos_b.amino_acid)
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

    cpdef track_glycosylation(self, long index, glycosylation):
        self.glycosylation_manager[self.index] = glycosylation
        self.modification_index.increment(glycosylation, 1)

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
                self.modification_index.increment(mod, 1)
        self.amino_acids_counter.increment(position.amino_acid, 1)
        composition = self.composition_of(position)
        self.running_mass += position.amino_acid.mass
        self.running_composition.add_from(composition)

    cpdef list _build_fragments(self):
        cdef:
            list fragments_from_site
            list partial_loss_fragments
            PeptideFragment fragment, f
            size_t i, n, j, m, k

        frag = PeptideFragment._create(
            self.series,
            self.name_index_of(),
            self.modification_index.copy(),
            self.running_mass,
            flanking_amino_acids=self.flanking_residues(),
            glycosylation=self.glycosylation_manager.copy() if self.glycosylation_manager else None,
            chemical_shift=None,
            composition=self.running_composition)
        shifts = self._get_viable_chemical_shift_combinations()
        partial_loss_fragments = self.partial_loss(frag)

        n = PyList_Size(partial_loss_fragments)
        m = PyList_Size(shifts)
        fragments_from_site = PyList_New(n * (m + 1))
        k = 0
        for i in range(n):
            fragment = <PeptideFragment>PyList_GetItem(partial_loss_fragments, i)
            Py_INCREF(fragment)
            PyList_SetItem(fragments_from_site, k, fragment)
            k += 1
            for j in range(m):
                shift = <ChemicalShiftBase>PyList_GetItem(shifts, j)
                f = fragment.clone()
                f.chemical_shift = shift
                Py_INCREF(f)
                PyList_SetItem(fragments_from_site, k, f)
                k += 1
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
    cdef ModificationConfiguration _create(CountTable modifications_of_interest, CountTable other_modifications,
                                           CComposition delta_composition, CountTable modification_set):
        cdef ModificationConfiguration inst = ModificationConfiguration.__new__(ModificationConfiguration)
        inst.modifications_of_interest = modifications_of_interest
        inst.other_modifications = other_modifications
        inst.delta_composition = delta_composition
        inst.modification_set = modification_set
        return inst

    def __init__(self, modifications_of_interest, other_modifications, delta_composition, modification_set):
        self.modifications_of_interest = modifications_of_interest
        self.other_modifications = other_modifications
        self.delta_composition = delta_composition
        self.modification_set = modification_set

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
        if not self.modifications_of_interest.equal_to(other.modifications_of_interest):
            return False
        if not self.other_modifications.equal_to(other.other_modifications):
            return False
        return True


cdef class HCDFragmentationStrategy(PeptideFragmentationStrategyBase):
    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1, include_neutral_losses=False):
        super(HCDFragmentationStrategy, self).__init__(
            peptide, series, chemical_shifts, max_chemical_shifts, include_neutral_losses)
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

    cpdef track_glycosylation(self, long index, glycosylation):
        # HCD strategy does not track intact topologies
        try:
            glycosylation = self._get_core_for(glycosylation)
        except KeyError:
            raise ValueError("Cannot determine which core to use for {}".format(
                glycosylation.rule.glycosylation_type))
        self.glycosylation_manager[self.index] = glycosylation
        self.modification_index.increment(glycosylation, 1)

    cpdef ModificationConfiguration _get_modifications_of_interest(self, PeptideFragment fragment):
        cdef:
            ModificationConfiguration result
            CComposition delta_composition
            CountTable modifications
            CountTableIterator modification_iterator
            CountTable modifications_of_interest, other_modifications
            Py_ssize_t pos, j
            PyObject *key
            PyObject *value
            ModificationBase mod
            long count
        if self._last_modification_set is not None:
            if self._last_modification_set.modification_set == fragment.modification_dict:
                return self._last_modification_set
        modifications = (fragment.modification_dict)
        delta_composition = CComposition._create(None)
        other_modifications = CountTable._create()
        modifications_of_interest = CountTable._create()
        pos = 0
        modification_iterator = CountTableIterator._create(modifications)
        while modification_iterator.has_more():
            pos = modification_iterator.get_next_value(&key, &count)
            mod = <ModificationBase>key
            if mod.name in hcd_modifications_of_interest:
                modifications_of_interest.setitem(mod, count)
                delta_composition.add_from(mod.composition * count)
            else:
                other_modifications.setitem(mod, count)
        result = ModificationConfiguration._create(
            modifications_of_interest, other_modifications,
            delta_composition, modifications)
        return result

    cpdef _replace_cores(self, CountTable modifications_of_interest):
        cdef:
            long n_cores, o_cores, gag_cores, hexnac_cores
        n_cores = modifications_of_interest.delitem(_n_glycosylation)
        o_cores = modifications_of_interest.delitem(_o_glycosylation)
        gag_cores = modifications_of_interest.delitem(_gag_linker_glycosylation)

        # Convert core glycosylation into the residual monosaccharides remaining
        # after dissociation
        hexnac_cores = n_cores * 2 + o_cores
        if hexnac_cores:
            modifications_of_interest.increment(_modification_hexnac, hexnac_cores)
        if gag_cores:
            modifications_of_interest.increment(_modification_xylose, gag_cores)
        return modifications_of_interest

    cpdef list _generate_modification_variants(self, CountTable interesting_modifications, CountTable other_modifications):
        cdef:
            list variants, variant_modification_list
            CountTable varied_modifications
            CountTable updated_modifications
            CountTableIterator it
            size_t i, n, j
            Py_ssize_t pos
            PyObject *key
            PyObject *value
            ModificationBase mod
            long count
        variants = []
        variant_modification_list = descending_combination_counter(interesting_modifications)
        n = PyList_Size(variant_modification_list)
        for i in range(n):
            varied_modifications = <CountTable>PyList_GetItem(variant_modification_list, i)
            updated_modifications = CountTable._create_from(other_modifications)
            extra_composition = CComposition._create(None)
            pos = 0
            j = 0
            it = CountTableIterator._create(varied_modifications)
            while it.has_more():
                pos = it.get_next_value(&key, &count)
                if pos != 0:
                    break
                j += 1
                mod = <ModificationBase>key
                if count != 0:
                    updated_modifications.increment(mod, count)
                    extra_composition.add_from(mod.composition * count)
            variants.append((updated_modifications, extra_composition))
        return variants

    cpdef list partial_loss(self, PeptideFragment fragment):
        cdef:
            CComposition base_composition, delta_composition, new_composition
            ModificationConfiguration mod_config
            list variants, fragments
            IonSeriesBase series
            tuple variant_pair
            CountTable updated_modifications
            size_t i, n
        mod_config = self._get_modifications_of_interest(fragment)
        base_composition = fragment.composition.clone()
        base_composition.subtract_from(mod_config.delta_composition)

        if mod_config is self._last_modification_set:
            variants = self._last_modification_variants
        else:
            self._replace_cores(mod_config.modifications_of_interest)
            self._last_modification_set = mod_config
            variants = self._last_modification_variants = self._generate_modification_variants(
                mod_config.modifications_of_interest,
                mod_config.other_modifications)

        series = self.series
        fragments = []
        n = PyList_Size(variants)
        for i in range(n):
            variant_pair = <tuple>PyList_GetItem(variants, i)
            updated_modifications = <CountTable>PyTuple_GetItem(variant_pair, 0)
            extra_composition = <CComposition>PyTuple_GetItem(variant_pair, 1)
            new_composition = base_composition.clone()
            new_composition.add_from(extra_composition)
            fragments.append(
                PeptideFragment._create(
                    series, fragment.position, updated_modifications.copy(), fragment.bare_mass,
                    flanking_amino_acids=fragment.flanking_amino_acids,
                    glycosylation=None,
                    chemical_shift=None,
                    composition=new_composition))
        return fragments
