# cython: embedsignature=True

from collections import defaultdict
from itertools import product, combinations

cimport cython

from cpython.ref cimport Py_INCREF
from cpython cimport PyObject
from cpython.tuple cimport PyTuple_GetItem
from cpython.list cimport PyList_GetItem, PyList_SetItem, PyList_Size, PyList_New, PyList_GET_ITEM, PyList_GET_SIZE
from cpython.dict cimport (PyDict_GetItem, PyDict_SetItem, PyDict_Next,
                           PyDict_Keys, PyDict_Update, PyDict_DelItem, PyDict_Size)

from cpython.int cimport PyInt_AsLong, PyInt_Check, PyInt_FromLong

from cython.view cimport array as cvarray

from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.collectiontools cimport descending_combination_counter
from glycopeptidepy._c.count_table cimport CountTable, CountTableIterator
from glycopeptidepy._c.structure.constants cimport Configuration
from glycopeptidepy._c.structure.fragment cimport IonSeriesBase, PeptideFragment, ChemicalShiftBase
from glycopeptidepy._c.structure.sequence_methods cimport _PeptideSequenceCore
from glycopeptidepy._c.structure.base cimport ModificationBase, AminoAcidResidueBase, SequencePosition
from glycopeptidepy._c.structure.modification.modification cimport ModificationInstanceBase

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

    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1, include_neutral_losses=False,
                 compute_compositions=False, **kwargs):
        if chemical_shifts is None:
            chemical_shifts = dict()
        elif not isinstance(chemical_shifts, dict):
            chemical_shifts = dict(chemical_shifts)
        super(PeptideFragmentationStrategyBase, self).__init__(peptide, compute_compositions, **kwargs)
        self.series = IonSeries(series)
        self.chemical_shift_rules = chemical_shifts
        self.max_chemical_shifts = max_chemical_shifts
        self.include_neutral_losses = include_neutral_losses
        self._initialize_fields()

    cpdef _initialize_fields(self):
        self.direction = self.series.direction

        # Null Values
        self.running_mass = 0
        self.running_delta_mass = 0
        if self.compute_compositions:
            self.running_composition = CComposition()
        else:
            self.running_composition = None
        self.index = -1
        self.size = self.peptide.get_size()

        self.modification_index = CountTable._create()
        self.glycosylation_manager = dict()
        self.amino_acids_counter = CountTable._create()

        self.running_mass += self.series.mass_shift
        if self.compute_compositions:
            self.running_composition.add_from(self.series.composition_shift)
        self._initialize_start_terminal()

    cpdef list _get_viable_chemical_shift_combinations(self):
        cdef:
            bint has_neutral_loss_rules
            list shifts
            PyObject* ptmp

        shifts = []
        if self.include_neutral_losses:
            shifts.append(ammonia_loss)
            shifts.append(water_loss)
        if not self.chemical_shift_rules:
            return shifts

        for residue, count in self.amino_acids_counter.items():
            ptmp = PyDict_GetItem(self.chemical_shift_rules, residue)
            if ptmp == NULL:
                pass
            else:
                loss_composition = <list?>ptmp
                if PyList_GET_SIZE(loss_composition):
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
            if self.compute_compositions:
                self.running_composition.add_from(self.peptide.get_n_term().get_composition())
            self.index = -1
        elif self.direction < 0:
            self.running_mass += self.peptide.get_c_term().mass
            if self.compute_compositions:
                self.running_composition.add_from(self.peptide.get_c_term().get_composition())
            self.index = self.peptide.get_size()
        else:
            raise ValueError("Unknown direction %r" % (self.series.direction,))

    cpdef CComposition composition_of(self, SequencePosition position):
        cdef:
            size_t i, n
            CComposition composition
            ModificationBase mod

        n = SequencePosition.get_modification_count(position)
        composition = position.amino_acid.composition.clone()
        for i in range(n):
            mod = SequencePosition.get_modification(position.modifications, i)
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
        self.running_delta_mass += glycosylation.mass

    cpdef _update_state(self):
        cdef:
            SequencePosition position
            ModificationBase mod
            size_t i, n
        self.index += self.direction
        position = self.peptide[self.index]
        n = SequencePosition.get_modification_count(position)
        for i in range(n):
            mod =  SequencePosition.get_modification(position, i)
            if mod.is_tracked_for(ModificationCategory_glycosylation):
                self.track_glycosylation(self.index, mod)
            else:
                self.modification_index.increment(mod, 1)
                self.running_delta_mass += mod.mass
        # self.amino_acids_counter.increment(position.amino_acid, 1)
        self.running_mass += position.amino_acid.mass
        if self.compute_compositions:
            composition = self.composition_of(position)
            self.running_composition.add_from(composition)

    cpdef list _build_fragments(self):
        cdef:
            list fragments_from_site
            list partial_loss_fragments
            list shifts
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
            composition=self.running_composition.clone() if self.compute_compositions else None,
            delta_mass=&self.running_delta_mass)

        shifts = None
        m = 0
        if PyDict_Size(self.chemical_shift_rules) or self.include_neutral_losses:
            shifts = self._get_viable_chemical_shift_combinations()
            m = PyList_GET_SIZE(shifts)

        partial_loss_fragments = self.partial_loss(frag)
        if not m:
            return partial_loss_fragments
        n = PyList_GET_SIZE(partial_loss_fragments)
        fragments_from_site = PyList_New(n * (m + 1))
        k = 0
        for i in range(n):
            fragment = <PeptideFragment>PyList_GET_ITEM(partial_loss_fragments, i)
            Py_INCREF(fragment)
            PyList_SetItem(fragments_from_site, k, fragment)
            k += 1
            for j in range(m):
                shift = <ChemicalShiftBase>PyList_GET_ITEM(shifts, j)
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

    cpdef reset(self):
        self._initialize_fields()


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


hcd_modifications_of_interest_to_variants_cache = dict()


@cython.freelist(100)
cdef class ModificationConfiguration(object):

    @staticmethod
    cdef ModificationConfiguration _create(CountTable modifications_of_interest, CountTable other_modifications,
                                           CComposition delta_composition, CountTable modification_set,
                                           double other_modifications_mass):
        cdef ModificationConfiguration inst = ModificationConfiguration.__new__(ModificationConfiguration)
        inst.modifications_of_interest = modifications_of_interest
        inst.other_modifications = other_modifications
        inst.delta_composition = delta_composition
        inst.modification_set = modification_set
        inst.other_modifications_mass = other_modifications_mass
        return inst

    def __init__(self, modifications_of_interest, other_modifications, delta_composition, modification_set, other_modifications_mass):
        self.modifications_of_interest = modifications_of_interest
        self.other_modifications = other_modifications
        self.delta_composition = delta_composition
        self.modification_set = modification_set
        self.other_modifications_mass = other_modifications_mass

    def __iter__(self):
        yield self.modifications_of_interest
        yield self.delta_composition
        yield self.other_modifications

    def __repr__(self):
        return "{self.__class__.__name__}({self.modifications_of_interest}, {self.other_modifications})".format(
            self=self)

    def __eq__(self, other):
        if not isinstance(other, ModificationConfiguration):
            return NotImplemented
        else:
            return self.equal_to(<ModificationConfiguration>other)

    cdef bint equal_to(self, ModificationConfiguration other):
        if other is None:
            return False
        if not self.modifications_of_interest.equal_to(other.modifications_of_interest):
            return False
        if not self.other_modifications.equal_to(other.other_modifications):
            return False
        return True


cdef class HCDFragmentationStrategy(PeptideFragmentationStrategyBase):
    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1, include_neutral_losses=False, compute_compositions=False, **kwargs):
        super(HCDFragmentationStrategy, self).__init__(
            peptide, series, chemical_shifts, max_chemical_shifts, include_neutral_losses,
            compute_compositions, **kwargs)
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

    @property
    def modification_variants_cache(self):
        return hcd_modifications_of_interest_to_variants_cache

    @modification_variants_cache.setter
    def modification_variants_cache(self, value):
        global hcd_modifications_of_interest_to_variants_cache
        hcd_modifications_of_interest_to_variants_cache = dict(value)

    cpdef CComposition composition_of(self, SequencePosition position):
        cdef:
            CComposition composition
            ModificationInstanceBase mod
            size_t i, n
        composition = position.amino_acid.composition.clone()
        n = SequencePosition.get_modification_count(position)
        for i in range(n):
            mod = <ModificationInstanceBase>SequencePosition.get_modification(position, i)
            if mod.is_tracked_for(ModificationCategory_glycosylation):
                mod = self._get_core_for(mod)
            composition.add_from(mod.composition)
        return composition

    cpdef track_glycosylation(self, long index, glycosylation):
        cdef:
            ModificationBase glycosylation_core

        # HCD strategy does not track intact topologies or compositions, just cores
        glycosylation_core = self._get_core_for(glycosylation)

        self.glycosylation_manager[self.index] = glycosylation_core
        self.modification_index.increment(glycosylation_core, 1)
        self.running_delta_mass += glycosylation_core.mass

    cpdef ModificationConfiguration _get_modifications_of_interest(self, PeptideFragment fragment):
        cdef:
            ModificationConfiguration result
            CComposition delta_composition
            CountTable modifications
            CountTableIterator modification_iterator
            CountTable modifications_of_interest, other_modifications
            Py_ssize_t pos, j
            PyObject *key
            ModificationBase mod
            long count
            double other_modifications_mass
        other_modifications_mass = 0
        if self._last_modification_set is not None:
            if self._last_modification_set.modification_set == fragment.modification_dict:
                return self._last_modification_set
        modifications = (fragment.modification_dict)
        if self.compute_compositions:
            delta_composition = CComposition._create(None)
        else:
            delta_composition = None
        other_modifications = CountTable._create()
        modifications_of_interest = CountTable._create()

        pos = 0
        modification_iterator = CountTableIterator._create(modifications)
        while modification_iterator.has_more():
            pos = modification_iterator.get_next_value(&key, &count)
            mod = <ModificationBase>key
            if mod.name in hcd_modifications_of_interest:
                modifications_of_interest.setitem(mod, count)
                if self.compute_compositions:
                    delta_composition.add_from(mod.composition * count)
            else:
                other_modifications.setitem(mod, count)
                other_modifications_mass += mod.mass * count
        result = ModificationConfiguration._create(
            modifications_of_interest, other_modifications,
            delta_composition, modifications, other_modifications_mass)
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
            double delta_mass
            long count
        variants = []
        variant_modification_list = descending_combination_counter(interesting_modifications)
        n = PyList_Size(variant_modification_list)
        for i in range(n):
            varied_modifications = <CountTable>PyList_GetItem(variant_modification_list, i)
            updated_modifications = CountTable._create_from(other_modifications)
            if self.compute_compositions:
                extra_composition = CComposition._create(None)
            else:
                extra_composition = None
            delta_mass = 0.0
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
                    delta_mass += mod.mass * count
                    if self.compute_compositions:
                        extra_composition.add_from(mod.composition * count)
            variants.append((updated_modifications, extra_composition, delta_mass))
        return variants

    cpdef list partial_loss(self, PeptideFragment fragment):
        cdef:
            CComposition base_composition, delta_composition, new_composition
            ModificationConfiguration mod_config
            list variants, fragments
            IonSeriesBase series
            tuple variant_pair
            CountTable updated_modifications
            PyObject* ptemp
            double delta_mass
            size_t i, n
        mod_config = self._get_modifications_of_interest(fragment)
        if self.compute_compositions:
            base_composition = fragment.composition.clone()
            # remove the composition shift associated with the interesting modifications
            base_composition.subtract_from(mod_config.delta_composition)
        else:
            base_composition = None

        # If the modification configuration hasn't changed, reuse the existing variants.
        #
        # If this is the first time entering :meth:`partial_loss`, :attr:`_last_modification_set`
        # is :const:`None`, and :meth:`_get_modifications_of_interest` can never return :const:`None`
        # so this block will only be entered after at least one pass through the else clause.
        if mod_config is self._last_modification_set:
            variants = self._last_modification_variants
        else:
            # Recalculate the modifications according to the HCD rules, and update the instance's cache
            # properties.
            #
            # Check the class-wide variant cache to see if these variants have already been calculated,
            # and if so, reuse them, otherwise calculate them and update the class cache.
            self._replace_cores(mod_config.modifications_of_interest)
            self._last_modification_set = mod_config
            ptemp = PyDict_GetItem(
                hcd_modifications_of_interest_to_variants_cache,
                self._last_modification_set.modifications_of_interest)
            # if the combination is new, generate them and update the cache
            if ptemp == NULL:
                variants = self._last_modification_variants = self._generate_modification_variants(
                    mod_config.modifications_of_interest,
                    # do not include extra modifications here so that caching of
                    # interesting modificatoin variants does not need to be recalculated
                    # when uninteresting modifications change
                    CountTable._create()
                )
                PyDict_SetItem(
                    hcd_modifications_of_interest_to_variants_cache,
                    self._last_modification_set.modifications_of_interest,
                    variants)
            # otherwise reuse the cached version
            else:
                variants = self._last_modification_variants = <list>ptemp

        series = self.series

        n = PyList_Size(variants)
        # If there is only one variant, it must be the unaltered form, so just return that
        if n == 1:
            return [fragment]
        fragments = []
        for i in range(n):
            variant_pair = <tuple>PyList_GetItem(variants, i)
            updated_modifications = CountTable._create_from(<CountTable>PyTuple_GetItem(variant_pair, 0))
            # add the current fragment's uninteresting modifications back on top of the
            # varied modifications here.
            updated_modifications._add_from(mod_config.other_modifications)
            if self.compute_compositions:
                extra_composition = <CComposition>PyTuple_GetItem(variant_pair, 1)
                new_composition = base_composition.clone()
                new_composition.add_from(extra_composition)
            else:
                new_composition = None
            delta_mass = <object>PyTuple_GetItem(variant_pair, 2) + mod_config.other_modifications_mass
            fragments.append(
                PeptideFragment._create(
                    series, fragment.position, updated_modifications, fragment.bare_mass,
                    flanking_amino_acids=fragment.flanking_amino_acids,
                    glycosylation=None,
                    chemical_shift=None,
                    composition=new_composition,
                    delta_mass=&delta_mass
                ))
        return fragments

    @cython.boundscheck(False)
    cpdef double[::1] mass_series(self):
        cdef:
            double[::1] masses
            size_t i

        masses = cvarray(shape=(self.size - 1,), itemsize=sizeof(double), format='d')

        i = 0
        while self.has_more():
            self._update_state()
            masses[i] = self.running_mass
            i += 1
        return masses
