cimport cython

from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.count_table cimport CountTable
from glycopeptidepy._c.structure.fragment cimport IonSeriesBase, PeptideFragment, ChemicalShiftBase
from glycopeptidepy._c.structure.sequence_methods cimport _PeptideSequenceCore
from glycopeptidepy._c.structure.base cimport ModificationBase, AminoAcidResidueBase, SequencePosition
from glycopeptidepy._c.structure.modification.modification cimport ModificationInstanceBase

from glycopeptidepy._c.structure.fragmentation_strategy.base cimport FragmentationStrategyBase

from glycopeptidepy.structure.modification import (
    Modification,
    NGlycanCoreGlycosylation,
    OGlycanCoreGlycosylation,
    GlycosaminoglycanLinkerGlycosylation)

from glycopeptidepy.structure.glycan import (GlycosylationType)


cdef:
    dict hcd_modifications_of_interest
    dict hcd_modification_compositions
    dict hcd_modifications_of_interest_to_variants_cache


cdef class PeptideFragmentationStrategyBase(FragmentationStrategyBase):
    cdef:
        public IonSeriesBase series
        public dict chemical_shift_rules
        public int max_chemical_shifts
        public int direction

        public bint include_neutral_losses

        public double running_mass
        public double running_delta_mass
        public CComposition running_composition
        public long index
        public long size
        public CountTable modification_index
        public dict glycosylation_manager
        public CountTable amino_acids_counter

    cpdef _initialize_fields(self)
    cpdef _initialize_start_terminal(self)
    cpdef list _get_viable_chemical_shift_combinations(self)

    cpdef CComposition composition_of(self, SequencePosition position)
    cpdef long name_index_of(self)
    cpdef list flanking_residues(self)
    cpdef track_glycosylation(self, long index, glycosylation)

    cpdef bint has_more(self)
    cpdef list _build_fragments(self)
    cpdef list partial_loss(self, PeptideFragment fragment)
    cpdef object step(self)

    cpdef _update_state(self)
    cpdef reset(self)


@cython.freelist(1000)
cdef class ModificationConfiguration(object):
    cdef:
        # the modifications of interest to HCDFragmentationStrategy
        public CountTable modifications_of_interest
        # the modifications not of interest to HCDFragmentationStrategy
        public CountTable other_modifications
        # the composition shift associated with the interesting modifications
        public CComposition delta_composition
        # the set of all modifications, both interesting and not interesting, used as
        # an identity key for this configuration
        public CountTable modification_set
        # the mass shift associated with the other modifications
        public double other_modifications_mass

    @staticmethod
    cdef ModificationConfiguration _create(CountTable modifications_of_interest, CountTable other_modifications,
                                           CComposition delta_composition, CountTable modification_set,
                                           double other_modifications_mass)

    cdef bint equal_to(self, ModificationConfiguration other)


cdef class HCDFragmentationStrategy(PeptideFragmentationStrategyBase):
    cdef:
        public ModificationConfiguration _last_modification_set
        public list _last_modification_variants

    cpdef _get_core_for(self, ModificationInstanceBase glycosylation)
    cpdef ModificationConfiguration _get_modifications_of_interest(self, PeptideFragment fragment)
    cpdef _replace_cores(self, CountTable modifications_of_interest)
    cpdef list _generate_modification_variants(self, CountTable interesting_modifications, CountTable other_modifications)

    cpdef double[::1] mass_series(self)
