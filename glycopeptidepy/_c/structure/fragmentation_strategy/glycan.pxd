from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.count_table cimport CountTable

from glycopeptidepy._c.structure.fragmentation_strategy.base cimport FragmentationStrategyBase


cdef class GlycanCompositionFragment(object):
    cdef:
        public double mass
        public CComposition composition
        public CountTable key
        public bint is_extended
        public Py_hash_t _hash_key

    cpdef GlycanCompositionFragment copy(self)
    cpdef _get_glycan_composition(self)

    @staticmethod
    cdef GlycanCompositionFragment _create(double mass, CComposition composition, CountTable key, bint is_extended)


cdef class GlycanCompositionFragmentStrategyBase(FragmentationStrategyBase):
    cdef:
        public bint _use_query
        public object _generator

    cpdef glycan_composition(self)
    cpdef bint _guess_query_mode(self, glycan_composition)
    cpdef long count_glycosylation_type(self, glycotype)


cdef class StubGlycopeptideStrategy(GlycanCompositionFragmentStrategyBase):
    cdef:
        public bint extended
        public bint extended_fucosylation

    cpdef GlycanCompositionFragment fucosylate_increment(self, GlycanCompositionFragment shift)
    cpdef GlycanCompositionFragment xylosylate_increment(self, GlycanCompositionFragment shift)
    cpdef GlycanCompositionFragment modified_increment(self, GlycanCompositionFragment base, object shift)
    cpdef list fucosylate_extended(self, GlycanCompositionFragment shift, long fucose_count)
    cpdef bint _validate_glycan_composition(self, CountTable aggregate_glycosylation, glycan)
    cpdef list _combinate_sites(self, list per_site_shifts, object glycan)

    cpdef list n_glycan_composition_fragments(self, object glycan, int core_count=*, int iteration_count=*)
    cpdef list o_glycan_composition_fragments(self, glycan, long core_count=*, long iteration_count=*)
    cpdef list gag_linker_composition_fragments(self, glycan, long core_count=*, long iteration_count=*)
