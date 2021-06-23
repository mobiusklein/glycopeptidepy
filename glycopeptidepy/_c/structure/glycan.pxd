from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.structure.base cimport PeptideSequenceBase

cdef class GlycanCompositionWithOffsetProxyBase(object):
    cdef:
        public object obj
        public CComposition composition_offset


cdef class GlycosylationManager(object):
    cdef:
        public dict mapping
        public object parent
        public object _aggregate
        public object _proxy
        public dict _type_track
        public int _total_glycosylation_size

    cpdef update(self, other)
    cpdef clear(self)
    cpdef keys(self)
    cpdef values(self)
    cpdef items(self)
    cpdef pop(self, key, default=*)

    cdef size_t get_size(self)

    cpdef GlycosylationManager copy(self)
    cpdef GlycosylationManager clone(self)

    cpdef object invalidate(self)

    cdef object get_aggregate(self)
    cdef int set_aggregate(self, object value)
    cpdef object _patch_aggregate(self)

    cdef dict get_glycosylation_types(self)
    cdef dict _update_track(self)

    cpdef CComposition total_composition(self)
    cpdef double mass(self)
    cpdef bint is_fully_specified_topologies(self)

    cdef object get_glycan_composition(self)
    cdef object _make_glycan_composition_proxy(self)
    cdef int total_glycosylation_size(self)

    cpdef int count_glycosylation_type(self, glycosylation_type)

    cdef void _init(self)

    @staticmethod
    cdef GlycosylationManager _create(PeptideSequenceBase parent, aggregate)