from glypy.composition.ccomposition cimport CComposition
from glypy._c.structure.glycan_composition cimport _CompositionBase

from glycopeptidepy._c.structure.base cimport PeptideSequenceBase

ctypedef fused glycan_composition_type:
    _CompositionBase
    GlycanCompositionProxy
    object


cdef class GlycanCompositionProxy(object):
    '''A mapping-like object that imitates the GlycanComposition interface in
    a read-only fashion.
    '''

    cdef:
        public _CompositionBase obj
        public str _serialized

    cpdef object _getitem_fast(self, key)

    cpdef str serialize(self)

    @staticmethod
    cdef GlycanCompositionProxy _create(glycan_composition_type obj)


cdef class GlycanCompositionWithOffsetProxy(GlycanCompositionProxy):
    cdef:
        public CComposition composition_offset


cdef class GlycosylationManager(object):
    cdef:
        public dict mapping
        public object parent
        public object _aggregate
        public GlycanCompositionProxy _proxy
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

    cdef GlycanCompositionProxy get_glycan_composition(self)
    cdef GlycanCompositionProxy _make_glycan_composition_proxy(self)
    cdef int get_total_glycosylation_size(self)
    cpdef int total_glycosylation_size(self)

    cpdef int count_glycosylation_type(self, glycosylation_type)

    cdef void _init(self)

    @staticmethod
    cdef GlycosylationManager _create(PeptideSequenceBase parent, aggregate)