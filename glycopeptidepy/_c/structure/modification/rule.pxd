from glypy.composition.ccomposition cimport CComposition
from glycopeptidepy._c.structure.base cimport ModificationBase
from glypy.utils.cenum cimport EnumValue

cdef class ModificationRuleBase(ModificationBase):
    cdef:
        public set names
        public list categories

        public basestring preferred_name
        public basestring title
        public basestring common_name
        public basestring unimod_name

        public Py_hash_t _hash

    cpdef basestring serialize(self)




cdef class NeutralLossBase(object):
    cdef:
        public CComposition composition
        public double mass
        public object label


cpdef bint is_tracked_for_glycosylation(self, object category)


cdef class GlycosylationBase(ModificationRuleBase):
    cdef:
        public bint _is_composition
        public bint _is_core
        public EnumValue _glycosylation_type

    cdef bint get_is_composition(self)
    cdef void set_is_composition(self, bint value)

    cdef bint get_is_core(self)
    cdef void set_is_core(self, bint value)

    cdef EnumValue get_glycosylation_type(self)
    cdef void set_glycosylation_type(self, EnumValue glycosylation_type)
