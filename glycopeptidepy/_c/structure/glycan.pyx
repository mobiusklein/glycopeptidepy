cimport cython
from cpython.object cimport PyObject
from cpython.list cimport PyList_Append, PyList_Size
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem, PyDict_Next, PyDict_Size

from glypy.composition.ccomposition cimport CComposition
from glypy._c.structure.glycan_composition cimport _CompositionBase
from glypy.structure.glycan_composition import HashableGlycanComposition

from glycopeptidepy._c.structure.modification.modification cimport ModificationInstanceBase
from glycopeptidepy._c.structure.modification.rule cimport ModificationRuleBase, GlycosylationBase
from glycopeptidepy._c.structure.base cimport PeptideSequenceBase

@cython.freelist(10000)
cdef class GlycanCompositionWithOffsetProxyBase(object):
    def __init__(self, obj, offset=None):
        if offset is None:
            offset = CComposition._create(None)
        self.obj = obj
        self.composition_offset = offset



cdef CComposition WATER_OFFSET = CComposition({"H": 2, "O": 1})


cdef object GlycanCompositionProxyType = None


def set_glycan_composition_proxy_type(tp):
    global GlycanCompositionProxyType
    GlycanCompositionProxyType = tp
    return GlycanCompositionProxyType


@cython.freelist(10000)
cdef class GlycosylationManager(object):
    def __init__(self, parent, aggregate=None):
        self._init()
        self.parent = parent
        if aggregate is not None:
            self.set_aggregate(aggregate)

    @staticmethod
    cdef GlycosylationManager _create(PeptideSequenceBase parent, object aggregate):
        cdef:
            GlycosylationManager self
        self = GlycosylationManager.__new__(GlycosylationManager)
        self._init()
        self.parent = parent
        if aggregate is not None:
            self.set_aggregate(aggregate)
        return self

    cdef void _init(self):
        self.mapping = {}
        self._aggregate = None
        self._proxy = None
        self._type_track = None
        self._total_glycosylation_size = -1

    def __getitem__(self, key):
        return self.mapping[key]

    def __setitem__(self, key, value):
        self.invalidate()
        self.mapping[key] = value

    def __delitem__(self, key):
        self.invalidate()
        del self.mapping[key]

    def __len__(self):
        return self.get_size()

    cdef size_t get_size(self):
        return PyDict_Size(self.mapping)

    cpdef keys(self):
        return self.mapping.keys()

    cpdef values(self):
        return self.mapping.values()

    cpdef items(self):
        return self.mapping.items()

    cpdef pop(self, key, default=None):
        self.invalidate()
        self.mapping.pop(key, default)

    cpdef clear(self):
        self.invalidate()
        self.set_aggregate(None)
        self.mapping.clear()

    cpdef update(self, other):
        self.invalidate()
        if isinstance(other, GlycosylationManager):
            self.mapping.update(other.mapping)
        else:
            self.mapping.update(other)

    cpdef GlycosylationManager copy(self):
        cdef:
            GlycosylationManager inst
            object aggregate
        aggregate = self.get_aggregate()

        inst = GlycosylationManager._create(
            self.parent,
            aggregate.clone() if aggregate is not None else None)
        inst.update(self)
        return inst

    cpdef GlycosylationManager clone(self):
        return self.copy()

    cpdef object _patch_aggregate(self):
        cdef:
            object aggregate
            CComposition offset
        aggregate = self.get_aggregate()
        if isinstance(aggregate, GlycanCompositionWithOffsetProxyBase):
            offset = (<GlycanCompositionWithOffsetProxyBase>aggregate).composition_offset
        elif isinstance(aggregate, _CompositionBase):
            offset = (<_CompositionBase>aggregate)._composition_offset
        else:
            offset = <CComposition>(aggregate).composition_offset
        offset.subtract_from(WATER_OFFSET)

    cpdef object invalidate(self):
        self._proxy = None
        self._type_track = None
        self._total_glycosylation_size = -1

    cdef object get_aggregate(self):
        return self._aggregate

    cdef int set_aggregate(self, object value):
        self._aggregate = value
        if self._aggregate is not None:
            self._patch_aggregate()
        self.invalidate()
        return 0

    @property
    def aggregate(self):
        return self.get_aggregate()

    @aggregate.setter
    def aggregate(self, value):
        self.set_aggregate(value)

    cdef dict get_glycosylation_types(self):
        if self._type_track is None:
            self._type_track = self._update_track()
        return self._type_track

    cdef dict _update_track(self):
        cdef:
            dict track
            Py_ssize_t pos
            PyObject* pkey
            PyObject* pval
            PyObject* tmp
            list bucket
            GlycosylationBase glycosylation

        track = dict()
        pos = 0
        while PyDict_Next(self.mapping, &pos, &pkey, &pval):
            glycosylation = <GlycosylationBase>(<ModificationInstanceBase>pval).rule
            glycosylation_type = glycosylation.get_glycosylation_type()
            tmp = PyDict_GetItem(track, glycosylation_type)
            if tmp == NULL:
                bucket = []
                PyDict_SetItem(track, glycosylation_type, bucket)
            else:
                bucket = <list>tmp
            PyList_Append(bucket, (<object>pkey, <object>pval))
        return  track

    cpdef CComposition total_composition(self):
        cdef:
            CComposition total
            bint has_aggregate
            PyObject* pkey
            PyObject* pval
            Py_ssize_t pos
            ModificationInstanceBase value
            object aggregate
        aggregate = self.get_aggregate()
        total = CComposition._create(None)
        has_aggregate = aggregate is not None
        pos = 0
        while PyDict_Next(self.mapping, &pos, &pkey, &pval):
            value = <ModificationInstanceBase>pval
            if has_aggregate and (<GlycosylationBase>(value.rule)).is_core:
                continue
            total.add_from(value.composition)
        if has_aggregate:
            total.add_from(aggregate.total_composition())
        return total

    @property
    def glycan_composition(self):
        return self.get_glycan_composition()

    cdef object get_glycan_composition(self):
        if self._proxy is None:
            self._proxy = self._make_glycan_composition_proxy()
        return self._proxy

    cpdef double mass(self):
        cdef:
            double total
            bint has_aggregate
            PyObject* pkey
            PyObject* pval
            Py_ssize_t pos
            ModificationInstanceBase value
            object aggregate
        total = 0.
        aggregate = self.get_aggregate()
        has_aggregate = aggregate is not None
        pos = 0
        while PyDict_Next(self.mapping, &pos, &pkey, &pval):
            value = <ModificationInstanceBase>pval
            if has_aggregate and (<GlycosylationBase>(value.rule)).is_core:
                continue
            total += value.mass
        if has_aggregate:
            total += aggregate.mass()
        return total

    cdef int total_glycosylation_size(self):
        if self._total_glycosylation_size == -1:
            self._total_glycosylation_size = sum(self.get_aggregate().values())
        return self._total_glycosylation_size

    cpdef bint is_fully_specified_topologies(self):
        cdef:
            bint is_fully_specified
            PyObject* pkey
            PyObject* pval
            Py_ssize_t pos
            ModificationInstanceBase value
            GlycosylationBase rule

        is_fully_specified = self.get_size() > 0
        pos = 0
        while PyDict_Next(self.mapping, &pos, &pkey, &pval):
            value = <ModificationInstanceBase>pval
            rule = <GlycosylationBase>value.rule
            if rule.is_core:
                is_fully_specified = False
                break
            elif rule.is_composition:
                is_fully_specified = False
                break
        return is_fully_specified

    cdef object _make_glycan_composition_proxy(self):
        cdef:
            object aggregate
            object base
            ModificationInstanceBase value
            GlycosylationBase rule
            Py_ssize_t pos
            PyObject* pkey
            PyObject* pval

        aggregate = self.get_aggregate()
        if aggregate is not None:
            base = aggregate.clone()
        else:
            base = HashableGlycanComposition()
            # Represent the initial amide bond between the peptide
            # and the first glycan. Subsequent glycans do not need
            # further chemical losses because of the dehyration built
            # directly into the Residue abstraction.
            base.composition_offset -= WATER_OFFSET
        pos = 0
        while PyDict_Next(self.mapping, &pos, &pkey, &pval):
            value = <ModificationInstanceBase>pval
            rule = <GlycosylationBase>(value.rule)
            if rule.is_core:
                continue
            elif rule.is_composition:
                base += rule.glycan
            else:
                # Convert Glycan object into a composition, using the original
                # detatched topology to omit the "aglycone" group which represents
                # the connection between the glycan and the peptide, which penalizes
                # the composition by H2O. This H2O is lost when that bond is formed,
                # but doesn't need to be explicitly included as the loss is tracked
                # when initializing the base above.
                gc = HashableGlycanComposition.from_glycan(rule._original)
                base += gc
        return GlycanCompositionProxyType(base)

    cpdef int count_glycosylation_type(self, glycosylation_type):
        cdef:
            dict tracks
            PyObject* tmp
        tracks = self.get_glycosylation_types()
        tmp = PyDict_GetItem(tracks, glycosylation_type)
        if tmp == NULL:
            return 0
        else:
            return PyList_Size(<list>tmp)

    def __repr__(self):
        return "{self.__class__.__name__}({self.mapping})".format(self=self)