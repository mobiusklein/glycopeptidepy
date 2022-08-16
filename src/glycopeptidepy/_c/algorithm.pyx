from cpython.object cimport PyObject
from cpython.ref cimport Py_INCREF, Py_DECREF
from cpython.iterator cimport PyIter_Next
from cpython.int cimport PyInt_AsLong
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem, PyDict_Next, PyDict_Keys, PyDict_Values
from cpython.list cimport PyList_Append, PyList_GetItem, PyList_Size, PyList_SetItem
from cpython.tuple cimport PyTuple_GetItem, PyTuple_Size
from cpython.set cimport PySet_Size

import itertools

cdef object product = itertools.product


cdef class LimitedCrossproduct(object):
    cdef:
        public list collections
        public size_t max_depth
        public size_t size
        public list precomputed
        public object iterator

    def __init__(self, collections, max_depth=4):
        self.collections = [list(v) for v in collections]
        self.max_depth = PyInt_AsLong(max_depth)
        self.size = PyList_Size(self.collections)

        self.precomputed = self.compose_iterative()
        self.iterator = iter(self.precomputed)

    def compose(self):
        i = 0
        current = []
        depth = 0
        for res in self.compose_inner(i, current, depth):
            yield res

    def compose_inner(self, layer_index, current, depth):
        if layer_index == self.size:
            yield current, depth
            return
        layer = self.collections[layer_index]
        for val in layer:
            new_depth = depth + (val is not None)
            if new_depth <= self.max_depth:
                for res in self.compose_inner(layer_index + 1, current + [val], new_depth):
                    yield res
            else:
                continue

    cdef list compose_iterative(self):
        cdef:
            list previous, current, layer, rec
            tuple rec_depth
            size_t depth, new_depth
            size_t layers_i, layers_n
            size_t previous_i, previous_n
            size_t current_layer_i, current_layer_n

        previous = [([], 0)]
        layers_i = 0
        layers_n = self.size
        for layers_i in range(layers_n):
            current = []
            layer = <list>PyList_GetItem(self.collections, layers_i)
            current_layer_n = PyList_Size(layer)
            previous_n = PyList_Size(previous)
            for previous_i in range(previous_n):
                rec_depth = <tuple>PyList_GetItem(previous, previous_i)
                rec = <list>PyTuple_GetItem(rec_depth, 0)
                depth = PyInt_AsLong(<object>PyTuple_GetItem(rec_depth, 1))
                for current_layer_i in range(current_layer_n):
                    val = <object>PyList_GetItem(layer, current_layer_i)
                    new_depth = depth + (val is not None)
                    if new_depth <= self.max_depth:
                        current.append((rec + [val], new_depth))
            previous = current
        return previous

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.iterator)


cdef class ModificationSiteAssignmentCombinator(object):
    cdef:
        public dict modification_to_site
        public dict site_to_modification
        public size_t max_modifications

    def __init__(self, variable_site_map, max_modifications=4):
        self.modification_to_site = variable_site_map
        self.site_to_modification = self.transpose_sites()
        self.max_modifications = PyInt_AsLong(max_modifications)

    cdef dict transpose_sites(self):
        """Given a dictionary mapping between modification names and
        an iterable of valid sites, create a dictionary mapping between
        modification names and a list of valid sites plus the constant `None`

        Returns
        -------
        dict
        """
        cdef:
            dict sites
            PyObject* pkey
            PyObject* pval
            PyObject* tmp
            Py_ssize_t pos
            object mod, site
            set varsites
            list bucket
            size_t i, n

        sites = dict()
        pos = 0
        while PyDict_Next(self.modification_to_site, &pos, &pkey, &pval):
            mod = <object>pkey
            varsites = <set>pval
            for site in varsites:
                tmp = PyDict_GetItem(sites, site)
                if tmp == NULL:
                    bucket = []
                    PyDict_SetItem(sites, site, bucket)
                else:
                    bucket = <list>tmp
                PyList_Append(bucket, mod)

        pos = 0
        while PyDict_Next(sites, &pos, &pkey, &pval):
            bucket = <list>pval
            PyList_Append(bucket, None)

        return sites


    cdef list _remove_empty_sites(self, list sites, list selected):
        cdef:
            size_t i, n
            list result
            tuple pair
            object p1, p2
        n = PyList_Size(selected)
        result = []
        for i in range(n):
            p2 = <object>PyList_GetItem(selected, i)
            if p2 is None:
                continue
            p1 = <object>PyList_GetItem(sites, i)
            pair = (p1, p2)
            PyList_Append(result, pair)
        return result

    def assign(self):
        cdef:
            list sites
            list choices
            list assigned
            object iterator
            tuple selected
            object tmp
        sites = PyDict_Keys(self.site_to_modification)
        choices = PyDict_Values(self.site_to_modification)

        combinator = LimitedCrossproduct(choices, self.max_modifications)
        n = PyList_Size(combinator.precomputed)
        for i in range(n):
            selected = <tuple>PyList_GetItem(combinator.precomputed, i)
            assigned = <list>PyTuple_GetItem(selected, 0)
            assigned = self._remove_empty_sites(sites, assigned)
            yield assigned

        # iterator = product(*choices)

        # while True:
        #     tmp = PyIter_Next(iterator)
        #     if (<PyObject*>tmp) == NULL:
        #         break
        #     else:
        #         selected = <tuple>tmp
        #         assigned = self._remove_empty_sites(sites, selected)
        #         if PyList_Size(assigned) > self.max_modifications:
        #             continue
        #         yield assigned

        # for selected in product(*choices):
        #     yield self._remove_empty_sites(sites, selected)

    def __iter__(self):
        return self.assign()