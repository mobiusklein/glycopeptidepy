from cpython.object cimport PyObject
from cpython.ref cimport Py_INCREF, Py_DECREF
from cpython.iterator cimport PyIter_Next
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem, PyDict_Next, PyDict_Keys, PyDict_Values
from cpython.list cimport PyList_Append, PyList_GetItem, PyList_Size, PyList_SetItem
from cpython.tuple cimport PyTuple_GetItem, PyTuple_Size
from cpython.set cimport PySet_Size

import itertools

cdef object product = itertools.product

cdef class ModificationSiteAssignmentCombinator(object):
    cdef:
        public dict modification_to_site
        public dict site_to_modification

    def __init__(self, variable_site_map):
        self.modification_to_site = variable_site_map
        self.site_to_modification = self.transpose_sites()

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


    cdef list _remove_empty_sites(self, list sites, tuple selected):
        cdef:
            size_t i, n
            list result
            tuple pair
            object p1, p2
        n = PyTuple_Size(selected)
        result = []
        for i in range(n):
            p2 = <object>PyTuple_GetItem(selected, i)
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
            object iterator
            tuple selected
            object tmp
        sites = PyDict_Keys(self.site_to_modification)
        choices = PyDict_Values(self.site_to_modification)
        iterator = product(*choices)

        while True:
            tmp = PyIter_Next(iterator)
            if (<PyObject*>tmp) == NULL:
                break
            else:
                selected = <tuple>tmp
                yield self._remove_empty_sites(sites, selected)

        # for selected in product(*choices):
        #     yield self._remove_empty_sites(sites, selected)

    def __iter__(self):
        return self.assign()