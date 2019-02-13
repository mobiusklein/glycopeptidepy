from cpython cimport Py_INCREF
from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_Append, PyList_GetItem, PyList_SetItem, PyList_New
from cpython.dict cimport PyDict_SetItem, PyDict_Keys, PyDict_Values
from cpython.int cimport PyInt_AsLong
from cpython.tuple cimport PyTuple_GetItem
from cpython.mapping cimport PyMapping_Keys, PyMapping_Values

from glycopeptidepy._c.count_table cimport CountTable

cdef object product

from itertools import product

cpdef list descending_combination_counter(mapping_t counter):
    cdef:
        list keys, values, results
        object interval, v
        list count_ranges, combos
        int i, j, n, k
        CountTable result
        tuple combo
    if mapping_t is dict:
        keys = list(PyDict_Keys(counter))
        values = list(PyDict_Values(counter))
    elif mapping_t is CountTable:
        keys = counter.keys()
        values = counter.values()
    else:
        keys = list(PyMapping_Keys(counter))
        values = list(PyMapping_Values(counter))
    n = PyList_GET_SIZE(keys)
    count_ranges = PyList_New(n)
    for i in range(n):
        v = (<object>PyList_GET_ITEM(values, i))
        interval = range(v + 1)
        Py_INCREF(interval)
        PyList_SetItem(count_ranges, i, interval)

    combos = list(product(*count_ranges))
    k = PyList_GET_SIZE(combos)
    results = PyList_New(k)
    for i in range(k):
        result = CountTable._create()
        combo = <tuple>PyList_GET_ITEM(combos, i)
        for j in range(n):
            result.setitem(
                <object>PyList_GET_ITEM(keys, j),
                PyInt_AsLong(<object>PyTuple_GetItem(combo, j)))
        Py_INCREF(result)
        PyList_SetItem(results, i, result)
    return results
