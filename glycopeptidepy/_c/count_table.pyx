import sys
from libc.stdlib cimport malloc, free, realloc
from cpython.ref cimport PyObject, Py_INCREF, Py_DECREF, Py_XDECREF
from cpython.int cimport PyInt_AsLong, PyInt_FromLong
from cpython.dict cimport PyDict_Next



cdef int initialize_count_table_bin(count_table_bin* bin, size_t size):
    bin.cells = <count_table_bin_cell*>malloc(sizeof(count_table_bin_cell) * size)
    if bin.cells == NULL:
        return 1
    for i in range(size):
        bin.cells[i].key = NULL
    bin.used = 0
    bin.size = size


cdef void free_count_table_bin(count_table_bin* bin):
    for i in range(bin.size):
        if bin.cells[i].key != NULL:
            Py_DECREF(<object>bin.cells[i].key)
            bin.cells[i].key = NULL
    free(bin.cells)


cdef int count_table_bin_append(count_table_bin* bin, PyObject* key, long value):
    if bin.used == bin.size - 1:
        bin.cells = <count_table_bin_cell*>realloc(bin.cells, sizeof(count_table_bin_cell) * bin.size * 2)
        if bin.cells == NULL:
            return 1
        bin.size *= 2
    bin.cells[bin.used].key = key
    bin.cells[bin.used].value = value
    bin.used += 1
    return 0


cdef int count_table_bin_find(count_table_bin* bin, PyObject* query, Py_ssize_t* cell_index):
    for i in range(bin.used):
        if bin.cells[i].key == NULL:
            continue
        Py_INCREF(<object>bin.cells[i].key)
        if (<object>bin.cells[i].key) == (<object>query):
            cell_index[0] = i
            return 0
    cell_index[0] = -1
    return 0


cdef count_table* make_count_table(size_t table_size, size_t bin_size):
    cdef:
        size_t i, j
        count_table* table
    table = <count_table*>malloc(sizeof(count_table))
    if table == NULL:
        raise MemoryError()    
    table.bins = <count_table_bin*>malloc(sizeof(count_table_bin) * table_size)
    if table.bins == NULL:
        raise MemoryError()
    for i in range(table_size):
        initialize_count_table_bin(&(table.bins[i]), bin_size)
    table.size = table_size
    return table


cdef void free_count_table(count_table* table):
    for i in range(table.size):
        free_count_table_bin(&table.bins[i])
    free(table.bins)
    free(table)


cdef int count_table_find_bin(count_table* table, PyObject* query, Py_ssize_t* bin_index):
    bin_index[0] = hash(<object>(query)) % table.size
    assert bin_index[0] < table.size
    return 0


cdef int count_table_put(count_table* table, PyObject* key, long value):
    cdef:
        Py_ssize_t bin_index, cell_index
        int status
    status = count_table_find_bin(table, key, &bin_index)
    assert status == 0
    status = count_table_bin_find(&table.bins[bin_index], key, &cell_index)
    assert status == 0
    if cell_index == -1:
        Py_INCREF(<object>key)
        status = count_table_bin_append(&table.bins[bin_index], key, value)
        if status != 0:
            return 1
    else:
        table.bins[bin_index].cells[cell_index].value = value
    return 0


cdef int count_table_del(count_table* table, PyObject* key, long* value):
    cdef:
        Py_ssize_t bin_index, cell_index
        int status
    status = count_table_find_bin(table, key, &bin_index)
    assert status == 0
    status = count_table_bin_find(&table.bins[bin_index], key, &cell_index)
    assert status == 0
    if cell_index == -1:
        value[0] = 0
        return 0
    else:
        value[0] = table.bins[bin_index].cells[cell_index].value
        Py_DECREF(<object>table.bins[bin_index].cells[cell_index].key)
        table.bins[bin_index].cells[cell_index].key = NULL
        table.bins[bin_index].cells[cell_index].value = 0
        return 0


cdef int count_table_get(count_table* table, PyObject* key, long* value):
    cdef:
        Py_ssize_t bin_index, cell_index
        int status
    status = count_table_find_bin(table, key, &bin_index)
    assert status == 0
    status = count_table_bin_find(&table.bins[bin_index], key, &cell_index)
    assert status == 0
    if cell_index == -1:
        value[0] = 0
        return 0
    else:
        value[0] = table.bins[bin_index].cells[cell_index].value
        return 0


cdef Py_ssize_t count_table_count(count_table* table):
    cdef:
        Py_ssize_t count
        size_t i, j
    count = 0
    for i in range(table.size):
        for j in range(table.bins[i].used):
            count += table.bins[i].cells[j].key != NULL
    return count


cdef list count_table_keys(count_table* table):
    cdef:
        list keys
        size_t i, j
    keys = []
    count = 0
    for i in range(table.size):
        for j in range(table.bins[i].used):
            if table.bins[i].cells[j].key != NULL:
                keys.append(<object>table.bins[i].cells[j].key)
    return keys


cdef list count_table_values(count_table* table):
    cdef:
        list values
        size_t i, j
    values = []
    count = 0
    for i in range(table.size):
        for j in range(table.bins[i].used):
            if table.bins[i].cells[j].key != NULL:
                values.append(PyInt_FromLong(table.bins[i].cells[j].value))
    return values


cdef list count_table_items(count_table* table):
    cdef:
        list items
        size_t i, j
    items = []
    for i in range(table.size):
        for j in range(table.bins[i].used):
            if table.bins[i].cells[j].key != NULL:
                items.append(
                    (<object>table.bins[i].cells[j].key,
                     PyInt_FromLong(table.bins[i].cells[j].value)))
    return items


cdef void count_table_add(count_table* table_a, count_table* table_b):
    cdef:
        size_t i, j
        long value, temp
    for i in range(table_b.size):
        for j in range(table_b.bins[i].used):
            if table_b.bins[i].cells[j].key != NULL:
                count_table_get(table_a, table_b.bins[i].cells[j].key, &temp)
                count_table_get(table_b, table_b.bins[i].cells[j].key, &value)
                count_table_put(table_a, table_b.bins[i].cells[j].key, value + temp)


cdef void count_table_subtract(count_table* table_a, count_table* table_b):
    cdef:
        size_t i, j
        long value, temp
    for i in range(table_b.size):
        for j in range(table_b.bins[i].used):
            if table_b.bins[i].cells[j].key != NULL:
                count_table_get(table_a, table_b.bins[i].cells[j].key, &temp)
                count_table_get(table_b, table_b.bins[i].cells[j].key, &value)
                count_table_put(table_a, table_b.bins[i].cells[j].key, value - temp)


cdef void count_table_update(count_table* table_a, count_table* table_b):
    cdef:
        size_t i, j
        long value
    for i in range(table_b.size):
        for j in range(table_b.bins[i].used):
            if table_b.bins[i].cells[j].key != NULL:
                count_table_get(table_b, table_b.bins[i].cells[j].key, &value)
                count_table_put(table_a, table_b.bins[i].cells[j].key, value)


cdef count_table* count_table_copy(count_table* table_a):
    cdef count_table* dup = make_count_table(table_a.size, 2)
    count_table_add(dup, table_a)
    return dup


cdef bint count_table_equals(count_table* table_a, count_table* table_b):
    cdef:
        size_t i, j
        long value
    for i in range(table_a.size):
        for j in range(table_a.bins[i].used):
            if table_a.bins[i].cells[j].key != NULL:
                count_table_get(table_b, table_a.bins[i].cells[j].key, &value)
                if table_a.bins[i].cells[j].value != value:
                    return False
    for i in range(table_b.size):
        for j in range(table_b.bins[i].used):
            if table_b.bins[i].cells[j].key != NULL:
                count_table_get(table_a, table_b.bins[i].cells[j].key, &value)
                if table_b.bins[i].cells[j].value != value:
                    return False
    return True


cdef class CountTable(object):

    @staticmethod
    cdef CountTable _create():
        cdef CountTable inst
        inst = CountTable.__new__(CountTable)
        return inst

    def __cinit__(self, *args, **kwargs):
        self.table = make_count_table(6, 2)

    def __init__(self, obj=None, **kwargs):
        if obj is not None:
            self.update(obj)
        if kwargs:
            self.update(kwargs)

    cpdef update(self, obj):
        if isinstance(obj, CountTable):
            self._update_from_count_table(<CountTable>obj)
        if isinstance(obj, dict):
            self._update_from_dict(<dict>obj)
        else:
            for k, v in obj.items():
                self[k] = v

    cdef void _update_from_dict(self, dict d):
        cdef:
            PyObject *key
            PyObject *value
            Py_ssize_t pos
        pos = 0
        while PyDict_Next(d, &pos, &key, &value):
            self.setitem(<object>key, PyInt_AsLong(<object>value))

    cdef void _update_from_count_table(self, CountTable other):
        count_table_update(self.table, other.table)

    cpdef CountTable copy(self):
        cdef CountTable inst = CountTable._create()
        inst._update_from_count_table(self)
        return inst

    def __reduce__(self):
        return self.__class__, (self._to_dict(), )

    def __dealloc__(self):
        free_count_table(self.table)

    def __getitem__(self, object key):
        cdef long value = self.getitem(key)
        return PyInt_FromLong(value)

    def __setitem__(self, object key, object value):
        self.setitem(key, PyInt_AsLong(value))

    def __delitem__(self, object key):
        self.pop(key)

    def __len__(self):
        return count_table_count(self.table)

    def __contains__(self, key):
        cdef long value = self.getitem(key)
        return value != 0

    def __iter__(self):
        return iter(self.keys())

    cpdef dict _to_dict(self):
        return dict(self.items())

    def __repr__(self):
        return "{self.__class__.__name__}({content})".format(
            self=self, content=self._to_dict())

    def __eq__(self, other):
        if isinstance(self, CountTable):
            if isinstance(other, CountTable):
                return count_table_equals(self.table, (<CountTable>other).table)
            else:
                return self._to_dict() == other
        else:
            return self == other._to_dict()

    def __ne__(self, other):
        return not (self == other)

    cpdef list keys(self):
        return count_table_keys(self.table)

    cpdef list values(self):
        return count_table_values(self.table)

    cpdef list items(self):
        return count_table_items(self.table)

    cpdef clear(self):
        for key in self.keys():
            self.setitem(key, 0)

    cpdef setdefault(self, key, value):
        if self.getitem(key) == 0:
            self.setitem(key, value)

    cdef long getitem(self, object key):
        cdef long value
        cdef PyObject* pkey = <PyObject*>key
        count_table_get(self.table, pkey, &value)
        return value        

    cdef void setitem(self, object key, long value):
        cdef PyObject* pkey = <PyObject*>key
        count_table_put(self.table, pkey, value)

    cdef long pop(self, object key):
        cdef PyObject* pkey = <PyObject*>key
        cdef long value
        count_table_del(self.table, pkey, &value)
        return value


def main():
    cdef long val
    cdef str string
    cdef object pint
    cdef PyObject* key
    string = "spam"
    pint = 12
    key = <PyObject*>string
    cdef count_table* table = make_count_table(6, 2)
    count_table_put(table, key, 42)
    count_table_get(table, key, &val)
    assert val == 42
    key = <PyObject*>pint
    count_table_put(table, key, 67)
    count_table_del(table, key, &val)
    assert val == 67
    val = 1
    count_table_get(table, key, &val)
    free_count_table(table)
    pytable = CountTable()
    Py_INCREF(<object>key)
    pytable[<object>key] = 2
    pytable[int(12)] = 11
    pytable[255.2] = 99
    assert len(pytable) == 2
    del pytable
    pytable = CountTable._create()
    print pytable
