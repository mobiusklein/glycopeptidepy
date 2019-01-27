from cpython cimport PyObject

cdef struct count_table_bin_cell:
    long value
    PyObject* key


cdef struct count_table_bin:
    count_table_bin_cell* cells
    size_t size
    size_t used


cdef struct count_table:
    count_table_bin* bins
    size_t size


cdef int initialize_count_table_bin(count_table_bin* bin, size_t size)
cdef void free_count_table_bin(count_table_bin* bin)
cdef int count_table_bin_append(count_table_bin* bin, PyObject* key, long value)
cdef int count_table_bin_find(count_table_bin* bin, PyObject* query, Py_ssize_t* cell_index)

cdef count_table* make_count_table(size_t table_size, size_t bin_size)
cdef void free_count_table(count_table* table)
cdef int count_table_find_bin(count_table* table, PyObject* query, Py_ssize_t* bin_index)
cdef int count_table_put(count_table* table, PyObject* key, long value)
cdef int count_table_del(count_table* table, PyObject* key, long* value)
cdef int count_table_get(count_table* table, PyObject* key, long* value)
cdef int count_table_increment(count_table* table, PyObject* key, long value)
cdef int count_table_decrement(count_table* table, PyObject* key, long value)

cdef Py_ssize_t count_table_count(count_table* table)

cdef list count_table_keys(count_table* table)
cdef list count_table_values(count_table* table)
cdef list count_table_items(count_table* table)

cdef void count_table_add(count_table* table_a, count_table* table_b)
cdef void count_table_subtract(count_table* table_a, count_table* table_b)
cdef void count_table_scale(count_table* table, long value)

cdef void count_table_update(count_table* table_a, count_table* table_b)
cdef count_table* count_table_copy(count_table* table_a)

cdef bint count_table_equals(count_table* table_a, count_table* table_b)


cdef class CountTable(object):
    cdef:
        count_table* table

    @staticmethod
    cdef CountTable _create()

    cdef void _update_from_dict(self, dict d)
    cdef void _update_from_count_table(self, CountTable other)

    cdef void _add_from(self, CountTable other)
    cdef void _add_from_dict(self, dict other)
    cdef void _subtract_from(self, CountTable other)
    cdef void _subtract_from_dict(self, dict other)
    cdef void _scale_by(self, long value)

    cpdef dict _to_dict(self)

    cpdef update(self, obj)
    cpdef CountTable copy(self)
    cpdef list items(self)
    cpdef list values(self)
    cpdef list keys(self)
    cpdef clear(self)
    cpdef setdefault(self, key, value)

    cdef void increment(self, object key, long value)
    cdef void decrement(self, object key, long value)
    cdef long getitem(self, object key)
    cdef void setitem(self, object key, long value)
    cdef long pop(self, object key)