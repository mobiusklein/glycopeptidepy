
cdef class Configuration(object):
    cdef:
        public str N_TERM_DEFAULT
        public str C_TERM_DEFAULT
        public int FRAG_OFFSET
        public bint ALLOW_MODIFIED_ASPARAGINE
        public bint PARTIAL_HEXNAC_LOSS
        public bint EXCLUDE_B1