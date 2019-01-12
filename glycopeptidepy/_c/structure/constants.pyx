
cdef class Configuration(object):
    def __init__(self):
        self.N_TERM_DEFAULT = "H"
        self.C_TERM_DEFAULT = "OH"
        self.FRAG_OFFSET = 1
        self.ALLOW_MODIFIED_ASPARAGINE = False
        self.PARTIAL_HEXNAC_LOSS = True
        self.EXCLUDE_B1 = True
