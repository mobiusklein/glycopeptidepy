ctypedef fused sequence_encoded_t:
    str
    object


cdef enum ParserState:
    start,
    n_term,
    aa,
    mod,
    c_term,


cdef object _cstring_sequence_tokenizer(str sequence, object implicit_n_term=*, object implicit_c_term=*, object glycan_parser_function=*)

cdef object _sequence_tokenizer(sequence_encoded_t sequence, object implicit_n_term=*, object implicit_c_term=*,
                                object glycan_parser_function=*)

cpdef sequence_tokenizer(object sequence, object implicit_n_term=*, object implicit_c_term=*,
                         object glycan_parser_function=*)

cpdef parse_simple(str sequence)