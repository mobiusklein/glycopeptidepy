from cpython cimport Py_INCREF, PY_MAJOR_VERSION
from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_Append, PyList_GetItem, PyList_SetItem, PyList_New, PyList_SET_ITEM
from cpython.dict cimport PyDict_SetItem, PyDict_Keys, PyDict_Values
from cpython.sequence cimport PySequence_GetSlice

from glycopeptidepy._c.structure.constants cimport Configuration

from glypy.structure.glycan_composition import FrozenGlycanComposition

cdef object glycan_parser = FrozenGlycanComposition.parse

from glycopeptidepy.structure import constants as _structure_constants
from glycopeptidepy._c.compat cimport PyStr_AsString, PyStr_FromStringAndSize

cdef Configuration structure_constants = _structure_constants


# this function should use some form of linked list for accumulating positions
# when the sequence is especially long. The repeated reallocation of the list
# storing most of the sequence repeatedly seems like a likely performance bottleneck
cdef object _cstring_sequence_tokenizer(str sequence, object implicit_n_term=None, object implicit_c_term=None, object glycan_parser_function=None):
    '''A simple stateful sequence parser implementing a formally context-free language
    describing components of a polypeptide sequence with N-, C- and internal modifications,
    as well as a glycan composition written at the end.

    Return
    ------
    chunks: list
    mods: list
    n_term: str
    c_term: str
    glycan: :class:`glypy.composition.glycan_composition.FrozenGlycanComposition`
    '''
    cdef:
        ParserState state
        str current_aa, current_mod
        char next_char
        char* csequence
        object glycan
        list chunks, current_mods
        int paren_level
        size_t i, n
        Py_ssize_t mod_start, mod_end


    if glycan_parser_function is None:
        glycan_parser_function = glycan_parser

    state = ParserState.start  # [start, n_term, aa, mod, c_term]
    if implicit_n_term is None:
        n_term = structure_constants.N_TERM_DEFAULT
    else:
        n_term = implicit_n_term
    if implicit_c_term is None:
        c_term = structure_constants.C_TERM_DEFAULT
    else:
        c_term = implicit_c_term
    chunks = []
    glycan = ""
    current_aa = ""
    current_mod = ""
    current_mods = None
    paren_level = 0
    i = 0
    mod_start = -1
    mod_end = -1
    n = len(sequence)
    csequence = PyStr_AsString(sequence)
    while i < n:
        next_char = csequence[i]

        # Transition from aa to mod when encountering the start of a modification
        # internal to the sequence
        if next_char == "(":
            if state == ParserState.aa:
                state = ParserState.mod
                paren_level += 1
                mod_start = i + 1


            # Transition to n_term when starting on an open parenthesis
            elif state == ParserState.start:
                state = ParserState.n_term
                paren_level += 1
                mod_start = i + 1

            else:
                paren_level += 1
                if mod_start == -1:
                    mod_start = i + 1

        elif next_char == ")":
            if state == ParserState.aa:
                raise Exception(
                    "Invalid Sequence. ) found outside of modification, Position {0}. {1}".format(i, sequence))
            else:
                paren_level -= 1
                if paren_level == 0:
                    mod_end = i
                    current_mod = PyStr_FromStringAndSize(&csequence[mod_start], mod_end - mod_start)
                    mod_start = -1
                    mod_end = -1
                    if current_mods is None:
                        current_mods = [current_mod]
                    else:
                        current_mods.append(current_mod)
                    if state == ParserState.mod:
                        state = ParserState.aa
                        # If we encounter multiple modifications in separate parentheses
                        # in a row, attach them to the previous position
                        if current_aa == "":
                            chunks[-1][1].extend(current_mods)
                        else:
                            chunks.append([current_aa, current_mods])

                    elif state == ParserState.n_term:
                        if csequence[i + 1] != "-":
                            raise Exception("Malformed N-terminus for " + sequence)
                        # Only one modification on termini
                        n_term = current_mod
                        state = ParserState.aa
                        # Jump ahead past - into the amino acid sequence
                        i += 1

                    elif state == ParserState.c_term:
                        # Only one modification on termini
                        c_term = current_mod

                    current_mods = None
                    current_mod = ""
                    current_aa = ""
                    mod_start = -1
                    mod_end = -1
                else:
                    mod_end = i

        elif next_char == "|":
            if state == ParserState.aa:
                raise Exception("Invalid Sequence. | found outside of modification")
            else:
                mod_end = i
                current_mod = PyStr_FromStringAndSize(&csequence[mod_start], mod_end - mod_start)
                if current_mods is None:
                    current_mods = [current_mod]
                else:
                    current_mods.append(current_mod)
                mod_start = i + 1
                mod_end = -1
                current_mod = ""

        elif next_char == "{":
            if (state == ParserState.aa or (state == ParserState.c_term and paren_level == 0)):
                glycan = PyStr_FromStringAndSize(&csequence[i], n - i)
                break
            elif (state in {ParserState.mod, ParserState.n_term, ParserState.c_term}) and paren_level > 0:
                mod_end = i

        elif next_char == "-":
            if state == ParserState.aa:
                state = ParserState.c_term
                if(current_aa != ""):
                    chunks.append([current_aa, current_mods])
                    current_mod = ""
                    current_mods = None
                    current_aa = ""
            else:
                pass

        elif state == ParserState.start:
            state = ParserState.aa
            current_aa = PyStr_FromStringAndSize(&next_char, 1)
        elif state == ParserState.aa:
            if(current_aa != ""):
                chunks.append([current_aa, current_mods])
                current_mod = ""
                current_aa = ""
                current_mods = None
            current_aa = PyStr_FromStringAndSize(&next_char, 1)
        elif state in {ParserState.n_term, ParserState.mod, ParserState.c_term}:
            mod_end = i
        else:
            raise Exception(
                "Unknown Tokenizer State", current_aa, current_mod, i, next_char)
        i += 1
    if current_aa != "":
        chunks.append([current_aa, current_mods])

    if glycan != "":
        try:
            glycan = glycan_parser_function(glycan)
        except Exception as e:
            print(e, glycan)
    else:
        glycan = None

    return chunks, None, glycan, n_term, c_term

cdef object _sequence_tokenizer(sequence_encoded_t sequence, object implicit_n_term=None, object implicit_c_term=None, object glycan_parser_function=None):
    '''A simple stateful sequence parser implementing a formally context-free language
    describing components of a polypeptide sequence with N-, C- and internal modifications,
    as well as a glycan composition written at the end.

    Return
    ------
    chunks: list
    mods: list
    n_term: str
    c_term: str
    glycan: :class:`glypy.composition.glycan_composition.FrozenGlycanComposition`
    '''
    cdef:
        ParserState state
        sequence_encoded_t current_aa, current_mod, next_char
        object glycan
        list chunks, current_mods
        int paren_level, i, n

    if glycan_parser_function is None:
        glycan_parser_function = glycan_parser

    if sequence_encoded_t is str:
        return _cstring_sequence_tokenizer(sequence, implicit_n_term, implicit_c_term, glycan_parser_function)

    state = ParserState.start  # [start, n_term, aa, mod, c_term]
    if implicit_n_term is None:
        n_term = structure_constants.N_TERM_DEFAULT
    else:
        n_term = implicit_n_term
    if implicit_c_term is None:
        c_term = structure_constants.C_TERM_DEFAULT
    else:
        c_term = implicit_c_term
    chunks = []
    glycan = ""
    current_aa = ""
    current_mod = ""
    current_mods = None
    paren_level = 0
    i = 0
    n = len(sequence)
    while i < n:
        next_char = sequence[i]

        # Transition from aa to mod when encountering the start of a modification
        # internal to the sequence
        if next_char == "(":
            if state == ParserState.aa:
                state = ParserState.mod
                assert paren_level == 0
                paren_level += 1

            # Transition to n_term when starting on an open parenthesis
            elif state == ParserState.start:
                state = ParserState.n_term
                paren_level += 1

            else:
                paren_level += 1
                if not (state in {ParserState.n_term, ParserState.c_term} and paren_level == 1):
                    current_mod += next_char

        elif next_char == ")":
            if state == ParserState.aa:
                raise Exception(
                    "Invalid Sequence. ) found outside of modification, Position {0}. {1}".format(i, sequence))
            else:
                paren_level -= 1
                if paren_level == 0:
                    if current_mods is None:
                        current_mods = [current_mod]
                    else:
                        current_mods.append(current_mod)
                    if state == ParserState.mod:
                        state = ParserState.aa
                        # If we encounter multiple modifications in separate parentheses
                        # in a row, attach them to the previous position
                        if current_aa == "":
                            chunks[-1][1].extend(current_mods)
                        else:
                            chunks.append([current_aa, current_mods])

                    elif state == ParserState.n_term:
                        if sequence[i + 1] != "-":
                            raise Exception("Malformed N-terminus for " + sequence)
                        # Only one modification on termini
                        n_term = current_mod
                        state = ParserState.aa
                        # Jump ahead past - into the amino acid sequence
                        i += 1

                    elif state == ParserState.c_term:
                        # Only one modification on termini
                        c_term = current_mod

                    current_mods = None
                    current_mod = ""
                    current_aa = ""
                else:
                    current_mod += next_char

        elif next_char == "|":
            if state == ParserState.aa:
                raise Exception("Invalid Sequence. | found outside of modification")
            else:
                if current_mods is None:
                    current_mods = [current_mod]
                else:
                    current_mods.append(current_mod)
                current_mod = ""

        elif next_char == "{":
            if (state == ParserState.aa or (state == ParserState.c_term and paren_level == 0)):
                glycan = sequence[i:]
                break
            elif (state in {ParserState.mod, ParserState.n_term, ParserState.c_term}) and paren_level > 0:
                current_mod += next_char

        elif next_char == "-":
            if state == ParserState.aa:
                state = ParserState.c_term
                if(current_aa != ""):
                    if current_mods is None:
                        current_mods = [current_mod]
                    else:
                        current_mods.append(current_mod)
                    chunks.append([current_aa, current_mods])
                    current_mod = ""
                    current_mods = None
                    current_aa = ""
            else:
                current_mod += next_char

        elif state == ParserState.start:
            state = ParserState.aa
            current_aa = next_char
        elif state == ParserState.aa:
            if(current_aa != ""):
                if current_mods is None:
                    current_mods = [current_mod]
                else:
                    current_mods.append(current_mod)
                chunks.append([current_aa, current_mods])
                current_mod = ""
                current_mods = None
                current_aa = ""
            current_aa = next_char
        elif state in {ParserState.n_term, ParserState.mod, ParserState.c_term}:
            current_mod += next_char
        else:
            raise Exception(
                "Unknown Tokenizer State", current_aa, current_mod, i, next_char)
        i += 1
    if current_aa != "":
        chunks.append([current_aa, current_mods])

    if glycan != "":
        try:
            glycan = glycan_parser_function(glycan)
        except Exception as e:
            print(e, glycan)
    else:
        glycan = None

    return chunks, None, glycan, n_term, c_term


cpdef sequence_tokenizer(object sequence, object implicit_n_term=None, object implicit_c_term=None, object glycan_parser_function=None):
    if isinstance(sequence, str):
        return _sequence_tokenizer[str](<str>sequence, implicit_n_term, implicit_c_term, glycan_parser_function)
    if PY_MAJOR_VERSION > 2:
        if isinstance(sequence, bytes):
            return _sequence_tokenizer[str](
                (<bytes>sequence).decode('utf8'), implicit_n_term, implicit_c_term, glycan_parser_function)
    else:
        if isinstance(sequence, unicode):
            return _sequence_tokenizer[str](
                (<unicode>sequence).encode('utf8'), implicit_n_term, implicit_c_term, glycan_parser_function)
    return _sequence_tokenizer[object](sequence, implicit_n_term, implicit_c_term, glycan_parser_function)


cpdef parse_simple(str sequence):
    cdef:
        size_t i, n
        list chunks, chunk, junk
        str current_aa
        char* csequence

    n = len(sequence)
    chunks = PyList_New(n)
    junk = []
    csequence = PyStr_AsString(sequence)
    for i in range(n):
        next_char = csequence[i]
        if next_char == '(' or next_char == '{':
            raise ValueError("Cannot parse a modified peptide sequence with parse_simple!")
        current_aa = PyStr_FromStringAndSize(&next_char, 1)
        chunk = [current_aa, junk]
        Py_INCREF(chunk)
        PyList_SET_ITEM(chunks, i, chunk)

    return chunks, [], None, structure_constants.N_TERM_DEFAULT, structure_constants.C_TERM_DEFAULT
