from glycresoft_sqlalchemy.structure import residue as pyresidue, modification as pymodification
from cpython.ref cimport PyObject
from cpython.string cimport PyString_AsString, PyString_FromString
from cpython.float cimport PyFloat_AsDouble
from libc.stdlib cimport abort, malloc, free, realloc, calloc
from libc.math cimport fabs
from libc.string cimport strcmp, strdup, strlen, strcat
from libc cimport *

cdef extern from * nogil:
    void qsort (void *base, unsigned short n, unsigned short w, int (*cmp_func)(void*, void*))
    int printf   (const char *template, ...)
    int sprintf  (char *s, const char *template, ...)


cdef AminoAcidStruct* get_residue_by_symbol(char* symbol, AminoAcidStructArray* residues) nogil:
    cdef:
        size_t i
        AminoAcidStruct* residue
    for i in range(residues.size):
        if strcmp(residues.residues[i].symbol, symbol) == 0:
            return &residues.residues[i]
    
cdef AminoAcidStruct* get_residue_by_name(char* name, AminoAcidStructArray* residues) nogil:
    cdef:
        size_t i
        AminoAcidStruct* residue
    for i in range(residues.size):
        if strcmp(residues.residues[i].name, name) == 0:
            return &residues.residues[i]

cdef ModificationStruct* get_modification(char* name, ModificationStructArray* modifications) nogil:
    cdef:
        size_t i
        ModificationStruct* mod
    for i in range(modifications.size):
        if strcmp(modifications.modifications[i].name, name) == 0:
            return &modifications.modifications[i]

cdef inline ModificationStruct* _get_modification(char* name) nogil:
    return get_modification(name, Modifications)

cdef inline AminoAcidStruct* _get_residue_by_symbol(char* symbol) nogil:
    return get_residue_by_symbol(symbol, AminoAcids)

cdef inline AminoAcidStruct* _get_residue_by_name(char* name) nogil:
    return get_residue_by_name(name, AminoAcids)

cdef inline double sequence_position_mass(PeptideSequencePositionStruct* position) nogil:
    cdef:
        size_t i
        double mass
    mass = 0
    mass += position.residue.mass
    for i in range(position.modifications.size):
        mass += position.modifications.modifications[i].mass
    return mass

cdef size_t _total_seqstring_length(PeptideSequenceStruct* seq) nogil:
    cdef:
        size_t size, i, j
    size = 0
    for i in range(seq.size):
        size += strlen(seq.sequence[i].residue.symbol)
        for j in range(seq.sequence[i].modifications.size):
            size += strlen(seq.sequence[i].modifications.modifications[j].name) + 2
    if seq.glycan != NULL:
        size += strlen(seq.glycan.name)
    return size

cdef char* sequence_to_string(PeptideSequenceStruct* seq) nogil:
    cdef:
        size_t size, i, j
        char* result
        char open_mod
        char close_mod
    open_mod = "("
    close_mod = ")"
    size = _total_seqstring_length(seq)

    result = <char*>malloc(size)
    result[0] = "\0"

    for i in range(seq.size):
        strcat(result, seq.sequence[i].residue.symbol)
        for j in range(seq.sequence[i].modifications.size):
            strcat(result, &open_mod)
            strcat(result, seq.sequence[i].modifications.modifications[j].name)
            strcat(result, &close_mod)
    if seq.glycan != NULL:
        strcat(result, seq.glycan.name)
    return result

cdef FragmentIonStructArray* fragment_series(PeptideSequenceStruct* sequence, char* series, int direction, double offset) nogil:
    cdef:
        size_t i
        double mass
        FragmentIonStructArray* result

    result = <FragmentIonStructArray*>malloc(sizeof(FragmentIonStructArray))
    result.fragments = <FragmentIonStruct*>malloc(sizeof(FragmentIonStruct) * sequence.size - 1)
    result.size = sequence.size - 1
    mass = offset
    if direction > 0:
        for i in range(sequence.size):
            mass += sequence.sequence[i].mass
            if i == 0:
                continue
            result.fragments[i - 1].series = series
            result.fragments[i - 1].mass = mass
            result.fragments[i - 1].index = i
    elif direction < 0:
        for i in range(sequence.size - 1, -1, -1):
            mass += sequence.sequence[i].mass
            if i == sequence.size:
                continue
            result.fragments[i - 1].series = series
            result.fragments[i - 1].mass = mass
            result.fragments[i - 1].index = sequence.size - i
    return result

cdef char* fragment_key(FragmentIonStruct* fragment) nogil:
    cdef char* buff
    buff = <char*>malloc(sizeof(char)*20)
    sprintf("%s%d", buff, fragment.series, fragment.index)
    return buff


# Convert Functions

cdef AminoAcidStructArray* create_amino_acid_table(dict residue_map):
    cdef:
        str key
        object value
        AminoAcidStructArray* residues
        size_t size, i
    size = len(residue_map)
    residues = <AminoAcidStructArray*>malloc(sizeof(AminoAcidStructArray))
    residues.residues = <AminoAcidStruct*>malloc(sizeof(AminoAcidStruct) * size)
    residues.size = size

    i = 0
    for key, value in residue_map.items():
        residues.residues[i].name = strdup(PyString_AsString(value.name))
        residues.residues[i].symbol = strdup(PyString_AsString(value.symbol))
        residues.residues[i].mass = PyFloat_AsDouble(value.mass)
        i += 1

    return residues

cdef ModificationStructArray* create_modification_table(dict modification_map):
    cdef:
        object key
        object value
        ModificationStructArray* modifications
        size_t size, i

    size = len(modification_map)
    modifications = <ModificationStructArray*>malloc(sizeof(ModificationStructArray))
    modifications.modifications = <ModificationStruct*>malloc(sizeof(ModificationStruct) * size)
    modifications.size = size

    i = 0
    for key, value in modification_map.items():
        modifications.modifications[i].name = strdup(PyString_AsString(value.preferred_name.encode('ascii', 'ignore')))
        modifications.modifications[i].mass = PyFloat_AsDouble(value.mass)
        i += 1
    return modifications

cdef inline ModificationStructArray* get_modification_list_from_sequence(list modlist):
    cdef:
        ModificationStructArray* array
        size_t size, i, j
    size = len(modlist)
    array = <ModificationStructArray*>malloc(sizeof(ModificationStructArray))
    array.modifications = <ModificationStruct*>malloc(sizeof(ModificationStruct) * size)
    array.size = size

    for i in range(size):
        array.modifications[i] = _get_modification(PyString_AsString(modlist[i].rule.preferred_name))[0]
    return array

cdef PeptideSequenceStruct* sequence_from_object(object obj):
    cdef:
        size_t size, i, j
        PeptideSequenceStruct* seq
        double mass_accum
        list obj_pos

    mass_accum = 0.
    size = len(obj)
    seq = <PeptideSequenceStruct*>malloc(sizeof(PeptideSequenceStruct))
    seq.sequence = <PeptideSequencePositionStruct*>malloc(sizeof(PeptideSequencePositionStruct)*size)
    seq.size = size

    for i in range(size):
        obj_pos = obj[i]
        seq.sequence[i].residue = _get_residue_by_symbol(obj_pos[0].symbol)
        seq.sequence[i].modifications = get_modification_list_from_sequence(obj_pos[1])
        seq.sequence[i].mass = sequence_position_mass(&seq.sequence[i])
    if not isinstance(obj.glycan, str) and obj.glycan is not None:
        seq.glycan = glycan_from_object(obj.glycan)
    else:
        seq.glycan = NULL
    seq.mass = PyFloat_AsDouble(obj.mass)
    return seq

cdef SimpleGlycanStruct* glycan_from_object(object obj):
    cdef:
        SimpleGlycanStruct* glycan

    glycan = <SimpleGlycanStruct*>malloc(sizeof(SimpleGlycanStruct))
    glycan.name = strdup(PyString_AsString(str(obj)))
    glycan.mass = PyFloat_AsDouble(obj.mass())
    return glycan


AminoAcids = create_amino_acid_table({k: pyresidue.Residue(k) for k in pyresidue.residue_table})
Modifications = create_modification_table(dict(**pymodification.Modification._table))


# Free Functions

cdef void free_sequence(PeptideSequenceStruct* seq) nogil:
    cdef:
        size_t i
    for i in range(seq.size):
        free_sequence_position(&seq.sequence[i])

    free_simple_glycan(seq.glycan)
    free(seq)

cdef void free_sequence_position(PeptideSequencePositionStruct* pos) nogil:
    free_modification_array(pos.modifications)
    free(pos)

cdef void free_modification_array(ModificationStructArray* mod_array) nogil:
    free(mod_array.modifications)
    free(mod_array)

cdef void free_amino_acid_array(AminoAcidStructArray* aa_array) nogil:
    free(aa_array.residues)
    free(aa_array)

cdef void free_fragment_ion_array(FragmentIonStructArray* frag_array) nogil:
    free(frag_array.fragments)
    free(frag_array)

cdef void free_simple_glycan(SimpleGlycanStruct* g) nogil:
    free(g)

#

def test(name):
    residue = _get_residue_by_name(name)
    return residue[0]

def test_residue_sym(symbol):
    residue = _get_residue_by_symbol(symbol)
    return residue[0]

def test_modifications(mapping, str name):
    cdef:
        ModificationStructArray* array
        ModificationStruct* modification
    array = create_modification_table(mapping)
    modification = get_modification(PyString_AsString(name), array)
    print modification.name
    return modification[0]

def test_sequence(obj):
    print(sequence_to_string(sequence_from_object(obj)))


# cpdef list parse(str sequence, str n_term="H", str c_term="OH"):
#     cdef:
#         list chunks = []
#         list current_chunk = []
#         list mods = []
#         str current_aa
#         str current_mod
#         str next_char
#         size_t i, length
#         # start: 0
#         # aa: 1
#         # mod: 2
#         # n_term: 3
#         # c_term: 4
#         # glycan: 5
#         str state = "start"
#         int paren_level = 0

#     while i < len(sequence):
#         next_char = sequence[i]

#         # Transition from aa to mod when encountering the start of a modification
#         # internal to the sequence
#         if next_char == "(":
#             if state == "aa":
#                 state = "mod"
#                 assert paren_level == 0
#                 paren_level += 1

#             # Transition to n_term when starting on an open parenthesis
#             elif state == "start":
#                 state = "n_term"
#                 paren_level += 1

#             else:
#                 paren_level += 1
#                 if not (state in {"n_term", "c_term"} and paren_level == 1):
#                     current_mod += next_char

#         elif next_char == ")":
#             if state == "aa":
#                 raise Exception(
#                     "Invalid Sequence. ) found outside of modification, Position {0}. {1}".format(i, sequence))
#             else:
#                 paren_level -= 1
#                 if paren_level == 0:
#                     mods.append(current_mod)
#                     current_mods.append(current_mod)
#                     if state == "mod":
#                         state = 'aa'
#                         # If we encounter multiple modifications in separate parentheses
#                         # in a row, attach them to the previous position
#                         if current_aa == "":
#                             chunks[-1][1].extend(current_mods)
#                         else:
#                             chunks.append([current_aa, current_mods])

#                     elif state == "n_term":
#                         if sequence[i+1] != "-":
#                             raise Exception("Malformed N-terminus for " + sequence)
#                         # Only one modification on termini
#                         n_term = current_mod
#                         state = "aa"
#                         # Jump ahead past - into the amino acid sequence
#                         i += 1

#                     elif state == "c_term":
#                         # Only one modification on termini
#                         c_term = current_mod

#                     current_mods = []
#                     current_mod = ""
#                     current_aa = ""
#                 else:
#                     current_mod += next_char

#         elif next_char == "|":
#             if state == "aa":
#                 raise Exception("Invalid Sequence. | found outside of modification")
#             else:
#                 current_mods.append(current_mod)
#                 mods.append(current_mod)
#                 current_mod = ""

#         elif next_char == "{":
#             if (state == 'aa' or (state == "c_term" and paren_level == 0)):
#                 glycan = sequence[i:]
#                 break
#             elif (state in {"mod", "n_term", "c_term"}) and paren_level > 0:
#                 current_mod += next_char

#         elif next_char == "-":
#             if state == "aa":
#                 state = "c_term"
#                 if(current_aa != ""):
#                     current_mods.append(current_mod)
#                     chunks.append([current_aa, current_mods])
#                     current_mod = ""
#                     current_mods = []
#                     current_aa = ""
#             else:
#                 current_mod += next_char

#         elif state == "start":
#             state = "aa"
#             current_aa = next_char
#         elif state == "aa":
#             if(current_aa != ""):
#                 current_mods.append(current_mod)
#                 chunks.append([current_aa, current_mods])
#                 current_mod = ""
#                 current_mods = []
#                 current_aa = ""
#             current_aa = next_char
#         elif state in {"n_term", "mod", "c_term"}:
#             current_mod += next_char
#         else:
#             raise Exception(
#                 "Unknown Tokenizer State", current_aa, current_mod, i, next_char)
#         i += 1
#     if current_aa != "":
#         current_mods.append("")
#         chunks.append([current_aa, current_mods])
#     if current_mod != "":
#         mods.append(current_mod)

#     #if glycan != "":
#     #    try:
#     #        glycan = glycan_parser(glycan)
#     #    except Exception, e:
#     #        logging.exception("Error in parser, %s and %s", glycan, sequence, exc_info=e)

#     # return chunks, mods, glycan, n_term, c_term

#     return chunks
