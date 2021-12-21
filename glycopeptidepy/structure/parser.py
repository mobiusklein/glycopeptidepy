import logging
import sys
from . import constants as structure_constants
from . import PeptideSequenceBase
from .modification import NGlycanCoreGlycosylation


from glypy.structure.glycan_composition import FrozenGlycanComposition

_n_glycan_core = NGlycanCoreGlycosylation()
glycan_parser = FrozenGlycanComposition.parse


def sequence_tokenizer(sequence, implicit_n_term=None, implicit_c_term=None, glycan_parser_function=None):
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
    if glycan_parser_function is None:
        glycan_parser_function = glycan_parser
    if sys.version_info[0] > 2 and isinstance(sequence, bytes):
        sequence = sequence.decode("utf8")

    state = "start"  # [start, n_term, aa, mod, c_term]
    n_term = implicit_n_term or structure_constants.N_TERM_DEFAULT
    c_term = implicit_c_term or structure_constants.C_TERM_DEFAULT
    mods = []
    chunks = []
    glycan = ""
    current_aa = ""
    current_mod = ""
    current_mods = []
    paren_level = 0
    i = 0
    while i < len(sequence):
        next_char = sequence[i]

        # Transition from aa to mod when encountering the start of a modification
        # internal to the sequence
        if next_char == "(":
            if state == "aa":
                state = "mod"
                assert paren_level == 0
                paren_level += 1

            # Transition to n_term when starting on an open parenthesis
            elif state == "start":
                state = "n_term"
                paren_level += 1

            else:
                paren_level += 1
                if not (state in {"n_term", "c_term"} and paren_level == 1):
                    current_mod += next_char

        elif next_char == ")":
            if state == "aa":
                raise Exception(
                    "Invalid Sequence. ) found outside of modification, Position {0}. {1}".format(i, sequence))
            else:
                paren_level -= 1
                if paren_level == 0:
                    mods.append(current_mod)
                    current_mods.append(current_mod)
                    if state == "mod":
                        state = 'aa'
                        # If we encounter multiple modifications in separate parentheses
                        # in a row, attach them to the previous position
                        if current_aa == "":
                            chunks[-1][1].extend(current_mods)
                        else:
                            chunks.append([current_aa, current_mods])

                    elif state == "n_term":
                        if sequence[i + 1] != "-":
                            raise Exception("Malformed N-terminus for " + sequence)
                        # Only one modification on termini
                        n_term = current_mod
                        state = "aa"
                        # Jump ahead past - into the amino acid sequence
                        i += 1

                    elif state == "c_term":
                        # Only one modification on termini
                        c_term = current_mod

                    current_mods = []
                    current_mod = ""
                    current_aa = ""
                else:
                    current_mod += next_char

        elif next_char == "|":
            if state == "aa":
                raise Exception("Invalid Sequence. | found outside of modification")
            else:
                current_mods.append(current_mod)
                mods.append(current_mod)
                current_mod = ""

        elif next_char == "{":
            if (state == 'aa' or (state == "c_term" and paren_level == 0)):
                glycan = sequence[i:]
                break
            elif (state in {"mod", "n_term", "c_term"}) and paren_level > 0:
                current_mod += next_char

        elif next_char == "-":
            if state == "aa":
                state = "c_term"
                if(current_aa != ""):
                    current_mods.append(current_mod)
                    chunks.append([current_aa, current_mods])
                    current_mod = ""
                    current_mods = []
                    current_aa = ""
            else:
                current_mod += next_char

        elif state == "start":
            state = "aa"
            current_aa = next_char
        elif state == "aa":
            if(current_aa != ""):
                current_mods.append(current_mod)
                chunks.append([current_aa, current_mods])
                current_mod = ""
                current_mods = []
                current_aa = ""
            current_aa = next_char
        elif state in {"n_term", "mod", "c_term"}:
            current_mod += next_char
        else:
            raise Exception(
                "Unknown Tokenizer State", current_aa, current_mod, i, next_char)
        i += 1
    if current_aa != "":
        current_mods.append("")
        chunks.append([current_aa, current_mods])
    if current_mod != "":
        mods.append(current_mod)

    if glycan != "":
        try:
            glycan = glycan_parser_function(glycan)
        except Exception as e:
            logging.exception("Error in parser, %s and %s", glycan, sequence, exc_info=e)
    else:
        glycan = None

    return chunks, mods, glycan, n_term, c_term


try:
    _sequence_tokenizer = sequence_tokenizer
    from glycopeptidepy._c.parser import sequence_tokenizer
except ImportError:
    pass


def prefix_to_postfix_modifications(sequence):
    '''Transform a chunk-list from modification-on-the-left to
    modification-on-the-right association rules.

    This code does not work with terminal modifications

    Return
    ------
    list
    '''
    sequence = "!" + sequence

    tokens, mods, glycan, n_term, c_term = sequence_tokenizer(sequence)
    last_mods = tokens.pop(0)[1]
    for token in tokens:
        mods = token[1]
        token[1] = last_mods
        last_mods = mods
    return tokens


def sequence_length(sequence):
    '''Compute the peptide length of a sequence string.
    This parses the sequence and counts the number of chunks.

    Return
    ------
    int
    '''
    sequence, modifications, glycan, n_term, c_term = sequence_tokenizer(sequence)
    return len(sequence)


def strip_modifications(sequence):
    '''A wrapper around `sequence_tokenizer` that discards all modifications

    Return
    ------
    str
    '''
    if isinstance(sequence, PeptideSequenceBase):
        sequence = str(sequence)
    chunks, _, glycan, n_term, c_term = sequence_tokenizer(sequence)
    return ''.join([residue for residue, mod in chunks])


def sequence_tokenizer_respect_sequons(sequence, known_sites=None):
    '''A wrapper around `sequence_tokenizer` that will treat an n-glycan sequon as a single unit

    Return
    ------
    list
    '''
    if known_sites is None:
        known_sites = tuple()
    chunks, _, glycan, n_term, c_term = sequence_tokenizer(sequence)
    positions = []
    i = 0
    sequence_length = len(chunks)
    while(i < sequence_length):
        cur_pos = chunks[i]
        if (cur_pos[0] == "N" and cur_pos[1] and _n_glycan_core in cur_pos[1]) or (i in known_sites):
            positions.append(chunks[i:(i + 3)])
            i += 2
        else:
            positions.append(cur_pos)
        i += 1
    return positions


def parse_simple(sequence):
    '''Parse a simple sequence string without *any* modifications
    into a chunked token list, matching the return signature of :func:`sequence_tokenizer`
    '''
    if sys.version_info[0] > 2 and isinstance(sequence, bytes):
        sequence = sequence.decode("utf8")

    chunks = [[a, []] for a in sequence]
    mods = None
    glycan = None
    n_term = structure_constants.N_TERM_DEFAULT
    c_term = structure_constants.C_TERM_DEFAULT
    return chunks, mods, glycan, n_term, c_term


try:
    from glycopeptidepy._c.parser import parse_simple
except ImportError:
    pass
