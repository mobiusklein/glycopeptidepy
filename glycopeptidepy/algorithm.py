from glycopeptidepy.enzyme import cleave
from glycopeptidepy.structure.sequence import PeptideSequence, list_to_sequence
from glycopeptidepy.structure.parser import sequence_tokenizer_respect_sequons


def pair_rotate(sequence):
    """Invert each token pair.

    ABCD -> BADC

    Parameters
    ----------
    sequence : iterable

    Returns
    -------
    list
    """
    chunks = []
    gen = iter(sequence)
    while True:
        chunk = []
        try:
            chunk = [gen.next()]
            chunk.append(gen.next())
            chunks.append(chunk)
        except StopIteration:
            if len(chunk) > 0:
                chunks.append(chunk)
            break
    rev_seq = []
    for chunk in reversed(chunks):
        rev_seq.extend(chunk)
    return rev_seq


def reverse_preserve_sequon(sequence, prefix_len=0, suffix_len=1, peptide_type=None):
    if peptide_type is None:
        peptide_type = PeptideSequence
    if isinstance(sequence, PeptideSequence):
        sequence = str(sequence)
    original = peptide_type(sequence)
    sequence_tokens = sequence_tokenizer_respect_sequons(sequence)
    pref = sequence_tokens[:prefix_len]
    if suffix_len == 0:
        suf = ""
        body = sequence_tokens[prefix_len:]
    else:
        suf = sequence_tokens[-suffix_len:]
        body = sequence_tokens[prefix_len:-suffix_len]
    body = body[::-1]
    rev_sequence = (list_to_sequence(pref + list(body) + suf))
    if str(list_to_sequence(sequence_tokens)) == str(rev_sequence):
        rot_body = pair_rotate(body)
        rev_sequence = (list_to_sequence(pref + list(rot_body) + suf))
    rev_sequence.n_term = original.n_term
    rev_sequence.c_term = original.c_term
    if original.glycan:
        rev_sequence.glycan = original.glycan.clone(propogate_composition_offset=False)
    return rev_sequence
