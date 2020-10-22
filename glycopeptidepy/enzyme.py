import warnings
import re
import itertools
import operator

from collections import deque

from .structure.sequence import PeptideSequence
from .structure.parser import sequence_length

_get1 = operator.itemgetter(1)
_get12 = operator.itemgetter(1, 2)


def cleave(sequence, protease, missed_cleavages=0, min_length=0, **kwargs):
    '''A reimplementation of pyteomics.parser.cleave which produces leaky cleavages
    of a peptide sequence by a regex rule. Includes the cut indices, not included in
    pyteomics.'''

    if isinstance(protease, str):
        protease = Protease(protease)
    return protease.cleave(sequence, missed_cleavages, min_length, **kwargs)


# Required information for generating peptides that can be cleaved by a given protease
enzymes = {
    "trypsin": {
        "cleavage_start": [""],
        "cleavage_end": ["K", "R"],
        "name": "trypsin"
    }
}


class Protease(object):
    def __init__(self, name, cleavage_start=None, cleavage_end=None):
        self.name = name
        self.cleavage_start = cleavage_start
        self.cleavage_end = cleavage_end
        try:
            self.regex = re.compile(enzyme_rules[name.lower()])
        except KeyError:
            self.regex = re.compile(name)

    def __eq__(self, other):
        try:
            return self.regex == other.regex
        except AttributeError:
            return False

    def __hash__(self):
        return hash(self.regex)

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return "Protease(%s, %s)" % (self.name, self.regex.pattern)

    def missed_cleavages(self, sequence):
        if isinstance(sequence, PeptideSequence):
            sequence = str(sequence)
        return len(self.regex.findall(sequence))

    def _find_sites(self, sequence):
        return itertools.chain(
            map(lambda x: x.end(), self.regex.finditer(sequence)), [None])

    def cleave(self, sequence, missed_cleavages=0, min_length=0, max_length=100,
               semispecific=False):
        if semispecific:
            return self.cleave_semispecific(sequence, missed_cleavages, min_length, max_length)
        peptides = []
        if isinstance(sequence, PeptideSequence):
            sequence = str(sequence)
        sequence_size = sequence_length(sequence)
        # Keep track of the previous indices a cut was made at, with
        # the start of the sequence (index 0) being treated as a possible
        # site. This will be used to backtrack and produce missed cleavages.
        #
        # The maxlen argument to deque is important, as when this value is
        # not None, the deque object automatically regulates its length. As
        # items are added to the end in excess of maxlen, the same number of
        # items are removed from the start of the container. This implicitly
        # enforces the missed cleavage count requirement. The +2 ensures room
        # will be available for at least one previous site to look up from,
        # and the current site being added during each step of the loop below.
        cleavage_sites = deque([0], maxlen=missed_cleavages + 2)
        for i in self._find_sites(sequence):
            # Add the new site to the backtrack list
            cleavage_sites.append(i)
            # Step back through previous cleavage sites before
            # this one
            for j in range(len(cleavage_sites) - 1):
                # Get the sequence between the the jth previous site
                # and the most recent cleavage site
                seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
                # Only deal with this if the sequence is non-empty
                if seq:
                    seq_size = sequence_length(seq)
                    if min_length is None or seq_size >= min_length:

                        # Get parent sequence coordinates
                        start_position = cleavage_sites[j]
                        end_position = cleavage_sites[-1]
                        if end_position is None:
                            end_position = sequence_size
                        # Store the sequence, start position, end position, and # of missed cleavages
                        peptides.append((seq, start_position, end_position, missed_cleavages - j))
        # Deduplicate and sort by starting postion
        return sorted(set(peptides), key=_get1)

    def cleave_semispecific(self, sequence, missed_cleavages=0, min_length=0, max_length=100):
        peptides = []
        if isinstance(sequence, PeptideSequence):
            sequence = str(sequence)

        cleavage_sites = deque([0], maxlen=missed_cleavages + 2)
        all_sites = list(self._find_sites(sequence))
        all_sites[-1] = sequence_length(sequence)
        n_sites = len(all_sites)

        # See the large comment in Protease.cleave for a descrption of this first
        # pair of for-loops
        for i in all_sites:
            # Add the new site to the backtrack list
            cleavage_sites.append(i)
            for j in range(len(cleavage_sites) - 1):
                # Get the sequence between the the jth previous site
                # and the most recent cleavage site
                seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
                # Only deal with this if the sequence is non-empty
                if seq:
                    seq_size = sequence_length(seq)
                    if (min_length is None or seq_size >= min_length) and (
                            max_length is None or seq_size <= max_length):
                        # Get parent sequence coordinates
                        start_position = cleavage_sites[j]
                        end_position = cleavage_sites[-1]
                        if end_position is None:
                            end_position = sequence_length(sequence)
                        # Store the sequence, start position, end position, and # of missed cleavages
                        peptides.append((seq, start_position, end_position, missed_cleavages - j))

            # C-terminus is enzyme-specific
            # Let the N-terminus be any position between the least recent proteolytic site
            # (cleavage_sites[0]) and the C-terminal site (i).
            j = cleavage_sites[0]
            if i is None:
                i = sequence_length(sequence)
            for k in range(j, i):
                if k in cleavage_sites:
                    continue
                seq = sequence[k:i]
                # Only deal with this if the sequence is non-empty
                if seq:
                    seq_size = sequence_length(seq)
                    if (min_length is None or seq_size >= min_length) and (
                            max_length is None or seq_size <= max_length):
                        # Get parent sequence coordinates
                        start_position = k
                        end_position = cleavage_sites[-1]
                        if end_position is None:
                            end_position = sequence_length(sequence)
                        n_missed = 0
                        for site in cleavage_sites:
                            if k < site < i:
                                n_missed += 1
                        # Store the sequence, start position, end position, and # of missed cleavages
                        peptides.append((seq, start_position, end_position, n_missed))

        # N-terminus is enzyme-specific
        for ix, j in enumerate(all_sites):
            terminal_ix = min(ix + missed_cleavages + 1, n_sites - 1)
            terminal_site = all_sites[terminal_ix]
            sites_between = all_sites[ix:terminal_ix]
            for k in range(j + 1, terminal_site):
                seq = sequence[j:k]
                if seq:
                    seq_size = sequence_length(seq)
                    if (min_length is None or seq_size >= min_length) and (
                            max_length is None or seq_size <= max_length):
                        # Get parent sequence coordinates
                        start_position = j
                        end_position = k
                        if end_position is None:
                            end_position = sequence_length(sequence)
                        n_missed = 0
                        for site in sites_between:
                            if j < site < k:
                                n_missed += 1
                        # Store the sequence, start position, end position, and # of missed cleavages
                        peptides.append((seq, start_position, end_position, n_missed))

        return sorted(set(peptides), key=_get12)

    @classmethod
    def combine(cls, *names):
        members = [Protease(name) if not isinstance(name, Protease) else name for name in names]
        patterns = merge_enzyme_rules([member.regex.pattern for member in members])
        name = ' + '.join([member.name for member in members])
        instance = cls(patterns)
        instance.name = name
        return instance

    @staticmethod
    def nonspecific_digest(sequence, min_length=0, max_length=100):
        if isinstance(sequence, PeptideSequence):
            sequence = str(sequence)
        n = len(sequence)
        out = []
        for i in range(n):
            for j in range(min_length, max_length + 1):
                if i + j > n:
                    break
                out.append((sequence[i:i + j], i, i + j, 0))
        return out


enzyme_rules = {
    'bnps-skatole': 'W',
    'caspase 1': '(?<=[FWYL]\\w[HAT])D(?=[^PEDQKR])',
    'caspase 10': '(?<=IEA)D',
    'caspase 2': '(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3': '(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4': '(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5': '(?<=[LW]EH)D',
    'caspase 6': '(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7': '(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8': '(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9': '(?<=LEH)D',
    'chymotrypsin': '(?<=[FYWL])(?!P)',
    'chymotrypsin high specificity': '([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity': '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain': 'R',
    'cnbr': 'M',
    'enterokinase': '(?<=[DE]{3})K',
    'factor xa': '(?<=[AFGILTVM][DE]G)R',
    'formic acid': 'D',
    'glutamyl endopeptidase': 'E',
    'glu-c': 'E',
    'granzyme b': '(?<=IEP)D',
    'hydroxylamine': 'N(?=G)',
    'iodosobenzoic acid': 'W',
    'lys-c': "(?<=K)(?!P)",
    'ntcb': '\\w(?=C)',
    'pepsin ph1.3': '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\\w[^P]))',
    'pepsin ph2.0': '((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\\w[^P]))',
    'proline endopeptidase': '(?<=[HKR])P(?=[^P])',
    'proteinase k': '[AEFILTVWY]',
    'staphylococcal peptidase i': '(?<=[^E])E',
    'thermolysin': '[^DE](?=[AFILMV])',
    'thrombin': '((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin': '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
    # No matches ever
    'none': '^&$',
    'trypchymo': '(?<=[FYWLKR])(?!P)',
    'trypsin/p': '(?<=[KR])',
    'lys-c/p': '(?<=K)',
    'pepsina': '(?<=[FL])',
    'v8-de': '(?<=[BDEZ])(?!P)',
    'v8-e': '(?<=[EZ])(?!P)',
    '2-iodobenzoate': '(?<=W)',
    'arg-c': '(?<=R)(?!P)',
    'asp-n_ambic': '(?=[DE])',
    'asp-n': '(?=[BD])',
    'leukocyte elastase': '(?<=[ALIV])(?!P)',
    'alpha-lytic protease': '[TASV]',
}


expasy_rules = enzyme_rules


def register_enzyme(name, pattern):
    if name in enzyme_rules and enzyme_rules[name] != pattern:
        warnings.warn("Overriding {0} with rule {1} with {2}".format(name, enzyme_rules[name], pattern))
    enzyme_rules[name] = pattern


def merge_enzyme_rules(enzyme_patterns):
    rules = ["(" + pattern + ")" for pattern in enzyme_patterns]
    return "|".join(rules)


def get_enzyme(name):
    try:
        data = enzymes[name.lower()]
        return Protease(**data)
    except KeyError:
        warnings.warn("Could not identify protease {}".format(name))
        return Protease(name=name)


trypsin = Protease("trypsin")
gluc = Protease("glutamyl endopeptidase")
