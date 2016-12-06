import warnings
import re
import itertools
import operator

from collections import deque

from .structure.sequence import PeptideSequence
from .structure.parser import sequence_length

_get1 = operator.itemgetter(1)


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
            self.regex = re.compile(expasy_rules[name.lower()])
        except KeyError:
            self.regex = re.compile(name)

    def __repr__(self):
        return "Protease(%s, %s)" % (self.name, self.regex.pattern)

    def missed_cleavages(self, sequence):
        if isinstance(sequence, PeptideSequence):
            sequence = str(sequence)
        return len(self.regex.findall(sequence))

    def cleave(self, sequence, missed_cleavages=0, min_length=0):
        peptides = []
        if isinstance(sequence, PeptideSequence):
            sequence = str(sequence)
        cleavage_sites = deque([0], maxlen=missed_cleavages + 2)
        for i in itertools.chain(map(lambda x: x.end(), self.regex.finditer(sequence)), [None]):
            cleavage_sites.append(i)
            for j in range(len(cleavage_sites) - 1):
                seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
                if seq:
                    if min_length is None or sequence_length(seq) >= min_length:
                        peptides.append((seq, cleavage_sites[j], cleavage_sites[-1] if cleavage_sites[-1]
                                         is not None else sequence_length(sequence)))
        return sorted(set(peptides), key=_get1)

    @classmethod
    def combine(cls, *names):
        patterns = merge_enzyme_rules(names)
        name = ' + '.join(names)
        instance = cls(patterns)
        instance.name = name
        return instance


expasy_rules = {'arg-c': 'R',
                'asp-n': '\\w(?=D)',
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
                'chymotrypsin high specificity': '([FY](?=[^P]))|(W(?=[^MP]))',
                'chymotrypsin low specificity': '([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
                'clostripain': 'R',
                'cnbr': 'M',
                'enterokinase': '(?<=[DE]{3})K',
                'factor xa': '(?<=[AFGILTVM][DE]G)R',
                'formic acid': 'D',
                'glutamyl endopeptidase': 'E',
                'granzyme b': '(?<=IEP)D',
                'hydroxylamine': 'N(?=G)',
                'iodosobenzoic acid': 'W',
                'lysc': 'K',
                'ntcb': '\\w(?=C)',
                'pepsin ph1.3': '((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|((?<=[^HKR][^P])[FLWY](?=\\w[^P]))',
                'pepsin ph2.0': '((?<=[^HKR][^P])[^R](?=[FL][^P]))|((?<=[^HKR][^P])[FL](?=\\w[^P]))',
                'proline endopeptidase': '(?<=[HKR])P(?=[^P])',
                'proteinase k': '[AEFILTVWY]',
                'staphylococcal peptidase i': '(?<=[^E])E',
                'thermolysin': '[^DE](?=[AFILMV])',
                'thrombin': '((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
                'trypsin': '([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))'}


def merge_enzyme_rules(enzyme_names):
    rules = ["(" + expasy_rules[name] + ")" for name in enzyme_names]
    return "|".join(rules)


def get_enzyme(name):
    try:
        data = enzymes[name.lower()]
        return Protease(**data)
    except KeyError:
        warnings.warn("Could not identify protease {}".format(name))
        return Protease(name=name)
