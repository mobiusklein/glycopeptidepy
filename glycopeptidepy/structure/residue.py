from . import ResidueBase
from .composition import Composition
from ..utils.memoize import memoize
from glypy.utils.multimap import MultiMap
from six import add_metaclass


try:
    # check to see if the C base class exists so we
    # do not allocate additional slots
    _has_c = True
    from glycopeptidepy._c.structure import base
    del base
except ImportError:
    _has_c = False

symbol_to_residue = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'C': 'Cys',
    'E': 'Glu',
    'Q': 'Gln',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'L': 'Leu',
    'K': 'Lys',
    'M': 'Met',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'T': 'Thr',
    'W': 'Trp',
    'Y': 'Tyr',
    'V': 'Val',
    "U": "Sec",
    "O": "Pyl",
    "X": "Unk"
}


residue_to_symbol = {value: key for key, value in symbol_to_residue.items()}


symbol_to_long = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic Acid",
    "C": "Cysteine",
    "E": "Glutamic Acid",
    "Q": "Glutamine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine",
    "X": "Unknown",
    "U": "Selenocysteine",
    "O": "Pyrrolysine",
}


long_to_symbol = {value: key for key, value in symbol_to_long.items()}


residue_table = {
    'Ala': 'C3H5NO',
    'Arg': 'C6H12N4O1',
    'Asn': 'C4H6N2O2',
    'Asp': 'C4H5N1O3',
    'Cys': 'C3H5N1O1S1',
    'Glu': 'C5H7NO3',
    'Gln': 'C5H8N2O2',
    'Gly': 'C2H3N1O1',
    'His': 'C6H7N3O1',
    'Ile': 'C6H11N1O1',
    'Leu': 'C6H11N1O1',
    'Lys': 'C6H12N2O1',
    'Met': 'C5H9N1O1S1',
    'Phe': 'C9H9N1O1',
    'Pro': 'C5H7N1O1',
    'Ser': 'C3H5N1O2',
    'Thr': 'C4H7N1O2',
    'Trp': 'C11H10N2O1',
    'Tyr': 'C9H9N1O2',
    'Val': 'C5H9N1O1',
    "Sec": "C3H5NO1Se",
    "Pyl": "C12H19N3O2",
    "Unk": "",
}

residue_chemical_property_group = {
    'Ala': 'hydrophobic',
    'Arg': 'positive_side_chain',
    'Asn': 'polar_uncharged',
    'Asp': 'negative_side_chain',
    'Cys': 'special_case',
    'Glu': 'negative_side_chain',
    'Gln': 'polar_uncharged',
    'Gly': 'special_case',
    'His': 'positive_side_chain',
    'Ile': 'hydrophobic',
    'Leu': 'hydrophobic',
    'Lys': 'positive_side_chain',
    'Met': 'hydrophobic',
    'Phe': 'hydrophobic',
    'Pro': 'special_case',
    'Ser': 'polar_uncharged',
    'Thr': 'polar_uncharged',
    'Trp': 'hydrophobic',
    'Tyr': 'hydrophobic',
    'Val': 'hydrophobic',
}

residue_to_neutral_loss = MultiMap()
residue_to_neutral_loss.update({
    "Ser": [-Composition("H2O")],
    "Thr": [-Composition("H2O")],
    "Glu": [-Composition("H2O")],
    "Asp": [-Composition("H2O")],

    "Arg": [-Composition("NH3")],
    "Lys": [-Composition("NH3")],
    "Gln": [-Composition("NH3")],
    "Asn": [-Composition("NH3")],
})

residue_to_neutral_loss.update({
    "Arg": [-Composition("H2O")],
    "His": [-Composition("H2O")],
    "Lys": [-Composition("H2O")]
})


degeneracy_index = {
}


standard_residues = [
    'Ala',
    'Arg',
    'Asn',
    'Asp',
    'Cys',
    'Glu',
    'Gln',
    'Gly',
    'His',
    'Ile',
    'Leu',
    'Lys',
    'Met',
    'Phe',
    'Pro',
    'Ser',
    'Thr',
    'Trp',
    'Tyr',
    'Val',
]


class UnknownAminoAcidException(KeyError):
    pass


class MemoizedResidueMetaclass(type):
    '''
    A metaclass that memoizes Residues as they are constructed
    by overriding the class __call__ method. It will attempt to
    look up previously created residues by symbol, then by name.
    If a previous instance is not found, it will be created and
    saved.

    Attributes
    ----------
    _cache: dict

    '''

    def __call__(self, symbol=None, name=None, *args, **kwargs):
        try:
            cache = self._cache
        except AttributeError:
            cache = self._cache = dict()
        try:
            if symbol is not None:
                return cache[symbol]
            elif name is not None:
                return cache[name]
            else:
                raise Exception("Must provide a symbol or name parameter")
        except KeyError:
            if symbol is not None:
                inst = type.__call__(self, symbol=symbol, *args, **kwargs)
                cache[inst.symbol] = inst
                cache[inst.name] = inst
                return inst

            elif name is not None:
                inst = type.__call__(self, name=name, *args, **kwargs)
                cache[inst.symbol] = inst
                cache[inst.name] = inst
                return inst
            else:
                raise UnknownAminoAcidException(
                    "Cannot find a AminoAcidResidue for %r" % ((symbol, name),))

    def parse(self, symbol):
        try:
            cache = self._cache
            return cache[symbol]
        except (AttributeError, KeyError):
            return self(symbol)


def _AminoAcidResidue_reconstructor():
    return AminoAcidResidue.__new__(AminoAcidResidue)


@add_metaclass(MemoizedResidueMetaclass)
class AminoAcidResidue(ResidueBase):
    '''
    Represent a single Amino Acid residue which compose peptide sequences. The
    structure itself is intended to be immutable.

    Attributes
    ----------
    name: str
        The three letter abbreviation for the amino acid
    symbol: str
        The single letter abbreviation for the amino acid
    mass: float
        The neutral mass of the amino acid
    composition: :class:`glypy.Composition`
        The chemical composition of the amino acid
    chemical_class: str
        One of hydrophobic, polar_uncharged, positive_side_chain, negative_side_chain, or special_case,
        corresponding to the gross molecular properties of the amino acid
    is_degenerate: bool:
        Whether this amino acid may stand for other amino acids

    '''
    # If the C extension base class is available, do not declare redundant slot
    # descriptors.
    if _has_c:
        __slots__ = ["neutral_loss", ]
    else:
        __slots__ = ["name", "symbol", "mass", "composition", "_hash", "neutral_loss"]

    @staticmethod
    @memoize()
    def mass_by_name(sym):
        name = symbol_to_residue.get(sym, sym)
        formula = residue_table.get(name)
        return Composition(formula).mass

    def __init__(self, symbol=None, name=None):
        self.symbol = symbol
        self.name = name
        self.mass = 0.0
        try:
            if symbol is not None:
                self._by_symbol(symbol)
            elif name is not None:
                self._by_name(name)
        except KeyError:
            raise UnknownAminoAcidException(
                "No definition for Amino Acid %s" % (symbol if symbol is not None else name))
        self._hash = hash(self.name)
        self.neutral_loss = residue_to_neutral_loss[self.name]

    def _by_name(self, name):
        """Configure this instance by name information

        Parameters
        ----------
        name : str
            Amino Acid Name
        """
        self.composition = Composition(residue_table[name])
        self.name = name
        self.mass = self.composition.mass
        self.symbol = residue_to_symbol[name]

    def _by_symbol(self, symbol):
        """Configure this instance by symbol information,
        by going from symbol to name, and from name to data

        Parameters
        ----------
        symbol : str
            Amino Acid symbol
        """
        try:
            name = symbol_to_residue[symbol]
            self._by_name(name)
        except KeyError:
            self._by_name(symbol)

    def __repr__(self):
        return self.name

    if not _has_c:
        def __hash__(self):
            return self._hash

        def __eq__(self, other):
            if self is other:
                return True
            try:
                return self.name == other.name or self.symbol == other.symbol
            except AttributeError:
                return self.name == other or self.symbol == other

        def __ne__(self, other):
            if self is other:
                return False
            try:
                return self.name != other.name and self.symbol != other.symbol
            except AttributeError:
                return self.name != other and self.symbol != other

    def __reduce__(self):
        return _AminoAcidResidue_reconstructor, (), self.__getstate__()

    def __getstate__(self):
        return [self.name, self.symbol, self.mass, self.composition, self.neutral_loss]

    def __setstate__(self, state):
        self.name = state[0]
        self.symbol = state[1]
        self.mass = state[2]
        self.composition = state[3]
        if len(state) > 4:
            self.neutral_loss = state[4]
        self._hash = hash(self.name)

    @property
    def chemical_class(self):
        return residue_chemical_property_group[self.name]

    @property
    def is_degenerate(self):
        try:
            return degeneracy_index[self.name]
        except KeyError:
            return False

    @property
    def long_name(self):
        try:
            return symbol_to_long[self.symbol]
        except KeyError:
            return self.name


Residue = AminoAcidResidue


def register_residue(name, symbol, formula, chemical_class=None):
    assert symbol not in symbol_to_residue
    assert name not in residue_table
    residue_to_symbol[name] = symbol
    symbol_to_residue[symbol] = name
    residue_table[name] = formula
    residue_chemical_property_group[name] = chemical_class
    return AminoAcidResidue(symbol=symbol)


def register_degenerate(name, symbol, mappings):
    assert symbol not in symbol_to_residue
    assert name not in residue_table
    residue_to_symbol[name] = symbol
    symbol_to_residue[symbol] = name
    residue_table[name] = residue_table[mappings[0]]
    residue_chemical_property_group[
        name] = residue_chemical_property_group[mappings[0]]
    residue_to_neutral_loss[name] = residue_to_neutral_loss[mappings[0]]
    degeneracy_index[name] = frozenset(mappings)
    return AminoAcidResidue(symbol=symbol)


register_degenerate("Xle", "J", ["Leu", "Ile"])


def get_all_residues(omit_unknown=True, standard_only=True):
    '''Get a collection of all :class:`AminoAcidResidue` types defined.

    Parameters
    ----------
    omit_unknown: bool
        Whether to discard ``X`` or not.
    standard_only: bool
        Whether to include only the 20 standard amino acids or all
        defined amino acids.
    Returns
    -------
    set
    '''
    if standard_only:
        return set(map(AminoAcidResidue, standard_residues))
    residues = set(map(AminoAcidResidue, symbol_to_residue))
    if omit_unknown:
        residues.remove(AminoAcidResidue("X"))
    return residues


def get_all_sequencing_residues(omit_unknown=True, standard_only=True):
    residues = set(get_all_residues(omit_unknown, standard_only))
    residues.add(AminoAcidResidue("J"))

    for residue in list(residues):
        degenerate = residue.is_degenerate
        if degenerate:
            for degen in degenerate:
                residues.remove(degen)
    return residues


AminoAcidResidue.get_all_residues = staticmethod(get_all_residues)
AminoAcidResidue.get_all_sequencing_residues = staticmethod(get_all_sequencing_residues)
