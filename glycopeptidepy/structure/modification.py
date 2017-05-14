import csv
import json
from copy import deepcopy
import re
from pkg_resources import resource_stream

from six import string_types as basestring
from io import StringIO

from collections import defaultdict
from collections import Iterable
from glypy.composition.glycan_composition import (
    FrozenMonosaccharideResidue, FrozenGlycanComposition)
from glypy.utils import Enum

from .residue import Residue
from .composition import Composition
from . import ModificationBase


target_string_pattern = re.compile(
    r'(?P<amino_acid>[A-Z]+)?([@ ]?(?P<n_term>[Nn][-_][tT]erm)|(?P<c_term>[Cc][-_][tT]erm))?')
title_cleaner = re.compile(
    r'(?P<name>.*)\s\((?P<target>.+)\)$')

is_mass_delta = re.compile(r"^Delta:")


class ModificationIndex(dict):
    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            return 0


class SequenceLocation(Enum):
    anywhere = None
    n_term = -1
    c_term = -2
    protein_n_term = -3
    protein_c_term = -4


SequenceLocation.n_term.add_name("N-term")
SequenceLocation.n_term.add_name("N-Term")
SequenceLocation.n_term.add_name("n-term")
SequenceLocation.n_term.add_name("N_term")
SequenceLocation.c_term.add_name("C-term")
SequenceLocation.c_term.add_name("C-Term")
SequenceLocation.c_term.add_name("c-term")
SequenceLocation.c_term.add_name("C_term")
SequenceLocation.anywhere.add_name("Anywhere")
SequenceLocation.protein_n_term.add_name("Protein N-term")
SequenceLocation.protein_c_term.add_name("Protein C-term")
SequenceLocation.protein_n_term.add_name("Protein N-Term")
SequenceLocation.protein_c_term.add_name("Protein C-Term")


class ModificationCategory(Enum):
    unknown = None
    glycosylation = 1
    artefact = 2
    substitution = 3
    chemical_derivative = 4
    non_standard_residue = 5
    isotopic_label = 6
    post_translational = 7
    other = 8
    other_glycosylation = 9
    multiple = 10
    pre_translational = 11
    co_translational = 12
    synthetic_peptide_protect = 13


ModificationCategory.substitution.add_name("AA substitution")
ModificationCategory.other_glycosylation.add_name("Other glycosylation")
ModificationCategory.glycosylation.add_name("N-linked glycosylation")
ModificationCategory.glycosylation.add_name("O-linked glycosylation")
ModificationCategory.chemical_derivative.add_name("Chemical derivative")
ModificationCategory.post_translational.add_name("Post-translational")
ModificationCategory.multiple.add_name("Multiple")
ModificationCategory.artefact.add_name("Artefact")
ModificationCategory.isotopic_label.add_name("Isotopic label")
ModificationCategory.other.add_name("Other")
ModificationCategory.pre_translational.add_name("Pre-translational")
ModificationCategory.non_standard_residue.add_name("Non-standard residue")
ModificationCategory.co_translational.add_name("Co-translational")
ModificationCategory.synthetic_peptide_protect.add_name("Synth. pep. protect. gp.")


def composition_delta_parser(formula):
    '''Parses a Unimod composition offset formula where
    the isotope specification is the 'opposite' of the
    specification used by this library.

    Warning - Currently unused internally.

    Parameters
    ----------
    formula: str
        The formula parsed

    Return
    ------
    dict
    '''
    elemental_composition = {}
    components = re.findall(r"([A-Za-z0-9]+)\(([\-0-9]+)\)", formula)
    for element, amount in components:
        parts = re.search(r"(?P<isotope>\d*)(?P<element>.*)", element).groupdict()
        if parts['isotope'] != '':
            element = "{element}[{isotope}]".format(**parts)
        elemental_composition[element] = amount
    return elemental_composition


def extract_targets_from_string(target_string):
    '''Parses the Protein Prospector modification name string to
    extract the modification target specificity

    Format:

    `"Modification Name (Residues Targeted)"`

    Parameters
    ----------
    target_string: str

    Return
    ------
    :class:`ModificationTarget`

    '''
    position_only = re.search(r"^([NnCc][-_][tT]erm)", target_string)
    if position_only:
        return ModificationTarget(None, SequenceLocation[position_only.groups()[0]])
    else:
        amino_acid_and_position = re.search(
            r'(?P<amino_acid>[A-Z]+)([@ ]?(?P<n_term>[Nn][-_][tT]erm)|(?P<c_term>[Cc][-_][tT]erm))?', target_string)
        params = amino_acid_and_position.groupdict()
        aa = list(map(Residue, params["amino_acid"]))
        if params['n_term']:
            return ModificationTarget(aa, SequenceLocation.n_term)
        elif params['c_term']:
            return ModificationTarget(aa, SequenceLocation.c_term)
        else:
            return ModificationTarget(aa, SequenceLocation.anywhere)


def get_position_modifier_rules_dict(sequence):
    '''Labels the start and end indices of the sequence

    A convenience dictionary initializer for ease-of-access in
    :class:`ModificationTarget` site validation functions

    Parameters
    ----------
    sequence: PeptideSequence
        The sequence for which to project labels

    Return
    ------
    defaultdict
    '''
    return defaultdict(lambda: SequenceLocation.anywhere, {
        0: SequenceLocation.n_term,
        (len(sequence) - 1): SequenceLocation.c_term
    })


class ModificationTarget(object):
    '''Specifies the range of targets that a given modification can be found at.

    A single instance of :class:`ModificationTarget` describes one or more position/residue
    rules for a single type of modification, bundled together underneath a single :class:`ModificationRule`
    instance. :class:`ModificationTarget` may be extracted from Unimod data dictionaries using the
    :meth:`from_unimod_specificity` class method.

    Attributes
    ----------
    amino_acid_targets: list
        The list of amino acid residues that this rule
        permits.
    position_modifier: SequenceLocation
        The sequence-position specific filter
        this rule enforces.
    '''

    @classmethod
    def from_unimod_specificity(cls, specificity):
        amino_acid = (specificity["site"])
        if amino_acid in (["N-term", "C-term", None]):
            amino_acid = None
        else:
            amino_acid = {Residue(amino_acid)}
        position_modifier = specificity["position"]
        nterm = ("Protein N-term", "Any N-term", "N-term")
        cterm = ("Protein C-term", "Any C-term", "C-term")
        site = specificity["site"]
        classification = ModificationCategory[specificity.get("classification")]
        if position_modifier == "Anywhere":
            position_modifier = SequenceLocation.anywhere
        if position_modifier in nterm or site in nterm:
            position_modifier = SequenceLocation.n_term
        elif position_modifier in cterm or site in cterm:
            position_modifier = SequenceLocation.c_term
        # elif position_modifier :
        #     raise Exception("Undefined Position, " + str(position_modifier))
        return cls(amino_acid, position_modifier, [classification])

    def __init__(self, site_targets=None, position_modifier=None, classification=None):
        if classification is None:
            classification = []
        self.amino_acid_targets = frozenset(map(Residue, site_targets)) if site_targets is not None else None
        self.position_modifier = position_modifier
        self.classification = classification

    def valid_site(self, amino_acid=None, position_modifier=SequenceLocation.anywhere):
        '''Return if a given residue at a sequence position (N-term, C-term, None)
        is valid'''
        valid = False

        # Validate amino acid target target
        valid = (self.amino_acid_targets is None) or (
            amino_acid in self.amino_acid_targets)
        valid = valid and ((self.position_modifier is SequenceLocation.anywhere) or
                           (position_modifier == self.position_modifier))

        if valid and (position_modifier is not SequenceLocation.anywhere) and (
                self.position_modifier is not SequenceLocation.anywhere):
            return position_modifier

        return valid

    def valid_site_seq(self, sequence, position, position_modifier=SequenceLocation.anywhere):
        aa_mod_pair = sequence[position]
        if len(aa_mod_pair[1]) > 0:
            return False
        amino_acid = aa_mod_pair[0]

        return self.valid_site(amino_acid, position_modifier)

    def __repr__(self):
        amino_acid_components = "{%s}" % ', '.join(map(str, self.amino_acid_targets or ()))
        rep = "{amino_acid_targets}@{position_modifier}".format(
            amino_acid_targets=amino_acid_components,
            position_modifier=self.position_modifier)
        return rep

    def __eq__(self, other):
        return repr(self) == repr(other)

    def __hash__(self):
        return hash(repr(self))

    def __len__(self):
        position_modifiers_count = 0
        amino_acid_targets_count = 0
        if self.position_modifiers is not None:
            position_modifiers_count = len(self.position_modifiers)
        if self.amino_acid_targets is not None:
            amino_acid_targets_count = len(self.amino_acid_targets)
        return amino_acid_targets_count + position_modifiers_count

    def serialize(self):
        parts = []
        if self.amino_acid_targets is not None:
            parts.append(''.join(aa.symbol for aa in self.amino_acid_targets))
            if self.position_modifier == SequenceLocation.anywhere:
                pass
            else:
                parts.append("@")
        if self.position_modifier in (SequenceLocation.protein_n_term, SequenceLocation.n_term):
            parts.append("N-term")
        if self.position_modifier in (SequenceLocation.protein_c_term, SequenceLocation.c_term):
            parts.append("C-term")
        return " ".join(parts)


class ModificationNeutralLoss(object):
    def __init__(self, composition):
        self.composition = Composition(composition)
        self.mass = self.composition.mass

    def __repr__(self):
        return "ModificationNeutralLoss(%s)" % self.composition


class ModificationRule(object):
    '''Represent the name, mass, and position specifities associated with a given modification.
    Additionally stores information on metadata about this rule computed from these values and
    categorical information from expert sources.

    Attributes
    ----------
    categories : list
        The list of named categories of PTMs that this modification
        falls under. Examples are "artefact", "glycosylation", or "biological process"
    mass : float
        The mass of the modification defined by this rule
    name : str
        A common name for the modification defined by this rule
    names : set
        The set of all names used to refer to this modification, beyond those
        commonly used and those for display purposes.
    aliases : set
        The set of all names that are allowed to be used to display this modification
    options : dict
        Arbitrary information from external sources about this modification
        rule and its phenomena.
    preferred_name : str
        The name for this modification preferred for using in user-facing
        displays
    title : str
        A formal or "full" name for this modification that may be its chemical name
        or some other more verbose representation which is used to identify it but
        not convenient to read or write.
    targets : set
        All applicable :class:`.ModificationTarget` instances which apply to this
        rule. This set may differ from instance to instance of the same rule by
        configuration as per the :class:`RestrictedModificationTable`, for example
    '''

    @staticmethod
    def reduce(rules):
        collector = defaultdict(list)
        for rule in rules:
            collector[rule.name].append(rule)
        reductions = []
        for rule_name, rules in collector.items():
            accum = rules[0]
            for other in rules[1:]:
                accum = accum + other
            reductions.append(accum)
        return reductions

    @classmethod
    def from_protein_prospector(cls, amino_acid_specificity, modification_name,
                                title=None, monoisotopic_mass=None,
                                neutral_losses=None):
        if neutral_losses is None:
            neutral_losses = []
        # If the specificity is a string, parse it into rules
        targets = None
        if isinstance(amino_acid_specificity, basestring):
            try:
                targets = set(
                    [extract_targets_from_string(amino_acid_specificity)])
            except TypeError:
                print(modification_name, amino_acid_specificity)
                raise TypeError("Could not extract targets")

        # If the specificity is already a rule, store it as a list of one rules
        elif isinstance(amino_acid_specificity, ModificationTarget):
            targets = set([amino_acid_specificity])

        # If the specificity is a collection of rules, store it as is
        elif isinstance(amino_acid_specificity, Iterable) and\
                isinstance(iter(amino_acid_specificity).next(), ModificationTarget):
            targets = set(amino_acid_specificity)

        return cls(targets, modification_name, title, monoisotopic_mass,
                   neutral_losses=neutral_losses)

    @classmethod
    def from_unimod(cls, unimod_entry):
        specificity = unimod_entry["specificity"]
        monoisotopic_mass = unimod_entry["mono_mass"]
        full_name = unimod_entry["full_name"]
        title = unimod_entry["title"]
        alt_names = set(unimod_entry.get('alt_names', ()))
        specificity = list(map(ModificationTarget.from_unimod_specificity, specificity))
        return cls(
            specificity, full_name, title, monoisotopic_mass,
            composition=Composition({str(k): int(v) for k, v in unimod_entry['composition'].items()}),
            alt_names=alt_names)

    def __init__(self, amino_acid_specificity, modification_name,
                 title=None, monoisotopic_mass=None, composition=None,
                 categories=None, alt_names=None, neutral_losses=None,
                 aliases=None,
                 **kwargs):
        if neutral_losses is None:
            neutral_losses = []
        if categories is None:
            categories = []
        if aliases is None:
            aliases = set()
        # Attempt to parse the protein prospector name which contains the
        # target in
        try:
            self.common_name = title_cleaner.search(title).groupdict()['name']
        except Exception:
            self.common_name = modification_name

        if alt_names is None:
            alt_names = set()

        self.unimod_name = modification_name
        self.mass = float(monoisotopic_mass)
        self.title = title if title is not None else modification_name
        self.composition = composition
        self.names = {self.unimod_name, self.common_name, self.title} | alt_names
        for name in list(self.names):
            self.names.update(name.split(" or "))
        self.categories = set(categories)
        self.options = kwargs
        self._n_term_target = None
        self._c_term_target = None
        self.fragile = kwargs.get('fragile', False)
        self.preferred_name = min(self.names, key=len)
        self.aliases = aliases
        # The type of the parameter passed for amino_acid_specificity is variable
        # so select the method correct for the passed type

        # If the specificity is a string, parse it into rules
        if isinstance(amino_acid_specificity, basestring):
            try:
                self.targets = set([
                    extract_targets_from_string(amino_acid_specificity)])
            except TypeError:
                raise TypeError("Could not extract targets from %s" % amino_acid_specificity)

        # If the specificity is already a rule, store it as a list of one rules
        elif isinstance(amino_acid_specificity, ModificationTarget):
            self.targets = set([amino_acid_specificity])

        # If the specificity is a collection of rules, store it as is
        elif isinstance(amino_acid_specificity, Iterable):
            self.targets = set(amino_acid_specificity)
        else:
            raise ValueError("Could not interpret target specificity from %r" % (
                amino_acid_specificity,))
        self.neutral_losses = neutral_losses
        for target in self.targets:
            self.categories.update(
                target.classification)

    @property
    def is_standard(self):
        return True

    @property
    def name(self):
        return self.preferred_name

    def clone(self, propagated_targets=None):
        if propagated_targets is None:
            propagated_targets = (self.targets)

        dup = self.__class__(
            set(propagated_targets), self.name,
            self.title, self.mass, self.composition,
            self.categories, **self.options)
        dup.names.update(self.names)
        return dup

    def valid_site(self, amino_acid=None, position_modifiers=None):
        return max([target.valid_site(amino_acid, position_modifiers) for
                    target in self.targets])

    def why_valid(self, amino_acid=None, position_modifiers=None):
        possible_targets = [target for target in self.targets if target.valid_site(
            amino_acid, position_modifiers)]
        if len(possible_targets) == 0:
            raise ValueError("No valid targets. %s for %s" %
                             (str(amino_acid) + str(position_modifiers), self))
        minimized_target = min(possible_targets, key=len)
        return minimized_target

    def valid_site_seq(self, sequence, position, position_modifiers):
        return any([target.valid_site_seq(sequence, position, position_modifiers) for
                    target in self.targets])

    def find_valid_sites(self, sequence):
        valid_indices = []
        position_modifier_rules = get_position_modifier_rules_dict(sequence)
        for index in range(len(sequence)):
            position = sequence[index]
            if len(position[1]) > 0:
                continue
            amino_acid = position[0].name
            position_modifiers = position_modifier_rules.get(index)
            is_valid = self.valid_site(amino_acid, position_modifiers)
            if(is_valid):
                if (is_valid is SequenceLocation.n_term) or (is_valid is SequenceLocation.c_term):
                    valid_indices.append(is_valid)
                else:
                    valid_indices.append(index)

        return valid_indices

    def serialize(self):
        '''A string representation for inclusion in sequences'''
        return self.preferred_name

    def __repr__(self):
        rep = "{preferred_name}:{mass}".format(
            **self.__dict__)
        return rep

    def __add__(self, other):
        if(isinstance(other, ModificationRule)):
            if (other.targets) is (self.targets):
                return self
            dup = self.clone()
            dup.targets = (set(self.targets) | set(other.targets))
            if self.composition is None:
                new_composition = other.composition
            else:
                new_composition = self.composition
            dup.composition = new_composition
            dup.names = self.names | other.names
            dup.aliases = self.aliases | other.aliases
            dup.preferred_name = min(dup.names, key=len)
            dup.options.update(other.options)
            return dup
        else:
            raise TypeError(
                "Unsupported types. Can only add two ModificationRules together.")

    def __sub__(self, other):
        if isinstance(other, ModificationRule):
            dup = self.clone()
            dup.targets = set(self.targets) - set(other.targets)
            return dup
        else:
            raise TypeError(
                "Unsupported types. Can only subtract two ModificationRules together.")

    def __call__(self, **kwargs):
        return Modification(self, **kwargs)

    def is_a(self, category):
        return category in self.categories

    @property
    def n_term_targets(self):
        if self._n_term_target is None:
            solutions = []
            for target in self.targets:
                if target.position_modifier == SequenceLocation.n_term:
                    solutions.append(target)
            self._n_term_target = solutions
        return self._n_term_target

    @property
    def c_term_targets(self):
        if self._c_term_target is None:
            solutions = []
            for target in self.targets:
                if target.position_modifier == SequenceLocation.c_term:
                    solutions.append(target)
            self._c_term_target = solutions
        return self._c_term_target

    def __hash__(self):
        return hash(self.preferred_name)

    def __eq__(self, other):
        if self is other:
            return True
        idents = self.names
        try:
            other_idents = other.names
        except AttributeError:
            other_idents = {other}
        return len(idents & other_idents) > 0

    def __ne__(self, other):
        return not self == other

    def losses(self):
        yield self.preferred_name, self.composition.clone()

    def as_spec_strings(self):
        for target in self.targets:
            yield "%s (%s)" % (self.title, target.serialize())


class AnonymousModificationRule(ModificationRule):
    '''
    Construct a modification rule from a string that can be placed
    in a peptide sequence and be reconstructed without being stored
    in a :class:`ModificationTable` instance and backed by a modification
    database.

    "...(@{name}-{mass})..."

    '''
    parser = re.compile(r"@(?P<name>.+?)-(?P<mass>-?[0-9\.]+)")

    @classmethod
    def try_parse(cls, rule_string):
        anon = cls.parser.search(rule_string)
        if anon is not None:
            anon_match_data = anon.groupdict()
            name = anon_match_data['name']
            mass = float(anon_match_data['mass'])
            mod = AnonymousModificationRule(name, mass)
            return mod
        else:
            return None

    def __init__(self, name, mass):
        try:
            super(AnonymousModificationRule, self).__init__(
                "", name, name, mass)
        except AttributeError:
            pass

    @property
    def is_standard(self):
        return False

    def valid_site(self, *args, **kwargs):
        raise TypeError(
            "AnonymousModificationRule does not support site validation")

    def find_valid_sites(self, *args, **kwargs):
        raise TypeError(
            "AnonymousModificationRule does not support site validation")

    def serialize(self):
        return "@" + self.name + "-" + str(self.mass)

    def clone(self):
        dup = self.__class__(self.name, self.mass, **self.options)
        return dup

    def __repr__(self):
        rep = "{name}:{mass}".format(name=self.name, mass=self.mass)
        return rep


def _simple_sequence_mass(seq_list):
    '''Because AminoAcidSubstitution may contain sequence-like entries that
    need massing but we cannot import `sequence` here, this defines a simple
    sequence list to total mass without any terminal groups, but handles residues
    and modifications.

    For a full sequence_to_mass converter, see the sequence module'''
    mass = 0
    for res, mods in seq_list:
        mass += Residue.mass_by_name(res)
        for mod in mods:
            if mod != "":
                mass += Modification.mass_by_name(mod)
    return mass


class AminoAcidSubstitution(AnonymousModificationRule):
    parser = re.compile(r"@(?P<original_residue>\S+)->(?P<substitution_residue>\S+)")

    @classmethod
    def try_parse(cls, rule_string):
        aa_sub = cls.parser.search(rule_string)
        if aa_sub is not None:
            match_data = aa_sub.groupdict()
            mod = cls(match_data['original_residue'], match_data['substitution_residue'])
            return mod
        else:
            return None

    def __init__(self, original_residue, substitution_residue, **kwargs):
        from .parser import prefix_to_postfix_modifications
        self.common_name = "{0}->{1}".format(original_residue, substitution_residue)
        self.mass = 0.0
        self.original = prefix_to_postfix_modifications(original_residue)
        self.substitution = prefix_to_postfix_modifications(substitution_residue)
        self.delta = _simple_sequence_mass(self.substitution) - _simple_sequence_mass(self.original)
        self.preferred_name = self.common_name
        self.names = {self.common_name}
        self.options = kwargs
        self.categories = [ModificationCategory['substitution']]
        self.aliases = set()

    def serialize(self):
        return "@" + self.name

    def __repr__(self):
        return "{name}:{delta}".format(**self.__dict__)


_hexnac = FrozenMonosaccharideResidue.from_iupac_lite("HexNAc")
_hexose = FrozenMonosaccharideResidue.from_iupac_lite("Hex")
_xylose = FrozenMonosaccharideResidue.from_iupac_lite("Xyl")


hexnac_modification = ModificationRule.from_unimod({
    "title": "HexNAc",
    "composition": {
        "H": 13,
        "C": 8,
        "O": 5,
        "N": 1
    },
    "mono_mass": 203.079373,
    "full_name": "N-Acetylhexosamine",
    "specificity": [{
        "position": "Anywhere",
        "hidden": True,
        "site": "T",
        "classification": "Other glycosylation",
        "spec_group": 3
    }, {
        "position": "Anywhere",
        "hidden": True,
        "site": "S",
        "classification": "Other glycosylation",
        "spec_group": 2
    }, {
        "position": "Anywhere",
        "hidden": True,
        "site": "N",
        "classification": "N-linked glycosylation",
        "spec_group": 1
    }],
})


xylose_modification = ModificationRule.from_unimod({
    "title": "Xyl",
    "composition": dict(_xylose.total_composition()),
    "mono_mass": _xylose.mass(),
    "full_name": "Xylose",
    "specificity": [{
        "position": "Anywhere",
        "hidden": True,
        "site": "T",
        "classification": "Other glycosylation",
        "spec_group": 3
    }, {
        "position": "Anywhere",
        "hidden": True,
        "site": "S",
        "classification": "Other glycosylation",
        "spec_group": 2
    }],
})


class Glycosylation(ModificationRule):
    """
    Incubator Idea - Represent occupied glycosylation sites
    using the Modification interface.

    Attributes
    ----------
    categories : list
        Description
    mass : TYPE
        Description
    name : TYPE
        Description
    names : TYPE
        Description
    options : dict
        Description
    parser : TYPE
        Description
    preferred_name : TYPE
        Description
    title : TYPE
        Description
    """
    parser = re.compile(r"@(?P<glycan_composition>\{[^\}]+\})")

    @classmethod
    def try_parse(cls, rule_string):
        try:
            glycosylation = cls.parser.search(rule_string)
            glycosylation = FrozenGlycanComposition.parse(glycosylation.group(0))
            return cls(glycosylation)
        except Exception:
            return None

    def __init__(self, glycan_composition):
        if isinstance(glycan_composition, basestring):
            glycan_composition = FrozenGlycanComposition.parse(glycan_composition)

        glycan_composition.composition_offset = Composition()

        self.common_name = "@" + str(glycan_composition)
        self.mass = glycan_composition.mass()
        self.title = self.common_name
        self.preferred_name = self.common_name
        self.categories = [ModificationCategory.glycosylation]
        self.names = {self.common_name, self.title, self.preferred_name}
        self.aliases = set()
        self.options = {}
        self.fragile = True

    def losses(self):
        for i in []:
            yield i

    def clone(self):
        return self.__class__(FrozenGlycanComposition.parse(self.name[1:]))


class NGlycanCoreGlycosylation(Glycosylation):
    mass_ladder = {k: FrozenGlycanComposition.parse(k).total_composition() - Composition("H2O") for k in {
        "{HexNAc:1}",
        "{HexNAc:2}",
        "{HexNAc:2; Hex:1}",
        "{HexNAc:2; Hex:2}",
        "{HexNAc:2; Hex:3}",
    }}

    def __init__(self, base_mass=_hexnac.mass()):
        self.common_name = "NGlycanCoreGlycosylation"
        self.mass = base_mass
        self.title = "N-Glycan Core Glycosylation"
        self.unimod_name = "N-Glycosylation"
        self.preferred_name = self.unimod_name
        self.targets = [(ModificationTarget("N"))]
        self.composition = _hexnac.total_composition().clone()
        self.names = {self.unimod_name, self.title, self.preferred_name, self.common_name}
        self.options = {}
        self.aliases = set()
        self.categories = [ModificationCategory.glycosylation]

    def losses(self):
        for label_loss in self.mass_ladder.items():
            yield label_loss

    def clone(self):
        return self.__class__(self.mass)


class OGlycanCoreGlycosylation(Glycosylation):
    mass_ladder = {k: FrozenGlycanComposition.parse(k).total_composition() - Composition("H2O") for k in {
        "{HexNAc:1}",
        "{HexNAc:1; Hex:1}",
    }}

    def __init__(self, base_mass=_hexnac.mass()):
        self.common_name = "OGlycanCoreGlycosylation"
        self.mass = base_mass
        self.title = "O-Glycan Core Glycosylation"
        self.unimod_name = "O-Glycosylation"
        self.preferred_name = self.unimod_name
        self.targets = [ModificationTarget("S"), ModificationTarget("T")]
        self.composition = _hexnac.total_composition().clone()
        self.names = {self.unimod_name, self.title, self.preferred_name, self.common_name}
        self.options = {}
        self.aliases = set()
        self.categories = [ModificationCategory.glycosylation]

    def losses(self):
        for label_loss in self.mass_ladder.items():
            yield label_loss

    def clone(self):
        return self.__class__(self.mass)


class GlycosaminoglycanLinkerGlycosylation(Glycosylation):
    mass_ladder = {k: FrozenGlycanComposition.parse(k).total_composition() - Composition("H2O") for k in {
        "{Xyl:1}",
        "{Xyl:1; Hex:1}",
        "{Xyl:1; Hex:2}",
        "{Xyl:1; Hex:2; HexA:1}",
    }}

    def __init__(self, base_mass=_xylose.mass()):
        self.common_name = "GlycosaminoglycanLinkerGlycosylation"
        self.mass = base_mass
        self.title = "Glycosaminoglycan Linker Glycosylation"
        self.unimod_name = "GAG-Linker"
        self.preferred_name = self.unimod_name
        self.targets = [ModificationTarget("S")]
        self.composition = _hexnac.total_composition().clone()
        self.names = {self.unimod_name, self.title, self.preferred_name, self.common_name}
        self.options = {}
        self.aliases = set()
        self.categories = [ModificationCategory.glycosylation]

    def losses(self):
        for label_loss in self.mass_ladder.items():
            yield label_loss

    def clone(self):
        return self.__class__(self.mass)


class OGlcNAcylation(Glycosylation):
    mass_ladder = {
        "{GlcNAc:1}": _hexnac.total_composition()
    }

    def __init__(self, base_mass=_hexnac.mass()):
        self.common_name = "O-GlcNAc"
        self.preferred_name = self.common_name
        self.mass = base_mass
        self.title = "O-GlcNAc"
        self.unimod_name = "O-GlcNAc"
        self.targets = [ModificationTarget("S"), ModificationTarget("T")]
        self.composition = _hexnac.total_composition().clone()
        self.title = self.common_name
        self.names = {self.unimod_name, self.title, self.preferred_name, self.common_name}
        self.categories = [ModificationCategory.glycosylation]
        self.options = {}
        self.aliases = set()

    def losses(self):
        for label_loss in self.mass_ladder.items():
            yield label_loss

    def clone(self):
        return self.__class__(self.mass)


def load_from_csv(stream):
    modification_definitions = []
    parser = csv.DictReader(stream)
    modification_definitions = [rec for rec in parser]
    return modification_definitions


def load_from_json(stream):
    modification_definitions = []
    modification_definitions = list(map(ModificationRule.from_unimod, json.load(stream)))
    return modification_definitions


class ModificationSource(object):
    _table_definition_file = staticmethod(lambda: StringIO(
        resource_stream(__name__, "data/ProteinProspectorModifications-for_gly2.csv").read(
        ).decode("utf-8")))
    _unimod_definitions = staticmethod(lambda: StringIO(
        resource_stream(__name__, "data/unimod.json").read().decode('utf-8')))

    use_protein_prospector = True

    @classmethod
    def _definitions_from_stream(cls, stream, format="csv"):
        if format == "csv":
            definitions = load_from_csv(stream)
        elif format == "json":
            definitions = load_from_json(stream)
        return definitions

    @classmethod
    def _definitions_from_stream_default(cls):
        defs = []
        defs += load_from_json(cls._unimod_definitions())
        if cls.use_protein_prospector:
            defs += load_from_csv(cls._table_definition_file())
        return defs

    @classmethod
    def load_from_file(cls, stream=None, format="csv"):
        '''Load the rules definitions from a CSV or JSON file and instantiate a ModificationSource from it'''
        if stream is None:
            return cls.load_from_file_default()
        if format == "csv":
            definitions = load_from_csv(stream)
        elif format == "json":
            definitions = load_from_json(stream)
        return cls(definitions)

    @classmethod
    def load_from_file_default(cls):
        '''Load the rules definitions from the default package files'''
        defs = []
        defs += load_from_json(cls._unimod_definitions())
        if cls.use_protein_prospector:
            defs += load_from_csv(cls._table_definition_file())
        return cls(defs)

    bootstrapped = None

    @classmethod
    def bootstrap(cls, reuse=True):
        '''Load the bootstrapping file's definitions and instantiate a ModificationSource from it'''
        if cls.bootstrapped is not None and reuse:
            return cls.bootstrapped
        instance = cls.load_from_file()
        cls.bootstrapped = instance
        return instance


class ModificationTable(ModificationSource):

    other_modifications = {
        "HexNAc": hexnac_modification,
        "Xyl": xylose_modification,
        "N-Glycosylation": NGlycanCoreGlycosylation(),
        "O-Glycosylation": OGlycanCoreGlycosylation(),
        "GAG-Linker": GlycosaminoglycanLinkerGlycosylation(),
        "O-GlcNAc": OGlcNAcylation(),
        "H": ModificationRule(
            "N-term", "H", "H",
            Composition("H").mass, composition=Composition("H")),
        "OH": ModificationRule(
            "C-term", "OH", "OH",
            Composition("OH").mass, composition=Composition("OH")),
    }

    _custom_rules = {}

    @classmethod
    def register_new_rule(cls, rule):
        cls._custom_rules[rule.preferred_name] = rule

    @classmethod
    def remove_rule(cls, rule):
        try:
            cls._custom_rules.pop(rule.preferred_name)
        except AttributeError:
            cls._custom_rules.pop(rule)

    def __init__(self, rules=None):
        if rules is None:
            rules = self._definitions_from_stream_default()
        self.store = dict()

        for rule in rules:
            if not isinstance(rule, ModificationRule):
                rule = ModificationRule(**rule)
            self.add(rule)

        self._include_other_rules()

    def __len__(self):
        return len(self.store) + len(self._custom_rules) + len(self.other_modifications)

    def __iter__(self):
        for key in iter(self.store):
            yield key
        for key in iter(self.other_modifications):
            yield key
        for key in iter(self._custom_rules):
            yield key

    def _include_other_rules(self):
        for name, rule in self.other_modifications.items():
            self.add(rule)

    def register_alias(self, name, alias):
        rule = self[name]
        rule.aliases.add(alias)

    def rules(self):
        return set(self.store.values()) | set(self._custom_rules.values())

    def __getitem__(self, key):
        try:
            try:
                # First try the key in the custom mapping
                return self._custom_rules[key]
            except KeyError:
                # Next try the normal store
                return self.store[key]
        except KeyError:
            try:
                # Clean the key in case it is decorated with extra target information
                name = title_cleaner.search(key).groupdict()["name"]
            except (AttributeError):
                raise ModificationStringParseError("Could not parse {0}".format(key))
            try:
                return self.store[name]
            except KeyError:
                try:
                    return self._custom_rules[name]
                except KeyError:
                    raise ModificationStringParseError("Could not resolve {0}".format(name))

    def add(self, rule):
        if not isinstance(rule, ModificationRule):
            rule = ModificationRule(**rule)
        for name in rule.names:
            try:
                if self.store[name] is rule:
                    continue
                else:
                    self.store[name] += rule
            except KeyError:
                self.store[name] = rule

    def get_modification(self, name):
        return self[name]()


class RestrictedModificationTable(ModificationTable):
    def __init__(self, rules=None, constant_modifications=None, variable_modifications=None):
        super(RestrictedModificationTable, self).__init__(rules)
        if constant_modifications is None:
            constant_modifications = []
        self.constant_modifications = constant_modifications
        if variable_modifications is None:
            variable_modifications = []
        self.variable_modifications = variable_modifications
        self.collect_available_modifications()

    def collect_available_modifications(self):
        '''Clears the rules stored in the core dictionary (on self, but not its data members),
        and copies all information in the data member sub-tables into the core dictionary'''
        modifications = {}

        for name in self.constant_modifications + self.variable_modifications:
            mod = deepcopy(self[name])
            try:
                target = extract_targets_from_string(title_cleaner.search(name).groupdict()["target"])
                mod.targets = {target}
            except Exception:
                pass
            if mod.name in modifications:
                modifications[mod.name] += mod
            else:
                modifications[mod.name] = mod

        self.store.clear()
        for mod in modifications.values():
            self.add(mod)
        self._include_other_rules()


class ModificationStringParseError(Exception):
    pass


class ModificationNameResolutionError(KeyError):
    pass


class Modification(ModificationBase):

    """Represents a molecular modification, which may be bound at a given position,
    or may be present at multiple locations. This class pulls double duty,
    and when describing the total number of modifications of a type on a molecule, its
    position attributes are unused."""

    _table = ModificationTable()

    __slots__ = ["name", "mass", "position", "number", "rule", "composition"]

    @classmethod
    def register_new_rule(cls, rule):
        ModificationTable.register_new_rule(rule)

    @classmethod
    def mass_by_name(cls, name, mod_num=1):
        if len(name) == 0:
            mass = 0.0
        elif "@" == name[0]:
            mass = float(
                AnonymousModificationRule.parser.search(name).groupdict()['mass'])
        else:
            mass = cls._table[name].mass
        return mass

    def __init__(self, rule, mod_pos=-1, mod_num=1):
        name = None
        if(isinstance(rule, basestring)):
            try:
                _rule_string = rule
                rule = Modification._table[rule]
                if _rule_string in rule.aliases:
                    name = _rule_string
                else:
                    name = rule.preferred_name
            except ModificationStringParseError:
                anon = AnonymousModificationRule.try_parse(rule)
                if anon is None:
                    aa_sub = AminoAcidSubstitution.try_parse(rule)
                    if aa_sub is None:
                        glycosylation = Glycosylation.try_parse(rule)
                        if glycosylation is None:
                            raise ModificationNameResolutionError(rule)
                        else:
                            rule = glycosylation
                    else:
                        rule = aa_sub
                else:
                    rule = anon
                name = rule.name
        else:
            name = rule.preferred_name

        self.name = name
        self.mass = rule.mass
        self.position = mod_pos
        self.number = mod_num
        self.rule = rule

        try:
            self.composition = rule.composition
        except AttributeError:
            self.composition = None

    def serialize(self):
        if self.rule.is_standard:
            rep = str(self.name)
        else:
            rep = self.rule.serialize()
        if self.number > 1:
            # Numbers large than 1 will cause errors in lookup. Need
            # to incorporate a method for inferring number from string
            rep = "{0}{{{1}}}".format(rep, self.number)
        return rep

    def valid_site(self, amino_acid=None, position_modifiers=None):
        return self.rule.valid_site(amino_acid, position_modifiers)

    def why_valid(self, amino_acid=None, position_modifiers=None):
        return self.rule.why_valid(amino_acid, position_modifiers)

    def find_valid_sites(self, sequence):
        return self.rule.find_valid_sites(sequence)

    __repr__ = serialize

    __str__ = serialize

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        try:
            return (self.name == other.name) and (self.number == other.number)
        except AttributeError:
            return other in self.rule.names

    def __ne__(self, other):
        return not self == other

    def __getstate__(self):
        return [self.name, self.mass, self.position, self.number, self.rule, self.composition]

    def __setstate__(self, state):
        self.name, self.mass, self.position, self.number, self.rule, self.composition = state

    def clone(self):
        return self.__class__(self.rule)

    def is_a(self, category):
        return self.rule.is_a(category)
