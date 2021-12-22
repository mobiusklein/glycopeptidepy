'''Represent rules describing peptide modifications with
particular chemical properties, specificities, and other
metadata.

The core data structure is :class:`ModificationRule`, which
in addition to its chemical composition, mass, target specificities,
and categorical descriptors, carries a substantial amount of information
about what its names are.
'''

import re

from collections import defaultdict
try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable

from six import string_types as basestring

from ..base import ModificationBase
from ..composition import Composition, formula
from ..residue import AminoAcidResidue

from .utils import ModificationNameResolutionError
from .descriptors import ModificationCategory, SequenceLocation
from .target import (
    extract_targets_from_string, ModificationTarget,
    title_cleaner, get_position_modifier_rules_dict)


class _ModificationResolver(object):
    def __init__(self, rule_types=None):
        if rule_types is None:
            rule_types = []
        self.rule_types = list(rule_types) or []

    def register(self, rule_type):
        """Register a :class:`~.ModificationRule` subclass with
        the parsing cascade this object represents.

        This method may be used as a decorator for target classes.

        Parameters
        ----------
        rule_type: :class:`ModificationRule` subclass
            The subclass to register

        Returns
        -------
        :class:`ModificationRule` subclass
            The `rule_type` value, allowing this method to be used
            as a decorator.
        """
        self.rule_types.append(rule_type)
        return rule_type

    def __call__(self, name):
        for rule_type in self.rule_types:
            result = rule_type.try_parse(name)
            if result is not None:
                return result
        raise ModificationNameResolutionError(name)


try:
    from glycopeptidepy._c.structure.modification.rule import ModificationRuleBase
except ImportError:
    class ModificationRuleBase(ModificationBase):
        """Represent common properties of modification rules.
        """
        def is_a(self, category):
            '''Returns whether or not this :class:`ModificationRule` object belongs to
            the specified :class:`~.ModificationCategory`.

            Returns
            -------
            bool
            '''
            return category in self.categories

        def __hash__(self):
            return self._hash

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

        def is_tracked_for(self, category):
            """Determine if this :class:`ModificationRule` is tracked by a particular
            behavioral pattern associated with a :class:`~.ModificationCategory`.

            This relationship is distinct from :meth:`is_a` which merely observes that
            the semantic relationship holds, not that any actual behavior is available.

            Parameters
            ----------
            category : :class:`~.ModificationCategory`
                The category to check

            Returns
            -------
            bool
            """
            return False

        @property
        def is_standard(self):
            """Indicates whether the :class:`ModificationRule` describes a modification which
            has been specified from a reference source, and may be fully reconstructed from
            its name alone.

            This property is not strictly enforced, and does not cover all aspects of the object's
            state. It does however provide a general guide for whether or not a modification is able
            to be translated from :mod:`glycopeptidepy`'s internal representation into something that
            another program *should* be able to understand.

            Returns
            -------
            bool
            """
            return True

        def serialize(self):
            '''A string representation for inclusion in sequences'''
            return self.name


def _ModificationRule_reconstructor(tp):
    return tp.__new__(tp)


def _merge_dictlists(adict, bdict):
    result = dict()
    for key in set(adict) | set(bdict):
        result[key] = list(set(adict.get(key, [])) | set(bdict.get(key, [])))
    return result



class ModificationRule(ModificationRuleBase):
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
    title : str
        A formal or "full" name for this modification that may be its chemical name
        or some other more verbose representation which is used to identify it but
        not convenient to read or write.
    targets : set
        All applicable :class:`.ModificationTarget` instances which apply to this
        rule. This set may differ from instance to instance of the same rule by
        configuration as per the :class:`RestrictedModificationTable`, for example
    '''

    resolve = _ModificationResolver()

    @staticmethod
    def modification_tp(*args, **kwargs):
        from .modification import Modification
        return Modification(*args, **kwargs)

    @staticmethod
    def reduce(rules):
        """Given a collection of :class:`ModificationRule` instances, combine
        those with the same name into a single object using :meth:`ModificationRule.__add__`

        Parameters
        ----------
        rules : :class:`list` of :class:`ModificationRule`
            The rules to collapse

        Returns
        -------
        list
        """
        collector = defaultdict(list)
        for rule in rules:
            collector[rule.name].append(rule)
        reductions = []
        for _, rules in collector.items():
            accum = rules[0]
            for other in rules[1:]:
                accum = accum + other
            reductions.append(accum)
        return reductions

    @classmethod
    def from_protein_prospector(cls, amino_acid_specificity, modification_name,
                                title=None, monoisotopic_mass=None,
                                neutral_losses=None):
        """Build a :class:`ModificationRule` from data scraped from Protein Prospector

        Parameters
        ----------
        amino_acid_specificity : str
            The amino acid and position target specificity
        modification_name : str
            The name of the modification
        title : str, optional
            The display title for the modification, by default None
        monoisotopic_mass : float, optional
            The monoisotopic mass shift of the modification, by default None
        neutral_losses : list, optional
            Any neutral losses associated with the modification, by default None

        Returns
        -------
        :class:`ModificationRule`
        """
        if neutral_losses is None:
            neutral_losses = []
        # If the specificity is a string, parse it into rules
        targets = None
        if isinstance(amino_acid_specificity, basestring):
            try:
                targets = set(
                    [extract_targets_from_string(amino_acid_specificity)])
            except TypeError:
                raise TypeError("Could not extract targets for %s with specificity %r" % (
                    modification_name, amino_acid_specificity))

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
        """Createa a :class:`ModificationRule` instances from a :class:`dict` scraped
        from Unimod.

        Parameters
        ----------
        unimod_entry : dict
            The Unimod entry description for this modificatoin

        Returns
        -------
        :class:`ModificationRule`
        """
        specificity = unimod_entry["specificity"]
        monoisotopic_mass = unimod_entry["mono_mass"]
        full_name = unimod_entry["full_name"]
        record_id = unimod_entry.get("record_id")
        title = unimod_entry["title"]
        alt_names = set(unimod_entry.get('alt_names', ()))
        if record_id is not None:
            alt_names.add("UNIMOD:%d" % (record_id))
        neutral_losses = list()
        for losses in map(NeutralLoss.from_unimod_specificity, specificity):
            neutral_losses.append(losses)
        specificity = list(map(ModificationTarget.from_unimod_specificity, specificity))
        neutral_losses = dict(zip(specificity, neutral_losses))
        return cls(
            specificity, full_name, title, monoisotopic_mass,
            composition=Composition({str(k): int(v) for k, v in unimod_entry['composition'].items()}),
            alt_names=alt_names, neutral_losses=neutral_losses)

    @staticmethod
    def _get_preferred_name(names):
        return min([x for x in names if not x.startswith("UNIMOD:")], key=len)

    def __init__(self, amino_acid_specificity, modification_name,
                 title=None, monoisotopic_mass=None, composition=None,
                 categories=None, alt_names=None, neutral_losses=None,
                 aliases=None, **kwargs):
        if neutral_losses is None:
            neutral_losses = dict()
        if categories is None:
            categories = []
        if aliases is None:
            aliases = set()

        self.mass = float(monoisotopic_mass)
        self.composition = composition
        self.options = kwargs

        if alt_names is None:
            alt_names = set()

        self._configure_names(modification_name, title, alt_names, aliases)
        self._configure_targets(amino_acid_specificity)
        self._configure_categories(categories)

        self.neutral_losses = neutral_losses

    def clear_targets(self):
        self._n_term_targets = None
        self._c_term_targets = None
        self.targets = set()

    def _configure_targets(self, amino_acid_specificity):
        self._n_term_targets = None
        self._c_term_targets = None

        # The type of the parameter passed for amino_acid_specificity is variable
        # so select the method correct for the passed type

        # If the specificity is a string, parse it into rules
        if isinstance(amino_acid_specificity, basestring):
            try:
                self.targets = set([
                    extract_targets_from_string(amino_acid_specificity)])
            except TypeError:
                raise TypeError("Could not extract targets from %s" %
                                amino_acid_specificity)

        # If the specificity is already a rule, store it as a list of one rules
        elif isinstance(amino_acid_specificity, ModificationTarget):
            self.targets = set([amino_acid_specificity])

        # If the specificity is a collection of rules, store it as is
        elif isinstance(amino_acid_specificity, Iterable):
            self.targets = set(amino_acid_specificity)
        else:
            raise ValueError("Could not interpret target specificity from %r" % (
                amino_acid_specificity,))

    def _configure_names(self, modification_name, title, alt_names, aliases):

        # Attempt to parse the protein prospector name which contains the
        # target in
        try:
            self.common_name = title_cleaner.search(title).groupdict()['name']
        except Exception:
            self.common_name = modification_name

        self.unimod_name = modification_name
        self.title = title if title is not None else modification_name
        self.names = {self.unimod_name, self.common_name, self.title} | alt_names - {''}
        for name in list(self.names):
            self.names.update(name.split(" or "))
        self.aliases = aliases
        preferred_name = self.options.get("preferred_name")
        if preferred_name is not None:
            self.names.add(preferred_name)
            self.preferred_name = self.name = str(preferred_name)
        else:
            self.preferred_name = self.name = str(
                self._get_preferred_name(self.names))
        self._hash = hash(self.name)

    def _configure_categories(self, categories):
        self.categories = list(categories)
        categories = set(self.categories)
        for target in self.targets:
            categories.update(
                target.classification)
        self.categories = list(categories)

    def clone(self, propagated_targets=None):
        if propagated_targets is None:
            propagated_targets = (self.targets)

        dup = self.__class__(
            set(propagated_targets), self.name,
            self.title, self.mass, self.composition.copy(),
            self.categories, **self.options)
        dup.names.update(self.names)
        return dup

    def valid_site(self, amino_acid=None, position_modifiers=None):
        valid_sites = []
        for target in self.targets:
            result, site = target.valid_site(amino_acid, position_modifiers)
            if not result:
                continue
            # if SequenceLocation.anywhere is encountered, store 0 since None
            # can't be ordered with numbers
            if site.value is None:
                valid_sites.append(0)
            else:
                # otherwise store the location's numeric value
                valid_sites.append(site.value)
        if not valid_sites:
            return False
        # Find the minimum site type. Terminal sites are negative numbers
        # and prefer those from 0 which denotes SequenceLocation.anywhere
        site_type = min(valid_sites)
        if site_type == 0:
            return SequenceLocation.anywhere
        return SequenceLocation[site_type]

    def why_valid(self, amino_acid=None, position_modifiers=None):
        possible_targets = [target for target in self.targets if target.valid_site(
            amino_acid, position_modifiers)[0]]
        if not possible_targets:
            return None
        minimized_target = min(possible_targets, key=len)
        return minimized_target

    def find_valid_sites(self, sequence, protein_n_term=False, protein_c_term=False):
        valid_indices = []
        position_modifier_rules = get_position_modifier_rules_dict(
            sequence, protein_n_term=protein_n_term, protein_c_term=protein_c_term)
        for index in range(len(sequence)):
            position = sequence[index]
            if position.modifications:
                continue
            amino_acid = position[0].name
            position_modifiers = position_modifier_rules.get(index)
            is_valid = self.valid_site(amino_acid, position_modifiers)
            if is_valid:
                if (is_valid is SequenceLocation.n_term) or (is_valid is SequenceLocation.c_term):
                    valid_indices.append(is_valid)
                else:
                    valid_indices.append(index)

        return valid_indices


    def __repr__(self):
        rep = "{name}:{mass}".format(
            name=self.name, mass=self.mass)
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
            dup.name = self._get_preferred_name(dup.names)
            dup.neutral_losses = _merge_dictlists(self.neutral_losses, other.neutral_losses)
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
        return self.modification_tp(self, **kwargs)

    def __reduce__(self):
        return _ModificationRule_reconstructor, (self.__class__, ), self.__getstate__()

    def __getstate__(self):
        state = {}
        state['targets'] = self.targets
        state['mass'] = self.mass
        state['composition'] = self.composition
        state['aliases'] = self.aliases
        state['options'] = self.options
        state['categories'] = self.categories
        state['neutral_losses'] = self.neutral_losses
        state['names'] = self.names
        state['unimod_name'] = self.unimod_name
        state['title'] = self.title
        state['name'] = self.name
        state['common_name'] = self.common_name
        state['_n_term_targets'] = self.n_term_targets
        state['_c_term_targets'] = self.c_term_targets
        return state

    def __setstate__(self, state):
        self.targets = state['targets']
        self.mass = state['mass']
        self.composition = state['composition']
        self.aliases = state['aliases']
        self.options = state['options']
        self.categories = state['categories']
        self.names = state['names']
        self.neutral_losses = state['neutral_losses']
        self.title = state['title']
        # defend against outdated pickles here
        self.name = state.get(
            'name', state.get(
                'preferred_name', self._get_preferred_name(
                    state.get('names', []))))
        self.common_name = state['common_name']
        self.unimod_name = state['unimod_name']
        self._n_term_targets = state.get('_n_term_targets', state.get('_n_term_target', None))
        self._c_term_targets = state.get('_c_term_targets', state.get('_c_term_target', None))
        self._hash = hash(self.name)

    @property
    def n_term_targets(self):
        if self._n_term_targets is None:
            solutions = []
            for target in self.targets:
                if target.position_modifier == SequenceLocation.n_term or target.position_modifier == SequenceLocation.protein_n_term:
                    solutions.append(target)
            self._n_term_targets = solutions
        return self._n_term_targets

    @property
    def c_term_targets(self):
        if self._c_term_targets is None:
            solutions = []
            for target in self.targets:
                if target.position_modifier == SequenceLocation.c_term or target.position_modifier == SequenceLocation.protein_c_term:
                    solutions.append(target)
            self._c_term_targets = solutions
        return self._c_term_targets

    def as_spec_strings(self):
        for target in self.targets:
            yield "%s (%s)" % (self.name, target.serialize())

    @property
    def is_core(self):
        return True


@ModificationRule.resolve.register
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

    def valid_site(self, *args, **kwargs):  # pylint: disable=arguments-differ
        raise TypeError(
            "AnonymousModificationRule does not support site validation")

    def find_valid_sites(self, *args, **kwargs):  # pylint: disable=arguments-differ
        raise TypeError(
            "AnonymousModificationRule does not support site validation")

    def serialize(self):
        return "@" + self.name + "-" + str(self.mass)

    def clone(self, propagated_targets=None):
        dup = self.__class__(self.name, self.mass, **self.options)
        return dup

    def __repr__(self):
        rep = "{name}:{mass}".format(name=self.name, mass=self.mass)
        return rep


@ModificationRule.resolve.register
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

    @classmethod
    def _simple_sequence_mass(cls, seq_list):
        '''Because AminoAcidSubstitution may contain sequence-like entries that
        need massing but we cannot import `sequence` here, this defines a simple
        sequence list to total mass without any terminal groups, but handles residues
        and modifications.

        For a full sequence_to_mass converter, see the sequence module'''
        from .modification import Modification
        mass = 0.0
        for res, mods in seq_list:
            mass += AminoAcidResidue(res).mass
            for mod in mods:
                if mod != "":
                    mass += Modification(mod).mass
        return mass

    def __init__(self, original_residue, substitution_residue, **kwargs): # pylint: disable=super-init-not-called
        if not isinstance(original_residue, AminoAcidResidue):
            original_residue = AminoAcidResidue(original_residue)
        if not isinstance(substitution_residue, AminoAcidResidue):
            substitution_residue = AminoAcidResidue(substitution_residue)
        self.common_name = "{0}->{1}".format(original_residue, substitution_residue)
        self.original = (original_residue)
        self.substitution = (substitution_residue)
        self.mass = (self.substitution).mass - (self.original).mass
        self.name = self.common_name
        self.names = {self.common_name}
        self.options = kwargs
        self.categories = [ModificationCategory['substitution']]
        self.aliases = set()
        self._hash = hash(self.name)

    def serialize(self):
        return "@" + self.name

    @property
    def neutral_losses(self):
        return []

    def __repr__(self):
        return "{name}:{delta}".format(name=self.name, delta=self.mass)

    def __reduce__(self):
        return _ModificationRule_reconstructor, (self.__class__, ), self.__getstate__()

    def __getstate__(self):
        state = {}
        state['mass'] = self.mass
        state['composition'] = self.composition
        state['options'] = self.options
        state['categories'] = self.categories
        state['names'] = self.names
        state['name'] = self.name
        state['common_name'] = self.common_name
        state['original'] = self.original
        state['substitution'] = self.substitution
        return state

    def __setstate__(self, state):
        self.mass = state['mass']
        self.composition = state['composition']
        self.options = state['options']
        self.categories = state['categories']
        self.names = state['names']
        # defend against outdated pickles here
        self.name = state.get(
            'name', state.get(
                'preferred_name', self._get_preferred_name(
                    state.get('names', []))))
        self.common_name = state['common_name']
        self.original = state['original']
        self.substitution = state['substitution']
        self._hash = hash(self.name)



TMT10plex = [
    ('TMT10-126', 125.12044953),
    ('TMT10-127N', 126.1174845),
    ('TMT10-127C', 126.1238045),
    ('TMT10-128N', 127.1208395),
    ('TMT10-128C', 127.1271595),
    ('TMT10-129N', 128.1241945),
    ('TMT10-129C', 128.1305135),
    ('TMT10-130N', 129.1275485),
    ('TMT10-130C', 129.1338685),
    ('TMT10-131', 130.13090353)
]


try:
    from glycopeptidepy._c.structure.modification.rule import NeutralLossBase
except ImportError:
    class NeutralLossBase(object):
        def __eq__(self, other):
            return self.label == other.label and self.mass == other.mass

        def __ne__(self, other):
            return not self == other

        def __hash__(self):
            return hash(self.label)


class NeutralLoss(NeutralLossBase):

    def __init__(self, composition, mass, label=None):
        if label is None:
            label = formula(composition)
        self.composition = composition
        self.mass = mass
        self.label = label

    def __repr__(self):
        template = "{self.__class__.__name__}({self.composition}, {self.mass}, {self.label})"
        return template.format(self=self)

    @classmethod
    def from_unimod_specificity(cls, specificity):
        result = []
        for loss in specificity.get('neutral_losses', []):
            if loss['mono_mass'] == 0:
                continue
            result.append(cls(Composition(loss['composition']), loss['mono_mass']))
        return result
