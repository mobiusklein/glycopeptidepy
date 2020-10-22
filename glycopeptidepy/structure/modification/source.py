'''Represent queryable databases of peptide modifications from diverse
sources.
'''
import csv
import json
import bisect

from io import StringIO
from copy import deepcopy

from pkg_resources import resource_stream

try:
    from collections.abc import Mapping
except ImportError:
    from collections import Mapping

from ..composition import Composition

from .target import title_cleaner, extract_targets_from_string
from .rule import ModificationRule
from .utils import ModificationStringParseError
from .glycosylation import (
    hexnac_modification, xylose_modification,
    GlycosaminoglycanLinkerGlycosylation,
    OGlcNAcylation,
    OGlycanCoreGlycosylation,
    NGlycanCoreGlycosylation,)


def load_from_csv(stream):
    """Load a sequence of :class:`~.ModificationRule` objects from
    a CSV stream

    Parameters
    ----------
    stream : file-like
        A file-like object over a CSV

    Returns
    -------
    list of :class:`~.ModificationRule`
    """
    modification_definitions = []
    parser = csv.DictReader(stream)
    modification_definitions = [ModificationRule(**rec) for rec in parser]
    return modification_definitions


def load_from_json(stream):
    """Load a sequence of :class:`~.ModificationRule` objects from
    a JSON stream.

    :class:`~.ModificationRule` objects are constructed using the
    :meth:`~.ModificationRule.from_unimod` constructor instead of
    the default flat constructor.

    Parameters
    ----------
    stream : file-like
        A file-like object over a JSON

    Returns
    -------
    list of :class:`~.ModificationRule`
    """
    modification_definitions = []
    modification_definitions = list(map(ModificationRule.from_unimod, json.load(stream)))
    return modification_definitions


class ModificationSource(object):
    """An abstract collection of :class:`~.ModificationRule` objects.

    This class primarily deals with loading rules from file sources and
    managing global states. It should not be used instantiated directly,
    instead, use :class:`ModificationTable`.
    """

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
        '''Load the rules definitions from a CSV or JSON file and instantiate a
        :class:`ModificationSource` from it'''
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

    def resolve(self, key):
        """Locate a :class:`~.Modification` instance from
        `key`.

        May call :meth:`~.ModificationRule.resolve`.

        Parameters
        ----------
        key : str
            The name of the :class:`~.ModificationRule`

        Returns
        -------
        :class:`~.ModificationRule`
        """
        raise NotImplementedError()

    def __getitem__(self, key):
        raise NotImplementedError()

    def get_modification(self, name):
        """Create a :class:`~.Modification` from the specified name

        Parameters
        ----------
        name : :class:`str`
            The name of the :class:`~.ModificationRule` to base the :class:`~.Modification`.

        Returns
        -------
        :class:`~.Modification`
        """
        return self.resolve(name)()

    def __call__(self, key):
        """Create a :class:`~.Modification` from the specified name

        Parameters
        ----------
        name : :class:`str`
            The name of the :class:`~.ModificationRule` to base the :class:`~.Modification`.

        Returns
        -------
        :class:`~.Modification`

        See Also
        --------
        :meth:`get_modification`
        """
        return self.get_modification(key)


try:
    from glycopeptidepy._c.structure.modification.source import ModificationTableBase
except ImportError as err:

    class ModificationTableBase(object):
        '''A mixin class to provide implementation for :class:`ModificationSource`
        methods for :class:`ModificationTable`
        '''

        def __getitem__(self, key):
            try:
                # First try the key in the custom mapping
                if key in self._instance_custom_rules:
                    return self._instance_custom_rules[key]
                else:
                    # Next try the normal store
                    return self.store[key]
            except KeyError:
                try:
                    # Clean the key in case it is decorated with extra target information
                    name = title_cleaner.search(key).groupdict()["name"]
                except (AttributeError):
                    raise ModificationStringParseError(
                        "Could not parse {0}".format(key))
                try:
                    return self.store[name]
                except KeyError:
                    try:
                        return self._instance_custom_rules[name]
                    except KeyError:
                        raise ModificationStringParseError(
                            "Could not resolve {0}".format(name))

        def resolve(self, key):
            try:
                rule = self[key]
            except ModificationStringParseError:
                rule = ModificationRule.resolve(key)
            return rule


class ModificationTable(ModificationTableBase, ModificationSource):
    '''A :class:`~.Mapping` connecting modification names to :class:`~.ModificationRule`
    instances.

    '''
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
    def register_custom_rule(cls, rule):
        cls._custom_rules[rule.name] = rule

    register_new_rule = register_custom_rule

    @classmethod
    def remove_custom_rule(cls, rule):
        try:
            cls._custom_rules.pop(rule.name)
        except AttributeError:
            cls._custom_rules.pop(rule)

    def __init__(self, rules=None):
        if rules is None:
            rules = self._definitions_from_stream_default()
        self.store = dict()
        # For
        self._instance_custom_rules = self.__class__._custom_rules

        for rule in rules:
            if isinstance(rule, Mapping):
                rule = ModificationRule(**rule)
            if rule.name in self.other_modifications:
                continue
            self.add(rule)

        self._include_other_rules()

    def __len__(self):
        return len(self.store) + len(self._instance_custom_rules) + len(self.other_modifications)

    def __iter__(self):
        for key in iter(self.store):
            yield key
        for key in iter(self.other_modifications):
            yield key
        for key in iter(self._instance_custom_rules):
            yield key

    def _include_other_rules(self):
        for _name, rule in self.other_modifications.items():
            self.add(rule)

    def register_alias(self, name, alias):
        """Add a new name to an existing rule, and make that name
        resolvable.

        Parameters
        ----------
        name : :class:`str`
            The name of an existing :class:`~.ModificationRule`
        alias : :class:`str`
            The new name to add to the rule denoted by `name`.
        """
        rule = self[name]
        rule.aliases.add(alias)
        self.store[alias] = rule

    def rules(self):
        """The set of all :class:`~.ModificationRule` instances in this object

        Returns
        -------
        set
        """
        return set(self.store.values()) | set(self._custom_rules.values())

    def add(self, rule):
        """Add a new :class:`~.ModificationRule` to this object.

        If `rule` shares names with other rules and their masses match,
        they will be merged.

        Parameters
        ----------
        rule : :class:`~.ModificationRule` or :class:`Mapping`
            The rule to add. If is a :class:`Mapping`, it will be coerced
            into a :class:`~.ModificationRule`
        """
        if isinstance(rule, Mapping):
            rule = ModificationRule(**rule)
        for name in rule.names:
            try:
                if self.store[name] is rule:
                    continue
                else:
                    alt_rule = self.store[name]
                    if abs(alt_rule.mass - rule.mass) > 1e-3:
                        continue
                    self.store[name] += rule
            except KeyError:
                self.store[name] = rule

    def remove(self, rule):
        """Remove an existing :class:`~.ModificationRule` from this object.

        If `rule` is a :class:`~.ModificationRule`, all of its `names` will be
        removed. If `rule` is a :class:`str`, then only that name will be removed.

        Parameters
        ----------
        rule : :class:`~.ModificationRule` or :class:`str`
            Either a :class:`~.ModificationRule`, or a single :class:`str` denoting
            the name of the rule to remove.
        """
        try:
            names = rule.names
        except AttributeError:
            names = [rule]
        for name in names:
            self.store.pop(name, None)


class RestrictedModificationTable(ModificationTable):

    other_modifications = {}

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
                mod.clear_targets()
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


def rule_string_to_specialized_rule(rule_string):
    from .modification import Modification
    name, targets = title_cleaner.search(rule_string).groups()
    targets = {extract_targets_from_string(targets)}
    # look up the rule in the main name cache, and copy it.
    rule = Modification(name).rule.clone()
    rule.targets = targets
    return rule


class ModificationMassSearch(object):
    def __init__(self, rules):
        self.index = [(rule.mass, rule) for rule in rules]
        self.index.sort(key=lambda x: (x[0], str(x[1])))

    def find_lower_bound(self, mass, error_tolerance=1e-5):
        i = bisect.bisect_left(self.index, (mass, ''))
        while i > 0:
            rmass, _rule = self.index[i]
            if error_tolerance <= abs(rmass - mass) / mass:
                i -= 1
            else:
                break
        return i

    def search_mass(self, mass, error_tolerance=1e-5):
        i = self.find_lower_bound(mass, error_tolerance)
        out = []
        for j in range(i, len(self.index)):
            rmass, rule = self.index[j]
            if abs(rmass - mass) / mass > error_tolerance:
                break
            out.append(rule)
        return out
