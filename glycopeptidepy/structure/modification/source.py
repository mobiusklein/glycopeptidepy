import csv
import json

from pkg_resources import resource_stream
from io import StringIO
from copy import deepcopy

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
        cls._custom_rules[rule.name] = rule

    @classmethod
    def remove_rule(cls, rule):
        try:
            cls._custom_rules.pop(rule.name)
        except AttributeError:
            cls._custom_rules.pop(rule)

    def __init__(self, rules=None):
        if rules is None:
            rules = self._definitions_from_stream_default()
        self.store = dict()

        for rule in rules:
            if not isinstance(rule, ModificationRule):
                rule = ModificationRule(**rule)
            if rule.name in self.other_modifications:
                continue
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

    def resolve(self, key):
        try:
            rule = self[key]
        except ModificationStringParseError:
            rule = ModificationRule.resolve(key)
        return rule

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


def rule_string_to_specialized_rule(rule_string):
    from .modificatoin import Modification
    name, targets = title_cleaner.search(rule_string).groups()
    targets = {extract_targets_from_string(targets)}
    # look up the rule in the main name cache, and copy it.
    rule = Modification(name).rule.clone()
    rule.targets = targets
    return rule

