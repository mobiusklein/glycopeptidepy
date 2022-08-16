from cpython.object cimport PyObject
from cpython.dict cimport PyDict_GetItem
from glycopeptidepy._c.structure.modification.rule cimport ModificationRuleBase

from glycopeptidepy.structure.modification.target import title_cleaner
from glycopeptidepy.structure.modification.utils import ModificationStringParseError
from glycopeptidepy.structure.modification.rule import ModificationRule


cdef title_cleaner_search = title_cleaner.search
cdef object ModificationRuleImpl = ModificationRule


cdef basestring clean_name(basestring key):
    match = title_cleaner_search(key)
    if match is not None:
        return match.groupdict()["name"]
    raise ModificationStringParseError("Could not parse {}".format(key))


cdef class ModificationTableBase(object):
    '''A mixin class to provide implementation for :class:`ModificationSource`
    methods for :class:`ModificationTable`
    '''

    def __getitem__(self, basestring key):
        return self.resolve(key)

    cpdef ModificationRuleBase resolve(self, basestring key):
        try:
            rule = self._find_modification_rule(key)
        except ModificationStringParseError:
            rule = ModificationRuleImpl.resolve(key)
        return rule

    cpdef get_modification(self, name):
        return self.resolve(name)()

    cdef ModificationRuleBase _find_modification_rule(self, basestring key):
        cdef:
            PyObject* ptemp

        # First try the key in the custom mapping
        ptemp = PyDict_GetItem(self._instance_custom_rules, key)
        if ptemp != NULL:
            return <ModificationRuleBase>ptemp
        # Next try the normal store
        ptemp = PyDict_GetItem(self.store, key)
        if ptemp != NULL:
            return <ModificationRuleBase>ptemp
        # Otherwise try to look fuzzy
        return self._find_modification_rule_fuzzy(key)

    cdef ModificationRuleBase _find_modification_rule_fuzzy(self, basestring key):
        cdef:
            PyObject* ptemp
            basestring name
        # Clean the key in case it is decorated with extra target information
        # or throw an error
        name = clean_name(key)
        ptemp = PyDict_GetItem(self.store, name)
        if ptemp != NULL:
            return <ModificationRuleBase>ptemp
        ptemp = PyDict_GetItem(self._instance_custom_rules, name)
        if ptemp != NULL:
            return <ModificationRuleBase>ptemp
        raise ModificationStringParseError("Could not resolve {0}".format(name))