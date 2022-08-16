from glycopeptidepy._c.structure.modification.rule cimport ModificationRuleBase

cdef class ModificationTableBase(object):
    cdef:
        public dict store
        public dict _instance_custom_rules

    cpdef ModificationRuleBase resolve(self, basestring key)
    cpdef get_modification(self, name)

    cdef ModificationRuleBase _find_modification_rule(self, basestring key)
    cdef ModificationRuleBase _find_modification_rule_fuzzy(self, basestring key)