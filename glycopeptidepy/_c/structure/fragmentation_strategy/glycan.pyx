cimport cython

from cpython.ref cimport Py_INCREF
from cpython cimport PyObject
from cpython.float cimport PyFloat_AsDouble
from cpython.tuple cimport PyTuple_GetItem, PyTuple_Size
from cpython.list cimport PyList_GetItem, PyList_SetItem, PyList_Size, PyList_New
from cpython.dict cimport (PyDict_GetItem, PyDict_SetItem, PyDict_Next,
                           PyDict_Keys, PyDict_Update, PyDict_DelItem, PyDict_Size)

from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.structure.fragmentation_strategy.base cimport FragmentationStrategyBase

from glycopeptidepy._c.count_table cimport CountTable
from glycopeptidepy._c.structure.fragment cimport IonSeriesBase, StubFragment, _NameTree
from glycopeptidepy._c.structure.sequence_methods cimport _PeptideSequenceCore

from glycopeptidepy.structure.fragment import (
    ChemicalShift, IonSeries, format_negative_composition)

from glycopeptidepy.structure.glycan import (
    GlycosylationType,
    HashableGlycanComposition)

from itertools import product as _product


cdef object product = _product

cdef object GlycosylationType_n_linked = GlycosylationType.n_linked
cdef IonSeriesBase IonSeries_stub_glycopeptide = IonSeries.stub_glycopeptide


@cython.final
cdef class _CompositionTree(object):
    cdef:
        public _NameTree root

    def __init__(self, root=None):
        if root is None:
            root = _NameTree()
        self.root = root

    def __getitem__(self, key):
        return self.root[key]

    cpdef object build(self, pair_sequence):
        node = self.root.traverse(pair_sequence)
        if node.name is None:
            # modify the hashable glycan composition before calculating
            # its hash value
            node.name = HashableGlycanComposition()
            for k, v in pair_sequence:
                node.name._setitem_fast(k, v)
            # Force the population of the _mass and _str caches
            node.name.mass()
            str(node.name)
        return node.name


cdef _CompositionTree _composition_tree_root = _CompositionTree()
cdef dict _composition_name_cache = dict()


cpdef _prepare_glycan_composition_from_mapping(mapping):
    pair_sequence = mapping.items()
    return _composition_tree_root.build(pair_sequence)


def n_glycan_stub_fragments(FragmentationStrategyBase self):
    cdef:
        int core_count
        list per_site_shifts
        dict site
        CComposition base_composition, composition
        double base_mass, mass
        set seen
        size_t i, n_positions
        CountTable aggregate_glycosylation
        bint is_extended, is_glycosylated
        str name_key, name
        object glycosylation
        tuple positions
        PyObject* ptemp


    glycan = self.glycan_composition()
    self._use_query = self._guess_query_mode(glycan)
    core_count = self.count_glycosylation_type(GlycosylationType_n_linked)
    per_site_shifts = []
    base_composition = self.peptide_composition()
    base_mass = base_composition.calc_mass()
    seen = set()
    for i in range(core_count):
        core_shifts = self.n_glycan_composition_fragments(glycan, core_count, i)
        per_site_shifts.append(core_shifts)
    for _positions in product(*per_site_shifts):
        positions = <tuple>_positions
        aggregate_glycosylation = CountTable._create()
        mass = base_mass
        composition = base_composition.clone()
        is_extended = False
        n_positions = PyTuple_Size(positions)
        for i in range(n_positions):
            site = <dict>PyTuple_GetItem(positions, i)
            mass += PyFloat_AsDouble(<object>PyDict_GetItem(site, 'mass'))
            if n_positions > 1:
                aggregate_glycosylation._add_from(<CountTable?>PyDict_GetItem(site, 'key'))
            else:
                aggregate_glycosylation = <CountTable?>PyDict_GetItem(site, 'key')
            composition.add_from(<CComposition>PyDict_GetItem(site, 'composition'))
            if <object>PyDict_GetItem(site, 'is_extended'):
                is_extended = True
        is_glycosylated = (mass != base_mass)
        if n_positions > 1:
            invalid = self._validate_glycan_composition(aggregate_glycosylation, glycan)
            if invalid:
                continue
        glycosylation = _prepare_glycan_composition_from_mapping(aggregate_glycosylation)
        name_key = str(glycosylation)
        ptemp = PyDict_GetItem(_composition_name_cache, name_key)
        if ptemp == NULL:
            name = StubFragment.build_name_from_composition(aggregate_glycosylation)
            PyDict_SetItem(_composition_name_cache, name_key, name)
        else:
            name = <str>ptemp

        if name in seen:
            continue
        seen.add(name)
        yield StubFragment._create(
            name=name,
            mass=mass,
            composition=composition,
            is_glycosylated=is_glycosylated,
            kind=IonSeries_stub_glycopeptide,
            chemical_shift=None,
            glycosylation=glycosylation,
            is_extended=is_extended)
