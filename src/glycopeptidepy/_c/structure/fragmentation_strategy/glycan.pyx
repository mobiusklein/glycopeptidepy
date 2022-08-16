# cython: embedsignature=True
cimport cython

from cpython.ref cimport Py_INCREF
from cpython cimport PyObject
from cpython.float cimport PyFloat_AsDouble
from cpython.tuple cimport PyTuple_GetItem, PyTuple_Size, PyTuple_GET_ITEM, PyTuple_GET_SIZE
from cpython.list cimport PyList_GetItem, PyList_SetItem, PyList_Size, PyList_New, PyList_GET_ITEM, PyList_GET_SIZE, PyList_Append, PyList_SET_ITEM
from cpython.int cimport PyInt_AsLong
from cpython.dict cimport (PyDict_GetItem, PyDict_SetItem, PyDict_Next,
                           PyDict_Keys, PyDict_Update, PyDict_DelItem, PyDict_Size)

from glypy.composition.ccomposition cimport CComposition
from glypy._c.structure.glycan_composition cimport _CompositionBase

from glycopeptidepy._c.structure.fragmentation_strategy.base cimport FragmentationStrategyBase

from glycopeptidepy._c.count_table cimport CountTable, CountTableIterator
from glycopeptidepy._c.structure.fragment cimport (
    IonSeriesBase, StubFragment, _NameTree,
    build_name_from_composition)

from glycopeptidepy._c.structure.sequence_methods cimport _PeptideSequenceCore
from glycopeptidepy._c.structure.glycan cimport GlycosylationManager, GlycanCompositionProxy, glycan_composition_type

from glycopeptidepy.structure.fragment import (
    ChemicalShift, IonSeries, format_negative_composition)

from glycopeptidepy.structure.glycan import (
    GlycosylationType,
    HashableGlycanComposition)

from itertools import product as _product
from glypy.structure.glycan_composition import (
    FrozenMonosaccharideResidue)

cdef object product = _product

cdef:
    object GlycosylationType_n_linked = GlycosylationType.n_linked
    object GlycosylationType_o_linked = GlycosylationType.o_linked
    object GlycosylationType_glycosaminoglycan = GlycosylationType.glycosaminoglycan

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

    cpdef _NameTree build(self, list pair_sequence):
        cdef:
            tuple item
            int i, n
            _CompositionBase value
            _NameTree node
        node = self.root.traverse(pair_sequence)
        if node.name is None:
            # modify the hashable glycan composition before calculating
            # its hash value
            value = node.name = <_CompositionBase>HashableGlycanComposition()
            n = PyList_Size(pair_sequence)
            for i in range(n):
                item = <tuple>PyList_GET_ITEM(pair_sequence, i)
                k = <object>PyTuple_GET_ITEM(item, 0)
                v = <object>PyTuple_GET_ITEM(item, 1)
                value._setitem_fast(k, v)
            # Force the population of the _mass and _str caches
            node.name.mass()
            node.value = str(node.name)
        return node


cdef _CompositionTree _composition_tree_root = _CompositionTree()
cdef dict _composition_name_cache = dict()


cdef _NameTree _prepare_glycan_composition_from_mapping(CountTable mapping):
    pair_sequence = mapping.items()
    return _composition_tree_root.build(pair_sequence)

cdef:
    object _HEX = FrozenMonosaccharideResidue.from_iupac_lite("Hex")
    CComposition _HEX_composition = _HEX.total_composition()
    double _HEX_mass = _HEX.mass()
    object _HEXNAC = FrozenMonosaccharideResidue.from_iupac_lite("HexNAc")
    CComposition _HEXNAC_composition = _HEXNAC.total_composition()
    double _HEXNAC_mass = _HEXNAC.mass()
    object _FUC = FrozenMonosaccharideResidue.from_iupac_lite("Fuc")
    CComposition _FUC_composition = _FUC.total_composition()
    double _FUC_mass = _FUC.mass()
    object _DHEX = FrozenMonosaccharideResidue.from_iupac_lite("dHex")
    CComposition _DHEX_composition = _DHEX.total_composition()
    double _DHEX_mass = _DHEX.mass()
    object _XYL = FrozenMonosaccharideResidue.from_iupac_lite("Xyl")
    CComposition _XYL_composition = _XYL.total_composition()
    double _XYL_mass = _XYL.mass()
    object _AHEX = FrozenMonosaccharideResidue.from_iupac_lite("aHex")
    CComposition _AHEX_composition = _AHEX.total_composition()
    double _AHEX_mass = _AHEX.mass()

    object _Hex2NAc = FrozenMonosaccharideResidue.from_iupac_lite("Hex2NAc")
    object _Glc2NAc = FrozenMonosaccharideResidue.from_iupac_lite("Glc2NAc")
    object _Gal2NAc = FrozenMonosaccharideResidue.from_iupac_lite("Gal2NAc")


cdef class GlycanCompositionFragmentStrategyBase(FragmentationStrategyBase):
    def __init__(self, peptide, use_query=False, *args, **kwargs):
        super(GlycanCompositionFragmentStrategyBase, self).__init__(peptide, *args, **kwargs)
        # to be initialized closer to when the glycan composition will
        # be used. Determine whether to use the fast-path __getitem__ or
        # go through the slower but more general GlycanComposition.query
        # method to look up monosaccharides.
        self._use_query = use_query
        self._generator = None
        if self.peptide is not None:
            self.glycosylation_manager = self.peptide._glycosylation_manager
        else:
            self.glycosylation_manager = None

    cpdef GlycanCompositionProxy glycan_composition(self):
        return self.glycosylation_manager.get_glycan_composition()

    cdef bint _guess_query_mode(self, GlycanCompositionProxy glycan_composition):
        # these guesses will work for N-glycans and common types of mucin-type O-glycans
        # and GAG linkers
        flag = PyInt_AsLong(glycan_composition._getitem_fast(_Hex2NAc)) +\
               PyInt_AsLong(glycan_composition._getitem_fast(_Glc2NAc)) +\
               PyInt_AsLong(glycan_composition._getitem_fast(_Gal2NAc))
        return flag or self._use_query

    cpdef long count_glycosylation_type(self, glycotype):
        return GlycosylationManager.count_glycosylation_type(self.glycosylation_manager, glycotype)


@cython.freelist(10000)
cdef class GlycanCompositionFragment(object):

    def __init__(self, mass, composition, key, is_extended=False):
        self.mass = mass
        self.composition = composition
        self.key = key
        self.is_extended = is_extended
        self._hash_key = -1
        self._glycosylation_size = -1
        self.cache = None

    def __hash__(self):
        cdef _NameTree node
        if self._hash_key == -1:
            node = self._get_glycan_composition()
            self._hash_key = hash(node.name)
        return self._hash_key

    cpdef _NameTree _get_glycan_composition(self):
        return _prepare_glycan_composition_from_mapping(self.key)

    def __eq__(self, other):
        return self.key == other.key and self.is_extended == other.is_extended

    @staticmethod
    cdef GlycanCompositionFragment _create(double mass, CComposition composition, CountTable key, bint is_extended):
        cdef:
            GlycanCompositionFragment self
        self = GlycanCompositionFragment.__new__(GlycanCompositionFragment)
        self.mass = mass
        self.composition = composition
        self.key = key
        self.is_extended = is_extended
        self._hash_key = -1
        self.cache = None
        return self

    def __getitem__(self, key):
        if key == "mass":
            return self.mass
        elif key == "composition":
            return self.composition
        elif key == "key":
            return self.key
        elif key == 'is_extended':
            return self.is_extended
        else:
            raise KeyError(key)

    def __setitem__(self, key, value):
        if key == "mass":
            self.mass = value
        elif key == "composition":
            self.composition = value
        elif key == "key":
            self.key = value
        elif key == 'is_extended':
            self.is_extended = value
        else:
            raise KeyError(key)

    def __repr__(self):
        template = "{self.__class__.__name__}({self.mass}, {self.composition}, {self.key}, {self.is_extended})"
        return template.format(self=self)

    cpdef GlycanCompositionFragment copy(self):
        cdef:
            GlycanCompositionFragment inst
        inst = GlycanCompositionFragment._create(
            self.mass,
            self.composition.clone() if self.composition is not None else None, self.key.copy(),
            self.is_extended)
        inst._hash_key = self._hash_key
        inst.cache = self.cache
        return inst

    cdef int get_glycosylation_size(self):
        if self._glycosylation_size == -1:
            self._glycosylation_size = self.key.total()
        return self._glycosylation_size

    @property
    def glycosylation_size(self):
        return self.get_glycosylation_size()


cdef dict hexnac_hex_composition_cache = dict()

cdef CComposition hexnac_hex_composition(long hexnac, long hexose):
    cdef:
        tuple key
        CComposition result
        PyObject* ptemp

    key = (hexnac, hexose)
    ptemp = PyDict_GetItem(hexnac_hex_composition_cache, key)
    if ptemp == NULL:
        result = (_HEXNAC_composition * hexnac)
        if hexose > 0:
            result.add_from(<CComposition>(_HEX_composition * hexose))
        PyDict_SetItem(hexnac_hex_composition_cache, key, result.clone())
        return result
    else:
        result = <CComposition>ptemp
        return result.clone()


cdef class StubGlycopeptideStrategy(GlycanCompositionFragmentStrategyBase):

    """A fragmentation strategy that generates intact peptide + glycan Y fragments from
    glycopeptides where the glycan structure is not fully known. When the structure is
    known, :class:`CADFragmentationStrategy` should be used instead

    Attributes
    ----------
    extended : bool
        Whether to fragment beyond the conserved core
    extended_fucosylation : bool
        Whether to consider multiple ``Fuc`` residues per fragment rather than one per core
    """

    _composition_tree = _composition_tree_root

    def __init__(self, peptide, extended=True, use_query=False, extended_fucosylation=False, **kwargs):
        self.extended = extended
        self.extended_fucosylation = extended_fucosylation
        super(StubGlycopeptideStrategy, self).__init__(peptide, use_query)

    def __next__(self):
        if self._generator is None:
            self._generator = iter(self.stub_fragments())
        return next(self._generator)

    cpdef bint _validate_glycan_composition(self, CountTable aggregate_glycosylation, glycan):
        cdef:
            bint invalid = False
            int i, n
            long v, ov
            object key
            PyObject* pkey
            CountTableIterator iterator
        if self._use_query:
            for key, value in aggregate_glycosylation.items():
                if glycan.query(key) < value:
                    invalid = True
                    break
        else:
            iterator = CountTableIterator._create(aggregate_glycosylation)
            pkey = NULL
            getter = glycan._getitem_fast
            while iterator.has_more():
                iterator.get_next_value(&pkey, &v)
                key = <object>pkey
                ov = PyInt_AsLong(getter(key))
                if ov < v:
                    invalid = True
                    break
        return invalid

    cpdef GlycanCompositionFragment fucosylate_increment(self, GlycanCompositionFragment shift):
        cdef:
            GlycanCompositionFragment fucosylated
        fucosylated = shift.copy()
        fucosylated.mass += _FUC_mass
        if self.compute_compositions:
            fucosylated.composition.add_from(_FUC_composition)
        fucosylated.key.increment(_FUC, 1)
        return fucosylated

    cpdef GlycanCompositionFragment xylosylate_increment(self, GlycanCompositionFragment shift):
        cdef:
            GlycanCompositionFragment xylosylated
        xylosylated = shift.copy()
        xylosylated.mass += _XYL_mass
        if self.compute_compositions:
            xylosylated.composition.add_from(_XYL_composition)
        xylosylated.key.increment(_XYL, 1)
        return xylosylated

    cpdef list fucosylate_extended(self, GlycanCompositionFragment shift, long fucose_count):
        cdef:
            GlycanCompositionFragment fucosylated
            size_t i
            list result
        result = [None for i in range(fucose_count)]
        for i in range(1, fucose_count + 1):
            fucosylated = shift.copy()
            fucosylated.mass += _FUC_mass * i
            if self.compute_compositions:
                fucosylated.composition.add_from(<CComposition>(i * _FUC_composition))
            fucosylated.key.increment(_FUC, i)
            result[i - 1] = fucosylated
        return result

    cpdef GlycanCompositionFragment modified_increment(self, GlycanCompositionFragment base, object shift):
        cdef:
            GlycanCompositionFragment modified
        modified = base.copy()

        modified.mass += shift.mass()
        if self.compute_compositions:
            modified.composition.add_from(shift.total_composition())
        modified.key.increment(shift, 1)
        return modified

    cpdef list n_glycan_composition_fragments(self, object glycan, int core_count=1, int  iteration_count=0):
        """Generate theoretical N-glycan Y fragment compositions containing the core motif
        plus a portion of the extended branches if :attr:`extended` is used. Unless :attr:`extend

        Parameters
        ----------
        glycan : :class:`~.GlycanComposition`
            Description
        core_count : int, optional
            The total number of glycans attached to the peptide
        iteration_count : int, optional
            The core index. Unless :attr:`extended_fucosylation` is true, this will
            limit the number of Fucose residues to one per core at most, and if
            the iteration count exceeds the number of cores no Fucose will be added

        Returns
        -------
        list of dict
            The compositions corresponding to the theoretical fragments
        """

        cdef:
            long fucose_count, xylose_count
            long hexnac_count, hexose_count
            long hexnac_in_aggregate, hexose_in_aggregate
            long base_hexnac, extra_hexnac_count , extra_hexose_count
            CountTable key_ct
            list core_shifts

        if self._use_query:
            fucose_count = glycan.query('Fuc') + glycan.query('dHex')
            xylose_count = glycan.query('Xyl')
            hexnac_in_aggregate = glycan.query('HexNAc')
            hexose_in_aggregate = glycan.query('Hex')
        else:
            fucose_count = glycan[_FUC] + glycan[_DHEX]
            xylose_count = glycan[_XYL]
            hexnac_in_aggregate = glycan[_HEXNAC]
            hexose_in_aggregate = glycan[_HEX]

        core_shifts = []
        base_hexnac = min(hexnac_in_aggregate + 1, 3)
        for hexnac_count in range(base_hexnac):
            if hexnac_count == 0:
                shift = GlycanCompositionFragment._create(0,
                    CComposition() if self.compute_compositions else None,
                    CountTable._create(),
                    False)
                core_shifts.append(shift)
            elif hexnac_count == 1:
                key_ct = CountTable._create()
                key_ct.setitem(_HEXNAC, hexnac_count)
                shift = GlycanCompositionFragment._create(
                    (hexnac_count * _HEXNAC_mass),
                    hexnac_hex_composition(hexnac_count, 0) if self.compute_compositions else None,
                    key_ct,
                    False,
                )
                core_shifts.append(shift)
                if iteration_count < fucose_count:
                    fucosylated = self.fucosylate_increment(shift)
                    core_shifts.append(fucosylated)
            elif hexnac_count == 2:
                key_ct = CountTable._create()
                key_ct.setitem(_HEXNAC, hexnac_count)
                shift = GlycanCompositionFragment._create(
                    (hexnac_count * _HEXNAC_mass),
                    hexnac_hex_composition(hexnac_count, 0) if self.compute_compositions else None,
                    key_ct,
                    False,
                )
                core_shifts.append(shift)

                if not self.extended_fucosylation:
                    if iteration_count < fucose_count:
                        fucosylated = self.fucosylate_increment(shift)
                        core_shifts.append(fucosylated)
                        if iteration_count < xylose_count:
                            xylosylated = self.xylosylate_increment(fucosylated)
                            core_shifts.append(xylosylated)
                elif fucose_count > 0:
                    core_shifts.extend(self.fucosylate_extended(shift, fucose_count))
                    if iteration_count < xylose_count:
                        xylosylated = self.xylosylate_increment(fucosylated)
                        core_shifts.append(xylosylated)

                if iteration_count < xylose_count:
                    xylosylated = self.xylosylate_increment(shift)
                    core_shifts.append(xylosylated)

                for hexose_count in range(1, min(hexose_in_aggregate + 1, 4)):
                    key_ct = CountTable._create()
                    key_ct.setitem(_HEXNAC, hexnac_count)
                    key_ct.setitem(_HEX, hexose_count)
                    shift = GlycanCompositionFragment._create(
                        (hexnac_count * _HEXNAC_mass) + (hexose_count * _HEX_mass),
                        hexnac_hex_composition(hexnac_count, hexose_count) if self.compute_compositions else None,
                        key_ct,
                        False,
                    )
                    core_shifts.append(shift)

                    if not self.extended_fucosylation:
                        if iteration_count < fucose_count:
                            fucosylated = self.fucosylate_increment(shift)
                            core_shifts.append(fucosylated)
                            if iteration_count < xylose_count:
                                xylosylated = self.xylosylate_increment(fucosylated)
                                core_shifts.append(xylosylated)
                    elif fucose_count > 0:
                        core_shifts.extend(self.fucosylate_extended(shift, fucose_count))
                        if iteration_count < xylose_count:
                            xylosylated = self.xylosylate_increment(fucosylated)
                            core_shifts.append(xylosylated)

                    if iteration_count < xylose_count:
                        xylosylated = self.xylosylate_increment(shift)
                        core_shifts.append(xylosylated)

                    # After the core motif has been exhausted, speculatively add
                    # on the remaining core monosaccharides sequentially until
                    # exhausted. This doesn't allow for one missing alpha-mannose while still
                    # showing an extension on the other.
                    if hexose_count == 3 and hexnac_in_aggregate >= 2 * core_count and self.extended:
                        for extra_hexnac_count in range(0, hexnac_in_aggregate - hexnac_count + 1):
                            if extra_hexnac_count + hexnac_count > hexnac_in_aggregate:
                                continue
                            if extra_hexnac_count > 0:
                                key_ct = CountTable._create()
                                key_ct.setitem(_HEXNAC, hexnac_count + extra_hexnac_count)
                                key_ct.setitem(_HEX, hexose_count)
                                shift = GlycanCompositionFragment._create(
                                    ((hexnac_count + extra_hexnac_count) * _HEXNAC_mass) + (
                                        hexose_count * _HEX_mass),
                                    hexnac_hex_composition(hexnac_count + extra_hexnac_count, hexose_count) if self.compute_compositions else None,
                                    key_ct,
                                    True
                                )
                                core_shifts.append(shift)

                                if not self.extended_fucosylation:
                                    if iteration_count < fucose_count:
                                        fucosylated = self.fucosylate_increment(shift)
                                        core_shifts.append(fucosylated)
                                        if iteration_count < xylose_count:
                                            xylosylated = self.xylosylate_increment(fucosylated)
                                            core_shifts.append(xylosylated)
                                elif fucose_count > 0:
                                    core_shifts.extend(self.fucosylate_extended(shift, fucose_count))

                                if iteration_count < xylose_count:
                                    xylosylated = self.xylosylate_increment(shift)
                                    core_shifts.append(xylosylated)

                            for extra_hexose_count in range(1, hexose_in_aggregate - hexose_count + 1):
                                if extra_hexose_count + hexose_count > hexose_in_aggregate:
                                    continue
                                key_ct = CountTable._create()
                                key_ct.setitem(_HEXNAC, hexnac_count + extra_hexnac_count)
                                key_ct.setitem(_HEX, hexose_count + extra_hexose_count)
                                shift = GlycanCompositionFragment._create(
                                    ((hexnac_count + extra_hexnac_count) * _HEXNAC_mass) + (
                                        (hexose_count + extra_hexose_count) * _HEX_mass),
                                    hexnac_hex_composition(
                                        hexnac_count + extra_hexnac_count, hexose_count + extra_hexose_count) if self.compute_compositions else None,

                                    key_ct,
                                    True
                                )
                                core_shifts.append(shift)

                                if not self.extended_fucosylation:
                                    if iteration_count < fucose_count:
                                        fucosylated = self.fucosylate_increment(shift)
                                        core_shifts.append(fucosylated)
                                        if iteration_count < xylose_count:
                                            xylosylated = self.xylosylate_increment(fucosylated)
                                            core_shifts.append(xylosylated)
                                elif fucose_count > 0:
                                    core_shifts.extend(self.fucosylate_extended(shift, fucose_count))
                                    if iteration_count < xylose_count:
                                        xylosylated = self.xylosylate_increment(fucosylated)
                                        core_shifts.append(xylosylated)

                                if iteration_count < xylose_count:
                                    xylosylated = self.xylosylate_increment(shift)
                                    core_shifts.append(xylosylated)
        return core_shifts

    cdef list _combinate_sites_single(self, list per_site_shifts, object glycan):
        cdef:
            list result
            Py_ssize_t i, n
            GlycanCompositionFragment intermediate
            CountTable aggregate_glycosylation
            _NameTree glycosylation
            PyObject* ptemp
            tuple triple

        result = []
        n = PyList_GET_SIZE(per_site_shifts)
        for i in range(n):
            intermediate = <GlycanCompositionFragment>PyList_GET_ITEM(per_site_shifts, i)
            aggregate_glycosylation = intermediate.key
            glycosylation = intermediate.cache
            if glycosylation is None:
                glycosylation = _prepare_glycan_composition_from_mapping(
                    aggregate_glycosylation)
                intermediate.cache = glycosylation
            name_key = glycosylation.value

            ptemp = PyDict_GetItem(_composition_name_cache, name_key)
            if ptemp == NULL:
                name = build_name_from_composition(aggregate_glycosylation)
                PyDict_SetItem(_composition_name_cache, name_key, name)
            else:
                name = <str>ptemp
            triple = (intermediate, name, glycosylation.name)
            PyList_Append(result, triple)
        return result


    cpdef list _combinate_sites(self, list per_site_shifts, object glycan):
        '''Combine the glycan fragments from each site together to form an aggregate
        fragment glycan composition.
        '''
        cdef:
            list result
            list combos
            Py_ssize_t i, n, j, m
            tuple positions
            bint is_extended, invalid
            CComposition composition
            GlycanCompositionFragment intermediate, site
            CountTable aggregate_glycosylation
            _NameTree glycosylation
            set seen
            int core_count
            double total_mass
            double acc_mass

        core_count = PyList_Size(per_site_shifts)
        if core_count == 1:
            # Fast path that doesn't need to use itertools.product or any of the extra
            # unpacking and validation
            combos = <list>PyList_GET_ITEM(per_site_shifts, 0)
            return self._combinate_sites_single(combos, glycan)
        combos = list(product(*per_site_shifts))
        result = []
        seen = set()
        n = PyList_Size(combos)
        total_mass = PyFloat_AsDouble(glycan.mass())
        for i in range(n):
            positions = <tuple>PyList_GET_ITEM(combos, i)
            aggregate_glycosylation = None
            glycosylation = None
            composition = None

            intermediate = <GlycanCompositionFragment>PyTuple_GET_ITEM(positions, 0)
            acc_mass = intermediate.mass
            m = core_count

            # A quick pre-scan to skip all combinations that exceed the mass of the total
            for j in range(1, m):
                site = <GlycanCompositionFragment>PyTuple_GET_ITEM(positions, j)
                acc_mass += site.mass
            invalid = total_mass < acc_mass
            if invalid:
                continue

            # Now do the expensive accumulation
            aggregate_glycosylation = intermediate.key.copy()
            if self.compute_compositions:
                composition = intermediate.composition.copy()
            else:
                composition = None
            is_extended = intermediate.is_extended
            for j in range(1, m):
                site = <GlycanCompositionFragment>PyTuple_GET_ITEM(positions, j)
                aggregate_glycosylation._add_from(site.key)
                if self.compute_compositions:
                    composition.add_from(site.composition)
                is_extended |= site.is_extended
            invalid = self._validate_glycan_composition(
                aggregate_glycosylation, glycan)
            if invalid:
                continue

            intermediate = GlycanCompositionFragment._create(
                acc_mass, composition, aggregate_glycosylation, is_extended)

            if glycosylation is None:
                glycosylation = _prepare_glycan_composition_from_mapping(
                    aggregate_glycosylation)
                if core_count == 1:
                    intermediate.cache = glycosylation
            name_key = glycosylation.value
            # May not need to do this second round of building.
            ptemp = PyDict_GetItem(_composition_name_cache, name_key)
            if ptemp == NULL:
                name = build_name_from_composition(aggregate_glycosylation)
                PyDict_SetItem(_composition_name_cache, name_key, name)
            else:
                name = <str>ptemp
            if name in seen:
                continue
            seen.add(name)
            result.append((intermediate, name, glycosylation.name))
        return result

    cpdef list n_glycan_stub_fragments(self):
        cdef:
            int core_count
            list per_site_shifts
            GlycanCompositionFragment site
            CComposition base_composition, composition
            double base_mass
            size_t i
            size_t n_combos
            bint is_extended, is_glycosylated
            str name_key, name
            object glycosylation
            GlycanCompositionProxy glycan
            list accumulator
            list combos
            tuple result

        glycan = self.glycan_composition()
        self._use_query = self._guess_query_mode(glycan)
        core_count = self.count_glycosylation_type(GlycosylationType_n_linked)
        per_site_shifts = []
        base_composition = None
        if self.compute_compositions:
            base_composition = self.peptide_composition()
            base_mass = base_composition.calc_mass()
        else:
            base_mass = self.peptide.get_peptide_backbone_mass()
        for i in range(core_count):
            core_shifts = self.n_glycan_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)
        combos = self._combinate_sites(per_site_shifts, glycan)
        n_combos = PyList_Size(combos)
        accumulator = PyList_New(n_combos)
        for i in range(n_combos):
            result = <tuple>PyList_GetItem(combos, i)
            site = <GlycanCompositionFragment>PyTuple_GetItem(result, 0)
            name = <str>PyTuple_GetItem(result, 1)
            glycosylation = <object>PyTuple_GetItem(result, 2)
            is_glycosylated = (site.mass != 0)

            composition = None
            if self.compute_compositions:
                composition = base_composition.clone()
                composition.add_from(site.composition)

            frag = StubFragment._create(
                name=name,
                mass=base_mass + site.mass,
                composition=composition,
                is_glycosylated=is_glycosylated,
                kind=IonSeries_stub_glycopeptide,
                chemical_shift=None,
                glycosylation=glycosylation,
                is_extended=site.is_extended,
                glycosylation_size=site.get_glycosylation_size())
            Py_INCREF(frag)
            PyList_SET_ITEM(accumulator, i, frag)
        return accumulator

    cpdef list o_glycan_stub_fragments(self):
        cdef:
            list per_site_shifts
            CComposition base_composition, composition
            double base_mass, mass
            size_t core_count, i, n_positions
            tuple positions
            GlycanCompositionFragment site
            CountTable aggregate_glycosylation
            bint is_glycosylated, is_extended
            GlycanCompositionProxy glycan
            str name
            set seen
            list accumulator

        glycan = self.glycan_composition()
        self._use_query = self._guess_query_mode(glycan)
        core_count = self.count_glycosylation_type(GlycosylationType_o_linked)
        per_site_shifts = []

        base_composition = None
        if self.compute_compositions:
            base_composition = self.peptide_composition()
            base_mass = base_composition.calc_mass()
        else:
            base_mass = self.peptide.get_peptide_backbone_mass()

        for i in range(core_count):
            core_shifts = self.o_glycan_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)
        combos = self._combinate_sites(per_site_shifts, glycan)
        n_combos = PyList_Size(combos)
        accumulator = PyList_New(n_combos)
        for i in range(n_combos):
            result = <tuple>PyList_GetItem(combos, i)
            site = <GlycanCompositionFragment>PyTuple_GetItem(result, 0)
            name = <str>PyTuple_GetItem(result, 1)
            glycosylation = <object>PyTuple_GetItem(result, 2)
            is_glycosylated = (site.mass != 0)
            composition = None
            if self.compute_compositions:
                composition = base_composition.clone()
                composition.add_from(site.composition)

            frag = StubFragment._create(
                name=name,
                mass=base_mass + site.mass,
                composition=composition,
                is_glycosylated=is_glycosylated,
                kind=IonSeries_stub_glycopeptide,
                chemical_shift=None,
                glycosylation=glycosylation,
                is_extended=site.is_extended,
                glycosylation_size=site.get_glycosylation_size())
            Py_INCREF(frag)
            PyList_SET_ITEM(accumulator, i, frag)
        return accumulator

    cpdef list o_glycan_composition_fragments(self, glycan, long core_count=1, long iteration_count=0):
        cdef:
            list core_shifts
            long fucose_count, hexnac_in_aggregate, hexose_in_aggregate
            long hexnac_count, hexose_count
            long extra_hexnac_count, extra_hexose_count
            CountTable key_ct

        if self._use_query:
            fucose_count = glycan.query('Fuc') + glycan.query('dHex')
            hexnac_in_aggregate = glycan.query('HexNAc')
            hexose_in_aggregate = glycan.query('Hex')
        else:
            fucose_count = glycan[_FUC] + glycan[_DHEX]
            hexnac_in_aggregate = glycan[_HEXNAC]
            hexose_in_aggregate = glycan[_HEX]
        core_shifts = []
        for hexnac_count in range(3):
            if hexnac_in_aggregate < hexnac_count:
                continue
            if hexnac_count == 0:
                shift = GlycanCompositionFragment._create(
                    0,
                    CComposition(),
                    CountTable._create(),
                    False
                )
                core_shifts.append(shift)
            elif hexnac_count >= 1:
                key_ct = CountTable._create()
                key_ct.setitem(_HEXNAC, hexnac_count)
                shift = GlycanCompositionFragment._create(
                    (hexnac_count * _HEXNAC_mass),
                    hexnac_hex_composition(hexnac_count, 0),
                    key_ct,
                    False
                )
                core_shifts.append(shift)
                if iteration_count < fucose_count:
                    fucosylated = self.fucosylate_increment(shift)
                    core_shifts.append(fucosylated)
                for hexose_count in range(0, 2):
                    if hexose_in_aggregate < hexose_count:
                        continue
                    if hexose_count > 0:
                        key_ct = CountTable._create()
                        key_ct.setitem(_HEXNAC, hexnac_count)
                        key_ct.setitem(_HEX, hexose_count)
                        shift = GlycanCompositionFragment._create(
                            (
                                (hexnac_count) * _HEXNAC_mass) + (
                                (hexose_count) * _HEX_mass),
                            hexnac_hex_composition(hexnac_count, hexose_count),
                            key_ct,
                            False
                        )
                        core_shifts.append(shift)
                        if iteration_count < fucose_count:
                            fucosylated = self.fucosylate_increment(shift)
                            core_shifts.append(fucosylated)
                    # After the core motif has been exhausted, speculatively add
                    # on the remaining core monosaccharides sequentially until
                    # exhausted.
                    if self.extended and hexnac_in_aggregate - hexnac_count >= 0:
                        for extra_hexnac_count in range(hexnac_in_aggregate - hexnac_count + 1):
                            if extra_hexnac_count:
                                key_ct = CountTable._create()
                                key_ct.setitem(_HEXNAC, hexnac_count + extra_hexnac_count)
                                if hexose_count > 0:
                                    key_ct.setitem(_HEX, hexose_count)
                                shift = GlycanCompositionFragment._create(
                                    (
                                        (hexnac_count + extra_hexnac_count) * _HEXNAC_mass) + (
                                        (hexose_count) * _HEX_mass),
                                    hexnac_hex_composition(hexnac_count + extra_hexnac_count, hexose_count),
                                    key_ct,
                                    True,
                                )
                                core_shifts.append(shift)
                                if iteration_count < fucose_count:
                                    fucosylated = self.fucosylate_increment(shift)
                                    core_shifts.append(fucosylated)
                            if hexose_in_aggregate > hexose_count and hexose_count > 0:
                                for extra_hexose_count in range(hexose_in_aggregate - hexose_count):
                                    if extra_hexose_count > 0 and extra_hexose_count + hexose_count > 0:
                                        key_ct = CountTable._create()
                                        if hexnac_count > 0:
                                            key_ct.setitem(_HEXNAC, hexnac_count + extra_hexnac_count)
                                        key_ct.setitem(_HEX, hexose_count + extra_hexose_count)
                                        shift = GlycanCompositionFragment._create(
                                            (
                                                (hexnac_count + extra_hexnac_count) * _HEXNAC_mass) + (
                                                (hexose_count + extra_hexose_count) * _HEX_mass),
                                            hexnac_hex_composition(
                                                hexnac_count + extra_hexnac_count, hexose_count + extra_hexose_count),
                                            key_ct,
                                            True
                                        )
                                        core_shifts.append(shift)
                                        if iteration_count < fucose_count:
                                            fucosylated = self.fucosylate_increment(
                                                shift)
                                            core_shifts.append(fucosylated)

        return core_shifts

    cpdef list gag_linker_stub_fragments(self):
        cdef:
            list per_site_shifts
            CComposition base_composition, composition
            double base_mass, mass
            size_t core_count, i, n_positions
            tuple positions
            GlycanCompositionFragment site
            CountTable aggregate_glycosylation
            bint is_glycosylated, is_extended
            GlycanCompositionProxy glycan
            str name
            set seen
            list accumulator

        glycan = self.glycan_composition()
        self._use_query = self._guess_query_mode(glycan)
        core_count = self.count_glycosylation_type(GlycosylationType_glycosaminoglycan)
        per_site_shifts = []
        base_composition = None
        if self.compute_compositions:
            base_composition = self.peptide_composition()
            base_mass = base_composition.calc_mass()
        else:
            base_mass = self.peptide.get_peptide_backbone_mass()
        for i in range(core_count):
            core_shifts = self.gag_linker_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)
        combos = self._combinate_sites(per_site_shifts, glycan)
        n_combos = PyList_Size(combos)
        accumulator = []
        for i in range(n_combos):
            result = <tuple>PyList_GetItem(combos, i)
            site = <GlycanCompositionFragment>PyTuple_GetItem(result, 0)
            name = <str>PyTuple_GetItem(result, 1)
            glycosylation = <object>PyTuple_GetItem(result, 2)
            is_glycosylated = (site.mass != 0)
            composition = None
            if self.compute_compositions:
                composition = base_composition.clone()
                composition.add_from(site.composition)
            accumulator.append(
                StubFragment._create(
                    name=name,
                    mass=base_mass + site.mass,
                    composition=composition,
                    is_glycosylated=is_glycosylated,
                    kind=IonSeries_stub_glycopeptide,
                    chemical_shift=None,
                    glycosylation=glycosylation,
                    is_extended=site.is_extended,
                    glycosylation_size=site.get_glycosylation_size()))
        return accumulator

    cpdef list gag_linker_composition_fragments(self, glycan, long core_count=1, long iteration_count=0):
        cdef:
            long xyl_in_aggregate
            long xyl_count, hexose_count
            list core_shifts
            CountTable key_ct
        xyl_in_aggregate = glycan[_XYL]
        if xyl_in_aggregate == 0:
            xyl_in_aggregate = glycan.query("Pen", exact=False)
        core_shifts = []
        for xyl_count in range(max(0, xyl_in_aggregate) + 1):
            if xyl_count == 0:
                shift = GlycanCompositionFragment._create(
                    0,
                    CComposition(),
                    CountTable._create(),
                    False
                )
                core_shifts.append(shift)
            else:
                key_ct = CountTable._create()
                key_ct.setitem(_XYL, xyl_count)
                shift = GlycanCompositionFragment._create(
                    _XYL_mass * xyl_count,
                    _XYL_composition * xyl_count,
                    key_ct,
                    False
                )
                core_shifts.append(shift)
            if xyl_count > 0:
                # TODO: Handle modified Hexose residues here too.
                for hexose_count in range(1, 3):
                    key_ct = CountTable._create()
                    key_ct.setitem(_XYL, xyl_count)
                    key_ct.setitem(_HEX, hexose_count)
                    shift = GlycanCompositionFragment._create(
                        ((_XYL_mass * xyl_count) + (
                            _HEX_mass * hexose_count)),
                        (
                            (_XYL_composition * xyl_count) + (
                                _HEX_composition * hexose_count)),
                        key_ct,
                        False
                    )
                    core_shifts.append(shift)
                    if hexose_count == 2:
                        key_ct = CountTable._create()
                        key_ct.setitem(_XYL, xyl_count)
                        key_ct.setitem(_HEX, hexose_count)
                        key_ct.setitem(_AHEX, 1)
                        shift = GlycanCompositionFragment._create(
                            ((_XYL_mass * xyl_count) + (
                                _HEX_mass * hexose_count) + _AHEX_mass),
                            (
                                (_XYL_composition * xyl_count) + (
                                    _HEX_composition * hexose_count) + _AHEX_composition),
                            key_ct,
                            False
                        )
                        core_shifts.append(shift)
        return core_shifts

    cpdef list mixed_stub_fragments(self):

        cdef:
            int core_count
            int total_core_count
            list per_site_shifts
            GlycanCompositionFragment site
            CComposition base_composition, composition
            double base_mass
            size_t i
            size_t n_combos
            bint is_extended, is_glycosylated
            str name_key, name
            object glycosylation
            GlycanCompositionProxy glycan
            list accumulator
            list combos
            tuple result

        total_core_count = 0
        glycan = self.glycan_composition()
        self._use_query = self._guess_query_mode(glycan)
        per_site_shifts = []
        base_composition = None
        if self.compute_compositions:
            base_composition = self.peptide_composition()
            base_mass = base_composition.calc_mass()
        else:
            base_mass = self.peptide.get_peptide_backbone_mass()

        core_count = self.count_glycosylation_type(GlycosylationType_n_linked)
        for i in range(core_count):
            core_shifts = self.n_glycan_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)

        total_core_count += core_count
        core_count = self.count_glycosylation_type(GlycosylationType_o_linked)
        for i in range(core_count):
            core_shifts = self.o_glycan_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)

        total_core_count += core_count
        core_count = self.count_glycosylation_type(GlycosylationType_glycosaminoglycan)
        for i in range(core_count):
            core_shifts = self.gag_linker_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)

        core_count += total_core_count
        combos = self._combinate_sites(per_site_shifts, glycan)
        n_combos = PyList_Size(combos)
        accumulator = PyList_New(n_combos)
        for i in range(n_combos):
            result = <tuple>PyList_GetItem(combos, i)
            site = <GlycanCompositionFragment>PyTuple_GetItem(result, 0)
            name = <str>PyTuple_GetItem(result, 1)
            glycosylation = <object>PyTuple_GetItem(result, 2)
            is_glycosylated = (site.mass != 0)

            composition = None
            if self.compute_compositions:
                composition = base_composition.clone()
                composition.add_from(site.composition)

            frag = StubFragment._create(
                name=name,
                mass=base_mass + site.mass,
                composition=composition,
                is_glycosylated=is_glycosylated,
                kind=IonSeries_stub_glycopeptide,
                chemical_shift=None,
                glycosylation=glycosylation,
                is_extended=site.is_extended,
                glycosylation_size=site.get_glycosylation_size())
            Py_INCREF(frag)
            PyList_SET_ITEM(accumulator, i, frag)
        return accumulator

    cpdef list stub_fragments(self):
        n_glycan = self.count_glycosylation_type(GlycosylationType_n_linked) > 0
        o_glycan = self.count_glycosylation_type(GlycosylationType_o_linked) > 0
        gag_linker = self.count_glycosylation_type(GlycosylationType_glycosaminoglycan) > 0
        if (n_glycan + o_glycan + gag_linker) > 1:
            return self.mixed_stub_fragments()
        elif n_glycan:
            return self.n_glycan_stub_fragments()
        elif o_glycan:
            return self.o_glycan_stub_fragments()
        elif gag_linker:
            return self.gag_linker_stub_fragments()
        else:
            if self.glycosylation_manager.get_size() > 0:
                raise ValueError("Unknown Glycan Class Detected")
        return []