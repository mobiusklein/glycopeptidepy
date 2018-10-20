# -*- coding: utf-8 -*-
from collections import defaultdict
from itertools import product, combinations

from glypy.structure.glycan_composition import (
    FrozenGlycanComposition, FrozenMonosaccharideResidue,
    MonosaccharideResidue)

from glycopeptidepy.utils.collectiontools import descending_combination_counter
from glycopeptidepy.structure import constants as structure_constants
from glycopeptidepy.structure.composition import (
    Composition)
from glycopeptidepy.structure.modification import (
    ModificationCategory,
    Modification,
    NGlycanCoreGlycosylation,
    OGlycanCoreGlycosylation,
    GlycosaminoglycanLinkerGlycosylation,
    GlycanFragment,
    ModificationIndex)
from glycopeptidepy.structure.glycan import (
    GlycosylationManager,
    GlycosylationType,
    HashableGlycanComposition)
from glycopeptidepy.structure.fragment import (
    IonSeries,
    SimpleFragment,
    StubFragment,
    PeptideFragment,
    format_negative_composition,
    ChemicalShift)
from glycopeptidepy.structure.residue import (
    AminoAcidResidue as R,
    residue_to_neutral_loss as _residue_to_neutral_loss)
from ..utils.collectiontools import _AccumulatorBag


_n_glycosylation = NGlycanCoreGlycosylation()
_o_glycosylation = OGlycanCoreGlycosylation()
_gag_linker_glycosylation = GlycosaminoglycanLinkerGlycosylation()
_modification_hexnac = Modification("HexNAc").rule
_modification_xylose = Modification("Xyl").rule


glycosylation_type_to_core = {
    GlycosylationType.n_linked: _n_glycosylation,
    GlycosylationType.o_linked: _o_glycosylation,
    GlycosylationType.glycosaminoglycan: _gag_linker_glycosylation
}


class FragmentationStrategyBase(object):
    def __init__(self, peptide):
        self.peptide = peptide

    def __repr__(self):
        return "%s(%s)" % (
            self.__class__.__name__,
            self.peptide)


class CADFragmentationStrategy(FragmentationStrategyBase):
    def __init__(self, peptide, max_cleavages=2):
        super(CADFragmentationStrategy, self).__init__(peptide)
        self.max_cleavages = max_cleavages
        self.iterator = self._generator(self.max_cleavages)

    def _generator(self, max_cleavages=2):
        glycans = self.peptide.glycosylation_manager.items()
        n = len(glycans)
        base_composition = self.peptide.peptide_composition()
        for i in range(1, n + 1):
            for dissociated in combinations(glycans, i):
                remaining = dict(glycans)
                for position, glycan in dissociated:
                    remaining.pop(position)
                reference_composition = base_composition.copy()
                for position, glycan in remaining.items():
                    reference_composition += glycan.rule.total_composition()

                glycan_B_ions = defaultdict(list)
                glycan_Y_ions = defaultdict(list)
                for position, glycan in dissociated:
                    for f in glycan.rule.get_fragments("BY", max_cleavages=max_cleavages):
                        if 'B' in f.series:
                            glycan_B_ions[position].append(f)
                        else:
                            glycan_Y_ions[position].append(f)
                for b_ion_set in glycan_B_ions.values():
                    for f in b_ion_set:
                        yield f
                key_order = list(glycan_Y_ions.keys())
                for y_ion_set in product(*glycan_Y_ions.values()):
                    name = ",".join("%d:%s" % (k, f.name) for k, f in zip(key_order, y_ion_set))
                    fragment_composition = reference_composition.copy()
                    is_glycosylated = False
                    for f in y_ion_set:
                        fragment_composition += f.composition
                        is_glycosylated = True
                    full_name = "peptide"
                    if name:
                        full_name = "%s+%s" % (full_name, name)
                    f = SimpleFragment(
                        name=full_name,
                        mass=fragment_composition.mass,
                        composition=fragment_composition,
                        is_glycosylated=is_glycosylated,
                        kind=IonSeries.stub_glycopeptide)
                    yield f

    def __next__(self):
        return next(self.iterator)

    def next(self):
        return self.__next__()

    def __iter__(self):
        return self


class StubGlycopeptideStrategy(FragmentationStrategyBase):

    def __init__(self, peptide, extended=True):
        super(StubGlycopeptideStrategy, self).__init__(peptide)
        self.extended = extended
        self._hexose = None
        self._hexnac = None
        self._neuac = None
        self._neugc = None
        self._fucose = None
        self._dhex = None
        self._xylose = None
        self._hexa = None
        self._generator = None
        # to be initialized closer to when the glycan composition will
        # be used. Determine whether to use the fast-path __getitem__ or
        # go through the slower but more general GlycanComposition.query
        # method to look up monosaccharides.
        self._use_query = False

    def _guess_query_mode(self, glycan_composition):
        # these guesses will work for N-glycans and common types of mucin-type O-glycans
        # and GAG linkers
        flag = ("Hex2NAc" in glycan_composition) or ("Glc2NAc" in glycan_composition
                                                     ) or ("Gal2NAc" in glycan_composition)
        return flag

    def _prepare_monosaccharide(self, name):
        return FrozenMonosaccharideResidue.from_iupac_lite(name)

    @property
    def hexose(self):
        if self._hexose is None:
            self._hexose = self._prepare_monosaccharide("Hex")
        return self._hexose

    @property
    def hexnac(self):
        if self._hexnac is None:
            self._hexnac = self._prepare_monosaccharide("HexNAc")
        return self._hexnac

    @property
    def fucose(self):
        if self._fucose is None:
            self._fucose = self._prepare_monosaccharide("Fuc")
        return self._fucose

    @property
    def dhex(self):
        if self._dhex is None:
            self._dhex = self._prepare_monosaccharide("dHex")
        return self._dhex

    @property
    def neuac(self):
        if self._neuac is None:
            self._neuac = self._prepare_monosaccharide("Neu5Ac")
        return self._neuac

    @property
    def neugc(self):
        if self._neugc is None:
            self._neugc = self._prepare_monosaccharide("Neu5Gc")
        return self._neugc

    @property
    def xylose(self):
        if self._xylose is None:
            self._xylose = self._prepare_monosaccharide("Xyl")
        return self._xylose

    @property
    def hexa(self):
        if self._hexa is None:
            self._hexa = self._prepare_monosaccharide("HexA")
        return self._hexa

    def count_glycosylation_type(self, glycotype):
        return self.peptide.glycosylation_manager.count_glycosylation_type(glycotype)

    def glycan_composition(self):
        return self.peptide.glycan_composition

    def peptide_composition(self):
        return self.peptide.peptide_composition()

    def fucosylate_increment(self, shift):
        fucosylated = shift.copy()
        fucosylated['key'] = fucosylated['key'].copy()
        fucosylated['mass'] += self.fucose.mass()
        fucosylated['composition'] = fucosylated[
            'composition'] + self.fucose.total_composition()
        fucosylated['key']["Fuc"] = 1
        return fucosylated

    def xylosylate_increment(self, shift):
        xylosylated = shift.copy()
        xylosylated['key'] = xylosylated['key'].copy()
        xylosylated['mass'] += self.xylose.mass()
        xylosylated['composition'] = xylosylated[
            'composition'] + self.xylose.total_composition()
        xylosylated['key']["Fuc"] = 1
        return xylosylated

    def n_glycan_stub_fragments(self):
        glycan = self.glycan_composition()
        self._use_query = self._guess_query_mode(glycan)
        core_count = self.count_glycosylation_type(GlycosylationType.n_linked)
        per_site_shifts = []
        base_composition = self.peptide_composition()
        base_mass = base_composition.mass

        for i in range(core_count):
            core_shifts = self.n_glycan_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)
        for positions in product(*per_site_shifts):
            key_base = 'peptide'
            names = _AccumulatorBag()
            mass = base_mass
            composition = base_composition.clone()
            is_extended = False
            for site in positions:
                mass += site['mass']
                names += (site['key'])
                composition += site['composition']
                is_extended |= site['is_extended']
            is_glycosylated = (composition != base_composition)
            invalid = False
            if self._use_query:
                for key, value in names.items():
                    if glycan.query(key) < value:
                        invalid = True
                        break
            else:
                for key, value in names.items():
                    if glycan[key] < value:
                        invalid = True
                        break
            if invalid:
                continue
            extended_key = ''.join("%s%d" % kv for kv in sorted(names.items()))
            if len(extended_key) > 0:
                key_base = "%s+%s" % (key_base, extended_key)
            glycosylation = HashableGlycanComposition(names)
            yield StubFragment(
                name=key_base,
                mass=mass,
                composition=composition,
                is_glycosylated=is_glycosylated,
                kind=IonSeries.stub_glycopeptide,
                glycosylation=glycosylation,
                is_extended=is_extended)

    def n_glycan_composition_fragments(self, glycan, core_count=1, iteration_count=1):
        hexose = self.hexose
        hexnac = self.hexnac

        if self._use_query:
            fucose_count = glycan.query('Fuc') + glycan.query('dHex')
            xylose_count = glycan.query('Xyl')
            hexnac_in_aggregate = glycan.query('HexNAc')
            hexose_in_aggregate = glycan.query('Hex')
        else:
            fucose_count = glycan['Fuc'] + glycan['dHex']
            xylose_count = glycan['Xyl']
            hexnac_in_aggregate = glycan['HexNAc']
            hexose_in_aggregate = glycan["Hex"]

        core_shifts = []
        base_hexnac = min(hexnac_in_aggregate + 1, 3)
        print(base_hexnac)
        for hexnac_count in range(base_hexnac):
            if hexnac_count == 0:
                shift = {
                    "mass": 0,
                    "composition": Composition(),
                    "key": {},
                    'is_extended': False,
                }
                core_shifts.append(shift)
            elif hexnac_count == 1:
                shift = {
                    "mass": (hexnac_count * hexnac.mass()),
                    "composition": hexnac_count * hexnac.total_composition(),
                    "key": {"HexNAc": hexnac_count},
                    'is_extended': False,
                }
                core_shifts.append(shift)
                if iteration_count < fucose_count:
                    fucosylated = self.fucosylate_increment(shift)
                    core_shifts.append(fucosylated)
            elif hexnac_count == 2:
                shift = {
                    "mass": (hexnac_count * hexnac.mass()),
                    "composition": hexnac_count * hexnac.total_composition(),
                    "key": {"HexNAc": hexnac_count},
                    'is_extended': False,
                }
                core_shifts.append(shift)

                if iteration_count < fucose_count:
                    fucosylated = self.fucosylate_increment(shift)
                    core_shifts.append(fucosylated)
                    if iteration_count < xylose_count:
                        xylosylated = self.xylosylate_increment(fucosylated)
                        core_shifts.append(xylosylated)
                if iteration_count < xylose_count:
                    xylosylated = self.xylosylate_increment(shift)
                    core_shifts.append(xylosylated)

                for hexose_count in range(1, min(hexose_in_aggregate + 1, 4)):
                    shift = {
                        "mass": (hexnac_count * hexnac.mass()) + (hexose_count * hexose.mass()),
                        "composition": (hexnac_count * hexnac.total_composition()) + (
                            hexose_count * hexose.total_composition()),
                        "key": {"HexNAc": hexnac_count, "Hex": hexose_count},
                        'is_extended': False,
                    }
                    core_shifts.append(shift)
                    if iteration_count < fucose_count:
                        fucosylated = self.fucosylate_increment(shift)
                        core_shifts.append(fucosylated)
                        if iteration_count < xylose_count:
                            xylosylated = self.xylosylate_increment(fucosylated)
                            core_shifts.append(xylosylated)
                    if iteration_count < xylose_count:
                        xylosylated = self.xylosylate_increment(shift)
                        core_shifts.append(xylosylated)

                    # After the core motif has been exhausted, speculatively add
                    # on the remaining core monosaccharides sequentially until
                    # exhausted.
                    if hexose_count == 3 and hexnac_in_aggregate >= 2 * core_count and self.extended:
                        for extra_hexnac_count in range(0, 3):
                            if extra_hexnac_count + hexnac_count > hexnac_in_aggregate:
                                continue
                            shift = {
                                "mass": (
                                    (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                    hexose_count * hexose.mass()),
                                "composition": (
                                    (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                    hexose_count * hexose.total_composition()),
                                "key": {"HexNAc": hexnac_count + extra_hexnac_count, "Hex": hexose_count},
                                'is_extended': True
                            }
                            core_shifts.append(shift)
                            if iteration_count < fucose_count:
                                fucosylated = self.fucosylate_increment(shift)
                                core_shifts.append(fucosylated)
                                if iteration_count < xylose_count:
                                    xylosylated = self.xylosylate_increment(fucosylated)
                                    core_shifts.append(xylosylated)
                            if iteration_count < xylose_count:
                                xylosylated = self.xylosylate_increment(shift)
                                core_shifts.append(xylosylated)
                            for extra_hexose_count in range(1, 3):
                                if extra_hexose_count + hexose_count > hexose_in_aggregate:
                                    continue
                                shift = {
                                    "mass": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                        (hexose_count + extra_hexose_count) * hexose.mass()),
                                    "composition": (
                                        (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                        (hexose_count + extra_hexose_count) * hexose.total_composition()),
                                    "key": {"HexNAc": hexnac_count + extra_hexnac_count, "Hex": (
                                        hexose_count + extra_hexose_count)},
                                    'is_extended': True
                                }
                                core_shifts.append(shift)
                                if iteration_count < fucose_count:
                                    fucosylated = self.fucosylate_increment(shift)
                                    core_shifts.append(fucosylated)
                                    if iteration_count < xylose_count:
                                        xylosylated = self.xylosylate_increment(fucosylated)
                                        core_shifts.append(xylosylated)
                                if iteration_count < xylose_count:
                                    xylosylated = self.xylosylate_increment(shift)
                                    core_shifts.append(xylosylated)
        return core_shifts

    def o_glycan_stub_fragments(self):
        glycan = self.glycan_composition()
        self._use_query = self._guess_query_mode(glycan)
        core_count = self.count_glycosylation_type(GlycosylationType.o_linked)
        per_site_shifts = []

        base_composition = self.peptide_composition()
        base_mass = base_composition.mass

        for i in range(core_count):
            core_shifts = self.o_glycan_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)
        seen = set()
        for positions in product(*per_site_shifts):
            key_base = 'peptide'
            names = _AccumulatorBag()
            mass = base_mass
            composition = base_composition.clone()
            for site in positions:
                mass += site['mass']
                names += (site['key'])
                composition += site['composition']
            is_glycosylated = (composition != base_composition)
            invalid = False
            if self._use_query:
                for key, value in names.items():
                    if glycan.query(key) < value:
                        invalid = True
                        break
            else:
                for key, value in names.items():
                    if glycan[key] < value:
                        invalid = True
                        break
            if invalid:
                continue
            extended_key = ''.join(
                "%s%d" % kv for kv in sorted(names.items()) if kv[1] > 0)
            if len(extended_key) > 0:
                key_base = "%s+%s" % (key_base, extended_key)
            if key_base in seen:
                continue
            seen.add(key_base)
            glycosylation = HashableGlycanComposition(names)
            yield StubFragment(
                name=key_base,
                mass=mass,
                composition=composition,
                is_glycosylated=is_glycosylated,
                kind=IonSeries.stub_glycopeptide,
                glycosylation=glycosylation)

    def o_glycan_composition_fragments(self, glycan, core_count=1, iteration_count=1):
        hexose = self.hexose
        hexnac = self.hexnac

        if self._use_query:
            fucose_count = glycan.query('Fuc') + glycan.query('dHex')
            hexnac_in_aggregate = glycan.query('HexNAc')
            hexose_in_aggregate = glycan.query('Hex')
        else:
            fucose_count = glycan['Fuc'] + glycan['dHex']
            hexnac_in_aggregate = glycan['HexNAc']
            hexose_in_aggregate = glycan["Hex"]
        core_shifts = []
        for hexnac_count in range(3):
            if hexnac_in_aggregate < hexnac_count:
                continue
            if hexnac_count == 0:
                shift = {
                    "mass": 0,
                    "composition": Composition(),
                    "key": {}
                }
                core_shifts.append(shift)
            elif hexnac_count >= 1:
                shift = {
                    "mass": (hexnac_count * hexnac.mass()),
                    "composition": hexnac_count * hexnac.total_composition(),
                    "key": {"HexNAc": hexnac_count}
                }
                core_shifts.append(shift)
                if iteration_count < fucose_count:
                    fucosylated = self.fucosylate_increment(shift)
                    core_shifts.append(fucosylated)
                for hexose_count in range(0, 2):
                    if hexose_in_aggregate < hexose_count:
                        continue
                    shift = {
                        "mass": (
                            (hexnac_count) * hexnac.mass()) + (
                            (hexose_count) * hexose.mass()),
                        "composition": (
                            (hexnac_count) * hexnac.total_composition()) + (
                            (hexose_count) * hexose.total_composition()),
                        "key": {"HexNAc": hexnac_count, "Hex": (
                            hexose_count)}
                    }
                    if hexose_count > 0:
                        core_shifts.append(shift)
                        if iteration_count < fucose_count:
                            fucosylated = self.fucosylate_increment(shift)
                            core_shifts.append(fucosylated)
                    # After the core motif has been exhausted, speculatively add
                    # on the remaining core monosaccharides sequentially until
                    # exhausted.
                    if self.extended and hexnac_in_aggregate - hexnac_count >= 0:
                        for extra_hexnac_count in range(hexnac_in_aggregate - hexnac_count):
                            shift = {
                                "mass": (
                                    (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                    (hexose_count) * hexose.mass()),
                                "composition": (
                                    (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                    (hexose_count) * hexose.total_composition()),
                                "key": {"HexNAc": hexnac_count + extra_hexnac_count, "Hex": (
                                    hexose_count)}
                            }
                            if iteration_count < fucose_count:
                                fucosylated = self.fucosylate_increment(shift)
                                core_shifts.append(fucosylated)
                            if hexose_in_aggregate > hexose_count:
                                for extra_hexose_count in range(hexose_in_aggregate - hexose_count):
                                    shift = {
                                        "mass": (
                                            (hexnac_count + extra_hexnac_count) * hexnac.mass()) + (
                                            (hexose_count + extra_hexose_count) * hexose.mass()),
                                        "composition": (
                                            (hexnac_count + extra_hexnac_count) * hexnac.total_composition()) + (
                                            (hexose_count + extra_hexose_count) * hexose.total_composition()),
                                        "key": {"HexNAc": hexnac_count + extra_hexnac_count, "Hex": (
                                            hexose_count + extra_hexose_count)}
                                    }
                                    if iteration_count < fucose_count:
                                        fucosylated = self.fucosylate_increment(
                                            shift)
                                        core_shifts.append(fucosylated)

        return core_shifts

    def gag_linker_stub_fragments(self):
        glycan = self.glycan_composition()
        self._use_query = self._guess_query_mode(glycan)
        core_count = self.count_glycosylation_type(GlycosylationType.glycosaminoglycan)

        per_site_shifts = []

        base_composition = self.peptide_composition()
        base_mass = base_composition.mass
        for i in range(core_count):
            core_shifts = self.gag_linker_composition_fragments(glycan, core_count, i)
            per_site_shifts.append(core_shifts)
        seen = set()
        for positions in product(*per_site_shifts):
            key_base = 'peptide'
            names = _AccumulatorBag()
            mass = base_mass
            composition = base_composition.clone()
            for site in positions:
                mass += site['mass']
                names += (site['key'])
                composition += site['composition']
            is_glycosylated = (composition != base_composition)
            invalid = False
            if self._use_query:
                for key, value in names.items():
                    if glycan.query(key) < value:
                        invalid = True
                        break
            else:
                for key, value in names.items():
                    if glycan[key] < value:
                        invalid = True
                        break
            if invalid:
                continue
            extended_key = ''.join("%s%d" % kv for kv in sorted(names.items()))
            if len(extended_key) > 0:
                key_base = "%s+%s" % (key_base, extended_key)
            if key_base in seen:
                continue
            seen.add(key_base)
            glycosylation = HashableGlycanComposition(names)
            yield StubFragment(
                name=key_base,
                mass=mass,
                composition=composition,
                is_glycosylated=is_glycosylated,
                kind=IonSeries.stub_glycopeptide,
                glycosylation=glycosylation)

    def gag_linker_composition_fragments(self, glycan, core_count=1, iteration_count=1):
        hexose = self.hexose
        xyl = self.xylose
        hexa = self.hexa

        core_shifts = []
        for xyl_count in range(2):
            if xyl_count == 0:
                shift = {
                    "mass": 0,
                    "composition": Composition(),
                    "key": {}
                }
                core_shifts.append(shift)
            else:
                shift = {
                    "mass": xyl.mass() * xyl_count,
                    "composition": xyl.total_composition() * xyl_count,
                    "key": {
                        "Xyl": xyl_count
                    }
                }
                core_shifts.append(shift)
            if xyl_count > 0:
                # TODO: Handle modified Hexose residues here too.
                for hexose_count in range(1, 3):
                    shift = {
                        "mass": ((xyl.mass() * xyl_count) + (
                            hexose.mass() * hexose_count)),
                        "composition": (
                            (xyl.total_composition() * xyl_count) + (
                                hexose.total_composition() * hexose_count)),
                        "key": {
                            "Xyl": xyl_count,
                            "Hex": hexose_count
                        }
                    }
                    core_shifts.append(shift)
                    if hexose_count == 2:
                        shift = {
                            "mass": ((xyl.mass() * xyl_count) + (
                                hexose.mass() * hexose_count) + hexa.mass()),
                            "composition": (
                                (xyl.total_composition() * xyl_count) + (
                                    hexose.total_composition() * hexose_count) + hexa.total_composition()),
                            "key": {
                                "Xyl": xyl_count,
                                "Hex": hexose_count,
                                "aHex": 1
                            }
                        }
                        core_shifts.append(shift)

        return core_shifts

    def stub_fragments(self):
        n_glycan = self.count_glycosylation_type(GlycosylationType.n_linked) > 0
        o_glycan = self.count_glycosylation_type(GlycosylationType.o_linked) > 0
        gag_linker = self.count_glycosylation_type(GlycosylationType.glycosaminoglycan) > 0
        if (n_glycan + o_glycan + gag_linker) > 2:
            raise ValueError(
                "Does not support mixed-type glycan fragmentation (yet)")
        if n_glycan:
            return self.n_glycan_stub_fragments()
        elif o_glycan:
            return self.o_glycan_stub_fragments()
        elif gag_linker:
            return self.gag_linker_stub_fragments()
        else:
            if len(self.peptide.glycosylation_manager) > 0:
                raise ValueError("Unknown Glycan Class Detected")
        return (a for a in [])

    def __next__(self):
        if self._generator is None:
            self._generator = self.stub_fragments()
        return next(self._generator)

    def next(self):
        return self.__next__()

    def __iter__(self):
        return self


class PeptideFragmentationStrategyBase(FragmentationStrategyBase):

    def __init__(self, peptide, series, chemical_shifts=None, max_chemical_shifts=1):
        if chemical_shifts is None:
            chemical_shifts = dict()
        super(PeptideFragmentationStrategyBase, self).__init__(peptide)
        self.series = IonSeries(series)
        self.chemical_shift_rules = defaultdict(list, chemical_shifts)
        self.max_chemical_shifts = max_chemical_shifts

        # Null Values
        self.running_mass = 0
        self.running_composition = Composition()
        self.index = None
        self.size = len(self.peptide)
        self.modification_index = ModificationIndex()
        self.glycosylation_manager = GlycosylationManager(self.peptide)
        self.amino_acids_counter = _AccumulatorBag()

        self.running_mass += self.series.mass_shift
        self.running_composition += self.series.composition_shift
        self._initialize_start_terminal()

    @property
    def direction(self):
        return self.series.direction

    def _get_viable_chemical_shift_combinations(self):
        shifts = []
        if not self.chemical_shift_rules:
            return shifts
        for residue, count in self.amino_acids_counter.items():
            loss_composition = self.chemical_shift_rules[residue]
            if loss_composition:
                for i in range(count):
                    shifts.append(loss_composition + [Composition()])
        loss_composition_combinations = {}
        for comb in combinations(shifts, self.max_chemical_shifts):
            for prod in product(*comb):
                loss_composition = sum(prod, Composition())
                loss_composition_combinations[format_negative_composition(loss_composition)] = loss_composition
        return [ChemicalShift(k, v) for k, v in loss_composition_combinations.items() if v]

    def _initialize_start_terminal(self):
        if self.direction > 0:
            self.running_mass += self.peptide.n_term.mass
            self.running_composition += self.peptide.n_term.composition
            self.index = -1
        elif self.direction < 0:
            self.running_mass += self.peptide.c_term.mass
            self.running_composition += self.peptide.c_term.composition
            self.index = len(self.peptide)
        else:
            raise ValueError("Unknown direction %r" % (self.series.direction,))

    def composition_of(self, residue, modifications):
        composition = Composition(residue.composition)
        for mod in modifications:
            composition += mod.composition
        return composition

    def has_more(self):
        if self.direction > 0:
            return self.index < self.size - 2
        else:
            return self.index > 1

    def flanking_residues(self):
        residues = [self.peptide[self.index][0],
                    self.peptide[self.index + self.direction][0]]
        if self.direction < 0:
            residues = residues[::-1]
        return residues

    def name_index_of(self):
        if self.direction > 0:
            return self.index + structure_constants.FRAG_OFFSET
        else:
            return (self.size - self.index - 1) + structure_constants.FRAG_OFFSET

    def __next__(self):
        if self.has_more():
            return self.step()
        else:
            raise StopIteration()

    def next(self):
        return self.__next__()

    def __iter__(self):
        return self

    def track_glycosylation(self, index, glycosylation):
        self.glycosylation_manager[self.index] = glycosylation
        self.modification_index[glycosylation] += 1

    def _update_state(self):
        self.index += self.direction
        residue, modifications = self.peptide[self.index]
        for mod in modifications:
            if mod.is_tracked_for(ModificationCategory.glycosylation):
                self.track_glycosylation(self.index, mod)
            else:
                self.modification_index[mod] += 1
        self.amino_acids_counter[residue] += 1
        composition = self.composition_of(residue, modifications)
        self.running_mass += residue.mass
        self.running_composition += composition

    def _build_fragments(self):
        fragments_from_site = []
        frag = PeptideFragment(
            self.series,
            self.name_index_of(),
            dict(self.modification_index),
            self.running_mass,
            glycosylation=self.glycosylation_manager.copy(),
            flanking_amino_acids=self.flanking_residues(),
            composition=self.running_composition)
        shifts = self._get_viable_chemical_shift_combinations()
        for fragment in self.partial_loss(frag):
            fragments_from_site.append(fragment)
            for shift in shifts:
                f = fragment.clone()
                f.chemical_shift = shift
                fragments_from_site.append(f)
        return fragments_from_site

    def step(self):
        self._update_state()
        fragments_from_site = self._build_fragments()
        return fragments_from_site

    def __repr__(self):
        return "%s(%s, %r, %r, %0.4f, %d)" % (
            self.__class__.__name__,
            self.peptide, self.series, self.running_composition,
            self.running_mass, self.index)


class HCDFragmentationStrategy(PeptideFragmentationStrategyBase):
    modifications_of_interest = set([
        _n_glycosylation,
        _modification_hexnac,
        _o_glycosylation,
        _gag_linker_glycosylation,
        _modification_xylose
    ])

    modification_compositions = {
        k: k.composition for k in modifications_of_interest
    }

    def _get_core_for(self, glycosylation):
        try:
            if not glycosylation.rule.is_core:
                glycosylation = glycosylation_type_to_core[glycosylation.rule.glycosylation_type]
            return glycosylation
        except KeyError:
            raise ValueError("Cannot determine which core to use for {}".format(
                glycosylation.rule.glycosylation_type))

    def composition_of(self, residue, modifications):
        composition = Composition(residue.composition)
        for mod in modifications:
            if mod.is_tracked_for(ModificationCategory.glycosylation):
                mod = self._get_core_for(mod)
            composition += mod.composition
        return composition

    def track_glycosylation(self, index, glycosylation):
        # HCD strategy does not track intact topologies
        try:
            glycosylation = self._get_core_for(glycosylation)
        except KeyError:
            raise ValueError("Cannot determine which core to use for {}".format(
                glycosylation.rule.glycosylation_type))
        self.glycosylation_manager[self.index] = glycosylation
        self.modification_index[glycosylation] += 1

    def _get_modifications_of_interest(self, fragment):
        modifications = dict(fragment.modification_dict)
        modifications_of_interest = defaultdict(
            int, {k: v for k, v in modifications.items()
                  if k in self.modifications_of_interest})
        delta_composition = sum(
            (self.modification_compositions[k] * v
             for k, v in modifications_of_interest.items()),
            Composition())
        return modifications_of_interest, delta_composition

    def _replace_cores(self, modifications_of_interest):
        n_cores = modifications_of_interest.pop(_n_glycosylation, 0)
        o_cores = modifications_of_interest.pop(_o_glycosylation, 0)
        gag_cores = modifications_of_interest.pop(_gag_linker_glycosylation, 0)

        # Convert core glycosylation into the residual monosaccharides remaining
        # after dissociation
        modifications_of_interest[_modification_hexnac] += n_cores * 2 + o_cores
        modifications_of_interest[_modification_xylose] += gag_cores
        return modifications_of_interest

    def _generate_modification_variants(self, interesting_modifications, other_modifications):
        for varied_modifications in descending_combination_counter(interesting_modifications):
            updated_modifications = other_modifications.copy()
            updated_modifications.update(
                {k: v for k, v in varied_modifications.items() if v != 0})

            extra_composition = Composition()
            for mod, mod_count in varied_modifications.items():
                extra_composition += self.modification_compositions[mod] * mod_count
            yield updated_modifications, extra_composition

    def partial_loss(self, fragment):
        (modifications_of_interest,
         delta_composition) = self._get_modifications_of_interest(fragment)
        base_composition = fragment.composition - delta_composition
        self._replace_cores(modifications_of_interest)
        other_modifications = {
            k: v for k, v in fragment.modification_dict.items()
            if k not in self.modifications_of_interest
        }
        for updated_modifications, extra_composition in self._generate_modification_variants(
                modifications_of_interest, other_modifications):
            yield PeptideFragment(
                fragment.series, fragment.position, dict(updated_modifications), fragment.bare_mass,
                golden_pairs=fragment.golden_pairs,
                flanking_amino_acids=fragment.flanking_amino_acids,
                composition=base_composition + extra_composition)


# Cooper, H. J., Hudgins, R. R., Håkansson, K., & Marshall, A. G. (2002).
# Characterization of amino acid side chain losses in electron capture dissociation.
# Journal of the American Society for Mass Spectrometry, 13(3), 241–249.
# https://doi.org/10.1016/S1044-0305(01)00357-9
exd_sidechain_losses = defaultdict(list, {
    R("His"): [-Composition("H6C4N2")],
    R("Arg"): [-Composition("CH5N3"),
               -Composition("C4H11N3"),
               -Composition("C1N2H4"),
               -Composition("N2H6")],
    R("Asn"): [-Composition("CONH3")],
    R("Gln"): [-Composition("CONH3")],
    R("Met"): [-Composition("C3H6S")]
})

for key, value in _residue_to_neutral_loss.items():
    exd_sidechain_losses[R(key)].extend(value)


class EXDFragmentationStrategy(PeptideFragmentationStrategyBase):
    glycan_fragment_ladder = "Y"

    def __init__(self, peptide, series, neutral_losses=None, max_chemical_shifts=1, max_glycan_cleavages=2):
        super(EXDFragmentationStrategy, self).__init__(peptide, series, neutral_losses, max_chemical_shifts)
        self.max_glycan_cleavages = max_glycan_cleavages

    def _strip_glycosylation_from_modification_dict(self, modification_dict):
        glycosylations = []
        stripped_modifications = {}
        for position, glycosylation in sorted(self.glycosylation_manager.items()):
            glycosylations.append(glycosylation)
        for modification, count in modification_dict.items():
            if modification.is_tracked_for(ModificationCategory.glycosylation):
                continue
            stripped_modifications[modification] = count
        return glycosylations, stripped_modifications

    def partial_loss(self, fragment, glycan_fragment_ladder=None):
        if glycan_fragment_ladder is None:
            glycan_fragment_ladder = self.glycan_fragment_ladder
        # make a clone in case we mutate anything
        base = fragment.clone()

        (glycosylations,
         stripped_modifications) = self._strip_glycosylation_from_modification_dict(
            base.modification_dict)
        base_composition = base.composition
        for glycosylation in glycosylations:
            base_composition -= glycosylation.composition
        bare_fragment = PeptideFragment(
            fragment.series, fragment.position, dict(stripped_modifications),
            fragment.bare_mass,
            golden_pairs=fragment.golden_pairs,
            flanking_amino_acids=fragment.flanking_amino_acids,
            composition=base_composition.clone())

        fragments = [
            [None] + sorted(glycosylation.rule.get_fragments(
                glycan_fragment_ladder, max_cleavages=self.max_glycan_cleavages),
                key=lambda x: x.mass)
            for glycosylation in glycosylations if not glycosylation.rule.is_composition
        ]
        fragment_combinations = product(*fragments)
        for fragment_set in fragment_combinations:
            if fragment_set == ():
                continue
            delta_composition = Composition()
            new_modifications = ModificationIndex(stripped_modifications)
            for glycan_fragment in fragment_set:
                if glycan_fragment is None:
                    continue
                subfragment = GlycanFragment(glycan_fragment)
                new_modifications[subfragment] += 1
                delta_composition += subfragment.composition

            extended_fragment = PeptideFragment(
                bare_fragment.series, bare_fragment.position, new_modifications,
                bare_fragment.bare_mass,
                golden_pairs=bare_fragment.golden_pairs,
                flanking_amino_acids=bare_fragment.flanking_amino_acids,
                composition=bare_fragment.composition + delta_composition)
            yield extended_fragment
        # yield the intact fragment
        yield fragment
