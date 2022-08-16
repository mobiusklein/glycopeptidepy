from glypy.composition.ccomposition cimport CComposition

from glycopeptidepy._c.structure.fragment cimport IonSeriesBase
from glycopeptidepy._c.structure.sequence_methods cimport _PeptideSequenceCore
from glycopeptidepy._c.structure.base cimport ModificationBase

from glycopeptidepy.structure.modification import (
    Modification,
    NGlycanCoreGlycosylation,
    OGlycanCoreGlycosylation,
    GlycosaminoglycanLinkerGlycosylation)


from glycopeptidepy.structure.glycan import (GlycosylationType)


cdef:
    ModificationBase _n_glycosylation
    ModificationBase _o_glycosylation
    ModificationBase _gag_linker_glycosylation
    ModificationBase _modification_hexnac
    ModificationBase _modification_xylose

    dict glycosylation_type_to_core



cdef class FragmentationStrategyBase(object):
    cdef:
        public _PeptideSequenceCore peptide
        public bint compute_compositions

    cpdef CComposition peptide_composition(self)
    cpdef CComposition total_composition(self)
