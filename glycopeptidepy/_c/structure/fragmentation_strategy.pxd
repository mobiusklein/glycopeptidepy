from glycopeptidepy._c.structure.fragment cimport IonSeriesBase
from glycopeptidepy._c.structure.base cimport PeptideSequenceBase


cdef class FragmentationStrategyBase(object):
    cdef:
        public PeptideSequenceBase peptide
