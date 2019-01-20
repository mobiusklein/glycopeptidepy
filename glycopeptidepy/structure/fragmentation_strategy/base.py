from glycopeptidepy.structure.modification import (
    Modification,
    NGlycanCoreGlycosylation,
    OGlycanCoreGlycosylation,
    GlycosaminoglycanLinkerGlycosylation)


from glycopeptidepy.structure.glycan import (GlycosylationType)

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
    def __init__(self, peptide, *args, **kwargs):
        self.peptide = peptide
        super(FragmentationStrategyBase, self).__init__(*args, **kwargs)

    def __repr__(self):
        return "%s(%s)" % (
            self.__class__.__name__,
            self.peptide)

    def next(self):
        return self.__next__()

    def __iter__(self):
        return self

    def peptide_composition(self):
        return self.peptide.peptide_composition()

    def total_composition(self):
        return self.peptide.total_composition()
