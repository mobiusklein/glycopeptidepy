from glypy.utils import make_struct


class MoleculeBase(object):
    mass = None

    def __copy__(self):
        return self.clone()


# A few base types for doing type-based behavior changes
class PeptideSequenceBase(MoleculeBase):
    '''
    A base type for classes describing peptide sequences, with or without modifiations
    '''

    def _invalidate(self):
        pass


class ModificationBase(MoleculeBase):
    '''
    A base type for classes describing peptide sequence modifications
    '''

    def serialize(self):
        '''A string representation for inclusion in sequences'''
        return self.name


class ResidueBase(MoleculeBase):
    '''
    A base type for classes describing amino acid residues
    '''


class SequencePosition(make_struct('SequencePosition', ['amino_acid', 'modifications'])):
    __slots__ = ()

    def __init__(self, parts):
        return super(SequencePosition, self).__init__(*parts)

    def __repr__(self):
        return repr(list(self))


try:
    from glycopeptidepy._c.structure.base import (
        AminoAcidResidueBase as ResidueBase,
        ModificationBase, PeptideSequenceBase,
        SequencePosition)
except ImportError:
    pass
