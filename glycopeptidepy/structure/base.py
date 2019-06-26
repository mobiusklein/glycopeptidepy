from glypy.utils import make_struct


class MoleculeBase(object):
    __slots__ = ()
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
    __slots__ = ()

    def serialize(self):
        '''A string representation for inclusion in sequences'''
        return self.name


class ResidueBase(MoleculeBase):
    '''
    A base type for classes describing amino acid residues
    '''
    __slots__ = ()


class SequencePosition(make_struct('SequencePosition', ['amino_acid', 'modifications'])):
    __slots__ = ()

    def __init__(self, parts):
        super(SequencePosition, self).__init__(*parts)

    def __repr__(self):
        return repr(list(self))

    @property
    def mass(self):
        mass = self.amino_acid.mass
        for mod in self.modifications:
            mass += mod.mass
        return mass


try:
    from glycopeptidepy._c.structure.base import (
        AminoAcidResidueBase as ResidueBase,
        ModificationBase, PeptideSequenceBase,
        SequencePosition)
except ImportError:
    pass
