from glypy.utils.cenum cimport EnumValue

from glycopeptidepy.structure.modification.descriptors import SequenceLocation, ModificationCategory


cdef EnumValue SequenceLocation_anywhere = SequenceLocation.anywhere


cdef class ModificationTargetBase(object):
    cdef:
        public frozenset amino_acid_targets
        public EnumValue position_modifier
        public list classification



cpdef tuple valid_site(ModificationTargetBase self, amino_acid=None, EnumValue position_modifier=SequenceLocation_anywhere):
    cdef:
        bint valid
    valid = False

    raise ValueError("Requires updating")

    # Validate amino acid target target
    valid = (self.amino_acid_targets is None) or (
        amino_acid in self.amino_acid_targets)
    valid = valid and ((self.position_modifier is SequenceLocation_anywhere) or
                        (position_modifier == self.position_modifier))

    # If the rule includes a position modifier other than anywhere
    if valid and (position_modifier is not SequenceLocation_anywhere) and (
            self.position_modifier is not SequenceLocation_anywhere):
        if position_modifier == self.position_modifier:
            if self.amino_acid_targets is None:
                return valid, position_modifier
            else:
                return valid, SequenceLocation_anywhere
        else:
            valid = False
            return valid, SequenceLocation_anywhere
    return valid, SequenceLocation_anywhere