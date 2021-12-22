# pylint: disable=assigning-non-slot
from ..base import SequencePosition
from ..modification import ModificationCategory, Modification, SequenceLocation
from ..terminal_group import _make_terminal_group
from .. import constants as structure_constants


class MutableSequenceMixin(object):
    __slots__ = ()

    def drop_modification(self, position, modification_type):
        '''
        Drop a modification by name from a specific residue. If the
        position is the N-term or the C-term, the terminal modification will
        be reset to the default.

        Parameters
        ----------
        position: int
            The position of the modification to drop
        modification_type: str or Modification
            The modification to drop

        Raises
        ------
        ValueError:
            If the `modification_type` is not found at `position`, a :class:`ValueError` will
            be raised.
        '''
        dropped_index = None
        self._invalidate()
        if position is SequenceLocation.n_term:
            self.n_term = _make_terminal_group(structure_constants.N_TERM_DEFAULT)
            return
        elif position is SequenceLocation.c_term:
            self.c_term = _make_terminal_group(structure_constants.C_TERM_DEFAULT)
            return

        for i, mod in enumerate(self.sequence[position].modifications):
            if modification_type == mod.rule:
                dropped_index = i
                break
        if dropped_index is None:
            raise ValueError("Modification not found! %s @ %s" %
                             (modification_type, position))
        try:
            drop_mod = self.sequence[position].modifications.pop(dropped_index)
            self.mass -= drop_mod.mass
            if drop_mod.is_tracked_for(ModificationCategory.glycosylation):
                self._glycosylation_manager.pop(position)
        except (IndexError, ValueError):
            raise ValueError("Modification not found! %s @ %s" %
                             (modification_type, position))

    def add_modification(self, position, modification_type):
        '''
        Add a modification by name to a specific residue. If the
        position is the N-term or the C-term, the terminal modification will
        be replaced.

        Parameters
        ----------
        position: int
            The position of the modification to add
        modification_type: str or Modification
            The modification to add
        '''
        self._invalidate()
        if isinstance(modification_type, Modification):
            mod = modification_type
        else:
            mod = Modification(rule=modification_type)

        if position is SequenceLocation.n_term:
            self.n_term = self.n_term.modify(mod)
        elif position is SequenceLocation.c_term:
            self.c_term = self.c_term.modify(mod)
        else:
            if (position == -1) or (position >= len(self.sequence)):
                raise IndexError(
                    "Invalid modification position. %s, %s, %s" %
                    (position, str(self.sequence), modification_type))

            self.sequence[position][1].append(mod)
            self.mass += mod.mass
            if mod.is_tracked_for(ModificationCategory.glycosylation):
                self._glycosylation_manager[position] = mod

    def _retrack_sequence(self):
        self._glycosylation_manager.clear()
        for i, position in enumerate(self):
            if position.modifications:
                for mod in position.modifications:
                    if mod.is_tracked_for(ModificationCategory.glycosylation):
                        self._glycosylation_manager[i] = mod

    def insert(self, position, residue, modifications=None):
        if modifications is None:
            modifications = []
        self._invalidate()
        self.sequence.insert(
            position, SequencePosition([residue, modifications]))
        self.mass += residue.mass
        for mod in modifications:
            self.mass += mod.mass
        self._retrack_sequence()

    def delete(self, position):
        self._invalidate()
        residue, mods = self.sequence.pop(position)
        self.mass -= residue.mass
        if mods:
            for mod in mods:
                self.mass -= mod.mass
        self._retrack_sequence()

    def substitute(self, position, residue):
        old_residue = self.sequence[position].amino_acid
        self.mass -= old_residue.mass
        self.mass += residue.mass
        self.sequence[position].amino_acid = residue
        self._invalidate()
        self._retrack_sequence()

    def append(self, residue, modification=None):
        self._invalidate()
        self.mass += residue.mass
        next_pos = [residue]
        if modification is None:
            next_pos.append([])
        else:
            next_pos.append([modification])
            self.mass += modification.mass
        self.sequence.append(SequencePosition(next_pos))
        self._retrack_sequence()

    def extend(self, sequence):
        self._invalidate()
        self.sequence.extend(sequence.sequence)
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        self._retrack_sequence()

    def leading_extend(self, sequence):
        self._invalidate()
        self.sequence = sequence.sequence + self.sequence
        self.mass += sequence.mass - sequence.n_term.mass - sequence.c_term.mass
        self._retrack_sequence()


try:
    from glycopeptidepy._c.structure.sequence_methods import add_modification as _c_add_modification

    MutableSequenceMixin.add_modification = _c_add_modification
except ImportError:
    pass
