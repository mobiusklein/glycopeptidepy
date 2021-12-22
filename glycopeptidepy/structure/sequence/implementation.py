from six import string_types as basestring

from ..fragmentation_strategy import HCDFragmentationStrategy
from ..fragment import IonSeries
from ..base import SequencePosition
from ..modification import SequenceLocation

from .base import _PeptideSequenceCore
from .glycosylated_sequence import GlycosylatedSequenceMixin
from .mutable_sequence import MutableSequenceMixin


b_series = IonSeries.b
y_series = IonSeries.y
stub_glycopeptide_series = IonSeries.stub_glycopeptide


class PeptideSequence(_PeptideSequenceCore, GlycosylatedSequenceMixin, MutableSequenceMixin):
    @classmethod
    def from_iterable(cls, iterable, glycan_composition=None, n_term=None, c_term=None, text=None, **kwargs):
        """Construct a :class:`PeptideSequence` instance from an arbitrary iterable of (amino acid, [modifications...])
        pairs.

        Parameters
        ----------
        iterable: :class:`~.Iterable`
            An arbitrary :class:`~.Iterable` of two-element pairs where the first element denotes
            the amino acid at that position while the second is a list of modifications. These elements
            may be denoted by actual objects of the appropriate type or as strings (so long as all are are
            typed object or strings) and will be coerced accordingly.
        glycan_composition: :class:`str` or :class:`~.GlycanComposition`
            An optional glycan composition given either as a string or a typed object. Will be coerced
            as needed.
        n_term: :class:`~.TerminalGroup`
            The amino terminal group
        c_term: :class:`~.TerminalGroup`
            The carboxyl terminal group
        text: :class:`bool`
            Whether or not to enforce text parsing. Defaults to :const:`None` which results in automatic
            detection.

        Returns
        -------
        :class:`PeptideSequence`
        """
        seq = cls()
        if not isinstance(iterable, list):
            iterable = list(iterable)
        if text is None:
            try:
                if isinstance(iterable[0][0], basestring):
                    text = True
            except IndexError:
                text = True
        if text:
            seq._init_from_parsed_string(iterable, glycan_composition, n_term, c_term)
        else:
            try:
                if not isinstance(iterable[0], SequencePosition):
                    # Coerce all pairs in the input to :class:`SequencePosition` because :meth:`_init_from_components`
                    # assumes that its input is a strictly typed list.
                    iterable = list(map(SequencePosition, iterable))
            except (ValueError, AttributeError, TypeError):
                iterable = list(map(SequencePosition, iterable))
            seq._init_from_components(
                iterable,
                glycan_composition, n_term, c_term, **kwargs)
        return seq

    def modified_residues(self):
        """Find all the modified residue positions in the sequence

        Returns
        -------
        list of :class:`~.SequencePosition`
        """
        return [(i, position) for i, position in enumerate(self) if position.modifications]

    def has_modification(self, position, modification_type):
        '''Check if a modification is present on the sequence

        Parameters
        ----------
        position: int
            The position in the sequence to check for the modification
        modification_type: str or Modification
            The modification to check for.

        Returns
        -------
        bool:
            Whether or not `modification_typ` is present at `position`
        '''
        if position is SequenceLocation.n_term:
            return self.n_term == modification_type
        elif position is SequenceLocation.c_term:
            return self.c_term == modification_type
        return self.sequence[position].has_modification(modification_type)


    def strip_modifications(self):
        """Return a copy of this sequence with all modifications removed.

        Returns
        -------
        :class:`PeptideSequence`
        """
        parts = [[pos.amino_acid, []] for pos in self]
        return self.from_iterable(parts)

    def subsequence(self, slice_obj):
        sub = self[slice_obj]
        subseq = self.from_iterable(sub)
        if slice_obj.start == 0:
            subseq.n_term = self.n_term
        if slice_obj.stop == len(self):
            subseq.c_term = self.c_term
        return subseq

    def reverse(self):
        slc = slice(None, None, -1)
        seq = self.subsequence(slc)
        if seq.glycan_composition != self.glycan_composition:
            if seq.glycan_composition:
                raise ValueError("Cannot recalculate partially ambiguous glycan composition")
            else:
                seq.glycan = self.glycan_composition.clone()
        return seq

    def get_fragments(self, kind, chemical_shifts=None, strategy=None, include_neutral_losses=False, **kwargs):
        """Generate fragments from this structure according to the strategy specified
        by ``strategy``, returning an iterator over the sequence of theoretical fragments.

        There may be multiple fragments for a single position, resulting in the iterator
        yielding lists of fragment objects per step.

        Parameters
        ----------
        kind : :class:`~.IonSeries` or :class:`str`
            The name of the ion series to produce fragments from, either as a string or the
            :class:`~.IonSeries` object to use
        chemical_shifts : dict, optional
            A :class:`~.Mapping` between :class:`~.AminoAcidResidue` and a listof acceptable chemical shifts to apply
            to the produced fragments containing that :class:`~.AminoAcidResidue`, to be applied combinatorially.
        strategy : :class:`~.FragmentationStrategyBase` type, optional
            The strategy type to employ when producing fragments. Defaults to :class:`~.HCDFragmentationStrategy`.
        include_neutral_losses: :class:`bool`
            Whether to include generic neutral losses (-NH3 and -H2O) on all fragments
        **kwargs
            Passed to ``strategy``

        Returns
        -------
        :class:`~.FragmentationStrategyBase`
            The fragmentation iterator
        """
        if strategy is None:
            strategy = HCDFragmentationStrategy
        losses = kwargs.pop('neutral_losses', {})
        if losses and not chemical_shifts:
            chemical_shifts = losses
        if kind == stub_glycopeptide_series:
            return ([frag] for frag in self.stub_fragments(True))
        else:
            return strategy(
                self, kind, chemical_shifts=chemical_shifts,
                include_neutral_losses=include_neutral_losses, **kwargs)


def list_to_sequence(seq_list, wrap=True):
    flat_chunks = []
    for chunk in seq_list:
        if(isinstance(chunk[0], list)):
            flat_chunks.extend(chunk)
        else:
            flat_chunks.append(chunk)
    seq = PeptideSequence.from_iterable(flat_chunks) if wrap else flat_chunks
    return seq


Sequence = PeptideSequence
parse = Sequence


class NamedSequence(PeptideSequence):

    def __init__(self, name=None, sequence=None, parser_function=None, **kwargs):
        super(NamedSequence, self).__init__(
            sequence, parser_function, **kwargs)
        self.name = name

    def clone(self):
        dup = super(NamedSequence, self).clone()
        dup.name = self.name
        return dup

    def __repr__(self):
        string = super(NamedSequence, self).__str__()
        name = str(self.name)
        if name.startswith(">"):
            name = name[1:]
        return ">%s\n%s" % (name, string)


class AnnotatedSequence(NamedSequence):

    def __init__(self, name=None, sequence=None, parser_function=None, annotations=None,
                 **kwargs):
        super(AnnotatedSequence, self).__init__(
            name, sequence, parser_function=parser_function, **kwargs)
        self.annotations = self._prepare_annotations(annotations)

    @staticmethod
    def _prepare_annotations(annotation_data):
        if annotation_data:
            return dict(annotation_data)
        else:
            return dict()

    def clone(self):
        dup = super(AnnotatedSequence, self).clone()
        dup.annotations = self._prepare_annotations(self.annotations)
        return dup


class ProteinSequence(AnnotatedSequence):
    pass
