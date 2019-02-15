from six import string_types as basestring

from ..fragmentation_strategy import HCDFragmentationStrategy
from ..fragment import IonSeries

from .base import _PeptideSequenceCore
from .glycosylated_sequence import GlycosylatedSequenceMixin
from .mutable_sequence import MutableSequenceMixin


b_series = IonSeries.b
y_series = IonSeries.y
stub_glycopeptide_series = IonSeries.stub_glycopeptide


class PeptideSequence(_PeptideSequenceCore, GlycosylatedSequenceMixin, MutableSequenceMixin):
    @classmethod
    def from_iterable(cls, iterable, glycan_composition=None, n_term=None, c_term=None, text=None, **kwargs):
        seq = cls()
        if text is None:
            iterable = list(iterable)
            try:
                if isinstance(iterable[0][0], basestring):
                    text = True
            except IndexError:
                text = True
        if text:
            seq._init_from_parsed_string(iterable, glycan_composition, n_term, c_term)
        else:
            seq._init_from_components(iterable, glycan_composition, n_term, c_term, **kwargs)
        return seq

    def subsequence(self, slice_obj):
        sub = self[slice_obj]
        subseq = self.from_iterable(sub)
        if slice_obj.start == 0:
            subseq.n_term = self.n_term
        if slice_obj.stop == len(self):
            subseq.c_term = self.c_term
        return subseq

    def break_at(self, idx):
        if self._fragment_index is None:
            self._build_fragment_index()
        return self._fragment_index[idx]

    def _build_fragment_index(self, types=tuple('bycz')):
        self._fragment_index = [[] for i in range(len(self) + 1)]
        for series in types:
            series = IonSeries(series)
            if series.direction > 0:
                g = self.get_fragments(
                    series)
                for frags in g:
                    position = self._fragment_index[frags[0].position]
                    position.append(frags)
            else:
                g = self.get_fragments(
                    series)
                for frags in g:
                    position = self._fragment_index[
                        len(self) - frags[0].position]
                    position.append(frags)

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
    seq = Sequence.from_iterable(flat_chunks) if wrap else flat_chunks
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
