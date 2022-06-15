'''A generic implementation of sequence features
'''

import json
import io

from typing import Dict, Iterable, List, Sequence, Tuple, TYPE_CHECKING

from six import add_metaclass

from glycopeptidepy.structure import Modification

from . import uniprot

if TYPE_CHECKING:
    from glycopeptidepy.io.fasta import FastaHeader


class AnnotationMeta(type):
    _cache = {}

    def __new__(cls, name, parents, attrs):
        new_type = type.__new__(cls, name, parents, attrs)
        try:
            cls._cache[name] = new_type
        except AttributeError:
            pass
        return new_type

    def _type_for_name(cls, feature_type):
        return cls._cache[feature_type]

    def from_dict(cls, d):
        name = d.pop('__class__')
        impl = cls._type_for_name(name)  # pylint: disable=no-value-for-parameter
        return impl(**d)


@add_metaclass(AnnotationMeta)
class AnnotationBase(object):
    feature_type: str
    description: str

    def __init__(self, feature_type, description, **kwargs):
        self.feature_type = feature_type
        self.description = description

    def to_dict(self):
        d = {}
        d['feature_type'] = self.feature_type
        d['description'] = self.description
        d['__class__'] = self.__class__.__name__
        return d

    def __eq__(self, other):
        if other is None:
            return False
        return self.feature_type == other.feature_type and self.description == other.description

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.feature_type, self.description))


class AnnotatedResidue(AnnotationBase):
    position: int

    def __init__(self, position, feature_type, description, **kwargs):
        super(AnnotatedResidue, self).__init__(feature_type, description, **kwargs)
        self.position = position

    @property
    def start(self):
        return self.position

    @property
    def end(self):
        return self.position + 1

    def to_dict(self):
        d = super(AnnotatedResidue, self).to_dict()
        d['position'] = self.position
        return d

    def __eq__(self, other):
        res = super(AnnotatedResidue, self).__eq__(other)
        if not res:
            return res
        return self.position == other.position

    def __hash__(self):
        return hash((self.feature_type, self.description, self.position))


class AnnotatedInterval(AnnotationBase):
    start: int
    end: int

    def __init__(self, start, end, feature_type, description, **kwargs):
        super(AnnotatedInterval, self).__init__(
            feature_type, description, **kwargs)
        self.start = start
        self.end = end

    def to_dict(self):
        d = super(AnnotatedInterval, self).to_dict()
        d['start'] = self.start
        d['end'] = self.end
        return d

    def __eq__(self, other):
        res = super(AnnotatedInterval, self).__eq__(other)
        if not res:
            return res
        return self.start == other.start and self.end == other.end

    def __hash__(self):
        return hash((self.feature_type, self.description, self.start, self.end))


class Domain(AnnotatedInterval):
    feature_type = 'domain'

    def __init__(self, domain_type, start, end, **kwargs):
        super(Domain, self).__init__(
            start, end, self.feature_type, self.feature_type, **kwargs)
        self.domain_type = domain_type

    @property
    def name(self):
        return self.domain_type

    def __repr__(self):
        template = '{self.__class__.__name__}({self.domain_type!r}, {self.start}, {self.end})'
        return template.format(self=self)

    def to_dict(self):
        d = super(Domain, self).to_dict()
        d['domain_type'] = self.domain_type
        return d


class PeptideBase(AnnotatedInterval):
    feature_type = None

    def __init__(self, start, end, **kwargs):
        super(PeptideBase, self).__init__(
            start, end, self.feature_type, self.feature_type, **kwargs)

    def __repr__(self):
        template = '{self.__class__.__name__}({self.start}, {self.end})'
        return template.format(self=self)


class SignalPeptide(PeptideBase):
    feature_type = 'signal peptide'


class Propeptide(PeptideBase):
    feature_type = 'propeptide'


class TransitPeptide(PeptideBase):
    feature_type = "transit peptide"


class MaturationPeptide(PeptideBase):
    feature_type = "maturation peptide"


class Peptide(PeptideBase):
    feature_type = 'peptide'


class MatureProtein(PeptideBase):
    feature_type = 'mature protein'


class ProteolyticSite(AnnotatedResidue):
    feature_type = 'proteolytic site'

    def __init__(self, position, **kwargs):
        super(ProteolyticSite, self).__init__(
            self.position, self.feature_type, self.feature_type, **kwargs)


class ModifiedResidue(AnnotatedResidue):
    feature_type = 'modified residue'

    def __init__(self, position, modification, **kwargs):
        super(ModifiedResidue, self).__init__(
            position, self.feature_type, modification, **kwargs)

    @property
    def modification(self):
        return self.description

    def to_dict(self):
        d = super(ModifiedResidue, self).to_dict()
        d['modification'] = d.pop('description')
        return d

    def __repr__(self):
        return "{self.__class__.__name__}({self.position}, {self.modification!r})".format(self=self)


class SimpleVariant(AnnotatedResidue):
    feature_type = "simple variant"

    def __init__(self, position, substitution, **kwargs):
        super(SimpleVariant, self).__init__(position, self.feature_type, substitution, **kwargs)

    @property
    def substitution(self):
        return self.description

    def to_dict(self):
        d = super(SimpleVariant, self).to_dict()
        d['substitution'] = d.pop('description')
        return d

    def __repr__(self):
        return "{self.__class__.__name__}({self.position}, {self.description!r})".format(self=self)

class ComplexVariant(AnnotatedInterval):
    feature_type = 'complex variant'

    def __init__(self, start, end, substitution, **kwargs):
        super(ComplexVariant, self).__init__(
            start, end, self.feature_type, substitution, **kwargs)

    @property
    def substitution(self):
        return self.description

    def to_dict(self):
        d = super(ComplexVariant, self).to_dict()
        d['substitution'] = d.pop('description')
        return d


class AnnotationCollection(Sequence[AnnotationBase]):
    def __init__(self, items):
        self.items = list(items)

    def append(self, item: AnnotationBase):
        self.items.append(item)

    def extend(self, items: Iterable[AnnotationBase]):
        self.items.extend(items)

    def __getitem__(self, i):
        return self.items[i]

    def __setitem__(self, i, item: AnnotationBase):
        self.items[i] = item

    def __len__(self):
        return len(self.items)

    def __iter__(self):
        return iter(self.items)

    def __repr__(self):
        return "{self.__class__.__name__}({self.items})".format(self=self)

    def to_json(self) -> List[Dict]:
        return [a.to_dict() for a in self]

    @classmethod
    def from_json(cls, d: Dict):
        return cls([AnnotationBase.from_dict(di) for di in d])

    def dump(self, fp):
        json.dump(self.to_json(), fp)

    @classmethod
    def load(cls, fp: io.IOBase):
        d = json.load(fp)
        inst = cls.from_json(d)
        return inst

    def __eq__(self, other):
        return self.items == list(other)

    def __ne__(self, other):
        return not (self == other)


def from_uniprot(record):
    mapping = {
        uniprot.SignalPeptide.feature_type: SignalPeptide,
        uniprot.Propeptide.feature_type: Propeptide,
        uniprot.TransitPeptide.feature_type: TransitPeptide,
        uniprot.MatureProtein.feature_type: MatureProtein,
        uniprot.ModifiedResidue.feature_type: ModifiedResidue,
    }
    annotations = []
    for feature in record.features:
        if not feature.is_defined:
            continue
        try:
            annotation_tp = mapping[feature.feature_type]
            if issubclass(annotation_tp, PeptideBase):
                annotations.append(
                    annotation_tp(feature.start, feature.end))
            elif issubclass(annotation_tp, AnnotatedResidue):
                annotations.append(
                    annotation_tp(feature.position, feature.description))
        except KeyError:
            continue
    return AnnotationCollection(annotations)


def tag_to_annotation(name: str, values: List[Tuple]) -> AnnotationCollection:
    result = []
    if name == "Processed":
        for spec in values:
            start = spec[0]
            end = spec[1]
            names = spec[2:]
            if SignalPeptide.feature_type in names:
                result.append(SignalPeptide(start - 1, end))
            elif MatureProtein.feature_type in names:
                result.append(MatureProtein(start - 1, end))
            elif MaturationPeptide.feature_type in names:
                result.append(MaturationPeptide(start - 1, end))
            elif Propeptide.feature_type in names:
                result.append(Propeptide(start - 1, end))
            else:
                annot = PeptideBase(start - 1, end)
                annot.feature_type = names[0]
                result.append(annot)
    elif name == "ModResPsi" or name == "ModRes":
        for pos, acc, name in values:
            try:
                if acc:
                    mod = Modification(acc)
                elif name:
                    mod = Modification(name)
                ModifiedResidue(pos - 1, mod.name)
            except KeyError:
                continue
    elif name == "VariantSimple":
        for pos, sub in values:
            result.append(SimpleVariant(pos - 1, sub))

    return AnnotationCollection(result)


def defline_to_annotation(defline: 'FastaHeader') -> AnnotationCollection:
    annots = [tag_to_annotation(k, v) for k, v in defline.items()]
    if annots:
        base = annots[0]
        for a in annots[1:]:
            base.extend(a)
        return base
    return AnnotationCollection([])
