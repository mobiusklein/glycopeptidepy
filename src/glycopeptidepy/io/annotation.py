import io
import json
import zlib
import warnings

from typing import List, Sequence, Union, TextIO

from six import add_metaclass

from glycopeptidepy.structure.glycan import GlycosylationType

from . import uniprot
from . import fasta


class AnnotationMeta(type):
    _cache = {}
    _is_cleavable = set()
    _is_cleavable_tags = set()

    def __new__(mcs, name, parents, attrs):
        new_type = type.__new__(mcs, name, parents, attrs)
        try:
            mcs._cache[name] = new_type
        except AttributeError:
            pass
        if new_type.cleavable and getattr(new_type, 'feature_type', None):
            mcs._is_cleavable.add(new_type)
            mcs._is_cleavable_tags.add(new_type.feature_type)
        return new_type

    def _type_for_name(cls, feature_type: str):
        return cls._cache[feature_type]

    def from_dict(cls, d: dict) -> 'AnnotationBase':
        name = d.pop('__class__')
        impl = cls._type_for_name(name)  # pylint: disable=no-value-for-parameter
        d.pop("feature_type", None)
        return impl._coerce_state(d)

    def _coerce_state(cls, state: dict):
        return cls(**state)

    @property
    def cleavable_types(cls):
        return cls._is_cleavable | cls._is_cleavable_tags


@add_metaclass(AnnotationMeta)
class AnnotationBase(object):
    feature_type: str
    description: str
    cleavable: bool = False

    def __init__(self, feature_type, description, **kwargs):
        self.feature_type = feature_type
        self.description = description

    def to_dict(self):
        d = {}
        d['feature_type'] = self.feature_type
        d['description'] = self.description
        d['__class__'] = self.__class__.__name__
        return d

    def __eq__(self, other: 'AnnotationBase'):
        if other is None:
            return False
        return self.feature_type == other.feature_type and self.description == other.description

    def __ne__(self, other: 'AnnotationBase'):
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

    domain_type: str

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
    cleavable: bool = True

    def __init__(self, start, end, **kwargs):
        super(PeptideBase, self).__init__(
            start, end, self.feature_type, self.feature_type, **kwargs)

    def __repr__(self):
        template = '{self.__class__.__name__}({self.start}, {self.end})'
        return template.format(self=self)

    def to_dict(self):
        state = super().to_dict()
        state.pop("description")
        return state


class SignalPeptide(PeptideBase):
    feature_type = 'signal peptide'


class Propeptide(PeptideBase):
    feature_type = 'propeptide'


class TransitPeptide(PeptideBase):
    feature_type = "transit peptide"


class Peptide(PeptideBase):
    feature_type = 'peptide'


class MaturationPeptide(PeptideBase):
    feature_type = 'maturation peptide'


class MatureProtein(PeptideBase):
    feature_type = 'mature protein'


class ProteolyticSite(AnnotatedResidue):
    feature_type = 'proteolytic site'

    cleavable: bool = True

    def __init__(self, position, **kwargs):
        super(ProteolyticSite, self).__init__(
            position,
            self.feature_type,
            self.feature_type,
            **kwargs)


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


class GlycosylationSite(AnnotatedResidue):
    feature_type = 'glycosylation site'

    glycosylation_type: GlycosylationType

    def __init__(self, position, glycosylation_type, **kwargs):
        if isinstance(glycosylation_type, str):
            glycosylation_type = GlycosylationType[glycosylation_type]
        super().__init__(position, self.feature_type, glycosylation_type, **kwargs)
        self.glycosylation_type = glycosylation_type

    def to_dict(self):
        d = super().to_dict()
        d['glycosylation_type'] = d.pop('description').name
        return d

    def __repr__(self):
        return "{self.__class__.__name__}({self.position}, {self.glycosylation_type!r})".format(self=self)


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


class AnnotationCollection(Sequence[Union[AnnotationBase, AnnotatedInterval, AnnotatedResidue]]):
    items: List[Union[AnnotatedResidue, AnnotatedInterval, AnnotationBase]]

    def cleavable(self):
        return sorted({item for item in self if item.cleavable}, key=lambda x: (x.start, x.end))

    def modifications(self):
        return sorted({item for item in self if isinstance(item, (ModifiedResidue, GlycosylationSite))},
                      key=lambda x: x.start)

    def __init__(self, items):
        self.items = list(items)

    def append(self, item: AnnotationBase):
        self.items.append(item)

    def __getitem__(self, i) -> AnnotationBase:
        return self.items[i]

    def __setitem__(self, i, item):
        self.items[i] = item

    def __len__(self):
        return len(self.items)

    def __iter__(self):
        return iter(self.items)

    def __repr__(self):
        return "{self.__class__.__name__}({self.items})".format(self=self)

    def to_dict(self):
        return [a.to_dict() for a in self]

    @classmethod
    def from_dict(cls, d):
        return cls([AnnotationBase.from_dict(di) for di in d])

    def dump(self, fp: TextIO):
        json.dump(self.to_dict(), fp)

    @classmethod
    def load(cls, fp: TextIO):
        d = json.load(fp)
        inst = cls.from_dict(d)
        return inst

    def to_bytes(self) -> bytes:
        buffer = io.StringIO()
        self.dump(buffer)
        return zlib.compress(buffer.getvalue().encode("utf8"))

    @classmethod
    def from_bytes(cls, data: bytes) -> 'AnnotationCollection':
        data = zlib.decompress(data)
        return cls.load(io.StringIO(data.decode("utf8")))

    def __reduce__(self):
        return self.from_bytes, (self.to_bytes(), )

    def __eq__(self, other):
        return self.items == list(other)

    def __ne__(self, other):
        return not (self == other)


def from_uniprot(record: uniprot.UniProtProtein) -> AnnotationCollection:
    mapping = {
        uniprot.SignalPeptide.feature_type: SignalPeptide,
        uniprot.Propeptide.feature_type: Propeptide,
        uniprot.TransitPeptide.feature_type: TransitPeptide,
        uniprot.MatureProtein.feature_type: MatureProtein,
        uniprot.ModifiedResidue.feature_type: ModifiedResidue,
        uniprot.GlycosylationSite: GlycosylationSite,
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
            elif isinstance(annotation_tp, GlycosylationSite):
                annotations.append(
                    annotation_tp(feature.position, feature.glycosylation_type))
            elif issubclass(annotation_tp, AnnotatedResidue):
                annotations.append(
                    annotation_tp(feature.position, feature.description))
        except KeyError:
            continue
    return AnnotationCollection(annotations)


def from_peff(record: fasta.PEFFFastaHeader) -> AnnotationCollection:
    annotations = []
    start: int
    end: int
    name: str
    accession: str
    node: tuple

    for node in record.get("Processed", []):
        if len(node) == 3:
            (start, end, name) = node
        else:
            (start, end, accession, name) = node
        start = start - 1
        if name == "signal peptide":
            annotations.append(
                SignalPeptide(start, end)
            )
        elif name == "mature protein":
            annotations.append(
                MatureProtein(start, end)
            )
        elif name == "transit peptide":
            annotations.append(
                TransitPeptide(start, end)
            )
        elif name == "propeptide":
            annotations.append(
                Propeptide(start, end)
            )
        elif name == "maturation peptide":
            annotations.append(
                MaturationPeptide(start, end)
            )
        else:
            warnings.warn(f"Unknown processing annotation {start}-{end} {name}")

    mods = record.get("ModRes", []) + record.get("ModResUnimod", [])
    for position, accession, name in mods:
        position = position - 1
        if accession:
            annotations.append(
                ModifiedResidue(position, accession))
        else:
            tokens = name.split(" ")
            if tokens[0] in GlycosylationType:
                glycosylation_type = GlycosylationType[tokens[0]]
                if glycosylation_type == GlycosylationType.o_linked:
                    if "glycosaminoglycan" in name:
                        glycosylation_type = GlycosylationType.glycosaminoglycan
                annotations.append(
                    GlycosylationSite(position, glycosylation_type))
    return AnnotationCollection(annotations)

