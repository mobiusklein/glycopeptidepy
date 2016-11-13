import threading

from lxml import etree
from glycopeptidepy.structure import glycan
from glypy.utils import make_struct
from glycopeptidepy.utils import simple_repr


uri_template = "http://www.uniprot.org/uniprot/{accession}.xml"


class UniProtFeatureBase(object):
    __repr__ = simple_repr

    known = True

    @property
    def is_defined(self):
        return self.known


class SignalPeptide(UniProtFeatureBase):
    feature_type = 'signal peptide'

    def __init__(self, start, end, known=True):
        self.start = start
        self.end = end
        self.known = known

    @classmethod
    def fromxml(cls, feature):
        known = True
        try:
            begin = int(feature.find(".//{http://uniprot.org/uniprot}begin").attrib['position'])
        except KeyError:
            begin = 0
            known = False
        try:
            end = int(feature.find(".//{http://uniprot.org/uniprot}end").attrib['position'])
        except KeyError:
            end = float('inf')
            known = False
        return cls(begin, end, known)

    @property
    def name(self):
        return "signal peptide"


class Domain(UniProtFeatureBase):
    feature_type = 'domain'

    def __init__(self, name, start, end, known=True):
        self.name = name
        self.start = start
        self.end = end
        self.known = known

    @classmethod
    def fromxml(cls, feature):
        known = True
        try:
            begin = int(feature.find(".//{http://uniprot.org/uniprot}begin").attrib['position'])
        except KeyError:
            begin = 0
            known = False
        try:
            end = int(feature.find(".//{http://uniprot.org/uniprot}end").attrib['position'])
        except KeyError:
            end = float('inf')
            known = False
        description = feature.attrib.get("description")
        return cls(description, begin, end, known)


class RegionOfInterest(Domain):
    feature_type = 'region of interest'


class DisulfideBond(Domain):
    feature_type = "disulfide bond"

    @classmethod
    def fromxml(cls, feature):
        name = "disulfide bond"
        known = True
        try:
            try:
                begin = int(feature.find(".//{http://uniprot.org/uniprot}begin").attrib['position'])
            except KeyError:
                begin = 0
                known = False
            try:
                end = int(feature.find(".//{http://uniprot.org/uniprot}end").attrib['position'])
            except KeyError:
                end = float('inf')
                known = False
        except AttributeError:
            begin = end = int(feature.find(".//{http://uniprot.org/uniprot}position").attrib['position'])
        return cls(name, begin, end, known)


class ModifiedResidue(UniProtFeatureBase):
    feature_type = 'modified residue'

    def __init__(self, position, modification):
        self.position = position
        self.modification = modification

    @classmethod
    def fromxml(cls, feature):
        return cls(
            int(feature.find(".//{http://uniprot.org/uniprot}position").attrib['position']),
            feature.attrib["description"])


class GlycosylationSite(UniProtFeatureBase):
    feature_type = 'glycosylation site'

    def __init__(self, position, glycosylation_type):
        self.position = position
        self.glycosylation_type = glycosylation_type

    @property
    def name(self):
        return self.glycosylation_type

    @property
    def start(self):
        return self.position

    @property
    def end(self):
        return self.position + 1

    @classmethod
    def fromxml(cls, feature):
        position = int(feature.find(".//{http://uniprot.org/uniprot}position").attrib['position'])
        glycosylation_type = feature.attrib["description"].split(" ")[0]
        return cls(position, glycosylation_type)


UniProtProtein = make_struct("UniProtProtein", ("sequence", "features", "names"))


def get_etree_for(accession):
    tree = etree.parse(uri_template.format(accession=accession))
    return tree


def get_features_for(accession):
    tree = etree.parse(uri_template.format(accession=accession))
    seq = tree.find(".//{http://uniprot.org/uniprot}entry/{http://uniprot.org/uniprot}sequence").text.replace("\n", '')
    names = [el.text for el in tree.findall(
        ".//{http://uniprot.org/uniprot}protein/*/{http://uniprot.org/uniprot}fullName")]
    features = []
    for tag in tree.findall(".//{http://uniprot.org/uniprot}feature"):
        feature_type = tag.attrib['type']
        try:
            if feature_type == 'signal peptide':
                features.append(SignalPeptide.fromxml(tag))
            elif feature_type == "modified residue":
                features.append(ModifiedResidue.fromxml(tag))
            elif feature_type == "glycosylation site":
                features.append(GlycosylationSite.fromxml(tag))
            elif feature_type == 'domain':
                features.append(Domain.fromxml(tag))
            elif feature_type == 'region of interest':
                features.append(RegionOfInterest.fromxml(tag))
            elif feature_type == 'disulfide bond':
                features.append(DisulfideBond.fromxml(tag))
        except Exception as e:
            print(e, feature_type, accession, etree.tostring(tag))
    return UniProtProtein(seq, features, names)


class ProteinDownloader(object):
    @staticmethod
    def chunk(seq, n=5):
        i = 0
        while i < len(seq):
            yield seq[i:i + n]
            i += n

    @staticmethod
    def fetch(chunk):
        jobs = []
        results = []

        def task(name):
            try:
                feat = get_features_for(name)
                results.append((name, feat))
            except Exception as e:
                print(e, name)

        for c in chunk:
            th = threading.Thread(target=task, args=(c,))
            th.start()
            jobs.append(th)

        for th in jobs:
            th.join()

        return results

    @classmethod
    def download(cls, accession_iterable, chunk_size=15):
        chunks = cls.chunk(accession_iterable, chunk_size)
        proteome = []
        for chunk in chunks:
            proteome.extend(cls.fetch(chunk))
        return proteome

    def __call__(self, *args, **kwargs):
        return self.download(*args, **kwargs)
