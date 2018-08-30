import re
import csv
import warnings
import threading
import gzip
import io

try:
    from urllib import urlopen
except ImportError:
    from urllib.request import urlopen

from collections import defaultdict

import requests

from lxml import etree
from glycopeptidepy.utils import simple_repr
from glypy.utils import make_struct
from glypy.io.glyspace import UniprotRDFClient

from six import add_metaclass, string_types as basestring

uri_template = "https://www.uniprot.org/uniprot/{accession}.xml"
nsmap = {"up": "http://uniprot.org/uniprot"}


class UniProtFeatureMeta(type):
    _cache = {}

    def __new__(cls, name, parents, attrs):
        new_type = type.__new__(cls, name, parents, attrs)
        try:
            cls._cache[new_type.feature_type] = new_type
        except AttributeError:
            pass
        return new_type

    def feature_type(self, feature_type):
        return self._cache[feature_type]

    def handle_tag(self, tag):
        feature_type = tag.attrib['type']
        try:
            feature_type_t = self.feature_type(feature_type)
            return feature_type_t.fromxml(tag)
        except KeyError:
            return None


@add_metaclass(UniProtFeatureMeta)
class UniProtFeatureBase(object):
    __repr__ = simple_repr

    known = True

    @property
    def is_defined(self):
        return self.known

    def describe(self):
        return self.description.split("; ")

    @property
    def name(self):
        return self.feature_type


class Keyword(str):
    def __new__(cls, content, idstr):
        obj = str.__new__(cls, content)
        obj.id = idstr
        return obj


class PeptideBase(UniProtFeatureBase):

    def __init__(self, start, end, known=True):
        self.start = start
        self.end = end
        self.known = known

    @classmethod
    def fromxml(cls, feature):
        known = True
        try:
            begin = int(feature.find(".//{http://uniprot.org/uniprot}begin", nsmap).attrib['position']) - 1
        except KeyError:
            begin = 0
            known = False
        try:
            end = int(feature.find(".//{http://uniprot.org/uniprot}end", nsmap).attrib['position'])
        except KeyError:
            end = float('inf')
            known = False
        return cls(begin, end, known)


class SignalPeptide(PeptideBase):
    feature_type = 'signal peptide'


class Propeptide(PeptideBase):
    feature_type = 'propeptide'


class TransitPeptide(PeptideBase):
    feature_type = "transit peptide"


class Peptide(PeptideBase):
    feature_type = 'peptide'


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
            try:
                begin = int(feature.find(
                    ".//up:begin", nsmap).attrib['position']) - 1
            except KeyError:
                begin = 0
                known = False
            try:
                end = int(feature.find(
                    ".//up:end", nsmap).attrib['position'])
            except KeyError:
                end = float('inf')
                known = False
        except AttributeError:
            begin = int(feature.find(
                ".//up:position", nsmap).attrib['position']) - 1
            end = begin + 1

        description = feature.attrib.get("description")
        return cls(description, begin, end, known)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value


class RegionOfInterest(Domain):
    feature_type = 'region of interest'


class InitiatorMethionine(Domain):
    feature_type = "initiator methionine"


class DisulfideBond(Domain):
    feature_type = "disulfide bond"

    @classmethod
    def fromxml(cls, feature):
        name = "disulfide bond"
        known = True
        try:
            try:
                begin = int(feature.find(
                    ".//up:begin", nsmap).attrib['position']) - 1
            except KeyError:
                begin = 0
                known = False
            try:
                end = int(feature.find(
                    ".//up:end", nsmap).attrib['position'])
            except KeyError:
                end = float('inf')
                known = False
        except AttributeError:
            begin = end = int(feature.find(
                ".//up:position", nsmap).attrib['position'])
        return cls(name, begin, end, known)


class ModifiedResidue(UniProtFeatureBase):
    feature_type = 'modified residue'

    def __init__(self, position, modification):
        self.position = position
        self.modification = modification

    @classmethod
    def fromxml(cls, feature):
        return cls(
            int(feature.find(
                ".//up:position", nsmap).attrib['position']) - 1,
            feature.attrib["description"])


class Site(Domain):
    feature_type = 'site'


class GlycosylationSite(UniProtFeatureBase):
    feature_type = 'glycosylation site'

    def __init__(self, position, glycosylation_type, description=None):
        self.position = position
        self.glycosylation_type = glycosylation_type
        self.description = description

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
        position = int(feature.find(
            ".//up:position", nsmap).attrib['position']) - 1
        glycosylation_type = feature.attrib["description"].split(" ")[0]
        return cls(position, glycosylation_type, feature.attrib["description"])


UniProtProtein = make_struct("UniProtProtein", (
    "sequence", "features", "recommended_name", "gene_name", "names", "accessions",
    "keywords", "dbreferences"))


def get_etree_for(accession):
    tree = etree.parse(urlopen(uri_template.format(accession=accession)))
    return tree


def parse(tree, error=False):
    seq = tree.find(
        ".//up:entry/up:sequence", nsmap).text.replace(
        "\n", '')
    names = [el.text for el in tree.findall(
        ".//up:protein/*/up:fullName", nsmap)]
    recommended_name_tag = tree.find(
        ".//up:protein/*/up:recommendedName", nsmap)
    if recommended_name_tag is not None:
        if recommended_name_tag.text.strip():
            recommended_name = recommended_name_tag.text.strip()
        else:
            recommended_name = ' '.join(c.text for c in recommended_name_tag)
    else:
        try:
            recommended_name = names[0]
        except IndexError:
            recommended_name = ""
    gene_name_tag = tree.find(".//up:entry/up:name", nsmap)
    if gene_name_tag is not None:
        gene_name = gene_name_tag.text
    else:
        gene_name = ""
    accessions = [el.text for el in tree.findall(
        ".//up:accession", nsmap)]
    dbreferences = defaultdict(list)
    for tag in tree.findall(".//up:dbReference", nsmap):
        db = tag.attrib['type']
        entry = dict(id=tag.attrib['id'])
        for prop in tag.findall(".//up:property", nsmap):
            entry[prop.attrib['type']] = prop.attrib['value']
        dbreferences[db].append(entry)

    features = []
    exc_type = Exception
    for tag in tree.findall(".//up:feature", nsmap):
        feature_type = tag.attrib['type']
        try:
            feature_obj = UniProtFeatureBase.handle_tag(tag)
            if feature_obj is not None:
                features.append(feature_obj)
        except exc_type as e:
            if error:
                raise
            else:
                warnings.warn("An exception %s occurred while parsing feature type %s" % (e, feature_type))
    keywords = set()
    for kw in tree.findall(".//up:keyword", nsmap):
        keywords.add(Keyword(kw.text, kw.attrib['id']))
    return UniProtProtein(
        seq, features, recommended_name, gene_name, names, accessions,
        keywords, dbreferences)


def get_features_for(accession, error=False):
    tree = get_etree_for(accession)
    return parse(tree, error)


def get(accessions):
    if accessions is None:
        raise TypeError("accessions cannot be `None`")
    if isinstance(accessions, basestring):
        return get_features_for(accessions)
    elif len(accessions) < 5:
        return list(map(get_features_for, accessions))
    else:
        chunk_size = min(len(accessions) // 5, 15)
        return ProteinDownloader.download(accessions, chunk_size)


class ProteinDownloader(object):  # pragma: no cover
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
                results.append(feat)
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
    def download(cls, accession_iterable, chunk_size=15, ordered=True):
        accession_iterable = list(accession_iterable)
        chunks = cls.chunk(accession_iterable, chunk_size)
        proteome = []
        for chunk in chunks:
            proteome.extend(cls.fetch(chunk))
        if not ordered:
            return proteome
        acc_map = {}
        for prot in proteome:
            for acc in prot.accessions:
                acc_map[acc] = prot
        proteome = []
        for acc in accession_iterable:
            proteome.append(acc_map.get(acc))
        return proteome

    def __call__(self, *args, **kwargs):
        return self.download(*args, **kwargs)

    @classmethod
    def to_fasta_file(cls, accession_iterable, file_handle, chunk_size=15):
        from .fasta import FastaFileWriter
        proteins = cls.download(accession_iterable, chunk_size=chunk_size)
        writer = FastaFileWriter(file_handle)
        for protein in proteins:
            defline = ">sp|%s|%s %s" % (protein.accessions[0], protein.gene_name, protein.recommended_name)
            defline = defline.replace("\n", " ")
            seq = protein.sequence
            writer.write(defline, seq)
        return file_handle


rdf = UniprotRDFClient()


def search(query):
    url = (
        "https://www.uniprot.org/uniprot/?sort=score&"
        "desc=&compress=no&query={query}&fil=&&format=tab"
        "&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length")
    response = requests.get(url.format(query=query), stream=True)
    response.raise_for_status()
    response_buffer = response.raw
    data = response_buffer.read()
    if response.headers.get("content-encoding") == 'gzip':
        response_buffer = gzip.GzipFile(mode='rb', fileobj=io.BytesIO(data))
        data = response_buffer.read().decode('utf-8')
    else:
        data = data.decode("utf-8")
    reader = csv.reader(data.splitlines(), delimiter='\t')
    header = next(reader)
    return [dict(zip(header, line)) for line in reader]


def proteome(format='fasta', raw=False, **kwargs):
    query = "&".join(["query=%s:%s" % (k, v) for k, v in kwargs])
    uri = "https://www.uniprot.org/uniprot/{query}&force=no&format={format}"
    url = uri.format(query=query, format=format)
    response = requests.get(url, stream=True)
    response_buffer = response.raw
    data = response_buffer.read()
    if response.headers.get("content-encoding") == 'gzip':
        response_buffer = gzip.GzipFile(mode='rb', fileobj=io.BytesIO(data))
        data = response_buffer.read()
    if raw:
        return data
    if format == 'fasta':
        from .fasta import ProteinFastaFileReader
        return ProteinFastaFileReader(io.BytesIO(data))
    elif format == 'xml':
        return etree.fromstring(data)
    else:
        return io.BytesIO(data)
