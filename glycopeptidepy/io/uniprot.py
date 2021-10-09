'''A REST API for the UniProt protein database and a parser for their XML format.
'''
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
batch_uri = "https://www.uniprot.org/uploadlists/"
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
            feature_type_t = self.feature_type(feature_type)  # pylint: disable=no-value-for-parameter
            return feature_type_t.fromxml(tag)
        except KeyError:
            return None


@add_metaclass(UniProtFeatureMeta)
class UniProtFeatureBase(object):
    __repr__ = simple_repr

    known = True
    _name = None

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

    def __reduce__(self):
        return self.__class__, (str(self), self.id)


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
        except (KeyError, AttributeError):
            begin = 0
            known = False
        try:
            end = int(feature.find(".//{http://uniprot.org/uniprot}end", nsmap).attrib['position'])
        except (KeyError, AttributeError):
            end = float('inf')
            known = False
        if not known:
            if begin == 0 and end == float('inf'):
                begin = int(feature.find(
                    ".//{http://uniprot.org/uniprot}location", nsmap).attrib['position']) - 1
                end = begin + 1
                known = True
        return cls(begin, end, known)


class SignalPeptide(PeptideBase):
    feature_type = 'signal peptide'


class Propeptide(PeptideBase):
    feature_type = 'propeptide'


class TransitPeptide(PeptideBase):
    feature_type = "transit peptide"


class Peptide(PeptideBase):
    feature_type = 'peptide'


class MatureProtein(PeptideBase):
    feature_type = 'mature protein'


class SpanBase(UniProtFeatureBase):
    feature_type = "__spanbase__"

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

    @property
    def description(self):
        return self.name

    def __repr__(self):
        repr_str = super(SpanBase, self).__repr__()
        prefix, suffix = repr_str.split("(", 1)
        augmented = "%s(name=%r, %s" % (prefix, self.name, suffix)
        return augmented


class Domain(SpanBase):
    feature_type = 'domain'


class Chain(SpanBase):
    feature_type = 'chain'


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

    @property
    def name(self):
        return self.modification

    @property
    def description(self):
        return self.modification

    @property
    def start(self):
        return self.position

    @property
    def end(self):
        return self.position + 1

    @classmethod
    def fromxml(cls, feature):
        return cls(
            int(feature.find(
                ".//up:position", nsmap).attrib['position']) - 1,
            feature.attrib["description"])


class SequenceVariant(UniProtFeatureBase):
    feature_type = "sequence variant"

    def __init__(self, position, original, variation):
        self.position = position
        self.original = original
        self.variation = variation

    @property
    def start(self):
        return self.position

    @property
    def end(self):
        return self.position + 1

    @classmethod
    def fromxml(cls, feature):
        position = feature.find(".//up:position", nsmap)
        if position is None:
            # Deletion mutations not supported
            return None
        else:
            position = int(position.attrib['position']) - 1
            original = feature.find(".//up:original", nsmap)
            if original is not None:
                original = original.text
            variation = feature.find(".//up:variation", nsmap)
            if variation is not None:
                variation = variation.text
        return cls(position, original, variation)

    @property
    def description(self):
        return "%s->%s" % (self.original, self.variation)


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


_UniProtProtein = make_struct("UniProtProtein", (
    "sequence", "features", "recommended_name", "gene_name", "names", "accessions",
    "keywords", "dbreferences"))


class UniProtProtein(_UniProtProtein):
    __slots__ = ()

    def features_of_type(self, feature_type):
        if isinstance(feature_type, basestring):
            return [feature for feature in self.features if feature.feature_type == feature_type]
        elif isinstance(feature_type, type):
            return [feature for feature in self.features if isinstance(feature, feature_type)]
        else:
            raise TypeError("Can't use object of type %r, %r" % (type(feature_type), feature_type))

    @property
    def chains(self):
        return self.features_of_type(Chain)

    @property
    def domains(self):
        return self.features_of_type(Domain)

    @property
    def variants(self):
        return self.features_of_type(SequenceVariant)

    @property
    def peptides(self):
        return self.features_of_type(PeptideBase)


def get_etree_for(accession):
    tree = etree.parse(_open_url_for(accession))
    return tree


def _open_url_for(accession):
    return urlopen(uri_template.format(accession=accession))


def parse(tree, error=False):
    '''Parse an XML document into a :class:`UniProtProtein`.

    Parameters
    ----------
    tree: lxml.etree.Element
        The root XML element of the tree to parse.
    errors: bool
        Whether to treat errors during feature parsing
        as warnings or exceptions.

    Returns
    -------
    UniProtProtein
    '''
    entry = tree.find(".//up:entry", nsmap)
    if entry is None:
        if hasattr(tree, "tag") and tree.tag.endswith("entry"):
            entry = tree
        else:
            raise ValueError("Could not find root entry!")
    seq = entry.find(
        "./up:sequence", nsmap).text
    if seq is not None:
        seq = seq.replace("\n", '')
    names = [el.text for el in entry.findall(
        ".//up:protein/*/up:fullName", nsmap)] + [el.text for el in entry.findall(
            ".//up:protein/*/up:shortName", nsmap)]
    recommended_name_tags = entry.findall(
        ".//up:protein/up:recommendedName", nsmap) + entry.findall(".//up:protein/*/up:recommendedName", nsmap)
    if recommended_name_tags:
        recommended_name_tag = recommended_name_tags[0]
    else:
        recommended_name_tag = None
    if recommended_name_tag is not None:
        if recommended_name_tag.text.strip():
            recommended_name = recommended_name_tag.text.strip()
        else:
            recommended_name = (recommended_name_tag[0].text)
    else:
        try:
            recommended_name = names[0]
        except IndexError:
            recommended_name = ""
    gene_name_tag = entry.find(".//up:name", nsmap)
    if gene_name_tag is not None:
        gene_name = gene_name_tag.text
    else:
        gene_name = ""
    accessions = [el.text for el in entry.findall(
        ".//up:accession", nsmap)]
    dbreferences = defaultdict(list)
    for tag in entry.findall(".//up:dbReference", nsmap):
        db = tag.attrib['type']
        rec = dict(id=tag.attrib['id'])
        for prop in tag.findall(".//up:property", nsmap):
            rec[prop.attrib['type']] = prop.attrib['value']
        dbreferences[db].append(rec)

    features = []
    exc_type = Exception
    for tag in entry.findall(".//up:feature", nsmap):
        feature_type = tag.attrib['type']
        try:
            feature_obj = UniProtFeatureBase.handle_tag(tag)
            if feature_obj is not None:
                features.append(feature_obj)
        except exc_type as e:
            if error:
                raise
            else:
                warnings.warn("An exception %r occurred while parsing feature type %s for %s" % (
                    e, feature_type, accessions[0]))
    keywords = set()
    for kw in entry.findall(".//up:keyword", nsmap):
        keywords.add(Keyword(kw.text, kw.attrib['id']))
    return UniProtProtein(
        seq, features, recommended_name, gene_name, names, accessions,
        keywords, dbreferences)


def parse_all(tree, error=False):
    """Parse an XML document containing multiple <up:entry> elements into
    a list of :class:`UniProtProtein` instances.

    Parameters
    ----------
    tree : lxml.etree.Element
        The root element of the XML document to parse.
    errors: bool
        Whether to treat errors during feature parsing
        as warnings or exceptions.

    Returns
    -------
    list
    """
    entries = tree.findall(".//up:entry", nsmap)
    return [parse(e, error=error) for e in entries]


def get_features_for(accession, error=False):
    tree = get_etree_for(accession)
    return parse(tree, error)


def get_features_for_many(accessions, error=False):
    query_to_doc = dict()
    r = requests.post(batch_uri, params={
        "format": "xml",
        "from": "ACC+ID",
        "to": "ACC",
        "query": " ".join(accessions)
    })
    tree = etree.fromstring(r.content)
    for doc in parse_all(tree, error=error):
        for acc in doc.accessions:
            query_to_doc[acc] = doc
    results = []
    for acc in accessions:
        if acc in query_to_doc:
            results.append((acc, query_to_doc[acc]))
        else:
            results.append((acc, None))
    return results


def get(accessions):
    if accessions is None:
        raise TypeError("accessions cannot be `None`")
    if isinstance(accessions, basestring):
        return get_features_for(accessions)
    elif len(accessions) < 5:
        return [b for a, b in get_features_for_many(accessions) if b is not None]
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
        proteins = []
        for chunk in chunks:
            proteins.extend(cls.fetch(chunk))
        if not ordered:
            return proteins
        acc_map = {}
        for prot in proteins:
            for acc in prot.accessions:
                acc_map[acc] = prot
        proteins = []
        for acc in accession_iterable:
            proteins.append(acc_map.get(acc))
        return proteins

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


class _UniProtRestClientBase(object):

    def _get_response_buffer(self, url):
        response = requests.get(url, stream=True)
        response.raise_for_status()
        response_buffer = response.raw
        data = response_buffer.read()
        if response.headers.get("content-encoding") == 'gzip':
            response_buffer = gzip.GzipFile(
                mode='rb', fileobj=io.BytesIO(data))
            data = response_buffer.read()
        return data



class ProteomeSearchResult(dict):
    bind = None

    def download(self, raw=False, format='fasta'):
        return self.bind.get(self['Proteome ID'], raw=raw, format=format)


class ProteomeClient(_UniProtRestClientBase):

    SearchResultType = ProteomeSearchResult

    def search(self, format='tab', raw=False, **kwargs):
        query = "&".join(["query=%s:%s" % (k, v) for k, v in kwargs.items()])
        uri = "https://www.uniprot.org/proteomes/?{query}&sort=score&format={format}&sort=score"
        url = uri.format(query=query, format=format)
        data = self._get_response_buffer(url)
        if raw:
            return io.BytesIO(data)
        if format == 'tab':
            result = list(map(self.SearchResultType, csv.DictReader(io.BytesIO(data), delimiter='\t')))
            for item in result:
                item.bind = self
            return result
        elif format == 'list':
            result = []
            for line in data.splitlines():
                r = self.SearchResultType({
                    "Proteome ID": line.strip()
                })
                r.bind = self
                result.append(r)
            return result
        return io.BytesIO(data)

    def get(self, proteome_id, include_isoforms=False, raw=False, format='fasta'):
        uri = "https://www.uniprot.org/uniprot/?{query}&{params}&force=no&format={format}"
        query = {"proteome": proteome_id, }
        query = '&'.join(["query=%s:%s" % (k, v) for k, v in query.items()])
        params = {"include": str(bool(include_isoforms)).lower()}
        params = '&'.join(["%s=%s" % (k, v) for k, v in params.items()])
        url = uri.format(query=query, params=params, format=format)
        data = self._get_response_buffer(url)
        if raw:
            return io.BytesIO(data)
        if format == 'fasta':
            from .fasta import ProteinFastaFileReader
            return ProteinFastaFileReader(io.BytesIO(data))
        elif format == 'xml':
            return etree.fromstring(data)
        else:
            return io.BytesIO(data)

proteome = ProteomeClient()


class TaxonomyClient(_UniProtRestClientBase):

    def search(self, format='tab', raw=False, params=None, **kwargs):
        if params is None:
            params = dict()
        query = "&".join(["query=%s:%s" % (k, v) for k, v in kwargs.items()])
        params = '&'.join(["%s=%s" % (k, v) for k, v in params.items()])
        if params:
            params += '&'
        uri = "https://www.uniprot.org/taxonomy/?{query}&{params}sort=score&format={format}&sort=score"
        url = uri.format(query=query, params=params, format=format)
        data = self._get_response_buffer(url)
        if raw:
            return io.BytesIO(data)
        if format == 'tab':
            result = list(map(dict, csv.DictReader(
                io.BytesIO(data), delimiter='\t')))
            return result
        elif format == 'list':
            result = data.splitlines()
            return result
        return io.BytesIO(data)


taxonomy = TaxonomyClient()
