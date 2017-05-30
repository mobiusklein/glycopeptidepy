import re
import threading
import urllib2

from lxml import etree
from glycopeptidepy.structure import glycan
from glycopeptidepy.structure.modification import ModificationRule
from glycopeptidepy.utils import simple_repr
from glypy import Composition
from glypy.utils import make_struct


uri_template = "http://www.uniprot.org/uniprot/{accession}.xml"


class UniProtFeatureBase(object):
    __repr__ = simple_repr

    known = True

    @property
    def is_defined(self):
        return self.known

    def describe(self):
        return self.description.split("; ")


class Keyword(str):
    def __new__(cls, content, idstr):
        obj = str.__new__(cls, content)
        obj.id = idstr
        return obj


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
            begin = int(feature.find(".//{http://uniprot.org/uniprot}position").attrib['position'])
            end = begin + 1

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


class Site(Domain):
    feature_type = 'site'


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


UniProtProtein = make_struct("UniProtProtein", (
    "sequence", "features", "recommended_name", "gene_name", "names", "accessions",
    "keywords"))


def get_etree_for(accession):
    tree = etree.parse(uri_template.format(accession=accession))
    return tree


def get_features_for(accession, error=False):
    tree = etree.parse(uri_template.format(accession=accession))
    seq = tree.find(
        ".//{http://uniprot.org/uniprot}entry/{http://uniprot.org/uniprot}sequence").text.replace(
        "\n", '')
    names = [el.text for el in tree.findall(
        ".//{http://uniprot.org/uniprot}protein/*/{http://uniprot.org/uniprot}fullName")]
    recommended_name_tag = tree.find(
        ".//{http://uniprot.org/uniprot}protein/*/{http://uniprot.org/uniprot}recommendedName")
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
    gene_name_tag = tree.find(".//{http://uniprot.org/uniprot}entry/{http://uniprot.org/uniprot}name")
    if gene_name_tag is not None:
        gene_name = gene_name_tag.text
    else:
        gene_name = ""
    accessions = [el.text for el in tree.findall(
        ".//{http://uniprot.org/uniprot}accession")]
    features = []
    exc_type = Exception if not error else KeyboardInterrupt
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
            elif feature_type == "site":
                features.append(Site.fromxml(tag))
        except exc_type as e:
            print(e, feature_type, accession, etree.tostring(tag))
    keywords = set()
    for kw in tree.findall(".//{http://uniprot.org/uniprot}keyword"):
        keywords.add(Keyword(kw.text, kw.attrib['id']))
    return UniProtProtein(
        seq, features, recommended_name, gene_name, names, accessions, keywords)


get = get_features_for


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
    def download(cls, accession_iterable, chunk_size=15):
        chunks = cls.chunk(accession_iterable, chunk_size)
        proteome = []
        for chunk in chunks:
            proteome.extend(cls.fetch(chunk))
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


class _UniProtPTMListParser(object):
    def __init__(self, path):
        self.path = path
        self.handle = open(path)
        self._find_starting_point()

    def _find_starting_point(self):
        self.handle.seek(0)
        for line in self.handle:
            if line.startswith("___"):
                break

    def _next_line(self):
        line = next(self.handle)
        line = line.strip("\n")
        tokens = re.split("\s+", line, maxsplit=1)
        if len(tokens) == 1:
            return tokens[0], ""
        else:
            typecode, content = tokens
            return typecode, content

    def _formula_parser(self, formula):
        counts = dict()
        for symbol, count in re.findall(r"([A-Za-z]+)(-?\d+)", formula):
            count = int(count)
            counts[symbol] = count
        return Composition(counts)

    def parse_entry(self):
        ptm_id = None
        accession = None
        feature_key = None
        mass_difference = None
        correction_formula = None
        keywords = []
        crossref = []
        while True:
            typecode, line = self._next_line()
            if typecode == "ID":
                ptm_id = line
            elif typecode == "AC":
                accession = line
            elif typecode == "FT":
                feature_key = line
            elif typecode == "CF":
                correction_formula = self._formula_parser(line)
            elif typecode == "MM":
                mass_difference = float(line)
            elif typecode == "KW":
                keywords.append(line.strip("."))
            elif typecode == "DR":
                crossref.append(line.strip('.'))
            elif typecode == "//":
                break
        if mass_difference is None or correction_formula is None:
            return None
        rule = ModificationRule(
            [], ptm_id, monoisotopic_mass=mass_difference, composition=correction_formula,
            alt_names={accession}, categories=keywords)
        return rule

    def build_table(self):
        table = {}
        while True:
            try:
                entry = self.parse_entry()
                if entry is None:
                    continue
                table[entry.name] = entry
            except StopIteration:
                break
        return table
