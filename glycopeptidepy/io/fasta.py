import re
import textwrap
import warnings

from collections import OrderedDict

from glycopeptidepy.structure.sequence import ProteinSequence


def tryopen(path):
    if hasattr(path, 'read'):
        return path
    return open(path)


class FastaHeader(object):
    def __init__(self, mapping, original=None):
        self._mapping = mapping
        if original is None:
            self.defline = self._make_defline()
        else:
            self.defline = original

    def __getitem__(self, key):
        return self._mapping[key]

    def __iter__(self):
        return iter(self._mapping)

    def items(self):
        return self._mapping.items()

    def keys(self):
        return self._mapping.keys()

    def __contains__(self, key):
        return key in self._mapping

    def __getattr__(self, key):
        if key == "_mapping":
            raise AttributeError(key)
        try:
            return self._mapping[key]
        except KeyError:
            raise AttributeError(key)

    def __str__(self):
        return self.defline

    def __repr__(self):
        return "{self.__class__.__name__}({self._mapping})".format(self=self)

    def _make_defline(self):
        raise NotImplementedError()

    def __hash__(self):
        return hash(self.defline)

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return not (self == other)


class UnparsableDeflineError(Exception):
    pass


class DefLineParserBase(object):
    def __call__(self, defline):
        if defline.startswith(">"):
            defline = defline[1:]
        return FastaHeader(self.parse(defline), defline)

    def parse(self, defline):
        raise NotImplementedError()


class CallableDefLineParser(DefLineParserBase):
    def __init__(self, callable_obj):
        self.callable_obj = callable_obj

    def parse(self, defline):
        return OrderedDict(enumerate(self.callable_obj(defline)))


def split_on_whitespace(string):
    return re.split("\s+", string)


space_delim_parser = CallableDefLineParser(split_on_whitespace)


class RegexDefLineParser(DefLineParserBase):
    def __init__(self, pattern):
        self.pattern = re.compile(pattern)
        self._is_keyed = "?P<" in pattern

    def parse(self, defline):
        match = self.pattern.search(defline)
        if match:
            if not self._is_keyed:
                return OrderedDict(enumerate(match.groups()))
            else:
                return (match.groupdict())
        else:
            raise UnparsableDeflineError(defline)


class DispatchingDefLineParser(DefLineParserBase):
    def __init__(self, patterns, graceful=False, header_type=FastaHeader):
        self.parsers = patterns
        self.graceful = graceful

    def parse(self, defline):
        for parser in self.parsers:
            try:
                parts = parser.parse(defline)
                return parts
            except UnparsableDeflineError:
                continue
        else:
            if self.graceful:
                return defline
            else:
                raise UnparsableDeflineError(defline)


uniprot_regex = r"(?P<db>[a-z^\|]+)\|(?P<accession>[A-Z0-9]+)\|(?P<name>\S+)(?:\s(?P<description>.+))?"
uniprot_parser = RegexDefLineParser(uniprot_regex)


default_parser = DispatchingDefLineParser([uniprot_parser, space_delim_parser])


class FastaFileParser(object):
    def __init__(self, path, defline_parser=default_parser):
        self.state = "defline"
        self.handle = tryopen(path)
        self.defline = None
        self.sequence_chunks = []
        self.defline_parser = defline_parser
        self._generator = None

    def process_result(self, d):
        return d

    def __iter__(self):
        return self

    def __next__(self):
        if self._generator is None:
            self._generator = self._parse_lines()
        return next(self._generator)

    def next(self):
        return self.__next__()

    def _parse_lines(self):
        for line in self.handle:
            if self.state == 'defline':
                if line[0] == ">":
                    self.defline = re.sub(r"[\n\r]", "", line[1:])
                    self.state = "sequence"
                else:
                    continue
            else:
                if not re.match(r"^(\s+|>)", line):
                    self.sequence_chunks.append(re.sub(r"[\n\r]", "", line))
                else:
                    if self.defline is not None:
                        try:
                            yield self.process_result({
                                "name": self.defline_parser(self.defline),
                                "sequence": ''.join(self.sequence_chunks)
                            })
                        except KeyError as e:
                            print(e)
                            pass
                    self.sequence_chunks = []
                    self.defline = None
                    self.state = 'defline'
                    if line[0] == '>':
                        self.defline = re.sub(r"[\n\r]", "", line[1:])
                        self.state = "sequence"

        if len(self.sequence_chunks) > 0:
            try:
                yield self.process_result({"name": self.defline_parser(self.defline), "sequence": ''.join(self.sequence_chunks)})
            except Exception as e:
                pass


class ProteinFastaFileParser(FastaFileParser):

    def __init__(self, path, defline_parser=default_parser):
        super(ProteinFastaFileParser, self).__init__(path, defline_parser)

    def process_result(self, d):
        p = ProteinSequence(d['name'], d['sequence'])
        return p


class FastaFileWriter(object):

    def __init__(self, handle):
        self.handle = handle

    def write(self, defline, sequence):
        self.handle.write(defline)
        self.handle.write("\n")
        self.handle.write(sequence)
        self.handle.write("\n\n")

    def writelines(self, iterable):
        for defline, seq in iterable:
            self.write(defline, seq)


class ProteinFastaFileWriter(FastaFileWriter):

    def write(self, protein):
        defline = ">%s" % protein.name
        seq = '\n'.join(textwrap.wrap(protein.get_sequence(), 80))
        super(ProteinFastaFileWriter, self).write(defline, seq)

    def writelines(self, iterable):
        for protein in iterable:
            self.write(protein)


def read(f, defline_parser=default_parser, reader_type=ProteinFastaFileParser):
    return reader_type(f, defline_parser=defline_parser)
