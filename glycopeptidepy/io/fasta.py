import re
import textwrap

from collections import OrderedDict
try:
    from collections.abc import Mapping, Sequence as SequenceABC
except ImportError:
    from collections import Mapping, Sequence as SequenceABC

from glycopeptidepy.structure.sequence import ProteinSequence
from glypy.utils.base import opener


class FastaHeader(Mapping):
    """Hold parsed properties of a FASTA sequence's
    definition line.

    This object supports the :class:`Mapping` interface, and
    keys may be accessed by attribute access notation.

    Attributes
    ----------
    defline : str
        The textual representation of definition line which
        lead to the collection of attributes stored therein.
    """
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

    def values(self):
        return self._mapping.values()

    def __len__(self):
        return len(self._mapping)

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

    def __dir__(self):
        base = set(dir(super(FastaHeader, self)))
        keys = set(self._mapping.keys())
        return list(base | keys)


class UnparsableDeflineError(ValueError):
    """Indicate that a definition line could
    not be parsed by the current parser.
    """
    pass


class DefLineParserBase(object):
    """A base class for definition line parser, providing
    some common machinery in the :meth:`__call__` method.

    This class requires that the subclass provide a :meth:`parse`
    method to convert a string (without a ">" prefix). and return
    an instance of :class:`.OrderedDict`.
    """
    def __call__(self, defline):
        if defline.startswith(">"):
            defline = defline[1:]
        return self._make_header(self.parse(defline), defline)

    def _make_header(self, mapping, defline):
        return FastaHeader(mapping, defline)

    def parse(self, defline):
        raise NotImplementedError()


class CallableDefLineParser(DefLineParserBase):
    """Parse a definition line using an arbitrary
    provided :class:`.Callable`.

    Attributes
    ----------
    callable_obj : Callable
    """
    def __init__(self, callable_obj):
        if not callable(callable_obj):
            raise TypeError("Must provide a callable object.")
        self.callable_obj = callable_obj

    def parse(self, defline):
        result = self.callable_obj(defline)
        if isinstance(result, Mapping):
            return OrderedDict(result.items())
        elif isinstance(result, SequenceABC):
            return OrderedDict(enumerate(result))
        else:
            raise UnparsableDeflineError(
                "Failed to parse {!r} using {!r} (return value: {!r})".format(
                    defline, self, result))


def split_on_whitespace(string):
    return re.split("\s+", string)


space_delim_parser = CallableDefLineParser(split_on_whitespace)


class RegexDefLineParser(DefLineParserBase):
    """Parse a definition line using a provided
    regular expression.

    If :attr:`pattern` contains named groups, this will
    produce a :class:`OrderedDict` where the group names
    are keys, otherwise their index will be used as the
    key.

    Attributes
    ----------
    pattern : re.pattern
    """
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
            raise UnparsableDeflineError(
                "Failed to parse {!r} using {!r}".format(
                    defline, self))


class DispatchingDefLineParser(DefLineParserBase):
    """Attempt to parse a definition line using one of many
    :class:`DeflineParserBase` instances.

    Attributes
    ----------
    graceful : bool
        If :attr:`graceful`, failure to parse a defline
        with all of :attr:`parsers` results in a :class:`FastaHeader`
        no properties other than *defline*
    header_type : type
        The type to use to contain the result of parsing.
        Defaults to :class:`FastaHeader`.
    parsers : list
        A :class:`list` of :class:`DeflineParserBase` instances
    """
    def __init__(self, parsers, graceful=False, header_type=FastaHeader):
        self.parsers = list(parsers)
        self.graceful = graceful
        self.header_type = header_type

    def _make_header(self, mapping, defline):
        return self.header_type(mapping, defline)

    def parse(self, defline):
        for parser in self.parsers:
            try:
                parts = parser.parse(defline)
                return parts
            except UnparsableDeflineError:
                continue
        else:
            if self.graceful:
                return FastaHeader({}, defline)
            else:
                raise UnparsableDeflineError(
                    "Failed to parse {:r} using {:r}".format(
                        defline, self))


uniprot_regex = r"(?P<db>[a-z^\|]+)\|(?P<accession>[A-Z0-9\-]+)\|(?P<name>\S*)(?:\s(?P<description>.+))?"
uniprot_parser = RegexDefLineParser(uniprot_regex)

partial_uniprot_regex = r"(?P<accession>[A-Z0-9]+)\|(?P<name>\S+)"
partial_uniprot_parser = RegexDefLineParser(partial_uniprot_regex)

default_parser = DispatchingDefLineParser([uniprot_parser, partial_uniprot_parser, space_delim_parser])


class FastaFileParser(object):
    def __init__(self, path, defline_parser=default_parser):
        self.state = "defline"
        self.handle = opener(path)
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
            line = line.strip()
            if self.state == 'defline':
                if line.startswith(">"):
                    self.defline = line[1:]
                    self.state = "sequence"
                else:
                    continue
            else:
                if not re.match(r"^(\s+|>)", line):
                    self.sequence_chunks.append(line)
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
                yield self.process_result(
                    {
                        "name": self.defline_parser(self.defline),
                        "sequence": ''.join(self.sequence_chunks)
                    })
            except Exception as e:
                pass


class ProteinFastaFileParser(FastaFileParser):

    def __init__(self, path, defline_parser=default_parser):
        super(ProteinFastaFileParser, self).__init__(path, defline_parser)

    def process_result(self, d):
        p = ProteinSequence(d['name'], d['sequence'], annotations=dict(d['name']))
        return p


class FastaFileWriter(object):

    def __init__(self, handle):
        self.handle = handle

    def write(self, defline, sequence):
        defline = str(defline)
        if not defline.startswith(">"):
            defline = ">" + defline
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
