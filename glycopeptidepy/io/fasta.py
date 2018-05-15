import re
import textwrap

from collections import OrderedDict, defaultdict
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


class PEFFDeflineParser(DefLineParserBase):
    kv_pattern = re.compile(r"\\(?P<key>\S+)=(?P<value>.+?)(?:\s(?=\\)|$)")
    detect_pattern = re.compile(r"^>?\S+:\S+")

    def __init__(self, validate=True):
        self.validate = validate

    def extract_parenthesis_list(self, text):
        chunks = []
        chunk = []
        paren_level = 0
        i = 0
        n = len(text)
        while i < n:
            c = text[i]
            i += 1
            if c == "(":
                if paren_level > 0:
                    chunk.append(c)
                paren_level += 1
            elif c == ")":
                if paren_level > 1:
                    chunk.append(c)
                paren_level -= 1
                if paren_level == 0:
                    if chunk:
                        chunks.append(chunk)
                    chunk = []
            else:
                chunk.append(c)
        chunks = list(map(''.join, chunks))
        return chunks

    def split_pipe_separated_tuple(self, text):
        parts = text.split("|")
        return parts

    def coerce_types(self, key, value):
        if "|" in value:
            value = self.split_pipe_separated_tuple(value)
            result = []
            for i, v in enumerate(value):
                result.append(self._coerce_value(key, v, i))
            return tuple(result)
        else:
            return self._coerce_value(key, value, 0)

    def _coerce_value(self, key, value, index):
        try:
            return int(value)
        except ValueError:
            pass
        try:
            return float(value)
        except ValueError:
            pass
        return str(value)

    def parse(self, line):
        if self.validate:
            match = self.detect_pattern.match(line)
            if not match:
                raise UnparsableDeflineError(
                    "Failed to parse {!r} using {!r}".format(
                        line, self))
        storage = OrderedDict()
        prefix = None
        db_uid = None
        if line.startswith(">"):
            line = line[1:]
        prefix, line = line.split(":", 1)
        db_uid, line = line.split(" ", 1)
        storage['Prefix'] = prefix
        storage['DbUniqueId'] = db_uid
        kv_pattern = re.compile(r"\\(?P<key>\S+)=(?P<value>.+?)(?:\s(?=\\)|$)")
        for key, value in kv_pattern.findall(line):
            if not value.startswith("(") or " (" in value:
                storage[key] = self.coerce_types(key, value)
            else:
                # multi-value
                storage[key] = [self.coerce_types(key, v) for v in self.extract_parenthesis_list(value)]
        return storage


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


uniprot_regex = r"(?P<db>[a-z^\|]+)\|(?P<accession>[a-zA-Z0-9\-]+)\|(?P<name>\S*)(?:\s(?P<description>.+))?"
uniprot_parser = RegexDefLineParser(uniprot_regex)

partial_uniprot_regex = r"(?P<accession>[A-Z0-9\-]+)\|(?P<name>\S+)"
partial_uniprot_parser = RegexDefLineParser(partial_uniprot_regex)

peff_parser = PEFFDeflineParser()

default_parser = DispatchingDefLineParser([uniprot_parser, partial_uniprot_parser, peff_parser, space_delim_parser])


class FastaFileReader(object):
    def __init__(self, path, defline_parser=default_parser, encoding='utf8'):
        self.state = "defline"
        self.handle = opener(path, 'rb')
        self.defline = None
        self.sequence_chunks = []
        self.defline_parser = defline_parser
        self.encoding = encoding
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
            try:
                line = line.decode(self.encoding)
            except AttributeError:
                pass
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


FastaFileParser = FastaFileReader


class ProteinFastaFileReader(FastaFileReader):

    def __init__(self, path, defline_parser=default_parser, encoding='utf8'):
        super(ProteinFastaFileReader, self).__init__(path, defline_parser, encoding)

    def process_result(self, d):
        p = ProteinSequence(d['name'], d['sequence'], annotations=dict(d['name']))
        return p


ProteinFastaFileParser = ProteinFastaFileReader


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


class PEFFHeaderBlock(Mapping):
    def __init__(self, storage=None, block_type=None):
        if storage is None:
            storage = OrderedDict()
        self.block_type = block_type
        self._storage = storage

    def keys(self):
        return self._storage.keys()

    def values(self):
        return self._storage.values()

    def items(self):
        return self._storage.items()

    def __iter__(self):
        return iter(self._storage)

    def __getitem__(self, key):
        return self._storage[key]

    def __len__(self):
        return len(self._storage)

    def __contains__(self, key):
        return key in self._storage

    def __getattr__(self, key):
        if key == "_storage":
            raise AttributeError(key)
        try:
            return self._storage[key]
        except KeyError:
            raise AttributeError(key)

    def __dir__(self):
        base = set(dir(super(PEFFHeaderBlock, self)))
        keys = set(self._storage.keys())
        return list(base | keys)

    def __repr__(self):
        return "# //\n%s\n# //" % '\n'.join('# %s=%s' % kv for kv in self.items())


class PEFFReader(ProteinFastaFileReader):
    def __init__(self, path, defline_parser=PEFFDeflineParser(False)):
        super(PEFFReader, self).__init__(opener(path, 'rb'), defline_parser, encoding='ascii')
        try:
            if 'b' not in self.handle.mode:
                raise ValueError("PEFF files must be opened in binary mode! Make sure to open the file with 'rb'.")
        except AttributeError:
            pass
        self.version = (0, 0)
        self.blocks = []
        self.comments = []
        self.number_of_entries = 0
        self._parse_header()

    def _parse_header(self):
        offset = 0
        line = self.handle.readline()
        line = line.decode('ascii')
        offset += len(line)
        if not line.startswith("# PEFF"):
            raise ValueError("Not a PEFF File")
        self.version = tuple(map(int, line.strip()[7:].split(".")))
        in_header = True
        current_block = defaultdict(list)
        while in_header:
            line = self.handle.readline()
            line = line.decode('ascii')
            if not line.startswith("#"):
                in_header = False
                self.handle.seek(offset)
            offset += len(line)
            line = line.strip()[2:]
            if '=' in line:
                key, value = line.split("=", 1)
                if key == "GeneralComment":
                    self.comments.append(value)
                else:
                    current_block[key].append(value)
            if line.startswith("//"):
                if current_block:
                    self.blocks.append(PEFFHeaderBlock({k: v if len(v) > 1 else v[0]
                                                        for k, v in current_block.items()}))
                current_block = defaultdict(list)
        number_of_entries = 0
        for block in self.blocks:
            try:
                number_of_entries += int(block['NumberOfEntries'])
            except KeyError:
                pass
        self.number_of_entries = number_of_entries
