import re
import textwrap
import warnings

from io import BytesIO
from collections import OrderedDict, defaultdict
try:
    from collections.abc import Mapping, Sequence as SequenceABC
except ImportError:
    from collections import Mapping, Sequence as SequenceABC

from six import text_type

from glycopeptidepy.structure.residue import UnknownAminoAcidException, symbol_to_residue
from glycopeptidepy.structure.sequence import ProteinSequence
from glycopeptidepy.structure.parser import parse_simple
from glycopeptidepy.utils.sequence_tree import SuffixTree
from glypy.utils.base import opener

from .cv.peff import peff_cv_term


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
        Mapping.__init__(self)
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

    def startswith(self, prefix):
        '''Check if the definition line starts with `prefix`.

        If `prefix` starts with ">", trim it.

        Parameters
        ----------
        prefix: str
            The prefix string to test for

        Returns
        -------
        bool
        '''
        if self.defline is None:
            return False
        if prefix[0] == ">"  and self.defline[0] != ">":
            prefix = prefix[1:]
        return self.defline.startswith(prefix)

    def endswith(self, suffix):
        '''Check if the definition line ends with `suffix`.

        Parameters
        ----------
        suffix: str
            The suffix string to test for

        Returns
        -------
        bool
        '''
        if self.defline is None:
            return False
        return self.defline.endswith(suffix)


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
    has_feature_index = re.compile(r"^\(?(\d+):")

    def __init__(self, validate=True):
        self.validate = validate

    def _make_header(self, mapping, defline):
        return PEFFFastaHeader(mapping, defline)

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
        value = value.strip()
        feature_id_match = self.has_feature_index.search(value)
        if feature_id_match:
            feature_id = int(feature_id_match.group(1))
            value = self.has_feature_index.sub('', value)
        else:
            feature_id = None
        if "|" in value:
            value = self.split_pipe_separated_tuple(value)
            result = []
            for i, v in enumerate(value):
                result.append(self._coerce_value(key, v, i))
            return PEFFFeature(*result, feature_type=key, id=feature_id)
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

    def parse(self, defline):
        if self.validate:
            match = self.detect_pattern.match(defline)
            if not match:
                raise UnparsableDeflineError(
                    "Failed to parse {!r} using {!r}".format(
                        defline, self))
        storage = OrderedDict()
        prefix = None
        db_uid = None
        if defline.startswith(">"):
            defline = defline[1:]
        try:
            prefix, defline = defline.split(":", 1)
            db_uid, defline = defline.split(" ", 1)
        except ValueError:
            raise UnparsableDeflineError(defline)
        storage['Prefix'] = prefix
        storage['Tag'] = db_uid
        kv_pattern = re.compile(r"\\(?P<key>\S+)=(?P<value>.+?)(?:\s(?=\\)|$)")
        for key, value in kv_pattern.findall(defline):
            if self.validate:
                try:
                    peff_cv_term(key, strict=True)
                except KeyError:
                    warnings.warn("PEFF Key {} not recognized".format(key))
            if not (value.startswith("(") or " (" in value):
                storage[key] = self.coerce_types(key, value)
            else:
                # multi-value
                storage[key] = [
                    self.coerce_types(key, v) for v in self.extract_parenthesis_list(value)]
        return storage


class PEFFFeature(SequenceABC):
    def __init__(self, *fields, **kwargs):
        SequenceABC.__init__(self)
        self.fields = tuple(fields)
        self.id = kwargs.get('id')
        self.feature_type = kwargs.get("feature_type")

    def __eq__(self, other):
        return tuple(self) == tuple(other)

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(str(self))

    def __getitem__(self, i):
        return self.fields[i]

    def __len__(self):
        return len(self.fields)

    def __repr__(self):
        return repr(tuple(self))

    def __str__(self):
        return "(%s%s)" % (
            '%r:' % self.id if self.id is not None else '',
            '|'.join(map(str, self)), )


class PEFFFastaHeader(FastaHeader):
    def __init__(self, mapping, original=None):
        self._feature_map = {}
        FastaHeader.__init__(self, mapping, original=original)
        self._init_feature_id_map()

    def _init_feature_id_map(self):
        for _, values in self.items():
            if isinstance(values, list):
                for feature in values:
                    if hasattr(feature, 'id'):
                        feature_id = feature.id
                        if feature_id is not None:
                            self._feature_map[feature_id] = feature
            else:
                feature = values
                if hasattr(feature, 'id'):
                    feature_id = feature.id
                    if feature_id is not None:
                        self._feature_map[feature_id] = feature

    def _make_defline(self):
        template = text_type(
            ">{self.Prefix}:{self.Tag} {rest}")
        rest_parts = []
        for key, value in self.items():
            if key in ("Prefix", "Tag"):
                continue
            part = '\\{key}={fmt_value}'
            if isinstance(value, ((int, float, text_type, str))):
                fmt_value = str(value)
            else:
                fmt_value = ''.join(map(str, value))
            rest_parts.append(part.format(key=key, fmt_value=fmt_value))
        return template.format(self=self, rest=' '.join(rest_parts))


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
    return re.split(r"\s+", string)


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
        if isinstance(pattern, str):
            self.pattern = re.compile(pattern)
        else:
            self.pattern = pattern
        self._is_keyed = "?P<" in self.pattern.pattern

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

default_parser = DispatchingDefLineParser(
    [uniprot_parser, partial_uniprot_parser, peff_parser, space_delim_parser])


class FastaFileReader(object):
    '''A base class for parsing amino acid Fasta Files. Supports
    sequential access through iteration, and random access if :attr:`index`
    is :const:`True`.
    '''
    def __init__(self, path, defline_parser=default_parser, encoding='utf8', index=False, case_sensitive=False):
        self.state = "defline"
        self.handle = opener(path, 'rb')
        self.defline_parser = defline_parser
        self.encoding = encoding
        self.case_sensitive = case_sensitive
        self._generator = None
        self._index = None
        if index:
            self._build_index()

    def _build_index(self):
        try:
            pos = self.handle.tell()
            self.handle.seek(0)
        except (AttributeError, ValueError, TypeError) as ex:
            print(ex)
            pos = None
        indexer = FastaIndexer(encoding=self.encoding, defline_parser=self.defline_parser)
        index = indexer.build_index(self.handle)
        if pos is not None:
            self.handle.seek(pos)
        self._index = index

    @property
    def index(self):
        return self._index

    def __len__(self):
        if self.index is not None:
            return len(self.index)
        else:
            raise TypeError("Cannot get the length of an un-indexed FASTA")

    def __getitem__(self, key):
        if self.index is not None:
            try:
                pos = self.handle.tell()
                self.handle.seek(0)
            except (AttributeError, ValueError, TypeError):
                pos = None
            if isinstance(key, int):
                start, end = tuple(self.index.values())[key]
            elif isinstance(key, slice):
                return [self[i] for i in range(*key.indices(len(self.index)))]
            else:
                start, end = self.index[key]
            self.handle.seek(start)
            data = self.handle.read(end - start)
            value = next(self._parse_lines(BytesIO(data)))
            if pos is not None:
                self.handle.seek(pos)
            return value
        else:
            raise ValueError("No index was built")

    def process_result(self, d):
        return d

    def __iter__(self):
        return self

    def __next__(self):
        if self._generator is None:
            if self.index is None:
                self._generator = self._parse_lines(self.handle)
            else:
                self._generator = (self[k] for k in self.index)
        return next(self._generator)

    def reset(self):
        """Restart the iterator from the beginning of the sequence list.

        Attempts to seek to the beginning of the file. If the input file
        is not seekable, this will throw an exception.
        """
        self.handle.seek(0)
        self._generator = None

    def next(self):
        return self.__next__()

    def _parse_lines(self, stream):
        sequence_pattern = re.compile(r"^(\s+|>)")
        sequence_chunks = []
        defline = None
        state = 'defline'
        for line in stream:
            line = line.strip()
            try:
                line = line.decode(self.encoding)
            except AttributeError:
                pass
            if state == 'defline':
                if line.startswith(">"):
                    defline = line[1:]
                    state = "sequence"
                else:
                    continue
            else:
                if not sequence_pattern.match(line):
                    sequence_chunks.append(line)
                else:
                    if defline is not None:
                        try:
                            value = ''.join(sequence_chunks).replace(" ", "")
                            if not self.case_sensitive:
                                value = value.upper()
                            yield self.process_result({
                                "name": self.defline_parser(defline),
                                "sequence": value
                            })
                        except KeyError as e:
                            warnings.warn(str(e))
                    sequence_chunks = []
                    defline = None
                    state = 'defline'
                    if line[0] == '>':
                        defline = re.sub(r"[\n\r]", "", line[1:])
                        state = "sequence"

        if sequence_chunks:
            try:
                value = ''.join(sequence_chunks).replace(" ", "")
                if not self.case_sensitive:
                    value = value.upper()
                yield self.process_result(
                    {
                        "name": self.defline_parser(defline),
                        "sequence": value
                    })
            except KeyError as e:
                warnings.warn(str(e))


FastaFileParser = FastaFileReader


class ProteinFastaFileReader(FastaFileReader):
    '''Parses Fasta Files into :class:`~.ProteinSequence` instances. Supports
    sequential access, as well as random access when :attr:`index` is :const:`True`.

    Attributes
    ----------
    replace_unknown: :class:`bool`
        Whether to replace amino acids which are not recognized with an 'X'
        instead of throwing an error. If not set, the parser will replace
        the amino acids and emit a warning instead.
    case_sensitive: :class:`bool`
        Whether to convert all characters in the amino acid sequence to upper case
        prior to parsing.
    '''
    def __init__(self, path, defline_parser=default_parser, encoding='utf8', index=False,
                 replace_unknown=None, case_sensitive=False):
        super(ProteinFastaFileReader, self).__init__(
            path, defline_parser, encoding, index=index, case_sensitive=case_sensitive)
        self.replace_unknown = replace_unknown
        self._replace_amino_acid_pattern = None

    def _make_replace_amino_acid_pattern(self):
        pattern = re.compile(r"[^%s]" % (
            ''.join(symbol_to_residue.keys()), ))
        return pattern

    def _replace_unknown_amino_acids(self, sequence):
        if self._replace_amino_acid_pattern is None:
            self._replace_amino_acid_pattern = self._make_replace_amino_acid_pattern()
        pattern = self._replace_amino_acid_pattern
        return pattern.sub("X", sequence)

    def process_result(self, d):
        try:
            p = ProteinSequence(
                d['name'], d['sequence'], parser_funciton=parse_simple, annotations=dict(d['name']))
        except UnknownAminoAcidException:
            if self.replace_unknown is None or self.replace_unknown:
                if self.replace_unknown is None:
                    warnings.warn("Replacing unknown amino acids in %s" % d['name'])
                d['sequence'] = self._replace_unknown_amino_acids(d['sequence'])
                p = ProteinSequence(
                    d['name'], d['sequence'], parser_funciton=parse_simple,
                    annotations=dict(d['name']))
            else:
                raise
        return p


ProteinFastaFileParser = ProteinFastaFileReader


class FastaFileWriter(object):

    def __init__(self, handle, encoding='utf8'):
        self.handle = handle
        self.encoding = encoding
        try:
            if 'b' not in self.handle.mode:
                raise ValueError(
                    "FASTA files writer must be opened in binary mode! Make sure to open the file with 'wb'.")
        except (AttributeError, TypeError):
            pass

    def write(self, defline, sequence):
        if isinstance(defline, FastaHeader):
            defline = str(defline).encode(self.encoding)
        elif isinstance(defline, text_type):
            defline = defline.encode(self.encoding)
        if not defline.startswith(b">"):
            defline = b">" + defline
        if isinstance(sequence, text_type):
            sequence = sequence.encode(self.encoding)
        self.handle.write(defline)
        self.handle.write(b"\n")
        self.handle.write(sequence)
        self.handle.write(b"\n\n")

    def writelines(self, iterable):
        for defline, seq in iterable:
            self.write(defline, seq)


class ProteinFastaFileWriter(FastaFileWriter):
    def __init__(self, handle, encoding='utf8', line_length=80):
        super(ProteinFastaFileWriter, self).__init__(handle, encoding)
        self.line_length = line_length

    def write(self, protein): # pylint: disable=arguments-differ
        defline = str(protein.name)
        seq = '\n'.join(textwrap.wrap(protein.get_sequence(), self.line_length))
        super(ProteinFastaFileWriter, self).write(defline, seq)

    def writelines(self, iterable):
        for protein in iterable:
            self.write(protein)


def read(f, defline_parser=default_parser, reader_type=ProteinFastaFileParser, **kwargs):
    return reader_type(f, defline_parser=defline_parser, **kwargs)


class PEFFHeaderBlock(Mapping):
    def __init__(self, storage=None, block_type=None):
        Mapping.__init__(self)
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

    def __str__(self):
        return "%s" % '\n'.join('# %s=%s' % kv for kv in self.items())

    def __repr__(self):
        return "{self.__class__.__name__}({self._storage})".format(self=self)


class PEFFReader(ProteinFastaFileReader):
    def __init__(self, path, defline_parser=PEFFDeflineParser(False), index=False):
        super(PEFFReader, self).__init__(
            opener(path, 'rb'), defline_parser, encoding='ascii', index=index)
        try:
            if 'b' not in self.handle.mode:
                raise ValueError(
                    "PEFF files must be opened in binary mode! Make sure to open the file with 'rb'.")
        except (AttributeError, TypeError):
            pass
        self.version = (0, 0)
        self.blocks = []
        self.comments = []
        self.number_of_entries = 0
        self._parse_header()

    def _parse_header(self):
        offset = 0
        line = self.handle.readline()
        offset += len(line)
        line = line.decode('ascii')
        if not line.startswith("# PEFF"):
            raise ValueError("Not a PEFF File")
        self.version = tuple(map(int, line.strip()[7:].split(".")))
        in_header = True
        current_block = defaultdict(list)
        while in_header:
            line = self.handle.readline()
            offset += len(line)
            line = line.decode('ascii')
            if not line.startswith("#"):
                in_header = False
                self.handle.seek(offset)
            line = line.strip()[2:]
            if '=' in line:
                key, value = line.split("=", 1)
                if key == "GeneralComment":
                    self.comments.append(value)
                else:
                    current_block[key].append(value)
            if line.startswith("//"):
                if current_block:
                    self.blocks.append(PEFFHeaderBlock(OrderedDict((k, v if len(v) > 1 else v[0])
                                                                   for k, v in current_block.items())))
                current_block = defaultdict(list)
        number_of_entries = 0
        for block in self.blocks:
            try:
                number_of_entries += int(block['NumberOfEntries'])
            except KeyError:
                pass
        self.number_of_entries = number_of_entries


class PEFFWriter(ProteinFastaFileWriter):
    version = (1, 0)

    def __init__(self, handle, line_length=80):
        super(PEFFWriter, self).__init__(handle, 'ascii', line_length)

    def write_header(self, blocks, comments=None):
        if comments is None:
            comments = []
        if blocks is None:
            blocks = []
        self.handle.write("# PEFF %d.%d\n" % self.version)
        for comment in comments:
            self.handle.write("# GeneralComment=%s\n" % comment)
        self.handle.write("# //\n")
        for block in blocks:
            self.handle.write(str(block))
            self.handle.write("\n# //\n")


class FastaIndexer(object):
    '''Encapsulates the process of building an index of byte-level offsets
    and entry-level metadata for a Fasta file.
    '''
    def __init__(self, read_size=1000000, encoding='utf8', defline_parser=default_parser):
        self.read_size = int(read_size)
        self.encoding = encoding
        self.defline_parser = defline_parser

    def _chunk_iterator(self, stream):
        delim = b"\n>"
        read_size = self.read_size
        f = stream
        buff = f.read(read_size)
        started_with_with_delim = buff.startswith(delim)
        parts = buff.split(delim)
        tail = parts[-1]
        front = parts[:-1]
        i = 0
        for part in front:
            i += 1
            if part == b'':
                continue
            if i == 1:
                if started_with_with_delim:
                    yield delim + part
                else:
                    yield part
            else:
                yield delim + part
        running = True
        while running:
            buff = f.read(read_size)
            if not buff:
                running = False
                buff = tail
            else:
                buff = tail + buff
            parts = buff.split(delim)
            tail = parts[-1]
            front = parts[:-1]
            for part in front:
                yield delim + part
        yield delim + tail

    def _generate_offsets(self, stream):
        i = 0
        defline_pattern = re.compile(br"\s*>([^\n\r]+)[\n\r]+")
        for line in self._chunk_iterator(stream):
            match = defline_pattern.match(line)
            if match:
                yield i, match.group(1)
            i += len(line)
        yield i, None

    def _handle_defline(self, defline):
        defline = defline.decode(self.encoding)
        if self.defline_parser:
            try:
                defline = self.defline_parser(defline)
            except UnparsableDeflineError:
                warnings.warn("Could not parse defline %r" % (defline,))
        return defline

    def build_index(self, stream):
        '''Build an random access and metadata index over the stream
        of Fasta File entries.

        Parameters
        ----------
        stream: file-like

        Returns
        -------
        :class:`FastaIndex`
        '''
        index = OrderedDict()
        g = self._generate_offsets(stream)
        last_offset = 0
        last_defline = None
        for offset, defline in g:
            if last_defline is not None:
                index[self._handle_defline(last_defline)] = (last_offset, offset)
            last_defline = defline
            last_offset = offset
        assert last_defline is None
        return FastaIndex(index)


class FastaIndex(object):
    '''A :class:`~.Mapping`-like object which stores the byte offsets
    and entry metadata for a Fasta file.
    '''
    def __init__(self, mapping):
        self.mapping = mapping
        self._suffix = None
        self._index_map = dict()

    def __repr__(self):
        template = "{self.__class__.__name__}({size:d})"
        return template.format(self=self, size=len(self))

    def __getitem__(self, key):
        return self.mapping[key]

    def __iter__(self):
        return iter(self.mapping)

    def __contains__(self, key):
        return key in self.mapping

    def keys(self):
        return self.mapping.keys()

    def values(self):
        return self.mapping.values()

    def items(self):
        return self.mapping.items()

    def __len__(self):
        return len(self.mapping)

    def _build_suffix_tree(self):
        sfx_tree = SuffixTree()
        for k in self:
            sfx_tree.add_ngram(str(k), k)
        return sfx_tree

    def suffix(self, key):
        '''Search for headers where ``key`` is a suffix of
        the header.

        This method builds a suffix tree which may be very
        expensive in memory

        Parameters
        ----------
        key: :class:`str`
            The suffix to search for

        Returns
        -------
        :class:`list`
        '''
        if self._suffix is None:
            self._suffix = self._build_suffix_tree()
        return self._suffix.subsequences_of(key)

    def _build_index_for_key(self, key):
        index = defaultdict(set)
        for header in self.keys():
            try:
                value = header[key]
                index[value].add(header)
            except KeyError:
                continue
        self._index_map[key] = index

    def index_by(self, key):
        '''Build an extra index over the values of ``key``
        '''
        self._build_index_for_key(key)

    def clear_indices(self):
        '''Discard all the extra indices.
        '''
        self._index_map.clear()
        self._suffix = None

    def query(self, queries=None, **kwargs):
        '''Query the index for headers whose keys match the values passed
        through ``kwargs``

        Parameters
        ----------
        queries: :class:`dict`, optional
            A mapping of :class:`FastaHeader` keys to values to match for
            inclusion in the result set. All pairs must match for a header
            to be included in the results.
        **kwargs:
            Additional key-value pairs to search for

        Results
        -------
        :class:`list`:
            The :class:`FastaHeader` instances that matched all of the constraints
            the query.
        '''
        if not queries:
            queries = dict()
        else:
            queries = dict(queries)
        queries.update(kwargs)
        queries = queries.items()
        # The indices fully cover this query, so make it fast.
        if set(kwargs) <= set(self._index_map):
            index_matches = []
            for k, v in queries:
                try:
                    index_matches.append(self._index_map[k][v])
                except KeyError:
                    continue
            if not index_matches:
                return []
            return list(index_matches[0].intersection(*index_matches[1:]))
        else:
            # Employ partial index coverage to shrink the set of elements to
            # explore exhaustively
            index_keys = set(kwargs) & set(self._index_map)
            if index_keys:
                query_set = self.query({k: queries[k] for k in index_keys})
                for k in index_keys:
                    queries.pop(k)
            else:
                # No helpful index, so iterate over all headers.
                query_set = self.keys()

            matches = []
            # Check each entry for key-value equality, keeping those
            # cases which match on *all* keys. This isn't SQL, so we
            # cannot express complex relationships or conditionals
            for header in query_set:
                for k, v in queries:
                    try:
                        if header[k] != v:
                            break
                    except KeyError:
                        break
                else:
                    matches.append(header)
        return matches

    def prefix_filter(self, prefix):
        return [entry for entry in self if entry.startswith(prefix)]

    def suffix_filter(self, suffix):
        return [entry for entry in self if entry.endswith(suffix)]
