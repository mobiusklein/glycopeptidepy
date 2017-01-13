import re
import textwrap
import warnings
from glycopeptidepy.structure.sequence import ProteinSequence


def tryopen(path):
    if hasattr(path, 'read'):
        return path
    return open(path)


class FastaFileParser(object):
    handled_exceptions = (Exception)

    def __init__(self, path):
        self.state = "defline"
        self.handle = tryopen(path)
        self.defline = None
        self.sequence_chunks = []

    def process_result(self, d):
        return d

    def __iter__(self):
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
                                "name": self.defline,
                                "sequence": ''.join(self.sequence_chunks)
                            })
                        except self.handled_exceptions as e:
                            pass
                    self.sequence_chunks = []
                    self.defline = None
                    self.state = 'defline'
                    if line[0] == '>':
                        self.defline = re.sub(r"[\n\r]", "", line[1:])
                        self.state = "sequence"

        if len(self.sequence_chunks) > 0:
            try:
                yield self.process_result({"name": self.defline, "sequence": ''.join(self.sequence_chunks)})
            except self.handled_exceptions as e:
                pass


class ProteinFastaFileParser(FastaFileParser):

    def __init__(self, path):
        super(ProteinFastaFileParser, self).__init__(path)

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
