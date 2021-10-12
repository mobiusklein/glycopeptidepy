import itertools
from collections import defaultdict

from glycopeptidepy.enzyme import Protease
from glycopeptidepy.structure.residue import UnknownAminoAcidException
from glycopeptidepy.structure.modification import SequenceLocation
from glycopeptidepy.structure.sequence import PeptideSequence, list_to_sequence
from glycopeptidepy.structure.parser import sequence_tokenizer_respect_sequons, sequence_tokenizer


def pair_rotate(sequence):
    """Invert each token pair.

    ABCD -> BADC

    Parameters
    ----------
    sequence : iterable

    Returns
    -------
    list
    """
    chunks = []
    gen = iter(sequence)
    while True:
        chunk = []
        try:
            chunk = [next(gen)]
            chunk.append(next(gen))
            chunks.append(chunk)
        except StopIteration:
            if len(chunk) > 0:
                chunks.append(chunk)
            break
    rev_seq = []
    for chunk in reversed(chunks):
        rev_seq.extend(chunk)
    return rev_seq


def edit_distance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2 + 1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]


def longest_common_substring(seq1, seq2):
    n = len(seq1)
    m = len(seq2)
    solution_matrix = [[0 for i in range(m)] for j in range(n)]
    # track the length of the current longest common subsequence
    z = 0
    result = []
    for i in range(n):
        for j in range(m):
            if seq1[i] == seq2[j]:
                # match the start of one string
                if i == 0 or j == 0:
                    solution_matrix[i][j] = 1
                else:
                    solution_matrix[i][j] = solution_matrix[i - 1][j - 1] + 1
                # new longest substring
                if solution_matrix[i][j] > z:
                    z = solution_matrix[i][j]
                    result = list(seq1[i - z + 1:i + 1])
                # extending existing solution (may be invalid)
                elif solution_matrix[i][j] == z:
                    result.extend(list(seq1[i - z + 1:i + 1]))
            else:
                solution_matrix[i][j] = 0
    return result


def reverse_preserve_sequon(sequence, prefix_len=0, suffix_len=1, peptide_type=None, known_sites=None):
    if peptide_type is None:
        peptide_type = PeptideSequence
    if isinstance(sequence, PeptideSequence):
        sequence = str(sequence)
    original = peptide_type(sequence)
    sequence_tokens = sequence_tokenizer_respect_sequons(sequence, known_sites=known_sites)
    pref = sequence_tokens[:prefix_len]
    if suffix_len == 0:
        suf = []
        body = sequence_tokens[prefix_len:]
    else:
        suf = sequence_tokens[-suffix_len:]
        body = sequence_tokens[prefix_len:-suffix_len]
    body = body[::-1]
    rev_sequence = (list_to_sequence(list(pref) + list(body) + list(suf)))
    if str(list_to_sequence(sequence_tokens)) == str(rev_sequence):
        rot_body = pair_rotate(body)
        rev_sequence = (list_to_sequence(pref + list(rot_body) + suf))
    rev_sequence.n_term = original.n_term.clone()
    rev_sequence.c_term = original.c_term.clone()
    if original.glycan:
        rev_sequence.glycan = original.glycan.clone(propogate_composition_offset=False)
    return rev_sequence


def reverse_sequence(sequence, prefix_len=0, suffix_len=1, peptide_type=None, known_sites=None):
    if peptide_type is None:
        peptide_type = PeptideSequence
    if isinstance(sequence, PeptideSequence):
        sequence = str(sequence)
    original = peptide_type(sequence)
    sequence_tokens = sequence_tokenizer(sequence)[0]
    pref = sequence_tokens[:prefix_len]
    if suffix_len == 0:
        suf = []
        body = sequence_tokens[prefix_len:]
    else:
        suf = sequence_tokens[-suffix_len:]
        body = sequence_tokens[prefix_len:-suffix_len]
    body = body[::-1]
    rev_sequence = (list_to_sequence(list(pref) + list(body) + list(suf)))
    if str(list_to_sequence(sequence_tokens)) == str(rev_sequence):
        rot_body = pair_rotate(body)
        rev_sequence = (list_to_sequence(pref + list(rot_body) + suf))
    rev_sequence.n_term = original.n_term.clone()
    rev_sequence.c_term = original.c_term.clone()
    if original.glycan:
        rev_sequence.glycan = original.glycan.clone(propogate_composition_offset=False)
    return rev_sequence


class LimitedCrossproduct(object):
    def __init__(self, collections, max_depth=4):
        self.collections = collections
        self.max_depth = max_depth
        self.size = len(collections)

        self.precomputed = self.compose_iterative()
        self.iterator = iter(self.precomputed)

    def compose(self):
        i = 0
        current = []
        depth = 0
        for res in self.compose_inner(i, current, depth):
            yield res

    def compose_inner(self, layer_index, current, depth):
        if layer_index == self.size:
            yield current, depth
            return
        layer = self.collections[layer_index]
        for val in layer:
            new_depth = depth + (val is not None)
            if new_depth <= self.max_depth:
                for res in self.compose_inner(layer_index + 1, current + [val], new_depth):
                    yield res
            else:
                continue

    def compose_iterative(self):
        previous = [([], 0)]
        layer_index = 0
        for layer in self.collections:
            current = []
            for rec, depth in previous:
                for val in layer:
                    new_depth = depth + (val is not None)
                    if new_depth <= self.max_depth:
                        current.append((rec + [val], new_depth))
            previous = current
        return previous

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.iterator)




class ModificationSiteAssignmentCombinator(object):
    def __init__(self, variable_site_map, max_modifications=4):
        self.modification_to_site = variable_site_map
        self.site_to_modification = self.transpose_sites()
        self.max_modifications = max_modifications

    def transpose_sites(self):
        """Given a dictionary mapping between modification names and
        an iterable of valid sites, create a dictionary mapping between
        modification names and a list of valid sites plus the constant `None`

        Returns
        -------
        dict
        """
        sites = defaultdict(list)
        for mod, varsites in self.modification_to_site.items():
            for site in varsites:
                sites[site].append(mod)
        for site in list(sites):
            sites[site].append(None)
        return sites

    def _remove_empty_sites(self, site_mod_pairs):
        return [sm for sm in site_mod_pairs if sm[1] is not None]

    def assign(self):
        sites = list(self.site_to_modification.keys())
        choices = list(self.site_to_modification.values())
        for selected in LimitedCrossproduct(choices, max_depth=self.max_modifications):
            site_mod_pairs = zip(sites, selected)
            assignment = self._remove_empty_sites(site_mod_pairs)
            if assignment > self.max_modifications:
                continue
            yield assignment

    def __iter__(self):
        return self.assign()


class PeptidoformGenerator(object):
    def __init__(self, constant_modifications, variable_modifications, max_variable_modifications=None):
        if max_variable_modifications is None:
            max_variable_modifications = 4
        self.constant_modifications = list(constant_modifications)
        self.variable_modifications = list(variable_modifications)
        self.max_variable_modifications = max_variable_modifications

        (self.n_term_modifications,
         self.c_term_modifications,
         self.variable_modifications) = self.split_terminal_modifications(self.variable_modifications)

    @staticmethod
    def split_terminal_modifications(modifications):
        """Group modification rules into three classes, N-terminal,
        C-terminal, and Internal modifications.

        A modification rule can be assigned to multiple groups if it
        is valid at multiple sites.

        Parameters
        ----------
        modifications : Iterable of ModificationRule
            The modification rules

        Returns
        -------
        n_terminal: list
            list of N-terminal modification rules
        c_terminal: list
            list of C-terminal modification rules
        internal: list
            list of Internal modification rules
        """
        n_terminal = []
        c_terminal = []
        internal = []

        for mod in modifications:
            n_term = mod.n_term_targets
            if n_term and all([t.amino_acid_targets is None for t in n_term]):
                n_term_rule = mod.clone(n_term)
                mod = mod - n_term_rule
                n_terminal.append(n_term_rule)
            c_term = mod.c_term_targets
            if c_term and all([t.amino_acid_targets is None for t in c_term]):
                c_term_rule = mod.clone(c_term)
                mod = mod - c_term_rule
                c_terminal.append(c_term_rule)
            if (mod.targets):
                internal.append(mod)

        return n_terminal, c_terminal, internal

    def prepare_peptide(self, sequence):
        if not isinstance(sequence, PeptideSequence):
            return PeptideSequence(str(sequence))
        return sequence

    def terminal_modifications(self, sequence, protein_n_term=False, protein_c_term=False):
        n_term_modifications = [
            mod for mod in self.n_term_modifications if mod.find_valid_sites(
                sequence, protein_n_term=protein_n_term)]
        c_term_modifications = [
            mod for mod in self.c_term_modifications if mod.find_valid_sites(
                sequence, protein_c_term=protein_c_term)]
        # the case with unmodified termini
        n_term_modifications.append(None)
        c_term_modifications.append(None)
        return (n_term_modifications, c_term_modifications)

    def apply_fixed_modifications(self, sequence, protein_n_term=False, protein_c_term=False):
        has_fixed_n_term = False
        has_fixed_c_term = False

        for mod in self.constant_modifications:
            for site in mod.find_valid_sites(sequence, protein_n_term=protein_n_term,
                                             protein_c_term=protein_c_term):
                if site == SequenceLocation.n_term or site == SequenceLocation.protein_n_term:
                    has_fixed_n_term = True
                elif site == SequenceLocation.c_term or site == SequenceLocation.protein_c_term:
                    has_fixed_c_term = True
                sequence.add_modification(site, mod.name)
        return has_fixed_n_term, has_fixed_c_term

    def modification_sites(self, sequence):
        variable_sites = {
            mod.name: set(
                mod.find_valid_sites(sequence)) for mod in self.variable_modifications}
        modification_sites = ModificationSiteAssignmentCombinator(variable_sites)
        return modification_sites

    def apply_variable_modifications(self, sequence, assignments, n_term, c_term):
        n_variable = 0
        result = sequence.clone()
        if n_term is not None:
            result.n_term = n_term
            n_variable += 1
        if c_term is not None:
            result.c_term = c_term
            n_variable += 1
        for site, mod in assignments:
            if mod is not None:
                result.add_modification(site, mod)
                n_variable += 1
        return result, n_variable

    def generate_peptidoforms(self, sequence, protein_n_term=False, protein_c_term=False):
        try:
            sequence = self.prepare_peptide(sequence)
        except UnknownAminoAcidException:
            return
        (n_term_modifications,
         c_term_modifications) = self.terminal_modifications(
             sequence, protein_n_term=protein_n_term, protein_c_term=protein_c_term)

        (has_fixed_n_term,
         has_fixed_c_term) = self.apply_fixed_modifications(
             sequence, protein_n_term=protein_n_term, protein_c_term=protein_c_term)

        if has_fixed_n_term:
            n_term_modifications = [None]
        if has_fixed_c_term:
            c_term_modifications = [None]

        modification_sites = self.modification_sites(sequence)

        for n_term, c_term in itertools.product(n_term_modifications, c_term_modifications):
            for assignments in modification_sites:
                if len(assignments) > self.max_variable_modifications:
                    continue
                yield self.apply_variable_modifications(
                    sequence, assignments, n_term, c_term)

    def __call__(self, peptide, protein_n_term=False, protein_c_term=False):
        return self.generate_peptidoforms(
            peptide, protein_n_term=protein_n_term, protein_c_term=protein_c_term)

    @classmethod
    def peptidoforms(cls, sequence, constant_modifications, variable_modifications, max_variable_modifications=4):
        inst = cls(constant_modifications, variable_modifications,
                   max_variable_modifications)
        return inst(sequence)


class ProteinDigestor(object):

    def __init__(self, protease, constant_modifications=None, variable_modifications=None,
                 max_missed_cleavages=2, min_length=6, max_length=60, semispecific=False,
                 max_variable_modifications=None):
        if constant_modifications is None:
            constant_modifications = []
        if variable_modifications is None:
            variable_modifications = []
        self.protease = self._prepare_protease(protease)
        self.constant_modifications = constant_modifications
        self.variable_modifications = variable_modifications
        self.peptidoform_generator = PeptidoformGenerator(
            self.constant_modifications,
            self.variable_modifications,
            max_variable_modifications=max_variable_modifications)
        self.max_missed_cleavages = max_missed_cleavages
        self.min_length = min_length
        self.max_length = max_length
        self.semispecific = semispecific
        self.max_variable_modifications = max_variable_modifications

    def _prepare_protease(self, protease):
        if isinstance(protease, Protease):
            pass
        elif isinstance(protease, basestring):
            protease = Protease(protease)
        elif isinstance(protease, (list, tuple)):
            protease = Protease.combine(*protease)
        return protease

    def cleave(self, sequence):
        seqs = self.protease.cleave(
            sequence, self.max_missed_cleavages,
            min_length=self.min_length, max_length=self.max_length,
            semispecific=self.semispecific)
        for peptide, start, end, missed in seqs:
            if "X" in peptide:
                continue
            yield peptide, start, end, missed

    def _prepare_protein(self, protein):
        if isinstance(protein, PeptideSequence):
            return protein
        return PeptideSequence(str(protein))

    def digest(self, protein):
        sequence = self._prepare_protein(protein)
        for peptide, start, end, n_missed_cleavages in self.cleave(sequence):
            if end - start > self.max_length:
                continue
            for inst in self.modify_string(peptide):
                inst.count_missed_cleavages = n_missed_cleavages
                inst.start_position = start
                inst.end_position = end
                yield inst

    def modify_string(self, peptide):
        for modified_peptide, n_variable_modifications in self.peptidoform_generator(peptide):
            modified_peptide.count_variable_modifications = n_variable_modifications
            yield modified_peptide


try:
    _ModificationSiteAssignmentCombinator = ModificationSiteAssignmentCombinator
    has_c = True
    from glycopeptidepy._c.algorithm import ModificationSiteAssignmentCombinator, LimitedCrossproduct
except ImportError:
    has_c = False
