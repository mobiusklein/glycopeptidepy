import re
import string
import warnings

from collections import defaultdict

from .descriptors import SequenceLocation, ModificationCategory

from ..residue import AminoAcidResidue


target_string_pattern = re.compile(
    r'''(?P<amino_acid>[A-Z]+)?
        (@
            (?P<n_term>[Nn][-_][tT]erm)|
            (?P<c_term>[Cc][-_][tT]erm)
        )?''',
    re.VERBOSE)
title_cleaner = re.compile(
    r'(?P<name>.*)\s\((?P<target>.+)\)$')

is_mass_delta = re.compile(r"^Delta:")


def extract_targets_from_string(target_string):
    position_only = re.search(r"^((?:Protein )?[NnCc][-_][tT]erm)", target_string)
    if position_only:
        return ModificationTarget(None, SequenceLocation[position_only.groups()[0]])
    else:
        state = 'aa'
        amino_acids = []
        position = []

        protein_terminal_only = False
        i = 0
        n = len(target_string)
        while i < n:
            c = target_string[i]

            if c == " ":
                if state == 'aa':
                    state = "aa-space"
                elif state == 'position-begin':
                    pass
                elif state == 'protein':
                    state = 'position-begin'
                    protein_terminal_only = True
                else:
                    raise ValueError(target_string)
            elif c == "@":
                if state == 'aa':
                    state = 'position-begin'
                elif state == 'aa-space':
                    state = 'position-begin'
                else:
                    raise ValueError("Found @ not between amino acids and position")
            elif c in string.ascii_uppercase:
                if state == 'aa':
                    amino_acids.append(c)
                elif state == 'position-begin':
                    if c in "NC":
                        position.append(c)
                        state = 'position'
                    elif c == "P":
                        state = 'protein'
                elif state == 'position':
                    if c == 'T':
                        position.append(c)
                    else:
                        raise ValueError("Invalid position specification")
            elif c in string.ascii_lowercase:
                if state == 'aa':
                    amino_acids.append(c.upper())
                elif state == 'position-begin':
                    if c in 'nc':
                        position.append(c)
                        state = 'position'
                    if c == 'p':
                        state = 'protein'
                elif state == "protein":
                    if c in 'protein':
                        pass
                    else:
                        raise ValueError("Invalid terminal speicification")
                elif state == 'position':
                    position.append(c)
            elif c in '-_':
                if state == 'position':
                    position.append(c)
                else:
                    raise ValueError("Found - not in position specification")
            i += 1

        position = ''.join(position)
        position = position.lower().replace("-", "_").strip()

        amino_acids = list(map(AminoAcidResidue, amino_acids))
        if position == "":
            position = SequenceLocation.anywhere
        elif position == "n_term":
            if not protein_terminal_only:
                position = SequenceLocation.n_term
            else:
                position = SequenceLocation.protein_n_term
        elif position == "c_term":
            if not protein_terminal_only:
                position = SequenceLocation.c_term
            else:
                position = SequenceLocation.protein_c_term
        else:
            raise ValueError("Unrecognized position specification {}".format(position))
        return ModificationTarget(amino_acids, position)


def extract_targets_from_rule_string(name):
    return extract_targets_from_string(title_cleaner.search(name).groupdict()["target"])


def get_position_modifier_rules_dict(sequence, protein_n_term=False, protein_c_term=False):
    '''Labels the start and end indices of the sequence

    A convenience dictionary initializer for ease-of-access in
    :class:`ModificationTarget` site validation functions

    Parameters
    ----------
    sequence : PeptideSequence
        The sequence for which to project labels
    protein_n_term : bool
        Whether the sequence includes a protein N-terminal
    protein_c_term : bool
        Whether the sequence includes a protein C-terminal

    Return
    ------
    defaultdict
    '''
    if protein_n_term:
        n_term = (SequenceLocation.n_term, SequenceLocation.protein_n_term)
    else:
        n_term = (SequenceLocation.n_term, )
    if protein_c_term:
        c_term = (SequenceLocation.c_term, SequenceLocation.protein_c_term)
    else:
        c_term = (SequenceLocation.c_term, )
    return defaultdict(lambda: SequenceLocation.anywhere, {
        0: n_term,
        (len(sequence) - 1): c_term
    })


class ModificationTarget(object):
    '''Specifies the range of targets that a given modification can be found at.

    A single instance of :class:`ModificationTarget` describes one or more position/residue
    rules for a single type of modification, bundled together underneath a single :class:`ModificationRule`
    instance. :class:`ModificationTarget` may be extracted from Unimod data dictionaries using the
    :meth:`from_unimod_specificity` class method.

    Attributes
    ----------
    amino_acid_targets: list
        The list of amino acid residues that this rule
        permits.
    position_modifier: SequenceLocation
        The sequence-position specific filter
        this rule enforces.
    '''

    @classmethod
    def from_unimod_specificity(cls, specificity):
        amino_acid = (specificity["site"])
        if amino_acid in (["N-term", "C-term", None]):
            amino_acid = None
        else:
            amino_acid = {AminoAcidResidue(amino_acid)}
        position_modifier = specificity["position"]
        nterm = ("Any N-term", "N-term")
        prot_nterm = ("Protein N-term", )
        cterm = ("Any C-term", "C-term")
        prot_cterm = ("Protein C-term", )
        site = specificity["site"]
        classification = ModificationCategory[specificity.get("classification")]
        if position_modifier == "Anywhere":
            position_modifier = SequenceLocation.anywhere
        elif position_modifier in nterm or site in nterm:
            position_modifier = SequenceLocation.n_term
        elif position_modifier in cterm or site in cterm:
            position_modifier = SequenceLocation.c_term
        elif position_modifier in prot_nterm:
            position_modifier = SequenceLocation.protein_n_term
        elif position_modifier in prot_cterm:
            position_modifier = SequenceLocation.protein_c_term
        else:
            warnings.warn("Could not find a location for %r" % (position_modifier,))

        return cls(amino_acid, position_modifier, [classification])

    def __init__(self, site_targets=None, position_modifier=None, classification=None):
        if classification is None:
            classification = []
        self.amino_acid_targets = frozenset(map(AminoAcidResidue, site_targets)) if site_targets is not None else None
        self.position_modifier = SequenceLocation[position_modifier]
        self.classification = classification

    def valid_site(self, amino_acid=None, position_modifier=SequenceLocation.anywhere):
        '''Return if a given residue at a sequence position (N-term, C-term, None)
        is valid

        Arguments
        ---------
        amino_acid: AminoAcidResidue
            The amino acid to test
        position_modifier: SequenceLocation
            The position in the sequence to test

        Returns
        -------
        valid: bool
            Whether the match is valid
        location_specification: SequenceLocation
            If the site is valid, this will be the type of location matched,
            otherwise it will default to SequenceLocation.anywhere
        '''
        valid = False
        if not isinstance(position_modifier, (tuple, list)):
            position_modifiers = (position_modifier, )
        else:
            position_modifiers = position_modifier

        for position_modifier in position_modifiers:
            # Validate amino acid target
            valid = (self.amino_acid_targets is None) or (
                amino_acid in self.amino_acid_targets)
            valid = valid and ((self.position_modifier == SequenceLocation.anywhere) or
                            (position_modifier == self.position_modifier))

            # If the rule includes a position modifier other than anywhere
            if valid and (position_modifier != SequenceLocation.anywhere) and (
                    self.position_modifier != SequenceLocation.anywhere):
                if position_modifier == self.position_modifier:
                    if self.amino_acid_targets is None:
                        return valid, position_modifier
                    else:
                        return valid, SequenceLocation.anywhere
                else:
                    valid = False
        return valid, SequenceLocation.anywhere

    def __repr__(self):
        amino_acid_components = "{%s}" % ', '.join(map(str, self.amino_acid_targets or ()))
        rep = "{amino_acid_targets}@{position_modifier}".format(
            amino_acid_targets=amino_acid_components,
            position_modifier=self.position_modifier)
        return rep

    def __eq__(self, other):
        if self.amino_acid_targets != other.amino_acid_targets:
            return False
        elif self.position_modifier != other.position_modifier:
            return False
        return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.amino_acid_targets)

    def __len__(self):
        position_modifier_count = 0
        amino_acid_targets_count = 0
        if self.amino_acid_targets is not None:
            amino_acid_targets_count = len(self.amino_acid_targets)
        return amino_acid_targets_count + position_modifier_count

    def serialize(self):
        parts = []
        if self.amino_acid_targets is not None:
            parts.append(''.join(aa.symbol for aa in self.amino_acid_targets))
            if self.position_modifier == SequenceLocation.anywhere:
                pass
            else:
                parts.append("@")
        if self.position_modifier == SequenceLocation.n_term:
            parts.append("N-term")
        elif self.position_modifier == SequenceLocation.protein_n_term:
            parts.append("Protein N-term")
        elif self.position_modifier == SequenceLocation.c_term:
            parts.append("C-term")
        elif self.position_modifier == SequenceLocation.protein_c_term:
            parts.append("Protein C-term")
        return " ".join(parts)
