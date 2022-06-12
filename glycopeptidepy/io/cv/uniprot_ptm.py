import re

import pkg_resources

from glycopeptidepy.structure.residue import long_to_symbol

from glycopeptidepy.structure.modification import (
    ModificationRule,
    ModificationTarget,
    SequenceLocation)

from glypy import Composition


class UniProtPTMListParser(object):  # pragma: no cover

    def __init__(self, path=None):
        if path is None:
            path = pkg_resources.resource_filename("glycopeptidepy.io.cv.data", "uniprot_ptmlist.txt")
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
        tokens = re.split(r"\s+", line, maxsplit=1)
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

    def _translate_position(self, position):
        t = {
            'Anywhere': SequenceLocation.anywhere,
            'N-terminal': SequenceLocation.n_term,
            'C-terminal': SequenceLocation.c_term
        }
        return t.get(position, 'Anywhere')

    def _translate_amino_acid(self, target):
        try:
            symbol = long_to_symbol[target]
        except KeyError:
            return None
        return symbol

    def parse_entry(self):
        ptm_id = None
        accession = None
        # feature_key = None
        mass_difference = None
        correction_formula = None
        target = None
        position = None

        keywords = []
        crossref = []
        psimod_id = None
        while True:
            typecode, line = self._next_line()
            if typecode == "ID":
                ptm_id = line
            elif typecode == "AC":
                accession = line
            elif typecode == "FT":
                # feature_key = line
                pass
            elif typecode == 'TG':
                target = self._translate_amino_acid(line.strip("."))
            elif typecode == 'PP':
                position = self._translate_position(line.strip('.'))
            elif typecode == "CF":
                correction_formula = self._formula_parser(line)
            elif typecode == "MM":
                mass_difference = float(line)
            elif typecode == "KW":
                keywords.append(line.strip("."))
            elif typecode == "DR":
                crossref.append(':'.join(line.strip('.').split("; ")))
                if crossref[-1].startswith("PSI-MOD"):
                    crossref[-1] = crossref[-1].replace("PSI-", "")
                    psimod_id = crossref[-1]
            elif typecode == "//":
                break
        if mass_difference is None or correction_formula is None:
            return None
        mod_target = ModificationTarget([target] if target else set(), position)
        rule = ModificationRule(
            [mod_target], ptm_id, monoisotopic_mass=mass_difference, composition=correction_formula,
            alt_names={accession} | set(crossref), categories=keywords)
        rule.preferred_name = ptm_id if psimod_id is None else psimod_id
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


def target_to_source(target: ModificationTarget) -> str:
    aas = {a.symbol for a in target.amino_acid_targets}
    pos = target.position_modifier
    return f"ModificationTarget({aas}, SequenceLocation[{pos.name!r}])"


def rule_to_source(rule: ModificationRule) -> str:
    targets = f"[{', '.join([target_to_source(t) for t in rule.targets])}]"
    name = rule.preferred_name
    mass = rule.mass
    formula = dict(rule.composition)
    return f"ModificationRule({targets}, {name!r}, monoisotopic_mass={mass}, composition=Composition({formula}), alt_names={rule.names})"
