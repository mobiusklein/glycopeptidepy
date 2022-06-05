from glypy import Composition

from glycopeptidepy.structure.modification import ModificationRule, ModificationTarget, SequenceLocation


def parse_formula(formula):
    tokens = formula.split(" ")
    token_stream = iter(tokens)
    n_items = len(tokens) // 2
    comp = Composition()
    for a, b in ((next(token_stream), next(token_stream)) for i in range(n_items)):
        b = int(b)
        comp[a] += b
    return comp


def psimod_to_rule(term):
    try:
        origin = term.Origin
    except KeyError:
        origin = None
    try:
        mass = term.DiffMono
    except KeyError:
        return None
    if origin is None or ',' in origin or ":" in origin or "X" in origin:
        return None
    if term.DiffFormula is None or mass is None:
        return None
    if mass == 0:
        return None
    comp = parse_formula(term.DiffFormula)
    try:
        loc = SequenceLocation[term.TermSpec] # pylint: disable=unsubscriptable-object
    except KeyError:
        loc = SequenceLocation.anywhere
    target = ModificationTarget([origin] if origin else [], loc)
    name = term.id
    return ModificationRule(target, name, name, mass, comp)


def target_to_source(target: ModificationTarget) -> str:
    aas = {a.symbol for a in target.amino_acid_targets}
    pos = target.position_modifier
    return f"ModificationTarget({aas}, SequenceLocation[{pos.name!r}])"


def rule_to_source(rule: ModificationRule) -> str:
    targets = f"[{', '.join([target_to_source(t) for t in rule.targets])}]"
    name = rule.name
    mass = rule.mass
    formula = dict(rule.composition)
    return f"ModificationRule({targets}, {name!r}, monoisotopic_mass={mass}, composition=Composition({formula}))"
