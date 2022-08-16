from glycopeptidepy.structure import (PeptideSequence, AminoAcidResidue,
                                      Modification, ModificationNameResolutionError, AnonymousModificationRule)


class MzIdentMLPeptideParser(object):
    def __init__(self, modification_translation_table=None):
        if modification_translation_table is None:
            modification_translation_table = {}
        self.modification_translation_table = modification_translation_table

    def process(self, peptide_dict):
        peptide_sequence = PeptideSequence(peptide_dict["PeptideSequence"])
        peptide_sequence.id = peptide_dict.get('id')
        self.handle_substitutions(peptide_dict, peptide_sequence)
        self.handle_modifications(peptide_dict, peptide_sequence)
        return peptide_sequence

    def __call__(self, peptide_dict):
        return self.process(peptide_dict)

    def handle_substitutions(self, peptide_dict, peptide_sequence):
        if "SubstitutionModification" in peptide_dict:
            subs = peptide_dict["SubstitutionModification"]
            for sub in subs:
                pos = sub['location'] - 1
                replace = AminoAcidResidue(sub["replacementResidue"])
                peptide_sequence.substitute(pos, replace)

    def add_modification(self, peptide_sequence, modification, position):
        if position == -1:
            targets = modification.rule.n_term_targets
            for t in targets:
                if t.position_modifier is not None and t.amino_acid_targets is None:
                    break
            else:
                position += 1
        if position == len(peptide_sequence):
            targets = modification.rule.c_term_targets
            for t in targets:
                if t.position_modifier is not None and t.amino_acid_targets is None:
                    break
            else:
                position -= 1

        if position == -1:
            peptide_sequence.n_term = modification
        elif position == len(peptide_sequence):
            peptide_sequence.c_term = modification
        else:
            peptide_sequence.add_modification(position, modification)

    def handle_modifications(self, peptide_dict, peptide_sequence):
        if "Modification" in peptide_dict:
            mods = peptide_dict["Modification"]
            for mod in mods:
                pos = mod["location"] - 1
                accession = None
                try:
                    if "unknown modification" in mod:
                        try:
                            _name = mod['unknown modification']
                            if _name in self.modification_translation_table:
                                modification = self.modification_translation_table[_name](
                                )
                            else:
                                modification = Modification(str(_name))
                        except ModificationNameResolutionError:
                            raise KeyError("Cannot find key in %r" % (mod,))
                    else:
                        try:
                            _name = mod["name"]
                            accession = getattr(_name, "accession", None)
                            _name = str(_name)
                            if accession is not None:
                                accession = str(accession)
                                try:
                                    modification = Modification(accession)
                                except ModificationNameResolutionError:
                                    modification = Modification(_name)
                            else:
                                modification = Modification(_name)

                        except (KeyError, ModificationNameResolutionError) as e:
                            raise KeyError(
                                "Cannot find key %s in %r" % (e, mod))

                    self.add_modification(peptide_sequence, modification, pos)
                except KeyError:
                    if "unknown modification" in mod:
                        if 'monoisotopicMassDelta' in mod:
                            mass = float(mod['monoisotopicMassDelta'])
                            modification = AnonymousModificationRule(
                                str(_name), mass)()
                            self.add_modification(
                                peptide_sequence, modification, pos)
                        else:
                            raise
                    else:
                        raise


parse = MzIdentMLPeptideParser()
