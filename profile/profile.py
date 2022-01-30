import time
from typing import List

from glycopeptidepy import PeptideSequence

from glycopeptidepy.io import fasta
from glycopeptidepy.structure import modification
from glycopeptidepy.structure.modification import RestrictedModificationTable
from glycopeptidepy.enzyme import Protease
from glycopeptidepy.algorithm import ProteinDigestor
from glycopeptidepy.test.common import datafile

from glypy.structure.glycan_composition import HashableGlycanComposition


constant_modifications = ["Carbamidomethyl (C)"]
variable_modifications = ["Deamidation (N)", "Pyro-glu (Q@N-term)"]

modification_table = RestrictedModificationTable(
    constant_modifications=constant_modifications,
    variable_modifications=variable_modifications)

constant_modifications: List[modification.ModificationRule] = [
    modification_table[name] for name in constant_modifications]
variable_modifications: List[modification.ModificationRule] = [
    modification_table[name] for name in variable_modifications]


digestor = ProteinDigestor(
    Protease("trypsin"),
    constant_modifications=constant_modifications,
    variable_modifications=variable_modifications)

with open("./glycopeptidepy/test/test_data/glycan_compositions.txt", 'rt') as fh:
    glycan_compositions: List[HashableGlycanComposition] = list(
        map(lambda x: HashableGlycanComposition.parse(x.strip()), fh))

protein_iterator = fasta.ProteinFastaFileReader(
    "./glycopeptidepy/test/test_data/proteins.fa")

peptide: PeptideSequence

n_peptides = 0
n_fragments = 0
start_time = time.monotonic()
end_time = None
for i, protein in enumerate(protein_iterator):
    print(f"Digesting {protein.name} {i}/{n_peptides}/{n_fragments} {time.monotonic() - start_time:0.3f} sec")
    for peptide in digestor.digest(protein):
        glycosylation_sites: List[int] = peptide.n_glycan_sequon_sites
        if glycosylation_sites:
            for glycan_comp in glycan_compositions:
                for site in glycosylation_sites:
                    template = peptide.clone()
                    template.add_modification(site, "N-Glycosylation")
                    template.glycan = glycan_comp
                    n_peptides += 1

                    for frags in template.get_fragments('b'):
                        for frag in frags:
                            frag.name
                            frag.mass
                            n_fragments += 1

                    for frags in template.get_fragments('y'):
                        for frag in frags:
                            frag.name
                            frag.mass
                            n_fragments += 1

                    for frag in template.stub_fragments(True, True):
                        frag.name
                        frag.mass
                        n_fragments += 1
    if i == 10:
        break
end_time = time.monotonic()
print(f"{end_time - start_time:0.3f} seconds elapsed")
