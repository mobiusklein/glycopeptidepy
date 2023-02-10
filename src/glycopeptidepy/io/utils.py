from typing import Dict, DefaultDict, List, Union


from glycopeptidepy.structure import ProteinSequence, ModificationTable, ModificationRule

from .fasta import PEFFFeature, PEFFFastaHeader
from . import uniprot
from .cv.uniprot_ptm import load as load_uniprot_ptms


class UniprotToPeffConverter(object):

    feature_to_peff_term = {
        uniprot.SignalPeptide.feature_type: ('Processed', 'PEFF:0001021'),
        uniprot.MatureProtein.feature_type: ('Processed', 'PEFF:0001020'),
        uniprot.TransitPeptide.feature_type: ('Processed', 'PEFF:0001022'),
        uniprot.Propeptide.feature_type: ('Processed', 'PEFF:0001034'),
        uniprot.SequenceVariant.feature_type: ('VariantSimple', 'PEFF:0001028'),
        uniprot.GlycosylationSite.feature_type: ('ModRes', ''),
        uniprot.ModifiedResidue.feature_type: ('ModResPsi', '')
    }

    uniprot_ptms: ModificationTable
    psimod_rules: Dict[str, ModificationRule]

    def __init__(self):
        self.uniprot_ptms = load_uniprot_ptms()
        self.psimod_rules = {}
        self.prepare_psimod_conversion_table()

    def prepare_psimod_conversion_table(self):
        psimod_mapping = {}
        for mod in self.uniprot_ptms.rules():

            mod_names = [name for name in mod.names if name.startswith('MOD:')]
            if mod_names:
                psimod_mapping[mod.name] = mod_names[0]
        self.psimod_rules = psimod_mapping

    def handle_Processed(self, feature: uniprot.PeptideBase, accession) -> PEFFFeature:
        peff_feature = PEFFFeature(feature.start + 1, feature.end, accession, feature.feature_type)
        return peff_feature

    def handle_VariantSimple(self, feature, accession) -> PEFFFeature:
        peff_feature = PEFFFeature(feature.position + 1, feature.original, feature.variation)
        return peff_feature

    def handle_ModRes(self, feature: uniprot.ModifiedResidue, accession) -> PEFFFeature:
        if accession:
            peff_feature = PEFFFeature(
                feature.position + 1, str(feature.description), accession)
        else:
            peff_feature = PEFFFeature(feature.position + 1, str(feature.description))
        return peff_feature

    def handle_ModResPsi(self, feature: uniprot.ModifiedResidue, accession) -> PEFFFeature:
        if feature.description in self.psimod_rules:
            accession = self.psimod_rules[feature.description]
            return PEFFFeature(feature.position + 1, accession, feature.description)

    def base_features(self, record: uniprot.UniProtProtein) -> DefaultDict[str, Union[str, List[PEFFFeature]]]:
        features = DefaultDict(list)
        features['Prefix'] = 'sp'
        features['Tag'] = record.accessions[0]
        features['PName'] = record.recommended_name
        features['GName'] = record.gene_name
        return features

    def __call__(self, record: uniprot.UniProtProtein):
        features = self.base_features(record)
        for feature in record.features:
            try:
                key = self.feature_to_peff_term[feature.feature_type]
                handler = getattr(self, 'handle_' + key[0], None)
                if handler is None:
                    continue
                f = handler(feature, key[1])
                if f is None:
                    continue
                features[key[0]].append(f)
            except KeyError:
                continue
        header = PEFFFastaHeader(dict(features))
        return ProteinSequence(header, record.sequence)
