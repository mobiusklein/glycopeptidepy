from collections import defaultdict


from ..structure.sequence import ProteinSequence

from .fasta import PEFFFeature, PEFFFastaHeader
from . import uniprot


class UniprotToPeffConverter(object):

    feature_to_peff_term = {
        uniprot.SignalPeptide.feature_type: ('Processed', 'PEFF:0001021'),
        uniprot.MatureProtein.feature_type: ('Processed', 'PEFF:0001020'),
        uniprot.TransitPeptide.feature_type: ('Processed', 'PEFF:0001022'),
        uniprot.Propeptide.feature_type: ('Processed', 'PEFF:0001034'),
        uniprot.SequenceVariant.feature_type: ('VariantSimple', 'PEFF:0001028'),
    }

    def handle_Processed(self, feature, accession):
        peff_feature = PEFFFeature(feature.start + 1, feature.end, feature.feature_type)
        return peff_feature

    def handle_VariantSimple(self, feature, accession):
        peff_feature = PEFFFeature(feature.position, feature.original, feature.variation)
        return peff_feature

    def __call__(self, record):
        features = defaultdict(list)
        features['Prefix'] = 'sp'
        features['Tag'] = record.accessions[0]
        features['PName'] = record.recommended_name
        features['GName'] = record.gene_name
        for feature in record.features:
            try:
                key = self.feature_to_peff_term[feature.feature_type]
                handler = getattr(self, 'handle_' + key[0], None)
                if handler is None:
                    continue
                features[key[0]].append(handler(feature, key[1]))
            except KeyError:
                continue
        header = PEFFFastaHeader(dict(features))
        return ProteinSequence(header, record.sequence)
