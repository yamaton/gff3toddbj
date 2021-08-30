import argparse
from . import formatter
from . import translators

TRANS_FEATURES = "./translate_features.yaml"
TRANS_QUALIFIERS = "./translate_qualifiers.yaml"
DDBJ_RULES = "./ddbj_feature-qualifier_lists.yaml"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3", help="Input GFF3 file")
    parser.add_argument("fasta", help="Input FASTA file")
    args = parser.parse_args()




    t_features = translators.TranslateFeatures(TRANS_FEATURES)
    t_qualifiers = translators.TranslateQualifiers(TRANS_QUALIFIERS)

