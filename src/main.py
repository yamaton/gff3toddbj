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

    records = translators.run(
        args.gff3, args.fasta, TRANS_FEATURES, TRANS_QUALIFIERS, is_joining=False
    )

    fmt = formatter.DDBJFormatter(DDBJ_RULES)
    for rec in records:
        gen = fmt.run(rec)
        for line in gen:
            print(line)
