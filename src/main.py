import argparse
import formatter
import translators
import logging

logging.basicConfig(level=logging.DEBUG)

TRANS_FEATURES = "src/translate_features.toml"
TRANS_QUALIFIERS = "src/translate_qualifiers.toml"
DDBJ_RULES = "src/ddbj_rules.toml"
COMMON = "samples/common.toml"

IGNORE_FEATURE_QUALIFIER_RULE = False
# JOINABLES = ("CDS", "exon", "intron")
JOINABLES = None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3", help="Input GFF3 file")
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("common", help="Input COMMON file in TOML (or TSV ... be be implemented)", nargs='?', default=COMMON)
    args = parser.parse_args()

    logging.info("Input GFF   : {}".format(args.gff3))
    logging.info("Input FASTA : {}".format(args.fasta))
    logging.info("Input COMMON: {}".format(args.common))

    records = translators.run(
        args.gff3, args.fasta, TRANS_FEATURES, TRANS_QUALIFIERS, joinables=JOINABLES
    )

    logging.info("Records: {}".format(records))

    fmt = formatter.DDBJFormatter(args.common, DDBJ_RULES)
    gen = fmt.run(records, ignore_rules=IGNORE_FEATURE_QUALIFIER_RULE)
    for line in gen:
        print(line)


if __name__ == "__main__":
    main()