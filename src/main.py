import argparse
import formatter
import translators
import logging

logging.basicConfig(level=logging.DEBUG)

TRANS_FEATURES = "src/translate_features.yaml"
TRANS_QUALIFIERS = "src/translate_qualifiers.yaml"
DDBJ_RULES = "src/ddbj_rules.yaml"
COMMON = "samples/common.yaml"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3", help="Input GFF3 file")
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("common", help="Input COMMON file in YAML or TSV", nargs='?', default=COMMON)
    args = parser.parse_args()

    logging.info("Input GFF   : {}".format(args.gff3))
    logging.info("Input FASTA : {}".format(args.fasta))
    logging.info("Input COMMON: {}".format(args.common))

    records = translators.run(
        args.gff3, args.fasta, TRANS_FEATURES, TRANS_QUALIFIERS, joinables=("CDS", "exon", "intron")
    )

    logging.info("Records: {}".format(records))

    fmt = formatter.DDBJFormatter(args.common, DDBJ_RULES)
    gen = fmt.run(records, ignore_rules=False)
    for line in gen:
        print(line)


if __name__ == "__main__":
    main()