import argparse
import formatter
import translators
import logging

import utils

logging.basicConfig(level=logging.DEBUG)

PATH_TRANS_FEATURES = "src/translate_features.toml"
PATH_TRANS_QUALIFIERS = "src/translate_qualifiers.toml"
PATH_DDBJ_RULES = "src/ddbj_rules.toml"
PATH_METADATA = "samples/common.toml"

LOCUS_TAG_PREFIX = "LOCUSTAGPREFIX_"

IGNORE_FEATURE_QUALIFIER_RULE = False
# JOINABLES = ("CDS", "exon", "intron")
JOINABLES = None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff3", help="Input GFF3 file")
    parser.add_argument("--fasta", help="Input FASTA file", required=True)
    parser.add_argument(
        "--metadata",
        help="Input COMMON file in TOML (or TSV ... be be implemented)",
        default=PATH_METADATA,
    )
    parser.add_argument(
        "-p",
        "--locus_tag_prefix",
        help="Prefix of locus_tag. See https://www.ddbj.nig.ac.jp/ddbj/locus_tag-e.html",
        default=LOCUS_TAG_PREFIX,
    )
    args = parser.parse_args()

    logging.info("Input GFF   : {}".format(args.gff3))
    logging.info("Input FASTA : {}".format(args.fasta))
    logging.info("Input meta info: {}".format(args.metadata))
    logging.info("Prefix of locus_tag: {}".format(args.locus_tag_prefix))

    metadata = utils.load_header_info(args.metadata)

    records = translators.run(
        args.gff3,
        args.fasta,
        PATH_TRANS_FEATURES,
        PATH_TRANS_QUALIFIERS,
        metadata,
        args.locus_tag_prefix,
        joinables=JOINABLES,
    )

    logging.info("Records: {}".format(records))

    fmt = formatter.DDBJFormatter(metadata, PATH_DDBJ_RULES)
    gen = fmt.run(records, ignore_rules=IGNORE_FEATURE_QUALIFIER_RULE)
    for line in gen:
        print(line)


if __name__ == "__main__":
    main()
