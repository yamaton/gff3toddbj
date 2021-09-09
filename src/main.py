import argparse
import formatter
import translators
import logging

import utils

FORMAT = '%(levelname)s: %(message)s'

logging.basicConfig(level=logging.INFO, format=FORMAT)

PATH_TRANS_FEATURES = "src/translate_features.toml"
PATH_TRANS_QUALIFIERS = "src/translate_qualifiers.toml"
PATH_DDBJ_RULES = "src/ddbj_rules.toml"
PATH_METADATA = "metadata.toml"

LOCUS_TAG_PREFIX = "LOCUSTAGPREFIX_"

IGNORE_FEATURE_QUALIFIER_RULE = False
# JOINABLES = ("CDS", "exon", "intron")
JOINABLES = None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff3", "--gff", help="Input GFF3 file")
    parser.add_argument("--fasta", help="Input FASTA file", required=True)
    parser.add_argument(
        "--metadata",
        help="Input metadata in TOML describing COMMON and other entires",
        default=PATH_METADATA,
    )
    parser.add_argument(
        "-p",
        "--prefix",
        "--locus_tag_prefix",
        help="Prefix of locus_tag. See https://www.ddbj.nig.ac.jp/ddbj/locus_tag-e.html",
        default=LOCUS_TAG_PREFIX,
    )
    parser.add_argument(
        "--transl_table",
        help="Genetic Code ID. 1 by default, and 11 for bacteria. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi",
        type=int,
        default=1,
    )
    parser.add_argument(
        "-o", "--output",
        help="Specify annotation file name as output",
    )
    args = parser.parse_args()
    output = args.output

    logging.info("Input GFF   : {}".format(args.gff3))
    logging.info("Input FASTA : {}".format(args.fasta))
    logging.info("Input meta info: {}".format(args.metadata))
    logging.info("Prefix of locus_tag: {}".format(args.prefix))
    logging.info("transl_table (The Genome Code): {}".format(args.transl_table))
    if output:
        logging.info("Output  : {}".format(output))

    metadata = utils.load_header_info(args.metadata)

    records = translators.run(
        args.gff3,
        args.fasta,
        PATH_TRANS_FEATURES,
        PATH_TRANS_QUALIFIERS,
        metadata,
        args.prefix,
        args.transl_table,
        joinables=JOINABLES,
    )

    logging.debug("Records: {}".format(records))

    fmt = formatter.DDBJFormatter(metadata, PATH_DDBJ_RULES)
    gen = fmt.run(records, ignore_rules=IGNORE_FEATURE_QUALIFIER_RULE)

    if output:
        with open(output, "w") as f:
            for line in gen:
                print(line, file=f)
    else:
        for line in gen:
            print(line)


if __name__ == "__main__":
    main()
