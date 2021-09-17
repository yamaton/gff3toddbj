#!/usr/bin/env python
import argparse
import logging
import pathlib

from . import formatter
from . import transforms
from . import utils

FORMAT = '%(levelname)s: %(message)s'

logging.basicConfig(level=logging.INFO, format=FORMAT)

_DIR = pathlib.Path(__file__).parent
PATH_TRANS_FEATURES = _DIR / "translate_features.toml"
PATH_TRANS_QUALIFIERS = _DIR / "translate_qualifiers.toml"
PATH_DDBJ_RULES = _DIR / "ddbj_rules.toml"
PATH_METADATA_DEFAULT = _DIR / "metadata.toml"

LOCUS_TAG_PREFIX = "LOCUSTAGPREFIX_"

IGNORE_FEATURE_QUALIFIER_RULE = False
JOINABLES = ("mRNA", "CDS")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff3", "--gff", help="Input GFF3 file")
    parser.add_argument("--fasta", help="Input FASTA file", required=True)
    parser.add_argument(
        "--metadata",
        help="Input metadata in TOML describing COMMON and other entires",
        default=PATH_METADATA_DEFAULT,
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
        "--translate_features",
        help="Translation table for features",
        default=PATH_TRANS_FEATURES,
    )
    parser.add_argument(
        "--translate_qualifiers",
        help="Translation table for qualifiers",
        default=PATH_TRANS_QUALIFIERS,
    )
    parser.add_argument(
        "-o", "--output",
        help="Specify annotation file name as output",
    )
    args = parser.parse_args()
    output = args.output

    logging.info("Input GFF   : {}".format(args.gff3))
    logging.info("Input FASTA : {}".format(args.fasta))
    logging.info("Input metadata: {}".format(args.metadata))
    logging.info("Prefix of locus_tag: {}".format(args.prefix))
    logging.info("transl_table (The Genome Code): {}".format(args.transl_table))
    if args.translate_features != PATH_TRANS_FEATURES:
        logging.info(f"feature translation: {args.translate_features}")
    if args.translate_qualifiers != PATH_TRANS_QUALIFIERS:
        logging.info(f"qualifier translation: {args.translate_qualifiers}")
    if output:
        logging.info("Output  : {}".format(output))

    metadata = utils.load_metadata_info(args.metadata)

    records = transforms.run(
        args.gff3,
        args.fasta,
        args.translate_features,
        args.translate_qualifiers,
        metadata,
        args.prefix,
        args.transl_table,
        joinables=JOINABLES,
    )

    logging.debug("Records: {}".format(records))

    fmt = formatter.DDBJFormatter(metadata, PATH_DDBJ_RULES)
    gen = fmt.run(records, ignore_rules=IGNORE_FEATURE_QUALIFIER_RULE)

    if output:
        parent_dir = pathlib.Path(output).parent
        parent_dir.mkdir(exist_ok=True)
        with open(output, "w") as f:
            for line in gen:
                print(line, file=f)
    else:
        for line in gen:
            print(line)


if __name__ == "__main__":
    main()
