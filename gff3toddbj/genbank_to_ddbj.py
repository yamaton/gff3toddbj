#!/usr/bin/env python
import argparse
import collections
import logging
import pathlib

from . import formatter
from . import io
from . import utils

_EXEC_NAME = "genbank-to-ddbj"
__version__ = "0.0.3"

FORMAT = '%(levelname)s: %(message)s'
_DIR = pathlib.Path(__file__).parent
PATH_TRANS_FEATURES_QUALIFIERS = _DIR / "translate_features_qualifiers.toml"

PATH_DDBJ_FILTER = _DIR / "ddbj_filter.toml"
LOCUS_TAG_PREFIX = "LOCUSTAGPREFIX_"

IGNORE_FILTERING_RULES = False


def main():
    parser = argparse.ArgumentParser(prog=_EXEC_NAME)
    parser.add_argument("--gbk", "--gbff", "--genbank", metavar="FILE", help="Input GenBank file")
    parser.add_argument(
        "--metadata",
        help="Input metadata in TOML describing COMMON and other entires",
        metavar="FILE",
        default=None,
    )
    parser.add_argument(
        "-p",
        "--prefix",
        "--locus_tag_prefix",
        help="Prefix of locus_tag. See https://www.ddbj.nig.ac.jp/ddbj/locus_tag-e.html",
        metavar="STR",
        default=LOCUS_TAG_PREFIX,
    )
    parser.add_argument(
        "--transl_table",
        help="Genetic Code ID. 1 by default, and 11 for bacteria. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi",
        type=int,
        metavar="INT",
        default=1,
    )
    parser.add_argument(
        "--config_filter",
        help="A set of Feature-Qualifier pairs allowed in the output. See https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-e.pdf",
        metavar="FILE",
        default=PATH_DDBJ_FILTER,
    )
    parser.add_argument(
        "-o", "--out", "--output",
        metavar="FILE",
        help="Specify annotation file name as output",
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s {}".format(__version__),
        help="Show version",
    )
    parser.add_argument(
        "--log",
        default="INFO",
        metavar="STR",
        help="[debug] Choose log level from (DEBUG, INFO, WARNING, ERROR) (default: INFO).",
    )
    args = parser.parse_args()

    logging.basicConfig(level=utils.to_loglevel(args.log), format=FORMAT)

    logging.info("Input GenBank       : {}".format(args.gbk))
    logging.info("Input metadata      : {}".format(args.metadata))
    logging.info("Prefix of locus_tag : {}".format(args.prefix))
    logging.info("transl_table        : {}".format(args.transl_table))
    if args.config_filter != PATH_DDBJ_FILTER:
        logging.info(f"Config-Filter       : {args.config_filter}")
    output = args.out
    if output:
        logging.info("Output              : {}".format(output))

    metadata = collections.OrderedDict() if args.metadata is None else utils.load_metadata_info(args.metadata)
    config_filter = utils.load_rules(args.config_filter)

    records = io.load_flatfile_as_seqrecords(args.gbk)
    fmt = formatter.DDBJFormatter(metadata, config_filter)
    gen = fmt.run(records, ignore_rules=IGNORE_FILTERING_RULES)

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
