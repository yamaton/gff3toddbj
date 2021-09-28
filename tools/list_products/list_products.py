#!/usr/bin/env python
"""
Extract and list up products from GFF3

DDBJ allows a single /product for a feature while annotation software often produce multiple.
This script helps choosing /product by listing up. See XXX for the detail.
"""
import argparse
import logging
import collections
import itertools
from typing import Generator, Optional, Tuple, Iterable
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from BCBio import GFF

__version__ = "0.1.0"
_PROG_NAME = "list-products"

class OrderedCounter(collections.Counter, collections.OrderedDict):
    pass


def gen_products(rec: SeqRecord, focused_types: Optional[Tuple]=None) -> Generator[str, None, None]:
    """Yield products


    """
    def _helper(features: Iterable[SeqFeature]):
        for f in features:
            if (focused_types is None) or (f.type in focused_types):
                p = f.qualifiers.get("product", None)
                if p:
                    yield from p
                yield from _helper(f.sub_features)

    yield from _helper(rec.features)



def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(prog=_PROG_NAME)
    parser.add_argument("gff3", help="Input GFF3 file")
    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s {}".format(__version__),
        help="Show version",
    )
    args = parser.parse_args()
    logging.info("Input GFF   : {}".format(args.gff3))

    records = GFF.parse(args.gff3)
    generators = [gen_products(rec) for rec in records]
    print(generators)
    gen = itertools.chain.from_iterable(generators)
    d = OrderedCounter(gen)

    for prod, cnt in d.items():
        logging.debug("{}\t\t\t\t(count {})".format(prod, cnt))

    for prod in d:
        print(prod)

if __name__ == "__main__":
    main()
