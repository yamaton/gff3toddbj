#!/usr/bin/env python
"""
Normalize the 1st column values ("Entries") in DDBJ annotation.

"""
from typing import Dict, Iterable, Set
import argparse
import re
import csv
import pathlib
import logging

__version__ = "0.1.0"
_PROG_NAME = "normalize-entry-names"

INVALID_LETTERS = r'[=|>" \[\]]'  # =|>" [] in regular expression
NEW_DELIMITER = ":"
INVALID_PATTERN = re.compile(INVALID_LETTERS)
FASTA_HEADER_PATTRN = re.compile(">(.*)")


def get_entry_names(path: str) -> Set[str]:
    """Get Entry values(1st column values) in the annotation file"""
    with open(path, "r") as f:
        spam = csv.reader(f, delimiter="\t")
        names = {row[0].strip() for row in spam if row and (row[0].strip() != "COMMON")}
    return names


def save_renamed_annotation(infile: str, outfile: str, table: Dict[str, str]):
    with open(infile, "r") as f, open(outfile, "w") as outfile:
        spamreader = csv.reader(f, delimiter="\t")
        spamwriter = csv.writer(outfile, delimiter="\t")
        for row in spamreader:
            if row[0].strip() not in table:
                spamwriter.writerow(row)
            else:
                renamed = table[row[0]]
                spamwriter.writerow([renamed, *row[1:]])


def get_rename_dictionary(names: Iterable[str]) -> Dict[str, str]:
    """Create renaming dictionary s.t. {before: after}"""
    result = dict()
    for name in names:
        if INVALID_PATTERN.search(name):
            xs = [x for x in INVALID_PATTERN.split(name) if x]
            new_name = NEW_DELIMITER.join(xs)
            result[name] = new_name
    return result


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(prog=_PROG_NAME)
    parser.add_argument("file", help="Input annotation file")
    parser.add_argument(
        "--suffix", help="Suffix to output filenames", default="_renamed"
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s {}".format(__version__),
        help="Show version",
    )
    args = parser.parse_args()

    path_input = args.file
    out_filename_suffix = args.suffix
    logging.info("Input annotation: {}".format(path_input))

    seqids = get_entry_names(path_input)
    table = get_rename_dictionary(seqids)

    if table:
        logging.info("Renaming ...")
        for key, value in table.items():
            logging.info("  {} \t ---> \t {}".format(key, value))

        p = pathlib.Path(path_input)
        path_output = p.parent / (p.stem + out_filename_suffix + p.suffix)
        logging.info("Output GFF   : {}".format(str(path_output)))
        save_renamed_annotation(path_input, path_output, table)
    else:
        logging.info("=====================================================")
        logging.info("      IDs are fine: No need to regularize them.      ")
        logging.info("=====================================================")


if __name__ == "__main__":
    main()
