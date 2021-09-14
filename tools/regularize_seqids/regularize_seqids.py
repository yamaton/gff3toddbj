#!/usr/bin/env python
"""
Regularize SeqIDs and FASTA headers

DDBJ annotation disallow SeqIDs containing characters  =|>" []

1st column in GFF3 file (= SeqID) and headers in FASTA file (= headers) must match.

"""
from typing import Dict, Iterable, Set
import argparse
import re
import csv
import pathlib
import logging

INVALID_LETTERS = r'[=|>" \[\]]'  # =|>" [] in regular expression
NEW_DELIMITER = ":"
INVALID_PATTERN = re.compile(INVALID_LETTERS)
FASTA_HEADER_PATTRN = re.compile(">(.*)")


def get_gff3_seqids(path: str) -> Set[str]:
    """Get SeqIDs (1st column) from a GFF3 file"""
    with open(path, "r") as f:
        spam = csv.reader(
            f,
            delimiter="\t",
        )
        seqids = {row[0].strip() for row in spam if not row[0].startswith("#")}
    return seqids


def save_renamed_gff3(infile: str, outfile: str, table: Dict[str, str]):
    with open(infile, "r") as f, open(outfile, "w") as outfile:
        spamreader = csv.reader(f, delimiter="\t")
        spamwriter = csv.writer(outfile, delimiter="\t")
        for row in spamreader:
            if row[0].strip() not in table:
                spamwriter.writerow(row)
            else:
                renamed = table[row[0]]
                spamwriter.writerow([renamed, *row[1:]])


def get_fasta_headers(path: str) -> Set[str]:
    """Get headers from a FASTA file"""
    result = set()
    with open(path, "r") as f:
        s = f.readlines()
        for line in s:
            m = FASTA_HEADER_PATTRN.match(line)
            if m:
                header = m.group(1).strip()
                name = header.split()[0]
                result.add(name)
    return result


def save_renamed_fasta(infile: str, outfile: str, table: Dict[str, str]):
    with open(infile, "r") as f, open(outfile, "w") as outfile:
        s = f.readlines()
        for line in s:
            if line.startswith(">"):
                header = line[1:].strip()
                name = header.split()[0]
                rest = header[len(name):]
                if name in table:
                    renamed = table[name]
                    s = renamed + rest
                    print(">{}".format(s), file=outfile)
                else:
                    print(line, file=outfile, end="")
            else:
                print(line, file=outfile, end="")


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

    parser = argparse.ArgumentParser()
    parser.add_argument("--gff3", "--gff", help="Input GFF3 file", required=True)
    parser.add_argument("--fasta", "--fa", help="Input FASTA file", required=True)
    parser.add_argument("--suffix", help="Suffix to output filenames", default="_renamed_ids")
    args = parser.parse_args()

    path_input_gff3 = args.gff3
    path_input_fasta = args.fasta
    out_filename_suffix = args.suffix
    logging.info("Input GFF   : {}".format(path_input_gff3))
    logging.info("Input FASTA : {}".format(path_input_fasta))

    seqids = get_gff3_seqids(args.gff3)
    headers = get_fasta_headers(args.fasta)
    seqids_only = seqids - headers

    if seqids_only:
        msg = "------> Following SeqID is not found in the FASTA file: {} <------"
        for name in seqids_only:
            logging.warning(msg.format(name))
    else:
        logging.info("GFF3 and FASTA have matching IDs!")

    table = get_rename_dictionary(seqids)
    if table:
        logging.info("Renaming ...")
        for key, value in table.items():
            logging.info("  {} \t ---> \t {}".format(key, value))

        p_gff3 = pathlib.Path(args.gff3)
        path_output_gff3 = p_gff3.parent / (p_gff3.stem + out_filename_suffix + p_gff3.suffix)
        logging.info("Output GFF   : {}".format(str(path_output_gff3)))
        save_renamed_gff3(path_input_gff3, path_output_gff3, table)

        p_fasta = pathlib.Path(args.fasta)
        path_output_fasta = p_fasta.parent / (p_fasta.stem + out_filename_suffix + p_fasta.suffix)
        logging.info("Output FASTA : {}".format(str(path_output_fasta)))
        save_renamed_fasta(path_input_fasta, path_output_fasta, table)

    else:
        logging.info("=====================================================")
        logging.info("      IDs are fine: No need to regularize them.      ")
        logging.info("=====================================================")


if __name__ == "__main__":
    main()
