#!/usr/bin/env python
"""
Split FASTA from GFF3 and save both FASTA and slimmed-out GFF3.

FASTA within GFF3 is detected via ##FASTA directive, i.e.
non-tab-separated lines following the directive are considered as FASTA.
"""
from typing import Tuple
import argparse
import pathlib
import logging

__version__ = "0.1.0"
_PROG_NAME = "split-fasta"

def _split_to_parent_stem_ext(p: pathlib.Path) -> Tuple[pathlib.Path, str, str]:
    """Split input file path into (parent, filename_without_ext, ext)"""
    if not isinstance(p, pathlib.Path):
        p = pathlib.Path(p)
    parent = p.parent
    ext = "." + "".join(suf[1:] for suf in p.suffixes)
    filename_wo_ext = p.name[: -len(ext)]
    return (parent, filename_wo_ext, ext)


def get_output_path(input_path: pathlib.Path, suffix: str) -> pathlib.Path:
    """
    Get output filename path by adding suffix.
    Also checks existing file, and add a digit to avoid a conflict.
    Raises an error if all digits 0..9 are exhausted.
    """
    if not get_fasta_output_path:
        input_path = pathlib.Path(input_path)

    parent, filename_wo_ext, ext = _split_to_parent_stem_ext(input_path)
    filename_stem = filename_wo_ext + suffix

    for i in range(10):
        if i == 0:
            out_path = parent / (filename_stem + ext)
        else:
            out_path = parent / (filename_stem + "__{}".format(i) + ext)
        if not out_path.exists():
            return out_path
    else:
        raise FileExistsError("Output file already exists: Clean up the folder.")


def get_fasta_output_path(input_gff3: pathlib.Path, suffix: str):
    parent, filename_wo_ext, ext = _split_to_parent_stem_ext(input_gff3)
    fasta_ext = ext.split(".")
    fasta_ext[1] = "fa"
    fasta_ext = ".".join(fasta_ext)
    p = parent / (filename_wo_ext + fasta_ext)
    return get_output_path(p, suffix)


def split(input_gff3, output_gff3, output_fasta) -> None:
    """Switch output stream line by line
    """
    is_reading_gff3 = True
    with open(input_gff3, "r") as fin, open(output_gff3, "w") as fout_gff3, open(output_fasta, "w") as fout_fasta:
        for line in fin:
            # check transition between GFF3 and FASTA
            if is_reading_gff3:
                if line.startswith("##FASTA"):
                    is_reading_gff3 = False
                    continue
            else:
                if (not line.startswith(">")) and ("\t" in line):
                    is_reading_gff3 = True

            # write to appropriate file
            if is_reading_gff3:
                print(line, end="", file=fout_gff3)
            else:
                print(line, end="", file=fout_fasta)


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        description="Split FASTA from GFF3",
        prog=_PROG_NAME,
    )
    parser.add_argument("gff3", help="Input GFF3")
    parser.add_argument(
        "--suffix", help="Suffix added to the output filenames", default="_splitted"
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s {}".format(__version__),
        help="Show version",
    )
    args = parser.parse_args()

    input_gff3 = args.gff3
    suffix = args.suffix
    logging.info("Input GFF3  : {}".format(input_gff3))

    output_gff3 = get_output_path(input_gff3, suffix)
    logging.info("Output GFF3 : {}".format(output_gff3))

    output_fasta = get_fasta_output_path(input_gff3, suffix)
    logging.info("Output FASTA: {}".format(output_fasta))

    split(input_gff3=input_gff3, output_gff3=output_gff3, output_fasta=output_fasta)

