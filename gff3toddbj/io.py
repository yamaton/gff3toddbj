import collections
import csv
import gzip
import logging
import pathlib
import pprint
import re
import shutil
import subprocess
import sys
import tempfile
import toml
from typing import Optional, OrderedDict, List, Any, Union
import urllib.parse
import uuid

import Bio.SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
import pysam

Faidx = pysam.libcfaidx.FastaFile

PATH_DATABASE_LOCAL_PREFIX = "temp_gff3-to-ddbj_"
PATH_DATABASE_LOCAL = PATH_DATABASE_LOCAL_PREFIX + str(uuid.uuid4()) + ".db"


class CommandNotFoundError(Exception):
    pass


def load_gff3_as_seqrecords(filepath, unquoting=False) -> List[SeqRecord]:
    """Load GFF3 as iterable of SeqRecord

    Args:
        filepath: path to load GFF3
        unquoting (bool): unquote characters.

    Unquoting means converesions like this:
        from %3B to ;
        from %2C to ,
    """
    ext = pathlib.Path(filepath).suffix

    logging.info("Loading GFF3. May take some time ...")
    if unquoting:
        # Save unquoted content to tempfile before parsing GFF3
        with tempfile.TemporaryFile(mode="w+") as ftemp:
            if ext == ".gz":
                with gzip.open(filepath, "rt") as fin:
                    ftemp.write(urllib.parse.unquote(fin.read()))
            else:
                with open(filepath, "r") as fin:
                    ftemp.write(urllib.parse.unquote(fin.read()))
            ftemp.seek(0)
            recs = _wrapper_gff_parse(ftemp)
    else:
        # If unquoting is unnecessary
        if ext == ".gz":
            with gzip.open(filepath, "rt") as fin:
                recs = _wrapper_gff_parse(fin)
        else:
            with open(filepath, "r") as fin:
                recs = _wrapper_gff_parse(fin)
    logging.info("... Done Loading GFF3")

    return recs


def _wrapper_gff_parse(fileobj) -> List[SeqRecord]:
    """
    Returns a list of SeqRecord from a GFF3 file object.
    Raises TypeError if the opening fails due to ##FASTA directive in GFF3.
    """
    try:
        records = list(GFF.parse(fileobj))
    except TypeError:
        if isinstance(fileobj, str) or isinstance(fileobj, pathlib.Path):
            fileobj = open(fileobj, "r")
        s = fileobj.read()
        patt = re.compile("##FASTA")
        matchobj = patt.search(s)
        if matchobj:
            msg = (
                "\n\n"
                "Directive ##FASTA exists in the input GFF3: Run \n\n"
                "    $ split-fasta <<this-gff3-file.gff3>>\n\n"
                "and split FASTA part from the file before running gff3-to-ddbj."
            )
            logging.error(msg)
            sys.exit(1)
        else:
            raise TypeError("Failed to read GFF3 for unknown reason.")

    # Revert sorting of IDs done by GFF.parse()
    if isinstance(fileobj, str) or isinstance(fileobj, pathlib.Path):
        fileobj = open(fileobj, "r")
    fileobj.seek(0)
    ids = _get_ids_gff3(fileobj)
    ids_score = {id_: (index, id_) for index, id_ in enumerate(ids)}
    records.sort(key=lambda rec: ids_score[rec.id])
    return records


def _get_ids_gff3(fileobj) -> List[str]:
    """
    Get SeqIDs from GFF3 file while keeping the order.

    [NOTE] BcBio.parse() returns a list of SeqRecords sorted by their IDs.
    I want to revert this operation.
    """
    seen = set()
    ids = []
    spamreader = csv.reader(fileobj, delimiter="\t")
    for row in spamreader:
        if not row:
            continue
        x = row[0]
        if (not x.startswith("#")) and (x not in seen):
            ids.append(x)
            seen.add(x)
    return ids


def load_toml_tables(filepath) -> OrderedDict[str, Any]:
    """
    Load TOML as python dictionary
    """
    with open(filepath) as fp:
        d = toml.load(fp, _dict=collections.OrderedDict)

    logging.debug("TOML table:\n{}".format(pprint.pformat(d)))
    return d


def load_fasta_as_faidx(filename: Union[str, pathlib.Path]) -> Faidx:
    """Load FASTA as indexed data with samtools faidx

    """
    logging.info("Loading FASTA...")
    try:
        rec = pysam.FastaFile(str(filename))
    except OSError:
        if pathlib.Path(filename).suffix == ".gz":
            create_bgzipped(filename)
            rec = pysam.FastaFile(str(filename))
        else:
            raise IOError("Failed to load FASTA")

    logging.info("             ... done loading FASTA")
    return rec


def create_bgzipped(path: Union[str, pathlib.Path]):
    """Create a file compressed in bgzip from

    """
    path = pathlib.Path(path)
    assert path.suffix == ".gz"
    if not path.exists():
        raise FileNotFoundError("No such file as {}".format(path))

    backup = path.parent / (path.name + ".bak")
    shutil.move(path, backup)

    logging.info("Re-compressing with bgzip (only once): {}".format(path))
    cmd = "gzip -c -d {} | bgzip --threads=4 > {}".format(str(backup), str(path))
    logging.info("   $ {}".format(cmd))
    ret = subprocess.run(cmd, shell=True)
    if ret.returncode == 127:
        msg = "bgzip is missing: bgzip is a part of samtools"
        logging.error(msg)
        raise CommandNotFoundError(msg)
    if ret.returncode != 0:
        msg = "Failed to re-compress with bgzip"
        logging.error(msg)
        raise ValueError(msg)
    logging.info("    ... done re-compressing FASTA with bgzip ")


def get_seqlens(faidx: Faidx) -> OrderedDict[str, int]:
    """
    Get sequence length as ordered dictionary with SeqID as the key.
    """
    assert len(faidx.references) == len(faidx.lengths)

    d = collections.OrderedDict()
    for (seqid, seqlen) in zip(faidx.references, faidx.lengths):
        d[seqid] = seqlen
    return d


def get_seq(faidx: Faidx, seqid: str, start:Optional[int]=None, end:Optional[int]=None) -> str:
    """
    Get sequence as str

    Args:
        faidx: samtool faidx object (pysam.libcfaidx.FastaFile)
        seqid: SeqID
        start: start position
        end:   end position

    [Note] `start` and `end` follow Python's indexing,
    i.e., 0-based, left-inclusive, right-exclusive.
    """
    return faidx.fetch(reference=seqid, start=start, end=end).upper()
