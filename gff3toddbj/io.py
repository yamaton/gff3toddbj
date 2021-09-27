import collections
import csv
import gzip
import logging
import pathlib
import pprint
import re
import sqlite3
import sys
import tempfile
import toml
from typing import OrderedDict, List, Any
import urllib.parse
import uuid

import Bio.SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF


PATH_DATABASE_LOCAL_PREFIX = "temp_gff3-to-ddbj_"
PATH_DATABASE_LOCAL = PATH_DATABASE_LOCAL_PREFIX + str(uuid.uuid4()) + ".db"


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


def load_fasta_as_database(filepath) -> sqlite3.Connection:
    """
    Load FASTA as sqlite3 database with (id, sequence) as keys,
    mainly for saving memory use.
    """
    # remove the database file if exists
    p = pathlib.Path(PATH_DATABASE_LOCAL)
    p.unlink(missing_ok=True)

    # create new database
    con = sqlite3.connect(p)
    cur = con.cursor()
    cur.execute(
        """CREATE TABLE fasta
            (SeqId VARCHAR(50) PRIMARY KEY, SeqLen INT, Sequence TEXT)"""
    )

    path_fasta = pathlib.Path(filepath)
    if path_fasta.suffix == ".gz":
        f = gzip.open(filepath, "rt")
    else:
        f = open(filepath, "r")

    logging.info("Loading FASTA data to a local database.")
    logging.info("    It may take minutes. Coffee break? ☕ ...")
    for rec in Bio.SeqIO.parse(f, "fasta"):
        seq = str(rec.upper().seq)
        seqlen = len(rec)
        cur.execute(
            """
            INSERT INTO fasta VALUES
            (?, ?, ?)""",
            (rec.id, seqlen, seq),
        )
    logging.info("... Done loading FASTA.")

    cur.close()
    f.close()
    return con


def close_and_remove_database(con: sqlite3.Connection) -> None:
    """
    Close the database connection, and REMOVE the sqlite3 file.
    """
    con.close()
    p = pathlib.Path(PATH_DATABASE_LOCAL)
    p.unlink()


def get_seqlens(cur: sqlite3.Cursor) -> OrderedDict[str, int]:
    """
    Get sequence length as ordered dictionary with SeqID as the key.
    """
    res = cur.execute(
        """
        SELECT SeqID, SeqLen FROM fasta"""
    )
    d = collections.OrderedDict()
    for (seqid, seqlen) in res.fetchall():
        d[seqid] = seqlen
    return d


def get_seq(cur: sqlite3.Cursor, seqid: str) -> str:
    """
    Get sequence as str

    Args:
        cur: database cursor
        seqid: SeqID
    """
    res = cur.execute(
        """
        SELECT Sequence FROM fasta
        WHERE SeqID = ?""",
        (seqid,),
    )
    tup = res.fetchone()
    if tup is None:
        logging.warning("Sequence NOT FOUND for SeqID: {}".format(seqid))
        return ""

    return tup[0]
