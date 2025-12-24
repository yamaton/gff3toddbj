"""
Sketch of DDBJ annotation parser providing `load_ddbj()`

"""
from typing import Generator, Union
import csv
import logging

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition, SeqFeature, SimpleLocation, CompoundLocation, BeforePosition, AfterPosition


FALLBACK_SEQLEN = 0

def _to_featureloc(s: str, is_positive_strand: bool) -> SimpleLocation:
    """Parse simple location of the form.

    Note that SimpleLocation takes 0-based both-inclusive indexing
    while texts from GFF3/DDBJ/flatfile have 1-based left-inclusive
    right-exclusive indexing.

    >>> _to_featureloc("123", False)
    SimpleLocation(ExactPosition(122), ExactPosition(123), strand=-1)

    >>> _to_featureloc("<23..35", True)
    SimpleLocation(BeforePosition(22), ExactPosition(35), strand=1)

    >>> _to_featureloc("123^124", True)
    SimpleLocation(ExactPosition(122), ExactPosition(123), strand=1)

    >>> _to_featureloc("<123", True)
    SimpleLocation(BeforePosition(122), ExactPosition(123), strand=1)

    >>> _to_featureloc(">123", True)
    SimpleLocation(ExactPosition(122), AfterPosition(123), strand=1)
    """
    if ".." in s:
        start, end = s.split("..")
        start = BeforePosition(int(start[1:]) - 1) if start.startswith("<") else int(start) - 1
        end = AfterPosition(int(end[1:])) if end.startswith(">") else int(end)
    elif "^" in s:
        # I'm not sure this is the right way to handle "between-position" notation
        # but this is the way RefSeq's GFF3 (GCF_902167145.1) shows as the counterpart
        # to "138683^138684"
        start, end = map(int, s.split("^"))
        start -= 1
        end -= 1
        assert start + 1 == end
    elif s.startswith("<"):
        n = int(s[1:]) - 1
        start = BeforePosition(n)
        end = ExactPosition(n + 1)
    elif s.startswith(">"):
        n = int(s[1:]) - 1
        start = ExactPosition(n)
        end = AfterPosition(n + 1)
    else:
        start = int(s) - 1
        end = start + 1

    strand = +1 if is_positive_strand else -1
    return SimpleLocation(start, end, strand=strand)


def _parse_loc(s: str) -> Union[SimpleLocation, CompoundLocation]:
    """Parse location item as either SimpleLocation or CompoundLocation

    >>> _parse_loc("10")
    SimpleLocation(ExactPosition(9), ExactPosition(10), strand=1)

    >>> _parse_loc("complement(333..>350)")
    SimpleLocation(ExactPosition(332), AfterPosition(350), strand=-1)

    >>> _parse_loc("join(1..5,100..>200)")
    CompoundLocation([SimpleLocation(ExactPosition(0), ExactPosition(5), strand=1), SimpleLocation(ExactPosition(99), AfterPosition(200), strand=1)], 'join')

    >>> _parse_loc("complement(join(1..5,100..>200))")
    CompoundLocation([SimpleLocation(ExactPosition(0), ExactPosition(5), strand=-1), SimpleLocation(ExactPosition(99), AfterPosition(200), strand=-1)], 'join')
    """
    is_positive_strand = True

    if s.startswith("complement("):
        size = len("complement(")
        s = s[size:-1]
        is_positive_strand = False

    if s.startswith("join("):
        size = len("join(")
        s = s[size:-1]
        parts = [_to_featureloc(x, is_positive_strand) for x in s.split(",")]
        res = CompoundLocation(parts)
    else:
        res = _to_featureloc(s, is_positive_strand)

    return res


def load_ddbj(path_ddbj) -> Generator[SeqRecord, None, None]:
    """Load DDBJ annotation file as a generator of SeqRecords

    [NOTE] COMMON part of annotation will be ignored because
    SeqRecord does not accomodate COMMON as its feature nicely.
    """
    skipping_entry = False
    record = None

    csvfile = open(path_ddbj)
    spamreader = csv.reader(csvfile, delimiter="\t")
    for row in spamreader:
        if row[0] == "COMMON":
            skipping_entry = True
            continue

        if skipping_entry and (not row[0]):
            continue

        if row[0]:
            skipping_entry = False
            if record is not None:
                yield record

            assert row[2].startswith("1..")
            try:
                length = int(row[2][3:])
            except:
                logging.error(f"{row} does not have location of the form 1..99")
                length = FALLBACK_SEQLEN
            record = SeqRecord(Seq(None, length), id=row[0])

        if row[1]:
            feature = SeqFeature(type=row[1])
            if row[2]:
                feature.location = _parse_loc(row[2])
            else:
                logging.warning(f"Location missing after a feature: {row}")
                continue
            if record is None:
                logging.error(f"Something is wrong: {row}")
                continue
            record.features.append(feature)

        if row[3]:
            if record is None:
                raise ValueError(f"Something is wrong {row}")
            feature = record.features[-1]
            qkey = row[3].strip()
            qval = row[4].strip()
            if qkey in feature.qualifiers:
                feature.qualifiers[qkey].append(qval)
            else:
                feature.qualifiers[qkey] = [qval]

    if record is None:
        logging.warning(f"No record in the data?")
    else:
        yield record

    csvfile.close()
