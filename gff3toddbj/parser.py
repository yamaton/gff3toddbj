"""
Sketch of DDBJ annotation parser

"""
from typing import Generator, Union
import csv
import logging

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, BeforePosition, AfterPosition


FALLBACK_SEQLEN = 100

def _to_featureloc(s: str, is_positive_strand: bool) -> FeatureLocation:
    """Parse simple location of the form like

    >>> _to_featureloc("123", False)
    FeatureLocation(ExactPosition(123), ExactPosition(123), strand=-1)

    >>> _to_featureloc("<23..35", True)
    FeatureLocation(BeforePosition(23), ExactPosition(35), strand=1)
    """
    if ".." in s:
        start, end = s.split("..")
        start = BeforePosition(int(start[1:])) if start.startswith("<") else int(start)
        end = AfterPosition(int(end[1:])) if end.startswith(">") else int(end)
    else:
        start = end = int(s)

    strand = +1 if is_positive_strand else -1
    return FeatureLocation(start, end, strand=strand)


def _parse_loc(s: str) -> Union[FeatureLocation, CompoundLocation]:
    """Parse location item as either FeatureLocation or CompoundLocation

    >>> _parse_loc("10")
    FeatureLocation(ExactPosition(10), ExactPosition(10), strand=1)

    >>> _parse_loc("complement(333..>350)")
    FeatureLocation(ExactPosition(333), AfterPosition(350), strand=-1)

    >>> _parse_loc("join(1..5,100..>200)")
    CompoundLocation([FeatureLocation(ExactPosition(1), ExactPosition(5), strand=1), FeatureLocation(ExactPosition(100), AfterPosition(200), strand=1)], 'join')

    >>> _parse_loc("complement(join(1..5,100..>200))")
    CompoundLocation([FeatureLocation(ExactPosition(1), ExactPosition(5), strand=-1), FeatureLocation(ExactPosition(100), AfterPosition(200), strand=-1)], 'join')
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

    """
    skipping_common = False
    record = None

    csvfile = open(path_ddbj)
    spamreader = csv.reader(csvfile, delimiter="\t")
    for row in spamreader:
        if row[0] == "COMMON":
            skipping_common = True
            continue

        if skipping_common and (not row[0]):
            continue

        if row[0]:
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
                logging.error(f"Something is wrong: {row}")
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
