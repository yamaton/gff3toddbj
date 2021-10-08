"""
For evaluation of DDBJ annotation

"""
from typing import Iterable, Generator, Set, OrderedDict, Dict, Counter, Tuple
import argparse
import collections

from . import parser

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


_EXEC_NAME = "compare-ddbj"
__version__ = "0.0.0"

Troika = Tuple[str, str, str]


def _get_feature_name(f: SeqFeature) -> str:
    """Return feature.type for all but ncRNA feature.
    /ncRNA_class value is appended in case of ncRNA.
    """
    if f.type == "ncRNA":
        res = "ncRNA_" + f.qualifiers["ncRNA_class"][0]
    else:
        res = f.type
    return res


def get_multiset_feature_qualifier(
    records: Iterable[SeqRecord],
) -> Counter[Troika]:
    """Get a multiset of (Entry, FeatureName, QualifierKey)
    """
    accum = collections.Counter(
        (rec.id, _get_feature_name(feat), qual_key)
        for rec in records
        for feat in rec.features
        for qual_key in feat.qualifiers
    )
    return accum


def _get_location_as_str(feature: SeqFeature, ignore_correction: bool) -> str:
    s = str(feature.location)
    if ignore_correction:
        s = s.replace("<", "").replace(">", "")
    return s


def _get_multiset_locations(
    records: Iterable[SeqRecord], ignore_loc_correction: bool = False
) -> Counter[Troika]:
    """Get an ordered dictionary of "Entry" --> "Feature Name" --> a set of locations

    Args:
        records: iterable of SeqRecords
        ignore_loc_correction: Ignore ">" and "<" in location if True
    """
    res = collections.Counter(
        (rec.id, _get_feature_name(feat), _get_location_as_str(feat, ignore_loc_correction))
        for rec in records
        for feat in rec.features
    )
    return res

def get_multiset_locations_wo_correction(records: Iterable[SeqRecord]) -> Counter[Troika]:
    return _get_multiset_locations(records, True)


def get_multiset_locations(records: Iterable[SeqRecord]) -> Counter[Troika]:
    return _get_multiset_locations(records, False)



def compare(path1: str, path2: str, func) -> Tuple[Counter[Troika], Counter[Troika], Counter[Troika]]:
    """

    """
    records1 = parser.load_ddbj(path1)
    records2 = parser.load_ddbj(path2)
    multiset1 = func(records1)
    multiset2 = func(records2)

    intersect = multiset1 & multiset2
    left_only = multiset1 - multiset2
    right_only = multiset2 - multiset1
    return intersect, left_only, right_only


def main():
    parser = argparse.ArgumentParser(prog=_EXEC_NAME)
    parser.add_argument("ddbj1", help="Input DDBJ annotation 1")
    parser.add_argument("ddbj2", help="Input DDBJ annotation 2")
    parser.add_argument(
        "--name1",
        metavar="STR",
        help="Specify name of the first annotation",
        default="ann1"
    )
    parser.add_argument(
        "--name2",
        metavar="STR",
        help="Specify name of the second annotation",
        default="ann2"
    )
    parser.add_argument(
        "--log",
        default="INFO",
        metavar="STR",
        help="[debug] Choose log level from (DEBUG, INFO, WARNING, ERROR) (default: INFO).",
    )

    args = parser.parse_args()

    intersect, left_only, right_only = compare(args.ddbj1, args.ddbj2, get_multiset_locations_wo_correction)
    cnt_total = sum(map(len, [intersect, left_only, right_only]))
    cnt_correct = len(intersect)
    print("Stat w/o  location correction: {}/{} ({:.1f} %)  ... (left-only: {}, right-only: {})".format(cnt_correct, cnt_total, 100 * cnt_correct / cnt_total, len(left_only),  len(right_only)))
    with open("loc_wo_correction_left-only.txt", "w") as fout:
        for tup in left_only:
            print("\t".join(tup), file=fout)
    with open("loc_wo_correction_right-only.txt", "w") as fout:
        for tup in right_only:
            print("\t".join(tup), file=fout)

    intersect, left_only, right_only = compare(args.ddbj1, args.ddbj2, get_multiset_locations)
    cnt_total = sum(map(len, [intersect, left_only, right_only]))
    cnt_correct = len(intersect)
    print("Stat with location correction: {}/{} ({:.1f} %)  ... (left-only: {}, right-only: {})".format(cnt_correct, cnt_total, 100 * cnt_correct / cnt_total, len(left_only),  len(right_only)))
    with open("loc_correction_left-only.txt", "w") as fout:
        for tup in left_only:
            print("\t".join(tup), file=fout)
    with open("loc_correction_right-only.txt", "w") as fout:
        for tup in right_only:
            print("\t".join(tup), file=fout)

    intersect, left_only, right_only = compare(args.ddbj1, args.ddbj2, get_multiset_feature_qualifier)
    cnt_total = sum(map(len, [intersect, left_only, right_only]))
    cnt_correct = len(intersect)
    print("Stat of feature-qualifier pairs: {}/{} ({:.1f} %)  ... (left-only: {}, right-only: {})".format(cnt_correct, cnt_total, 100 * cnt_correct / cnt_total, len(left_only),  len(right_only)))
    with open("quals_left-only.txt", "w") as fout:
        for tup in left_only:
            print("\t".join(tup), file=fout)
    with open("quals_right-only.txt", "w") as fout:
        for tup in right_only:
            print("\t".join(tup), file=fout)

