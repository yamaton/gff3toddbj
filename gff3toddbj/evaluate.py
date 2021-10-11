"""
For evaluation of DDBJ annotation

"""
from typing import Callable, Iterable, List, Counter, Tuple, Union
import argparse
import collections

from . import parser

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature


_EXEC_NAME = "compare-ddbj"
__version__ = "0.0.2"

Troika = Tuple[str, str, str]


def patch_up_short_introns(record: SeqRecord, gap_size_to: int=5) -> None:
    """Fix `CompundLocation` of "CDS" and "exon" features and sub_features of a record
    if intron gap is shorter than `min_gap_size`.
    """

    def _runner(features: List[SeqFeature]) -> None:
        for f in features:
            if isinstance(f.location, CompoundLocation):
                f.location = _fix(f.location)
            if hasattr(f, "sub_features"):
                _runner(f.sub_features)

    def _fix(loc: CompoundLocation) -> Union[CompoundLocation, FeatureLocation]:
        parts = sorted(loc.parts, key=lambda x: (x.start.position, x.end.position))
        acc = []
        for curr in parts:
            if not acc:
                acc.append(curr)
            elif curr.start.position - acc[-1].end.position - 1 > gap_size_to:
                acc.append(curr)
            else:
                prev = acc.pop()
                assert prev.strand == curr.strand, "Unmatched strand!"
                x = FeatureLocation(prev.start, curr.end, strand=curr.strand)
                acc.append(x)

        if len(acc) == 1:
            return acc[0]

        return CompoundLocation(acc)

    _runner(record.features)


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



def compare(records1: List[SeqRecord], records2: List[SeqRecord], func) -> Tuple[Counter[Troika], Counter[Troika], Counter[Troika]]:
    """

    """
    multiset1 = func(records1)
    multiset2 = func(records2)

    intersect = multiset1 & multiset2
    left_only = multiset1 - multiset2
    right_only = multiset2 - multiset1
    return intersect, left_only, right_only



def compare_and_report(records1: List[SeqRecord], records2: List[SeqRecord], f: Callable[[Iterable[SeqRecord]], Counter[Troika]], title: str, file_prefix: str):
    intersect, left_only, right_only = compare(records1, records2, f)
    cnt_match, cnt_leftonly, cnt_rightonly = map(len, [intersect, left_only, right_only])

    print(title)
    print("    Left  mismatching: {} / {} \t({:.2f} %)".format(cnt_leftonly, cnt_match + cnt_leftonly, 100 * cnt_leftonly / (cnt_match + cnt_leftonly)))
    print("    Right mismatching: {} / {} \t({:.2f} %)".format(cnt_rightonly, cnt_match + cnt_rightonly, 100 * cnt_rightonly / (cnt_match + cnt_rightonly)))

    with open(file_prefix + "left-only.txt", "w") as fout:
        for tup in left_only:
            print("\t".join(tup), file=fout)
    with open(file_prefix + "right-only.txt", "w") as fout:
        for tup in right_only:
            print("\t".join(tup), file=fout)


def main():
    argparser = argparse.ArgumentParser(prog=_EXEC_NAME)
    argparser.add_argument("ddbj1", help="Input DDBJ annotation 1")
    argparser.add_argument("ddbj2", help="Input DDBJ annotation 2")
    argparser.add_argument(
        "--name1",
        metavar="STR",
        help="Specify name of the first annotation",
        default="ann1"
    )
    argparser.add_argument(
        "--name2",
        metavar="STR",
        help="Specify name of the second annotation",
        default="ann2"
    )
    argparser.add_argument(
        "--log",
        default="INFO",
        metavar="STR",
        help="[debug] Choose log level from (DEBUG, INFO, WARNING, ERROR) (default: INFO).",
    )

    args = argparser.parse_args()

    records1 = list(parser.load_ddbj(args.ddbj1))
    records2 = list(parser.load_ddbj(args.ddbj2))
    for rec in records1:
        patch_up_short_introns(rec)

    for rec in records2:
        patch_up_short_introns(rec)

    compare_and_report(records1, records2, get_multiset_locations_wo_correction, "Stat w/o location correction:", "loc_wo_correction_")
    compare_and_report(records1, records2, get_multiset_locations, "Stat with location correction:", "loc_corrected_")
    compare_and_report(records1, records2, get_multiset_feature_qualifier, "Stat of feature-qualifier pairs:", "quals_")
