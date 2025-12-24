"""
For evaluation of DDBJ annotation

"""
from typing import Callable, Hashable, Iterable, List, Counter, Tuple, Union
import argparse
import collections

from . import parser

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, SimpleLocation, SeqFeature


_EXEC_NAME = "compare-ddbj"
__version__ = "0.0.4"

Troika = Tuple[str, str, str]


def remove_short_introns(record: SeqRecord, gap_size_to: int=5) -> None:
    """Patch parts of `CompundLocation` together
    if intron gap is shorter than `min_gap_size`.
    """

    def _runner(features: List[SeqFeature]) -> None:
        for f in features:
            if isinstance(f.location, CompoundLocation):
                f.location = _fix(f.location)
            if hasattr(f, "sub_features"):
                _runner(f.sub_features)

    def _fix(loc: CompoundLocation) -> Union[CompoundLocation, SimpleLocation]:
        parts = sorted(loc.parts, key=lambda x: (x.start, x.end))
        acc = []
        for curr in parts:
            if not acc:
                acc.append(curr)
            elif curr.start - acc[-1].end - 1 > gap_size_to:
                acc.append(curr)
            else:
                prev = acc.pop()
                assert prev.strand == curr.strand, "Unmatched strand!"
                x = SimpleLocation(prev.start, curr.end, strand=curr.strand)
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
    elif f.type == "regulatory":
        res = "regulatory_" + f.qualifiers["regulatory_class"][0]
    elif f.type == "repeat_region":
        qual_value = f.qualifiers.get("rpt_type", ["other"])[0]
        res = "repeat_region_" + qual_value
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

def get_multiset_entry_name(records: Iterable[SeqRecord]) -> Counter[Tuple[str]]:
    return collections.Counter((rec.id, ) for rec in records)


def compare(records1: List[SeqRecord], records2: List[SeqRecord], func) -> Tuple[Counter[Troika], Counter[Troika], Counter[Troika]]:
    """

    """
    multiset1 = func(records1)
    multiset2 = func(records2)

    intersect = multiset1 & multiset2
    left_only = multiset1 - multiset2
    right_only = multiset2 - multiset1
    return intersect, left_only, right_only


def compare_and_report(records1: List[SeqRecord], records2: List[SeqRecord], f: Callable[[Iterable[SeqRecord]], Counter[Hashable]], title: str, file_prefix: str):
    intersect, left_only, right_only = compare(records1, records2, f)
    cnt_match, cnt_leftonly, cnt_rightonly = map(lambda multiset: sum(multiset.values()), [intersect, left_only, right_only])

    print(title)
    print("    Left  mismatching: {} / {} \t({:.2f} %)".format(cnt_leftonly, cnt_match + cnt_leftonly, 100 * cnt_leftonly / (cnt_match + cnt_leftonly)))
    print("    Right mismatching: {} / {} \t({:.2f} %)".format(cnt_rightonly, cnt_match + cnt_rightonly, 100 * cnt_rightonly / (cnt_match + cnt_rightonly)))

    with open(file_prefix + "left-only.txt", "w") as fout:
        for tup, count in left_only.items():
            for _ in range(count):
                print("\t".join(tup), file=fout)
    with open(file_prefix + "right-only.txt", "w") as fout:
        for tup, count in right_only.items():
            for _ in range(count):
                print("\t".join(tup), file=fout)


def set_accession_as_entry_name(rec: SeqRecord):
    """Replace the ID with accession part if it contains '|'

    According to FASTA format, `dbj|accession|locus` is the one in DDBJ.
    https://en.wikipedia.org/wiki/FASTA_format
    """
    if "|" in rec.id:
        accession = rec.id.split("|")[1]
        rec.id = accession


def main():
    argparser = argparse.ArgumentParser(prog=_EXEC_NAME)
    argparser.add_argument("ddbj1", help="Input DDBJ annotation 1")
    argparser.add_argument("ddbj2", help="Input DDBJ annotation 2")
    argparser.add_argument(
        "--no-rename-entry",
        default=False,
        action="store_true",
        help="Disable renaming of entries by extracting accession part assuming dbj|accession|locus format",
    )
    argparser.add_argument(
        "--patch-features",
        default=False,
        action="store_true",
        help="Remove short (< 10bp) introns by patching feature gaps",
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
        if args.patch_features:
            remove_short_introns(rec)
        if not args.no_rename_entry:
            set_accession_as_entry_name(rec)

    for rec in records2:
        if args.patch_features:
            remove_short_introns(rec)
        if not args.no_rename_entry:
          set_accession_as_entry_name(rec)

    print("------------------------------------------------------------")
    print("!! Other stats are affected if entry names are different !!")
    compare_and_report(records1, records2, get_multiset_entry_name, "    Stat of entry names ...", "entry_names_")
    print("------------------------------------------------------------")
    print()
    compare_and_report(records1, records2, get_multiset_locations_wo_correction, "Stat w/o location correction:", "loc_wo_correction_")
    compare_and_report(records1, records2, get_multiset_locations, "Stat with location correction:", "loc_corrected_")
    compare_and_report(records1, records2, get_multiset_feature_qualifier, "Stat of feature-qualifier pairs:", "quals_")
