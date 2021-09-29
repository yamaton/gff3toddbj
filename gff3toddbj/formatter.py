import collections
import logging
from typing import Any, DefaultDict, FrozenSet, List, Optional, Tuple, Iterable, Union, Generator

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, ExactPosition, FeatureLocation, SeqFeature

from . import utils

Location = Union[FeatureLocation, CompoundLocation]

def format_location(loc: Location) -> str:
    """Format location

    Note that BioPython's indexing is 0-based, left-inclusive, right-exclusive
    while Annotation is 1-based, both-inclusive.
    """

    def _format_part(loc: FeatureLocation) -> str:
        PositionType = type(loc.start)
        start = str(PositionType(loc.start.position + 1))
        end = str(loc.end)
        if start == end:
            msg = "Something is wrong with position type!"
            assert all(isinstance(x, ExactPosition) for x in (loc.start, loc.end)), msg
            s = "{}".format(start)  # convert to 1-based
        else:
            s = "{}..{}".format(start, end)  # convert to 1-based, both-inclusive
        return s

    def _format_compound(loc: CompoundLocation) -> str:
        locs = sorted(
            loc.parts, key=lambda feature_loc: min(feature_loc.start, feature_loc.end)
        )
        parts = ",".join(_format_part(x) for x in locs)
        s = "join({})".format(parts)
        return s

    if isinstance(loc, FeatureLocation):
        s = _format_part(loc)
    elif isinstance(loc, CompoundLocation):
        s = _format_compound(loc)
    else:
        raise ValueError("Wrong type: type(loc) = {}".format(type(loc)))

    if loc.strand == -1:
        s = "complement({})".format(s)
    return s


def table_to_tsv(table: List[List[str]]) -> str:
    """Convert from table (list of list) to tab-separated variables (TSV)"""
    return "\n".join("\t".join(items) for items in table)




def get_common(header_info: DefaultDict[str, DefaultDict[str, Any]]) -> Optional[SeqRecord]:
    """Load header info as SeqRecord for COMMON entry

    [NOTE] COMMON entry can take "source" feature in DDBJ.
    https://www.ddbj.nig.ac.jp/ddbj/file-format.html#common
    But this "source" feature with id "source" does not carry
    location in current implementation. Need to insert 1..E
    as the location somewhere else..
    """
    features = [
        SeqFeature(type=key, qualifiers=xs)
        for (key, xs) in header_info.get("COMMON", collections.defaultdict()).items()
    ]

    record = None
    if features:
        record = SeqRecord("", id="COMMON", features=features)
    return record


class DDBJFormatter(object):
    """Format SeqRecord to DDBJ annotation table."""

    def __init__(self, header_info, ddbj_filter_path: str):
        self.common = get_common(header_info)
        self.rules = utils.load_rules(ddbj_filter_path)
        ## Counter of ignored feature keys
        self.ignored_feature_count = collections.defaultdict(int)
        ## Counter of ignored (feature key, qualifier key) pairs
        self.ignored_pair_count = collections.defaultdict(int)

    def _is_allowed_feature(self, feature_key: str) -> bool:
        """
        Check if feature_key is good for DDBJ annotation
            See DDBJ's feature-qualifier matrix for the detail.
            https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-e.pdf
        """
        return feature_key in self.rules

    def _is_allowed_pair(self, feature_key: str, qualifier_key: str) -> bool:
        """
        Check if (feature_key, qualifier_key) pair is good for DDBJ annotation
            See DDBJ's feature-qualifier matrix for the detail.
            https://www.ddbj.nig.ac.jp/assets/files/pdf/ddbj/fq-e.pdf
        """
        return self._is_allowed_feature(feature_key) and (
            qualifier_key in self.rules[feature_key]
        )

    def to_ddbj_table(self, rec: SeqRecord, ignore_rules=False) -> List[List[str]]:
        """Convert GFF.SeqRecord into DDBJ annotation table format

        Args:
            rec (SeqRecord): Annotated record to be converted
            ignore_rules (bool): Ignore DDBJ-recommended feature-qualifier choices if True

        DDBJ annotation table format is TSV (tab-separated variables) with 5 columns
            - Column 1: Sequence ID
            - Column 2: Feature Key
            - Column 3: Location
            - Column 4: Qualifier Key
            - Column 5: Qualifier Value

        Example:
            ["CLN01", "source", "1..12297"                         , "organism"   ,  "Mus musculus"  ],
            [       ,         ,                                    , "mol_type"   ,  "genomic DNA"   ],
            [       ,         ,                                    , "clone"      ,  "PC0110"        ],
        ...
        """
        # Collect discarded features by the filter
        if not ignore_rules:
            feature_keys = collections.Counter(f.type for f in utils.flatten_features(rec.features))
            for k, k_count in feature_keys.items():
                if not self._is_allowed_feature(k):
                    self.ignored_feature_count[k] += k_count  # update ignored feature count

        table = [
            row
            for feature in rec.features
            for row in self._gen_ddbj_table_feature_rows(feature, ignore_rules)
        ]

        if table:
            table[0][0] = rec.id

        return table

    def _gen_ddbj_table_feature_rows(
        self, feature: SeqFeature, ignore_rules=True
    ) -> Generator[List[str], None, None]:
        """Convert SeqFeature into DDBJ annotation table format"""
        if ignore_rules or self._is_allowed_feature(feature.type):
            is_first_line = True

            if feature.qualifiers:
                # normal case listing up qualifier values
                for (qualifier_key, values) in feature.qualifiers.items():
                    values = values if isinstance(values, list) else [values]
                    for qualifier_value in values:
                        xs = ["" for _ in range(5)]
                        feature_key = feature.type
                        is_keeping = self._is_allowed_pair(feature_key, qualifier_key)
                        if ignore_rules or is_keeping:
                            if is_first_line:
                                is_first_line = False
                                xs[1] = feature_key
                                if feature.location is not None:
                                    xs[2] = format_location(feature.location)
                            xs[3] = qualifier_key
                            xs[4] = str(qualifier_value)
                            yield xs
                        elif not ignore_rules and not is_keeping:
                            pair = (feature_key, qualifier_key)
                            self.ignored_pair_count[pair] += 1
            else:
                # special case when .qualifiers is empty
                xs = ["" for _ in range(5)]
                xs[1] = feature.type
                if feature.location is not None:
                    xs[2] = format_location(feature.location)
                yield xs

        if hasattr(feature, "sub_features"):
            for subfeature in feature.sub_features:
                yield from self._gen_ddbj_table_feature_rows(subfeature, ignore_rules)


    def _log_ignored_items(self):
        """Log tally of ignore feature keys, and feature-qualifier pairs
        """
        ## feature keys
        for fkey, cnt in self.ignored_feature_count.items():
            if fkey != utils.DUMMY_ORIGINALLY_EXON:
                msg = "[Discarded] feature ------->  {}  <------- \t (count: {})".format(fkey, cnt)
                logging.warning(msg)

        ## feature-qualifier pairs
        for (fkey, qkey), cnt in self.ignored_pair_count.items():
            msg = "[Discarded] (Feature, Qualifier) = ({}, {}) \t (count: {})".format(fkey, qkey, cnt)
            logging.warning(msg)


    def run(
        self, records: Iterable[SeqRecord], ignore_rules: bool
    ) -> Generator[str, None, None]:
        """Format records and generate string line by line."""
        # Format COMMON
        if self.common is not None:
            table_common = self.to_ddbj_table(self.common, ignore_rules=True)
            ## this is an ad-hoc handling to insert Location of "source" in "COMMON"
            for rows in table_common:
                if rows[1] == "source":
                    rows[2] = "1..E"
                yield "\t".join(rows)
        else:
            logging.warning("COMMON is unavailable")

        # Format main genome entries
        for rec in records:
            tbl = self.to_ddbj_table(rec, ignore_rules)
            for rows in tbl:
                yield "\t".join(rows)

        self._log_ignored_items()
