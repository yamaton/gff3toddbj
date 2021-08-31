from io import UnsupportedOperation
import logging
from typing import Any, Dict, FrozenSet, List, Tuple, Iterable, Union, Generator
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
import yaml

Location = Union[FeatureLocation, CompoundLocation]


def format_location(loc: Location) -> str:
    """Format location"""

    def _format_forward_segment(loc: FeatureLocation) -> str:
        start = min(loc.start, loc.end)
        end = max(loc.start, loc.end)
        return "{}..{}".format(start, end)

    def _format_forward_compound_segments(loc: CompoundLocation) -> str:
        locs = sorted(
            loc.parts, key=lambda feature_loc: min(feature_loc.start, feature_loc.end)
        )
        parts = ",".join(_format_forward_segment(x) for x in locs)
        s = "join({})".format(parts)
        return s

    if isinstance(loc, FeatureLocation):
        s = _format_forward_segment(loc)
    elif isinstance(loc, CompoundLocation):
        s = _format_forward_compound_segments(loc)
    else:
        raise ValueError("Wrong type: type(loc) = {}".format(type(loc)))

    if loc.strand == -1:
        s = "complement({})".format(s)
    return s


def table_to_tsv(table: List[List[str]]) -> str:
    """Convert from table (list of list) to tab-separated variables (TSV)"""
    return "\n".join("\t".join(items) for items in table)


def load_common(path) -> SeqRecord:
    """Create COMMON entry as SeqRecord"""
    try:
        with open(path, "r") as f:
            header_info = yaml.safe_load(f)
    except:
        raise UnsupportedOperation("COMMON input other than YAML is not implemented yet!")

    features = [
        SeqFeature(type=key, qualifiers=xs) for (key, xs) in header_info.items()
    ]
    record = SeqRecord("", id="COMMON", features=features)
    return record


def load_rules(path: str) -> Dict[str, FrozenSet[str]]:
    """Load DDBJ feature-qualifier rules in YAML"""
    with open(path, "r") as f:
        rules = yaml.safe_load(f)
    for feature_key in rules:
        rules[feature_key] = frozenset(rules[feature_key])
    return rules


class DDBJFormatter(object):
    """Format SeqRecord to DDBJ annotation table."""

    def __init__(self, common_header_path, ddbj_rule_path: str):
        self.common = load_common(common_header_path)
        self.rules = load_rules(ddbj_rule_path)

    def is_allowed(self, feature_key, qualifier_key):
        return (feature_key in self.rules) and (
            qualifier_key in self.rules[feature_key]
        )

    def to_ddbj_table(self, rec: SeqRecord, ignore_ddbj=False) -> List[List[str]]:
        """Convert GFF.SeqRecord into DDBJ annotation table format

        Args:
            rec (SeqRecord): Annotated record to be converted
            ignore_ddbj (bool): Ignore DDBJ-recommended feature-qualifier filtering if True

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
        table = [
            row
            for feature in rec.features
            for row in self._gen_ddbj_table_feature_rows(feature, ignore_ddbj)
        ]
        table[0][0] = rec.id
        return table

    def _gen_ddbj_table_feature_rows(
        self, feature: SeqFeature, ignore_ddbj=True
    ) -> Generator[List[str], None, None]:
        """Convert SeqFeature into DDBJ annotation table format"""
        is_first_line = True
        for (qualifier_key, values) in feature.qualifiers.items():
            values = values if isinstance(values, list) else [values]
            for qualifier_value in values:
                xs = ["" for _ in range(5)]
                feature_key = feature.type
                is_keeping = self.is_allowed(feature_key, qualifier_key)
                if ignore_ddbj or is_keeping:
                    if is_first_line:
                        is_first_line = False
                        xs[1] = feature_key
                        if feature.location is not None:
                            xs[2] = format_location(feature.location)
                    xs[3] = qualifier_key
                    xs[4] = str(qualifier_value)
                    yield xs
                elif not ignore_ddbj and not is_keeping:
                    troika = (feature_key, qualifier_key, qualifier_value)
                    logging.warn("skipping:  {}, {}, {}\n".format(*troika))

            if hasattr(feature, "sub_features"):
                for subfeature in feature.sub_features:
                    yield from self._gen_ddbj_table_feature_rows(
                        subfeature, ignore_ddbj
                    )

    def run(self, records: Iterable[SeqRecord], ignore_ddbj: bool) -> Generator[str, None, None]:
        """Format records and generate string line by line."""
        table_header = self.to_ddbj_table(self.common, ignore_ddbj=True)
        for rows in table_header:
            yield "\t".join(rows)

        for rec in records:
            tbl = self.to_ddbj_table(rec, ignore_ddbj)
            for rows in tbl:
                yield "\t".join(rows)
