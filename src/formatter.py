from typing import Any, Dict, List, Tuple, Iterable, Union, Generator
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
import yaml

Location = Union[FeatureLocation, CompoundLocation]


def format_location(loc: Location) -> str:
    """Format location"""

    def _format_forward_segment(loc: FeatureLocation) -> str:
        return "{start}..{end}".format(start=loc.start, end=loc.end)

    def _format_forward_compound_segments(loc: CompoundLocation) -> str:
        parts = ",".join(_format_forward_segment(x) for x in loc.parts)
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


def record_to_ddbj_table(rec: SeqRecord, limit_to_ddbj=True) -> List[List[str]]:
    """Convert GFF.SeqRecord into DDBJ annotation table format

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
        for row in _gen_ddbj_table_feature_rows(feature, limit_to_ddbj)
    ]
    table[0][0] = rec.id
    return table

def _gen_ddbj_table_feature_rows(feature: SeqFeature, limit_to_ddbj=True) -> Generator[List[str], None, None]:
    """Convert SeqFeature into DDBJ annotation table format
    """
    is_first_line = True
    for (qualifier_key, values) in feature.qualifiers.items():
        values = values if isinstance(values, list) else [values]
        for qualifier_value in values:
            xs = ["" for _ in range(5)]
            if is_first_line:
                is_first_line = False
                xs[1] = feature.type
                if feature.location is not None:
                    xs[2] = format_location(feature.location)
            xs[3] = qualifier_key
            xs[4] = str(qualifier_value)
            yield xs

        if hasattr(feature, "sub_features"):
            for subfeature in feature.sub_features:
                yield from _gen_ddbj_table_feature_rows(subfeature)


def table_to_tsv(table: List[List[str]]) -> str:
    """Convert from table (list of list) to tab-separated variables (TSV)"""
    return "\n".join("\t".join(items) for items in table)


def load_common(path) -> SeqRecord:
    """Create COMMON entry as SeqRecord"""
    with open(path, "r") as f:
        header_info = yaml.safe_load(f)

    features = [
        SeqFeature(type=key, qualifiers=xs) for (key, xs) in header_info.items()
    ]
    record = SeqRecord("", id="COMMON", features=features)
    return record

