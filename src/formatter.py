from typing import Any, Dict, List, Tuple, Iterable, Union
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
import yaml

Location = Union[FeatureLocation, CompoundLocation]

def format_location(loc: Location) -> str:
    """Format location
    """
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


def record_to_ddbj_table(rec: SeqRecord) -> List[List[str]]:
    """Convert GFF.SeqRecord into DDBJ annotation table format

    DDBJ annotation table format is TSV (tab-separated variables) with 5 columns
        - Column 1: Sequence ID
        - Column 2: Feature Key
        - Column 3: Location
        - Column 4: Qualifier Key
        - Column 5: Qualifier Value

    Many entries are blank.

    Example:
    [
        ["CLN01", "source", "1..12297"                         , "organism"   ,  "Mus musculus"  ],
        [       ,         ,                                    , "mol_type"   ,  "genomic DNA"   ],
        [       ,         ,                                    , "clone"      ,  "PC0110"        ],
        [       ,         ,                                    , "chromosome" ,  "8"             ],
        [       ,  "CDS"  , "join(<1..456,609..879,1070..1213)", "product"    ,  "protein kinase"],
        [       ,         ,                                    , "codon_start",  "2"             ],
    ]

    """
    num_lines = sum(
        len(vals) if isinstance(vals, list) else 1
        for feature in rec.features
        for vals in feature.qualifiers.values()
    )
    table = [["" for _ in range(5)] for _ in range(num_lines)]
    table[0][0] = str(rec.id)

    idx = 0
    for feature in rec.features:
        is_first_line = True
        for (qualifier_key, values) in feature.qualifiers.items():
            values = values if isinstance(values, list) else [values]
            for qualifier_value in values:
                if is_first_line:
                    is_first_line = False
                    table[idx][1] = feature.type
                    if feature.location is not None:
                        table[idx][2] = format_location(feature.location)
                table[idx][3] = qualifier_key
                table[idx][4] = str(qualifier_value)
                idx += 1

    return table


def table_to_tsv(table: List[List[str]]) -> str:
    """Convert from table (list of list) to tab-separated variables (TSV)
    """
    return "\n".join("\t".join(items) for items in table)


def load_common(path) -> SeqRecord:
    """Create COMMON entry as SeqRecord
    """
    with open(path, "r") as f:
        header_info = yaml.safe_load(f)

    features = [SeqFeature(type=key, qualifiers=xs) for (key, xs) in header_info.items()]
    record = SeqRecord("", id="COMMON", features=features)
    return record
