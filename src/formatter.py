from typing import Any, Dict, List, Tuple, Iterable, Union
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, FeatureLocation

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


def to_ddbj_table(rec: SeqRecord) -> List[List[str]]:
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
    lines = sum(
        len(vals)
        for feature in rec.features
        for vals in feature.qualifiers.values()
    )
    print("lines: {}".format(lines))
    table = [["" for _ in range(5)] for _ in range(lines)]

    idx = 0
    for feature in rec.features:
        for (qualifier_idx, (qualifier_key, values)) in enumerate(feature.qualifiers.items()):
            for qualifier_value in values:
                if qualifier_idx == 0:
                    table[idx][1] = feature.type
                    if feature.location is not None:
                        table[idx][2] = format_location(feature.location)
                table[idx][3] = qualifier_key
                table[idx][4] = qualifier_value
                idx += 1
    return table
