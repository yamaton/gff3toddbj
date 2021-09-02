from typing import Iterable, Generator, Dict, Any
from io import UnsupportedOperation
import toml
from Bio.SeqFeature import SeqFeature

# Supported metadata keys
METADATA_COMMON_KEYS = {
    "DIVISION",
    "DATATYPE",
    "KEYWORD",
    "DBLINK",
    "SUBMITTER",
    "REFERENCE",
    "DATE",
    "TOPOLOGY",
    "COMMENT",
    "ST_COMMENT",
    "source"
}

METADATA_KEYS = {
    "COMMON",
    "source",
    "assembly_gap",
}

def load_header_info(path) -> Dict[str, Dict[str, Any]]:
    """Load metadata as dictionary from TOML file."""
    try:
        with open(path, "r") as f:
            header_info = toml.load(f)
    except:
        raise UnsupportedOperation(
            "COMMON input other than TOML is not implemented yet!"
        )

    validate_metadata_keys(header_info)

    return header_info


def validate_metadata_keys(d: Dict[str, Dict[str, Any]]) -> None:
    """Check if metadata has proper table keys."""
    keys = set(d.keys())
    msg = "Some meta-info keys are invalid: {}".format(keys - METADATA_KEYS)
    assert keys.issubset(METADATA_KEYS), msg

    common_keys = set(d["COMMON"])
    msg2 = "Some COMMON keys are invalid: {}".format(common_keys - METADATA_COMMON_KEYS)
    assert common_keys.issubset(METADATA_COMMON_KEYS), msg2


def flatten_features(
    features: Iterable[SeqFeature],
) -> Generator[SeqFeature, None, None]:
    """
    Flatten features nested with .sub_features attribute.
    [NOTE] .sub_features attribute stays to keep the function pure.
    """
    for x in features:
        yield x
        if hasattr(x, "sub_features"):
            yield from flatten_features(x.sub_features)
