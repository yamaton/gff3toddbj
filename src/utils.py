from typing import Iterable, Generator, Dict, Any
from io import UnsupportedOperation
import toml
from Bio.SeqFeature import SeqFeature


def load_header_info(path) -> Dict[str, Dict[str, Any]]:
    """Create COMMON entry as SeqRecord"""
    try:
        with open(path, "r") as f:
            header_info = toml.load(f)
    except:
        raise UnsupportedOperation(
            "COMMON input other than TOML is not implemented yet!"
        )
    return header_info


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
