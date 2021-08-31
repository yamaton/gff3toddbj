from typing import Iterable, Generator
from Bio.SeqFeature import SeqFeature


def flatten_features(features: Iterable[SeqFeature]) -> Generator[SeqFeature, None, None]:
    """
    Flatten features nested with .sub_features attribute.
    [NOTE] .sub_features attribute stays to keep the function pure.
    """
    for x in features:
        yield x
        if hasattr(x, "sub_features"):
            yield from flatten_features(x.sub_features)
