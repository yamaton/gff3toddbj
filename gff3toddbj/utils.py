from typing import Iterable, Generator, Dict, List, FrozenSet, Optional, OrderedDict, Union, Any
import logging
import re
import collections
import toml

from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

# Dummy name of exon after joining them to suppress outputs
DUMMY_ORIGINALLY_EXON = "__exon"

# Supported COMMON Features
COMMON_FEATURES = {
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
}

METADATA_KEYS = {"COMMON", "source", "assembly_gap"}

INVALID_LETTERS_AS_SEQID = r'[=|>" \[\]]'
INVALID_PATTERN_AS_SEQID = re.compile(INVALID_LETTERS_AS_SEQID)

TARGET_KEYS_IN_TRANSLATION = {
    "feature_key",
    "qualifier_key",
    "qualifier_value",
    "qualifier_value_prefix",
    "attribute_value",
}


def load_rules(path: str) -> Dict[str, FrozenSet[str]]:
    """Load DDBJ feature-qualifier rules in TOML"""
    with open(path, "r") as f:
        rules = toml.load(f, _dict=collections.OrderedDict)
    for feature_key in rules:
        rules[feature_key] = frozenset(rules[feature_key])

    return rules


def load_metadata_info(path) -> OrderedDict[str, OrderedDict[str, Any]]:
    """Load metadata as dictionary from TOML file."""
    try:
        with open(path, "r") as f:
            header_info = toml.load(f, _dict=collections.OrderedDict)
    except:
        msg = "Failed to load metadata from {}:".format(path)
        raise FileNotFoundError(msg)

    # Change qualifier type to list
    for k in header_info:
        if k != "COMMON":
            for qkey, qval in header_info[k].items():
                if not isinstance(qval, list):
                    header_info[k][qkey] = [qval]

    # [FIXME] Disabled till DDBJ rule is counted in the validation
    # validate_metadata_keys(header_info)

    return header_info


def validate_metadata_keys(d: Dict[str, Dict[str, Any]]) -> None:
    """Check if metadata has proper table keys.

    [FIXME] Need to take care of DDBJ features used in main.py
    """
    keys = set(d.keys())
    msg = "Some meta-info keys are invalid: {}".format(keys - METADATA_KEYS)
    assert keys.issubset(METADATA_KEYS), msg

    common_keys = set(d["COMMON"])
    msg2 = "Some COMMON keys are invalid: {}".format(common_keys - COMMON_FEATURES)
    assert common_keys.issubset(COMMON_FEATURES), msg2


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


def is_invalid_as_seqid(name: str) -> bool:
    """Check if `name` has a invalid characters as SeqID.
    [NOTE] DDBJ annotation disallows names including =|>" []

    >>> is_invalid_as_seqid("hello|world")
    True

    >>> is_invalid_as_seqid("baba+keke")
    False
    """
    return bool(INVALID_PATTERN_AS_SEQID.search(name))


def has_start_codon(
    seq: Seq,
    location: Union[FeatureLocation, CompoundLocation],
    transl_table: int,
    phase: int = 0,
) -> bool:
    """Check if the codon starting at index_location in seq
    is start codon according to the Genetic Code transl_table.
    Location is shifted if phase is nonzero.

    >>> from Bio.Seq import Seq
    >>> seq = Seq("AATTCGAGGGG")
    >>> loc = FeatureLocation(1, 7, strand=1)
    >>> table = 11
    >>> phase = 0
    >>> has_start_codon(seq, loc, table, phase)
    True
    """
    strand = location.strand
    start_codons = CodonTable.unambiguous_dna_by_id[transl_table].start_codons

    if strand > 0:
        p = location.start.position + phase
        codon = seq[p : p + 3]
    elif strand < 0:
        p = location.end.position - phase
        codon = seq[p - 3 : p].reverse_complement()
    else:
        raise ValueError("Cannot determine as strand = {}".format(strand))
    return codon in start_codons


def has_stop_codon(
    seq: Seq,
    location: Union[FeatureLocation, CompoundLocation],
    transl_table: int,
) -> bool:
    """Check if stop codon exists in seq.
    Stop codons correspond to the Genetic Code transl_table.

    >>> from Bio.Seq import Seq
    >>> seq = Seq("GAATGCGAGGGTAG")
    >>> loc = FeatureLocation(2, 14, strand=1)
    >>> table = 1
    >>> has_stop_codon(seq, loc, table)
    True
    """
    strand = location.strand
    stop_codons = CodonTable.unambiguous_dna_by_id[transl_table].stop_codons

    if strand > 0:
        p = location.end.position
        codon = seq[p - 3 : p]
    elif strand < 0:
        p = location.start.position
        codon = seq[p : p + 3].reverse_complement()
    else:
        raise ValueError("Cannot determine as strand = {}".format(strand))
    return codon in stop_codons




def debug_checker(
    recs: List[SeqRecord], tag="", qual_key="note", qual_value="ID:cds-XP_008596228.1"
) -> None:
    def _helper(features: List[SeqFeature]) -> None:
        for f in features:
            if qual_value in f.qualifiers.get(qual_key, []):
                logging.info(
                    "{} f.qualifiers['{}'] = {}".format(
                        tag, qual_key, f.qualifiers[qual_key]
                    )
                )
            if hasattr(f, "sub_features"):
                _helper(f.sub_features)

    for rec in recs:
        _helper(rec.features)



def get_attribute_keys(subtree: OrderedDict[str, Any]) -> List[str]:
    """
    Get associated attribute keys in the translation subtable.
    Returns an empty list if the translation input is "type" only.
    """
    if set(subtree.keys()) < TARGET_KEYS_IN_TRANSLATION:
        return []
    return subtree.keys()

