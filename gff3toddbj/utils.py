from typing import Iterable, Generator, Dict, List, FrozenSet, Optional, OrderedDict, Union, Any
import logging
import re
import collections
import toml

from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, SimpleLocation, SeqFeature

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
    except FileNotFoundError:
        msg = "Failed to load metadata from {}:".format(path)
        raise FileNotFoundError(msg)

    # Change qualifier-value type to list
    for k in header_info:
        if k != "COMMON":
            for qkey, qval in header_info[k].items():
                if not isinstance(qval, list):
                    header_info[k][qkey] = [qval]

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
    location: Union[SimpleLocation, CompoundLocation],
    transl_table: int,
    phase: int = 0,
) -> bool:
    """Check if the codon starting at index_location in seq
    is start codon according to the Genetic Code transl_table.
    Location is shifted if phase is nonzero.

    >>> from Bio.Seq import Seq
    >>> seq = Seq("AATTCGAGGGG")
    >>> loc = SimpleLocation(1, 7, strand=1)
    >>> genetic_code = 11
    >>> phase = 0
    >>> has_start_codon(seq, loc, genetic_code, phase)
    True
    """
    # when phase is > 0 (i.e. codon_start > 1) at the head of joined featute,
    # there must be some problems or corrections about the start codon.
    if phase > 0:
        return False

    strand = location.strand
    start_codons = CodonTable.unambiguous_dna_by_id[transl_table].start_codons
    # Bio.Data.CodonTable says start codons are ATG, CTG, TTG,
    # but limit it to ATG because of following.
    # https://www.ddbj.nig.ac.jp/ddbj/geneticcode.html#1
    # https://www.ddbj.nig.ac.jp/faq/ja/how-to-describe-not-standard-genetic-code.html
    if transl_table == 1:
        start_codons = ["ATG"]

    codon = "XXX"
    if not isinstance(strand, int):
        raise ValueError("Cannot determine as strand = {}".format(strand))

    cds = _get_cds(seq, location)
    codon = cds[phase : phase + 3]

    if len(codon) != 3:
        msg = "Failed to access codon (loc={})".format(location)
        logging.error(msg)

    return str(codon) in start_codons


def has_stop_codon(
    seq: Seq,
    location: Union[SimpleLocation, CompoundLocation],
    transl_table: int,
    phase : int,
) -> bool:
    """Check if stop codon exists in seq.
    Stop codons correspond to the Genetic Code transl_table.

    >>> from Bio.Seq import Seq
    >>> seq = Seq("GAATGCGAGGGTAGT")
    >>> loc = SimpleLocation(2, 14, strand=1)
    >>> genetic_code = 1
    >>> has_stop_codon(seq, loc, genetic_code, phase=0)
    True
    """
    # if the length is not multiple of 3, stop codon does not exist.
    if (len(location) - phase) % 3 > 0:
        return False

    strand = location.strand
    stop_codons = CodonTable.unambiguous_dna_by_id[transl_table].stop_codons

    codon = "XXX"
    if not isinstance(strand, int):
        raise ValueError("Cannot determine as strand = {}".format(strand))

    cds = _get_cds(seq, location)
    codon = cds[-3:]

    if len(codon) != 3:
        msg = "Failed to access codon (loc={})".format(location)
        logging.error(msg)

    return str(codon) in stop_codons


def _get_cds(seq: Seq, location: Union[SimpleLocation, CompoundLocation]) -> Seq:
    """wrapper of SeqFeature.extract()

    CompoundLocation.part is ordered by (loc.start, loc.end)
    in this code base, hence loc.extract() gives wrong answer if loc is
    a CompoundLocation with strand (-1).

    Ref: https://github.com/biopython/biopython/issues/570
    """
    if isinstance(location, SimpleLocation) or location.strand > 0:
        f = SeqFeature(location, type="tmp")
    else:
        location = sorted(location.parts, key=lambda x: (x.end, x.start), reverse=True)
        x = CompoundLocation(location)
        f = SeqFeature(x, type="domain")
    return f.extract(seq)


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


def to_loglevel(s: str) -> int:
    """
    Convert string "INFO" to logging.INFO, "DEBUG" to logging.DEBUG, etc.

    >>> to_loglevel("info")
    20
    >>> to_loglevel("30")
    30
    """
    try:
        x = int(s)
        return x
    except ValueError:
        pass

    s = s.upper()
    valid = ("DEBUG", "INFO", "WARN", "WARNING", "ERROR", "FATAL")
    if s in valid:
        res = getattr(logging, s)
    else:
        res = logging.INFO
    return res


def get_attribute_keys(subtree: OrderedDict[str, Any]) -> List[str]:
    """
    Get associated attribute keys in the translation subtable.
    Returns an empty list if the translation input is "type" only.
    """
    if set(subtree.keys()) < TARGET_KEYS_IN_TRANSLATION:
        return []
    return list(subtree.keys())


get_attribute_values = get_attribute_keys


def to_lowercase_keys(table: OrderedDict[str, Any]) -> OrderedDict[str, Any]:
    """
    Convert all dictionary keys to lowercase. It can be nested.
    """
    res = collections.OrderedDict()

    for key in table:
        val = table[key]
        if ("feature_key" in val) or ("qualifier_key" in val):
            res[key.lower()] = val
        elif key == "__ANY__":
            res["__ANY__"] = to_lowercase_keys(val)
        elif isinstance(val, OrderedDict):
            res[key.lower()] = to_lowercase_keys(val)
        else:
            logging.error("Something is wrong with the dictionary!")
            return table

    return res
