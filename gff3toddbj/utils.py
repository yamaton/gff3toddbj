from typing import Iterable, Generator, Dict, List, FrozenSet, OrderedDict, Union, Any
import logging
import pprint
import re
import collections
import toml

from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature


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
}

METADATA_KEYS = {"COMMON", "source"}

INVALID_LETTERS_AS_SEQID = r'[=|>" \[\]]'
INVALID_PATTERN_AS_SEQID = re.compile(INVALID_LETTERS_AS_SEQID)


def load_rules(path: str) -> Dict[str, FrozenSet[str]]:
    """Load DDBJ feature-qualifier rules in TOML"""
    with open(path, "r") as f:
        rules = toml.load(f, _dict=collections.OrderedDict)
    for feature_key in rules:
        rules[feature_key] = frozenset(rules[feature_key])

    logging.debug("rules:\n{}".format(pprint.pformat(rules)))
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
    rec: SeqRecord,
    location: Union[FeatureLocation, CompoundLocation],
    transl_table: int,
    phase=0,
) -> bool:
    """Check if the codon starting at index_location in seq
    is start codon according to the Genetic Code transl_table.
    Location is shifted if phase is nonzero.

    >>> from Bio.Seq import Seq
    >>> seq = Seq("GAATTCGAGGGG")
    >>> rec = SeqRecord(seq)
    >>> loc = FeatureLocation(1, 7, strand=1)
    >>> has_start_codon(rec, loc, 11, phase=1)
    True
    """
    strand = location.strand
    start_codons = CodonTable.unambiguous_dna_by_id[transl_table].start_codons
    seq = rec.seq

    if strand > 0:
        p = location.start + phase
        codon = seq[p : p + 3]
    elif strand < 0:
        p = location.end - phase
        codon = seq[p - 3 : p].reverse_complement()
    else:
        raise ValueError("Cannot determine as strand = {}".format(strand))
    return codon in start_codons


def check_cds(
    rec: SeqRecord, fasta_record: Dict[str, SeqRecord], transl_table: int
) -> None:
    """
    Find an inconsistency in codon_start qualifier and displays suggestion as WARNING log.
    """
    if rec.id not in fasta_record:
        logging.error("Following SeqID from GFF3 not found in FASTA: {}".format(rec.id))
        return

    fasta_seq = fasta_record[rec.id]

    def helper(features: List[SeqFeature]):
        for f in features:
            if f.type == "CDS" and int(f.qualifiers["codon_start"][0]) != 1:
                cs_list = f.qualifiers["codon_start"]
                codon_start = cs_list[0]
                msg = "f.qualifiers['codon_start'] = {}".format(cs_list)
                assert len(cs_list) == 1, msg
                assert isinstance(cs_list[0], int)
                try:
                    shift = codon_start - 1  # 1-based indexing to shift
                    if not has_start_codon(fasta_seq, f.location, transl_table, shift):
                        msg = "A start codon NOTFOUND with shift={}, transl_table={}".format(
                            shift, transl_table
                        )
                        logging.warning(msg)
                        logging.warning("   CDS: location = {}".format(f.location))
                        if has_start_codon(fasta_seq, f.location, transl_table):
                            template = "Found a start codon with zero shift (while codon_start={} currently). Fix /codon_start?"
                            msg = template.format(codon_start)
                            logging.warning(msg)
                        else:
                            msg = "A start codon NOTFOUND. Fix the feature location?"
                            logging.warning(msg)
                    else:
                        temp = "Found a start codon as indicated by codon_start={}"
                        msg = temp.format(codon_start)
                        logging.info(msg)
                except:
                    raise ValueError(
                        "fix_cds...codon_start = {}".format(f.qualifiers["codon_start"])
                    )

            if hasattr(f, "sub_features"):
                helper(f.sub_features)

    helper(rec.features)


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
