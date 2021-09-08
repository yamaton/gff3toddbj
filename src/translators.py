import pprint
from typing import Any, Dict, Generator, List, Optional, Tuple, Iterable
import collections
import toml
import re
import pathlib
import gzip
import logging
import urllib
import tempfile
import sys
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation
from BCBio import GFF

import utils

Interval = Tuple[int, int]


def load_gff3_as_seqrecords(filepath, unquoting=False) -> List[SeqRecord]:
    """Load GFF3 as iterable of SeqRecord

    Args:
        filepath: path to load GFF3
        unquoting (bool): unquote characters.

    Unquoting means converesions like this:
        from %3B to ;
        from %2C to ,
    """
    ext = pathlib.Path(filepath).suffix

    if unquoting:
        # Save unquoted content to tempfile before feeding GFF.parse()
        with tempfile.TemporaryFile(mode="w+") as ftemp:
            if ext == ".gz":
                with gzip.open(filepath, "rt") as fin:
                    ftemp.write(urllib.parse.unquote(fin.read()))
            else:
                with open(filepath, "r") as fin:
                    ftemp.write(urllib.parse.unquote(fin.read()))
            ftemp.seek(0)
            recs = list(GFF.parse(ftemp))
    else:
        # If unquoting is unnecessary
        if ext == ".gz":
            with gzip.open(filepath, "rt") as fin:
                recs = list(GFF.parse(fin))
        else:
            recs = list(GFF.parse(filepath))

    return recs


def load_fasta_as_seq(filepath) -> Dict[str, SeqRecord]:
    """
    Load FASTA file as Seq
    """
    p = pathlib.Path(filepath)
    if p.suffix == ".gz":
        with gzip.open(filepath, "rt") as f:
            recs = {seq.id: seq.upper() for seq in Bio.SeqIO.parse(f, "fasta")}
    else:
        recs = {seq.id: seq.upper() for seq in Bio.SeqIO.parse(filepath, "fasta")}
    return recs


def load_toml_tables(filepath) -> Dict[str, Any]:
    """
    Load TOML as python dictionary
    """
    with open(filepath) as fp:
        d = toml.load(fp)

    logging.debug("TOML table:\n{}".format(pprint.pformat(d)))
    return d


def _merge_dicts(*dict_args) -> Dict:
    """
    Given any number of dictionaries, shallow copy and merge into a new dict,
    precedence goes to key-value pairs in latter dictionaries.

    >>> merge_dicts({'us', 1, 'canada, 1}, {'egypt': 20}, {'greece': 30, 'netherlands', 31})
    {'us', 1, 'canada, 1, 'egypt': 20, 'greece': 30, 'netherlands', 31}
    """
    res = {}
    for d in dict_args:
        res.update(d)
    return res


def _get_assembly_gap(seq: Seq, qualifiers: Dict[str, Any]) -> List[SeqFeature]:
    """
    Get assembly_gap features from seq.
    """
    locs = [
        FeatureLocation(start, end, strand=1)
        for (start, end) in _get_assembly_gap_locations(seq)
    ]

    features = []
    for loc in locs:
        f = SeqFeature(loc, type="assembly_gap", qualifiers=qualifiers)
        if f.qualifiers.get("estimated_length", None) == ["<COMPUTE>"]:
            length = loc.end - loc.start + 1
            f.qualifiers["estimated_length"] = [length]
        features.append(f)
    return features


def _get_assembly_gap_locations(seq: Seq) -> List[Interval]:
    """
    Get assembly_gap locations from seq as a list of intervals.

    [NOTE] This interval format [begin, end) is for Biopython
    such that 0-based, left-inclusive and right-exclusive.

    >>> s = Seq("atatnnngattacanccc")
    >>> _get_assembly_gap_locations(s)
    [(4, 7), (14, 15)]
    """
    s = str(seq)
    patt = re.compile("N+")
    matches = patt.finditer(s)

    segments = []
    for m in matches:
        tup = m.span()  # 0-based, left-inclusive, right-exclusive
        segments.append(tup)

    return segments


def _get_source(length: int, source_qualifiers: Dict[str, Any]) -> SeqFeature:
    """Create "source" feature"""
    loc = FeatureLocation(0, length, strand=1)  # 0-based and [0, length)
    return SeqFeature(loc, type="source", qualifiers=source_qualifiers)


class RenameQualifiers(object):
    """Translate annotation qualifiers according to table given in JSON.

    Meant for renaming GFF3 attributes to DDBJ annotation qualifier keys.
    JSON should have the form like following.

    ```
    {
        "ID": {
            "target": "note",
            "prefix": "ID:",
        },
        "Name": {
            "target": "",
        }
    }
    ```

    Empty "target" value means the name key is dropped.
    """

    def __init__(self, filepath: str, locus_tag_prefix: str):
        self.path = filepath
        self.locus_tag_prefix = locus_tag_prefix
        self.trans_table = load_toml_tables(filepath)

        # modify the table by inserting locus_tag prefix
        assert "locus_tag" in self.trans_table
        self.trans_table["locus_tag"]["prefix"] = locus_tag_prefix

    def run(self, record: SeqRecord) -> SeqRecord:
        """Modifies record according to the GFF3-attributes-to-qualifiers translation JSON data."""
        record.features = self._run_on_features(record.features)
        return record

    def _run_on_features(self, features: List[SeqFeature]) -> List[SeqFeature]:
        for feature in features:
            if hasattr(feature, "sub_features"):
                feature.sub_features = self._run_on_features(feature.sub_features)
            if hasattr(feature, "qualifiers"):
                feature.qualifiers = self._run_on_qualifiers(feature.qualifiers)
        return features

    def _run_on_qualifiers(self, qualifiers: Dict) -> Dict:
        res = collections.defaultdict(list)
        for name, vals in qualifiers.items():
            if name in self.trans_table:
                # Replace the item name
                new_name = self.trans_table[name]["target"]
                if new_name:
                    prefix = self.trans_table[name].get("prefix", "")
                    res[new_name] += [prefix + v for v in vals]
                else:
                    # Remove the item from qualifiers if "target" is not set, or emtpy
                    pass
            else:
                res[name] += vals
        return dict(res)


class RenameFeatures(object):
    """Translate annotation features according to table given in JSON.

    Meant for renaming GFF3 types to DDBJ annotation feature keys.
    JSON should have the form like following.

    ```
    {
        "five_prime_UTR": {
            "target": "5'UTR"
        },
        "three_prime_UTR": {
            "target": "3'UTR"
        }
    }
    ```

    Empty "target" value means the name key is dropped.
    Unlike qualifiers, "prefix" is not supported here.
    """

    def __init__(self, filepath: str):
        self.path = filepath
        self.d = load_toml_tables(filepath)

    def run(self, record: SeqRecord) -> SeqRecord:
        """Modifies record according to the GFF3-types-to-features translation JSON data."""
        record.features = self._run_on_features(record.features)
        return record

    def _run_on_features(self, features: List[SeqFeature]) -> List[SeqFeature]:
        for feature in features:
            if hasattr(feature, "sub_features"):
                feature.sub_features = self._run_on_features(feature.sub_features)
            name = feature.type
            if name in self.d:
                new_name = self.d[name]["target"]
                feature.type = new_name
        return features


def _join_features(record: SeqRecord, joinables: Optional[Tuple[str]]) -> SeqRecord:
    """
    Join features
    """
    joinables = [] if joinables is None else joinables

    def _join(features: List[SeqFeature]) -> SeqFeature:
        assert len(features) > 1
        locations = []
        qualifiers = collections.OrderedDict()
        sub_features = []
        for f in features:
            locations.append(f.location)
            qualifiers.update(f.qualifiers)
            if hasattr(f, "sub_features") and f.sub_features:
                sub_features.extend(f.sub_features)
        if not sub_features:
            sub_features = None

        compound_loc = CompoundLocation(locations)

        return SeqFeature(
            compound_loc,
            type=features[0].type,
            qualifiers=qualifiers,
            sub_features=sub_features,
        )

    def having_different_codon_starts(features: List[SeqFeature]) -> bool:
        dummy = 99
        codon_starts = {
            x for f in features for x in f.qualifiers.get("codon_start", [dummy])
        }
        return len(codon_starts) > 1 or dummy in codon_starts

    def _join_helper(features: List[SeqFeature]) -> List[SeqFeature]:
        # triples_of_features takes SeqFeature OR tuple (str, str, str) as its key
        triples_or_features = collections.defaultdict(list)
        for f in features:
            if f.type not in joinables:
                triples_or_features[f] = [True]  # dummy values
            else:
                if "product" in f.qualifiers:
                    prod = tuple(f.qualifiers["product"])
                else:
                    prod = None
                triple = (f.type, f.id, prod)
                triples_or_features[triple].append(f)

        res = []
        seen = set()
        for triple_or_f, fs in triples_or_features.items():
            if isinstance(triple_or_f, SeqFeature):
                res.append(triple_or_f)
            elif len(fs) == 1:
                res.extend(fs)
            # Don't join CDS features with different codon_start
            elif triple_or_f[0] == "CDS" and having_different_codon_starts(fs):
                logging.warning(
                    "Skip joining CDSs because of non-unique codon_start: {}".format(
                        triple_or_f
                    )
                )
                res.extend(fs)
            else:
                if triple_or_f not in seen:
                    seen.add(triple_or_f)
                    joined_feature = _join(fs)
                    res.append(joined_feature)

        # join after upper levels
        for f in res:
            if hasattr(f, "sub_features") and f.sub_features:
                f.sub_features = _join_helper(f.sub_features)

        return res

    record.features = _join_helper(record.features)
    return record


def _fix_codon_start_values(rec: SeqRecord) -> None:
    """Convert `codon_start` qualifier value
    from 0-based (in GFF3 'phase' column)
    to   1-based (in INSDC table definition)
    """

    def _fix_feature(feature: SeqFeature) -> None:
        trans = {"0": 1, "1": 2, "2": 3}

        if hasattr(feature, "sub_features"):
            for f in feature.sub_features:
                _fix_feature(f)

        if "codon_start" in feature.qualifiers:
            x = feature.qualifiers["codon_start"][0]
            y = trans[x]
            feature.qualifiers["codon_start"] = [y]

    for f in rec.features:
        _fix_feature(f)


def _add_transl_table(rec: SeqRecord, transl_table: int) -> None:
    """Add transl_table qualifier to all CDS"""

    def _apply(feature: SeqFeature) -> None:
        if hasattr(feature, "sub_features"):
            for f in feature.sub_features:
                _apply(f)
        if feature.type == "CDS":
            feature.qualifiers["transl_table"] = [transl_table]

    for f in rec.features:
        _apply(f)


def _regularize_qualifier_value_letters(rec: SeqRecord) -> None:
    """Fix qualifier values by removing backslash \ or double quote " """

    def _run(features: Iterable[SeqFeature]):
        for f in features:
            for (key, xs) in f.qualifiers.items():
                if isinstance(xs, list):
                    for i, x in enumerate(xs):
                        if isinstance(x, str):
                            if '"' in x:
                                msg = (
                                    "Rename by removing double quotes: ({}, {})".format(
                                        key, x
                                    )
                                )
                                logging.warning(msg)
                                f.qualifiers[key][i] = x.replace('"', "")
                            elif "\\" in x:
                                msg = "Rename by replacing backslash with space: ({}, {})".format(
                                    key, x
                                )
                                logging.warning(msg)
                                f.qualifiers[key][i] = x.replace("\\", " ")
                else:
                    logging.warning(
                        "WTF?? qualifier value type is not a list:  ({}, {}, {})".format(
                            f.type, key, xs
                        )
                    )

    _run(rec.features)


def _remove_duplicates_in_qualifiers(rec: SeqRecord) -> None:
    """Remove duplicate values within a qualifier"""

    def _run(features: SeqFeature) -> None:
        for f in features:
            f.qualifiers = {
                qkey: list(set(qval)) for (qkey, qval) in f.qualifiers.items()
            }
            if hasattr(f, "sub_features"):
                _run(f.sub_features)

    _run(rec.features)


def fix_location_or_start_codon(
    rec: SeqRecord, fasta_record: Dict[str, SeqRecord], transl_table: int
) -> None:
    """Fix location or start_codon according to the DDBJ FAQ
    https://www.ddbj.nig.ac.jp/faq/en/how-to-fix-error-msg-codon-start-e.html

    Args:
        seq may contain CDSs as features to be fixed.
        seq_dict is {SeqID: Record} dict containing sequence info.
        transl_table is the Genetic Code.
    """
    utils.fix_cds(rec, fasta_record, transl_table)


def run(
    path_gff3: Optional[str],
    path_fasta: str,
    path_trans_features: str,
    path_trans_qualifiers: str,
    meta_info: Dict[str, Dict[str, Any]],
    locus_tag_prefix: str,
    transl_table: int,
    joinables: Tuple[str, ...],
) -> List[SeqRecord]:
    """Create a list of `SeqRecord`s and apply various transformations"""

    # get sequence info
    fasta_records = load_fasta_as_seq(path_fasta)
    seq_lengths = dict()
    gaps: Dict[str, List[SeqFeature]] = dict()
    for rec_id, rec in fasta_records.items():
        seq_lengths[rec_id] = len(rec)
        gaps[rec_id] = _get_assembly_gap(rec.seq, meta_info["assembly_gap"])

    # Create record from GFF3 (or dummy if unavailable)
    if path_gff3 is not None:
        records = load_gff3_as_seqrecords(path_gff3)

        # Check if SeqID uses valid characters
        msg1 = "Found invalid letter(s) in the 1st column of the GFF3: {}"
        msg2 = "Run following script to generate corrected GFF3 and FASTA files:\n"
        msg3 = "  $ tools/regularize_seqids.py --gff3={} --fasta={}\n".format(path_gff3, path_fasta)
        for rec in records:
            if utils.is_invalid_as_seqid(rec.id):
                logging.error(msg1.format(rec.id))
                logging.error(msg2)
                logging.error(msg3)
                sys.exit(1)

        # Rename feature and qualifier keys
        f = RenameFeatures(path_trans_features).run
        g = RenameQualifiers(path_trans_qualifiers, locus_tag_prefix).run
        records = [g(f(rec)) for rec in records]

        # Fix codon_start value to 1-based indexing
        for rec in records:
            _fix_codon_start_values(rec)
    else:
        # Create dummy SeqRecords with IDs from FASTA
        records = [SeqRecord("", id=seq_id) for seq_id in fasta_records.keys()]

    # Add "assembly_gap" features
    for rec in records:
        if rec.id in gaps:
            rec.features.extend(gaps[rec.id])

    # Add the transl_table qualifier to CDS feature each
    for rec in records:
        _add_transl_table(rec, transl_table)

    # fix CDS codon_start or location (not implemented yet)
    fix_location_or_start_codon(rec, fasta_records, transl_table)

    # Add "source" feature if unavailable:
    #   [NOTE] GFF3's "region" type corresponds to annotation's "source" feature
    #   [NOTE] User-input metadata may contain "source" feature under the COMMON entry
    #   In either case, a "source" feature is NOT added for each entry.

    if ("source" in meta_info) and ("source" in meta_info["COMMON"]):
        msg = "[COMMON.source] overrides [source] items in metadata."
        logging.warning(msg)
    else:
        for rec in records:
            if not rec.features:
                if rec.features[0].type == "source":
                    msg = 'Ignore [source] in metadata as GFF3 already has "region" line at SeqID = {}'.format(
                        rec.id
                    )
                    logging.warning(msg)
                else:
                    src_length = seq_lengths[rec.id]
                    src_qualifiers = meta_info["source"]
                    src = _get_source(src_length, src_qualifiers)
                    rec.features.insert(0, src)

    # Remove invalid characters from qualifier values
    for rec in records:
        _regularize_qualifier_value_letters(rec)

    # Join features in `joinables` tuple
    if joinables:
        records = [_join_features(rec, joinables) for rec in records]

    # Remove duplicates within a qualifier
    for rec in records:
        _remove_duplicates_in_qualifiers(rec)

    return records
