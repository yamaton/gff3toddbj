import pprint
from typing import Any, Dict, Generator, List, Optional, Tuple, Iterable
import collections
import toml
import re
import pathlib
import gzip
import logging
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation
from BCBio import GFF

Interval = Tuple[int, int]


def load_gff3_as_seqrecords(filepath) -> List[SeqRecord]:
    """
    Load GFF3 as iterable of SeqRecord
    """
    p = pathlib.Path(filepath)
    if p.suffix == ".gz":
        with gzip.open(filepath, "rt") as f:
            recs = list(GFF.parse(f))
    else:
        recs = list(GFF.parse(filepath))
    return recs


def load_fasta_as_seq(filepath) -> List[SeqRecord]:
    """
    Load FASTA file as Seq
    """
    p = pathlib.Path(filepath)
    if p.suffix == ".gz":
        with gzip.open(filepath, "rt") as f:
            recs = list(Bio.SeqIO.parse(f, "fasta"))
    else:
        recs = list(Bio.SeqIO.parse(filepath, "fasta"))
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

    [NOTE] A location is of format (begin, end)
    with 1-based indexing, and inclusive in both sides.

    >>> s = Seq("atatnnngattacanccc")
    >>> get_assembly_gap(s)
    [(5, 7), (15, 15)]
    """
    locs = [
        FeatureLocation(start, end, strand=1)
        for (start, end) in _get_assembly_gap_locations(seq)
    ]

    features = []
    for loc in locs:
        f = SeqFeature(loc, type="assembly_gap", qualifiers=qualifiers)
        if f.qualifiers.get("estimated_length", None) == "<COMPUTE>":
            length = loc.end - loc.start + 1
            f.qualifiers["estimated_length"] = length
        features.append(f)
    return features


def _get_assembly_gap_locations(seq: Seq) -> List[Interval]:
    """
    Get assembly_gap locations from seq.

    [NOTE] A location is of format (begin, end)
    with 1-based indexing, and inclusive in both sides.

    >>> s = Seq("atatnnngattacanccc")
    >>> _get_assembly_gap_locations(s)
    [(5, 7), (15, 15)]
    """
    s = str(seq)
    patt = re.compile("n+")
    matches = patt.finditer(s)

    segments = []
    for m in matches:
        a, b = m.span()  # 0-based, left-inclusive, right-exclusive
        tup = (a + 1, b)  # 1-based, both-inclusive
        segments.append(tup)

    return segments


def _get_source(length: int, source_qualifiers: Dict[str, Any]) -> SeqFeature:
    """Create "source" feature"""
    loc = FeatureLocation(1, length, strand=1)
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
        codon_starts = {f.qualifiers.get("codon_start", dummy) for f in features}
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
                logging.warn("Skip joining CDSs because of non-unique codon_start: {}".format(triple_or_f))
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
        trans = {"0": "1", "1": "2", "2": "3"}

        if hasattr(feature, "sub_features"):
            for f in feature.sub_features:
                _fix_feature(f)

        if "codon_start" in feature.qualifiers:
            x = feature.qualifiers["codon_start"][0]
            y = trans[x]
            feature.qualifiers["codon_start"] = [y]

    for f in rec.features:
        _fix_feature(f)


def run(
    path_gff3: Optional[str],
    path_fasta: str,
    path_trans_features: str,
    path_trans_qualifiers: str,
    meta_info: Dict[str, Dict[str, Any]],
    locus_tag_prefix: str,
    joinables: Tuple[str, ...],
) -> List[SeqRecord]:
    """Create a list of `SeqRecord`s and apply various transformations
    """
    # get sequence info
    fasta_records = load_fasta_as_seq(path_fasta)
    seq_lengths = dict()
    gaps = dict()
    seq_ids = []
    for rec in fasta_records:
        seq_lengths[rec.id] = len(rec)
        gaps[rec.id] = _get_assembly_gap(rec.seq, meta_info["assembly_gap"])
        seq_ids.append(rec.id)

    # Create record from GFF3 (or dummy if unavailable)
    if path_gff3 is not None:
        records = load_gff3_as_seqrecords(path_gff3)
        f = RenameFeatures(path_trans_features).run
        g = RenameQualifiers(path_trans_qualifiers, locus_tag_prefix).run

        # Rename feature and qualifier keys
        records = [g(f(rec)) for rec in records]

        # Fix codon_start value to 1-based indexing
        for rec in records:
            _fix_codon_start_values(rec)
    else:
        # Create dummy SeqRecords with IDs from FASTA
        records = [SeqRecord("", id=seq_id) for seq_id in seq_ids]

    # Add "assembly_gap" feature
    for rec in records:
        if rec.id in gaps:
            rec.features.extend(gaps[rec.id])

    # Add "source" feature if unavailable
    # [NOTE] GFF3's "region" type corresponds to annotation's "source" feature
    for rec in records:
        if (not rec.features) or (rec.features[0].type != "region"):
            src_length = seq_lengths[rec.id]
            src_qualifiers = meta_info["source"]
            src = _get_source(src_length, src_qualifiers)
            rec.features.insert(0, src)
        else:
            logging.warn('Ignore [source] in metadata as GFF3 already has "region" line at SeqID = {}'.format(rec.id))

    # Join features in `joinables` tuple
    if joinables:
        records = [_join_features(rec, joinables) for rec in records]

    return records
