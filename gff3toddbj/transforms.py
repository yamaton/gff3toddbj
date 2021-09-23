from typing import Any, Dict, List, Optional, Tuple, Iterable, OrderedDict
import collections
import re
import logging
import sqlite3

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation

from . import utils
from . import io

Interval = Tuple[int, int]


def _set_assembly_gap(records: List[SeqRecord], cur: sqlite3.Cursor, metadata: OrderedDict[str, OrderedDict[str, Any]]) -> None:
    """
    Set assembly gaps feature to records.
    """
    gap_qualifiers = metadata.get("assembly_gap", collections.OrderedDict())
    for rec in records:
        seq = io.get_seq(cur, rec.id)
        if seq:
            gap = _get_assembly_gap(seq, gap_qualifiers)
            rec.features.extend(gap)


def _get_assembly_gap(seq: Seq, qualifiers: OrderedDict[str, Any]) -> List[SeqFeature]:
    """
    Get assembly_gap features from seq.

    Args:
        seq (Seq): sequence to check assembly gaps
        qualifiers (dict): qualifier values
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

    >>> s = Seq("ATATNNNGATTACANCCC")
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
        self.trans_table = io.load_toml_tables(filepath)

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

    def _run_on_qualifiers(self, qualifiers: OrderedDict) -> OrderedDict:
        res = collections.OrderedDict()
        for name, vals in qualifiers.items():
            if name in self.trans_table:
                # Replace the item name
                new_name = self.trans_table[name].get("target", "")
                if new_name:
                    prefix = self.trans_table[name].get("prefix", "")
                    if prefix:
                        vals = [prefix + v for v in vals]

                    if new_name not in res:
                        res[new_name] = []
                    res[new_name] += vals
                else:
                    # Remove the item from qualifiers if "target" is not set, or emtpy
                    pass
            else:
                if name not in res:
                    res[name] = []
                res[name] += vals
        return res


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
        self.d = io.load_toml_tables(filepath)

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


def _join_features(record: SeqRecord, joinables: Optional[Tuple[str, ...]]) -> SeqRecord:
    """
    Join features of the same hierarchy if the type is found in `joinables`.
    """
    joinables_ = [] if joinables is None else joinables

    def _join(features: List[SeqFeature]) -> SeqFeature:
        """Join features into a single feature assuming the list already has right members
        """
        assert len(features) > 1
        locations = []
        qualifiers = collections.OrderedDict()
        sub_features = []
        for f in features:
            locations.append(f.location)
            if hasattr(f, "sub_features") and f.sub_features:
                sub_features.extend(f.sub_features)

        # this is how to set qualifiers of the joined feature
        # and this matters in setting /codon_start right.
        qualifiers = features[0].qualifiers

        if not sub_features:
            sub_features = None

        compound_loc = CompoundLocation(locations)

        return SeqFeature(
            compound_loc,
            type=features[0].type,
            qualifiers=qualifiers,
            sub_features=sub_features,
        )

    def _join_helper(features: List[SeqFeature]) -> List[SeqFeature]:
        """
        """
        # `triples_or_features` has either `SeqFeature` or (type, id, product) as its key.
        # An item with a `SeqFeature` key has dummy while an item with a tuple key has
        # List[SeqFeature] as its value.
        triples_or_features = collections.defaultdict(list)
        for f in features:
            if f.type not in joinables_:
                triples_or_features[f] = [True]  # dummy values
            else:
                if "product" in f.qualifiers:
                    prod = tuple(f.qualifiers["product"])
                else:
                    prod = None
                triple = (f.type, prod)
                triples_or_features[triple].append(f)

        res = []
        seen = set()
        for triple_or_f, fs in triples_or_features.items():
            if isinstance(triple_or_f, SeqFeature):
                res.append(triple_or_f)
            elif len(fs) == 1:
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


def _convert_codon_start_to_1_based(rec: SeqRecord) -> None:
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
    """Fix qualifier values by removing backslash \\ or double quote \" """

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
    """Remove duplicates in each qualifier's values"""

    def _unique(xs: List) -> List:
        seen = set()
        result = []
        for x in xs:
            if x not in seen:
                result.append(x)
                seen.add(x)
        return result

    def _run(features: List[SeqFeature]) -> None:
        for f in features:
            f.qualifiers = {
                qkey: _unique(qval) for (qkey, qval) in f.qualifiers.items()
            }
            if hasattr(f, "sub_features"):
                _run(f.sub_features)

    _run(rec.features)


def _check_start_codons(
    rec: SeqRecord, fasta_record: Dict[str, SeqRecord], transl_table: int
) -> None:
    """Fix location or start_codon according to the DDBJ FAQ
    https://www.ddbj.nig.ac.jp/faq/en/how-to-fix-error-msg-codon-start-e.html

    Args:
        seq may contain CDSs as features to be fixed.
        seq_dict is {SeqID: Record} dict containing sequence info.
        transl_table is the Genetic Code.
    """
    utils.check_cds(rec, fasta_record, transl_table)


def _merge_mrna_qualifiers(rec: SeqRecord) -> None:
    """Set qualifiers in __mRNA (orignally mRNA type)
    as the qualifiers of the merged mRNA feature.
    """
    def _helper(features: List[SeqFeature]) -> None:
        for f in features:
            if f.type == utils.DUMMY_ORIGINAL_MRNA:
                mrnas = [subf for subf in f.sub_features if subf.type == "mRNA"]
                if len(mrnas) != 1:
                    logging.warning("Something is wrong with mRNA and exons: {}".format(f))
                    return
                for subf in f.sub_features:
                    if subf.type == "mRNA":
                        for key, val in f.qualifiers.items():
                            subf.qualifiers[key] = val

            if hasattr(f, "sub_features"):
                _helper(f.sub_features)

    # just call the _helper
    _helper(rec.features)


def _assign_single_product(rec: SeqRecord) -> None:
    """Take the first item in /product values as the value of /product
    and put the rest as /note values
    """
    DEFAULT = ["hypothetical protein"]

    def _helper(features: List[SeqFeature]) -> None:
        for f in features:
            if f.type == "CDS":
                if ("product" not in f.qualifiers) or (not f.qualifiers["product"]):
                    f.qualifiers["product"] = DEFAULT
                elif len(f.qualifiers["product"]) > 1:
                    values = f.qualifiers["product"]
                    head = values[0]
                    f.qualifiers["product"] = [head]

                    prefix = "product:"
                    rest = [prefix + x for x in values[1:]]
                    if "note" in f.qualifiers:
                        f.qualifiers["note"].extend(rest)
                    else:
                        f.qualifiers["note"] = rest

            if hasattr(f, "sub_features"):
                _helper(f.sub_features)

    # just call the _helper
    _helper(rec.features)


def _sort_features(features: List[SeqFeature]) -> None:
    """Sort features
    """
    def type_priority(type_name: str) -> float:
        d = {
            "source": -1,
            "gene": 0,
            "mRNA": 1,
            "rRNA": 1,
            "ncRNA": 1,
            "tmRNA": 1,
            "tRNA": 1,
            "misc_RNA": 1,
            "V_region": 1,
            "C_region": 1.5,
            "D_segment": 1.5,
            "J_segment": 1.5,
            "V_segment": 1.5,
            "5'UTR": 2,
            "3'UTR": 2,
            "CDS": 2,
            "mat_peptide": 2.5,
            "exon": 3,
            "intron": 3,
            "assembly_gap": 8,
        }
        return d[type_name] if (type_name in d) else 10

    def keyfunc(f: SeqFeature) -> Tuple[int, float, int, str]:
        return (f.location.start, type_priority(f.type), f.location.end, f.id)

    features.sort(key=keyfunc)
    for f in features:
        if hasattr(f, "sub_features"):
            _sort_features(f.sub_features)


def run(
    path_gff3: Optional[str],
    path_fasta: str,
    path_trans_features: str,
    path_trans_qualifiers: str,
    metadata: OrderedDict[str, OrderedDict[str, Any]],
    locus_tag_prefix: str,
    transl_table: int,
    joinables: Tuple[str, ...],
) -> List[SeqRecord]:
    """Create a list of `SeqRecord`s and apply various transformations"""

    # Load FASTA to database
    con = io.load_fasta_as_database(path_fasta)
    cur = con.cursor()

    # Get SeqIDs and Sequence lengths
    id_to_seqlen = io.get_seqlens(cur)
    fasta_ids = list(id_to_seqlen.keys())

    # Create record from GFF3 (or dummy if unavailable)
    if path_gff3 is not None:
        records = io.load_gff3_as_seqrecords(path_gff3)

        # Rename feature and qualifier keys
        f = RenameFeatures(path_trans_features).run
        g = RenameQualifiers(path_trans_qualifiers, locus_tag_prefix).run
        records = [g(f(rec)) for rec in records]
    else:
        # Create dummy SeqRecords with IDs from FASTA
        records = [SeqRecord("", id=seq_id) for seq_id in fasta_ids]

    # Convert codon_start value to 1-based indexing
    for rec in records:
        _convert_codon_start_to_1_based(rec)

    # Add "assembly_gap" features
    _set_assembly_gap(records, cur, metadata)

    # Add the transl_table qualifier to CDS feature each
    for rec in records:
        _add_transl_table(rec, transl_table)

    # Add "source" features in certain cases:
    #   [source] in metadata will inserts "source" to each entry
    #   UNLESS either of the following applies.
    #   [NOTE] GFF3's "region" type corresponds to annotation's "source" feature.
    #   [NOTE] User-input metadata may contain "[COMMON.source]" items.
    if ("source" in metadata) and ("source" in metadata["COMMON"]):
        msg = "[COMMON.source] overrides [source] items in metadata."
        logging.warning(msg)
    elif "source" in metadata:
        for rec in records:
            if rec.features:
                if rec.features[0].type == "source":
                    msg = 'Skip [source] in metadata as GFF3 already has "region" line at SeqID = {}'.format(
                        rec.id
                    )
                    logging.warning(msg)
                else:
                    src_length = id_to_seqlen[rec.id]
                    src_qualifiers = metadata["source"]
                    src = _get_source(src_length, src_qualifiers)
                    rec.features.insert(0, src)

    # Regularize characters in qualifier values
    for rec in records:
        _regularize_qualifier_value_letters(rec)

    # Join features in `joinables` tuple
    if joinables:
        records = [_join_features(rec, joinables) for rec in records]

    # Merge __mRNA with mRNAs (originally mRNA and exons)
    for rec in records:
        _merge_mrna_qualifiers(rec)

    # # check start codons in CDSs
    # for rec in records:
    #     _check_start_codons(rec, fasta_records, transl_table)

    # assign single value to /product and put the rest to /inference
    for rec in records:
        _assign_single_product(rec)

    # Remove duplicates within a qualifier
    for rec in records:
        _remove_duplicates_in_qualifiers(rec)

    # Sort features
    for rec in records:
        _sort_features(rec.features)

    # Check if SeqIDs have valid characters
    msg = (
        "\n\n"
        "Found invalid letter(s) in the 1st column of the GFF3: {}\n"
        "Run following script to generate corrected GFF3 and FASTA files:\n\n"
    ) + "  $ rename-ids --gff3={} --fasta={}\n".format(path_gff3, path_fasta)
    for rec in records:
        if utils.is_invalid_as_seqid(rec.id):
            logging.warning(msg.format(rec.id))

    cur.close()
    io.close_and_remove_database(con)
    return records
