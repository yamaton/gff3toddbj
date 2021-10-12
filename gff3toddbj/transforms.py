from typing import Any, Dict, Generator, List, Optional, Tuple, Iterable, OrderedDict
import collections
import re
import logging

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import (
    FeatureLocation,
    SeqFeature,
    CompoundLocation,
    BeforePosition,
    AfterPosition,
)

from . import utils
from . import io

Interval = Tuple[int, int]


def _set_assembly_gap(
    record: SeqRecord,
    metadata: OrderedDict[str, OrderedDict[str, Any]],
    faidx: Optional[io.Faidx]=None,
) -> None:
    """
    Set assembly gaps feature to records.
    """
    gap_qualifiers = metadata.get("assembly_gap", collections.OrderedDict())
    seq = str(record.seq) if faidx is None else io.get_seq(faidx, record.id)
    if seq:
        gap = _get_assembly_gap(seq, gap_qualifiers)
        record.features.extend(gap)


def _get_assembly_gap(seq: str, qualifiers: OrderedDict[str, Any]) -> List[SeqFeature]:
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
            length = loc.end.position - loc.start.position + 1
            f.qualifiers["estimated_length"] = [length]
        features.append(f)
    return features


def _get_assembly_gap_locations(seq: str) -> List[Interval]:
    """
    Get assembly_gap locations from seq as a list of intervals.

    [NOTE] This interval format [begin, end) is for Biopython
    such that 0-based, left-inclusive and right-exclusive.

    >>> s = "ATATNNNGATTACANCCC"
    >>> _get_assembly_gap_locations(s)
    [(4, 7), (14, 15)]
    """
    patt = re.compile("N+")
    matches = patt.finditer(seq)

    segments = []
    for m in matches:
        tup = m.span()  # 0-based, left-inclusive, right-exclusive
        segments.append(tup)

    return segments


def _get_source(length: int, source_qualifiers: Dict[str, Any]) -> SeqFeature:
    """Create "source" feature"""
    loc = FeatureLocation(0, length, strand=1)  # 0-based and [0, length)
    return SeqFeature(loc, type="source", qualifiers=source_qualifiers)


class RenameHandler(object):
    """Rename feature and qualifiers according to the structure in TOML.

    By default, "type" value in GFF3 becomes "feature" key in annotation,
    and "attribute" keys and values in GFF3 becomes "qualifier" keys and values
    in annotation. This handler renames/transforms the names and values.

    [TODO] Add description of Case 1, 2, 3, 4.
    """

    _DUMMY_PREFIX = "__tmpname__"

    def __init__(self, filepath: str, locus_tag_prefix: str):
        self.d = io.load_toml_tables(filepath)
        # overwrite with locus_tag_prefix
        self.d["__ANY__"]["locus_tag"]["qualifier_value_prefix"] = locus_tag_prefix

    def run(self, record: SeqRecord) -> SeqRecord:
        """Modifies record according to the GFF3-types-to-features translation JSON data."""
        #
        # Some GFFs contain the exon that is a child of multiple mRNAs,
        # which breaks the tree structure, and my traversal visits
        # some nodes more than once.
        #
        # OTOH, the renaming scheme can contain {A -> B, B -> C}
        # hence the application twice leads unexpected results like A -> C.
        #
        # A dummy constant defined as `RenameHandler._DUMMY_PREFIX = "__tmpname__"`
        # prevents collision of keys and values, hence avoiding multiple applications by
        # {A -> __tmpname__B, B -> __tmpname__C}
        #
        # The dummy prefix is removed in the second scan.
        #
        record.features = self._run_on_features(record.features)
        record.features = self._remove_dummy_prefix_features(record.features)
        return record

    def _run_on_features(self, features: List[SeqFeature]) -> List[SeqFeature]:
        for feature in features:
            if hasattr(feature, "sub_features"):
                feature.sub_features = self._run_on_features(feature.sub_features)

            if "__ANY__" in self.d:
                # (Case 1) Attribute --> Qualifier case
                feature.qualifiers = self._run_on_qualifiers(feature.qualifiers, self.d["__ANY__"])

            type_ = feature.type
            if type_ in self.d:
                below_type = self.d[type_]
                attribute_keys = utils.get_attribute_keys(below_type)
                new_type = below_type.get("feature_key", "")
                if new_type:
                    new_type = RenameHandler._DUMMY_PREFIX + new_type
                    if not attribute_keys:
                        if "qualifier_key" not in below_type.keys():
                            # (Case 2) Type --> Feature key
                            feature.type = new_type
                        else:
                            # (Case 3) Type --> (Feature key, Qualifier key-value)
                            key = below_type["qualifier_key"]
                            value = below_type.get("qualifier_value", None)
                            if value is None:
                                msg = "[{}] /{} ... qualifier_value is missing in the metadata (case 3)".format(type_, key)
                                logging.error(msg)
                                continue
                            feature.type = new_type
                            if key in feature.qualifiers:
                                if not isinstance(feature.qualifiers[key], list):
                                    feature.qualifiers[key] = [feature.qualifiers[key]]
                                if isinstance(value, list):
                                    feature.qualifiers[key].extend(value)
                                else:
                                    feature.qualifiers[key].append(value)
                            else:
                                feature.qualifiers[key] = [value]

                else:
                    # (Case 4) (Type, Attribute) --> Feature key
                    for (key, vals) in feature.qualifiers.items():
                        if key in attribute_keys:
                            below_attribute = below_type[key]
                            new_type = below_attribute.get("feature_key", "")
                            new_type = RenameHandler._DUMMY_PREFIX + new_type
                            if "attribute_value" in below_attribute:
                                attribute_value = below_attribute["attribute_value"]
                                if attribute_value in vals:
                                    feature.type = new_type
                                    feature.qualifiers[key].remove(attribute_value)

        return features


    def _run_on_qualifiers(self, qualifiers: OrderedDict, subtree: OrderedDict) -> OrderedDict:
        """Handle attribute-qualifier renamining
        """
        res = collections.OrderedDict()
        for name, vals in qualifiers.items():
            if name in subtree:
                # Replace the item name
                new_name = subtree[name].get("qualifier_key", "")
                if new_name:
                    if "qualifier_value" in subtree[name]:
                        # overwrite qualifier value with the setting
                        vals = [subtree[name]["qualifier_value"]]
                    else:
                        prefix = subtree[name].get("qualifier_value_prefix", "")
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

    def _remove_dummy_prefix_features(self, features: List[SeqFeature]) -> List[SeqFeature]:
        for f in features:
            if hasattr(f, "sub_features"):
                f.sub_features = self._remove_dummy_prefix_features(f.sub_features)
            if f.type.startswith(RenameHandler._DUMMY_PREFIX):
                size = len(RenameHandler._DUMMY_PREFIX)
                f.type = f.type[size:]
        return features


def _join_features(record: SeqRecord, joinables: Optional[Tuple[str, ...]]) -> SeqRecord:
    """
    Join features of the same hierarchy if the type is found in `joinables`.
    """
    joinables_ = [] if joinables is None else joinables

    def sortkey(f: SeqFeature) -> Tuple[int, int]:
        return (f.location.start.position, f.location.end.position)

    def _join(features: List[SeqFeature]) -> SeqFeature:
        """Join features into a single feature assuming the list already has right members"""
        assert len(features) > 1
        locations = []
        qualifiers = collections.OrderedDict()
        sub_features = []
        features.sort(key=sortkey)
        for f in features:
            locations.append(f.location)
            if hasattr(f, "sub_features") and f.sub_features:
                sub_features.extend(f.sub_features)

        compound_loc = CompoundLocation(locations)
        if compound_loc.strand is None:
            ids = [f.id for f in features]
            logging.error("Something is wrong in joining features:\n    ids = {}".format(ids))

        # set qualifiers of the joined feature from
        # head of the parts. Here "head" is meant by
        #   5' terminal for (+) strand
        #   3' terminal for (-) strand
        idx = 0 if compound_loc.strand > 0 else -1
        qualifiers = features[idx].qualifiers
        type_ = features[idx].type

        result = SeqFeature(
            location=compound_loc,
            type=type_,
            qualifiers=qualifiers,
        )
        if sub_features:
            result.sub_features = sub_features

        return result


    def _runner(features: List[SeqFeature], parent: Optional[str]) -> List[SeqFeature]:
        """Scan features and apply _join """

        # Don't join features if the parent is None or "gene"
        if (parent is None) or (parent == "gene"):
            for f in features:
                if hasattr(f, "sub_features") and f.sub_features:
                    f.sub_features = _runner(f.sub_features, f.type)
            return features

        # Otherwise, group features that are joined.
        # `groups_or_features` is a bag both to-be-joined and not-joined features.
        # Type of `groups_or_features` key is either `SeqFeature` or a group == (type, product).
        # Values for `SeqFeature` keys are dummy, while the type of values for the tuple keys are
        # List[SeqFeature] that are to be joined.
        groups_or_features = collections.defaultdict(list)
        for f in features:
            if f.type not in joinables_:
                groups_or_features[f] = [True]  # dummy values
            else:
                if "product" in f.qualifiers:
                    prod = tuple(f.qualifiers["product"])
                else:
                    prod = None
                group = (f.type, prod)
                groups_or_features[group].append(f)

        res = []
        seen = set()
        for group_or_f, fs in groups_or_features.items():
            if isinstance(group_or_f, SeqFeature):
                res.append(group_or_f)
            elif len(fs) == 1 or parent == "gene":
                res.extend(fs)
            else:
                if group_or_f not in seen:
                    seen.add(group_or_f)
                    joined_feature = _join(fs)
                    res.append(joined_feature)

        # join sublevels after the current level
        for f in res:
            if hasattr(f, "sub_features") and f.sub_features:
                f.sub_features = _runner(f.sub_features, f.type)

        return res


    # --- body of _join_features ---
    # [NOTE] Do not join features at the top-level (So set the parent as None)
    record.features = _runner(record.features, parent=None)

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


def _add_transl_table_to_cds(rec: SeqRecord, transl_table: int) -> None:
    """Add transl_table qualifier to all CDS"""

    def _apply(feature: SeqFeature) -> None:
        if hasattr(feature, "sub_features"):
            for f in feature.sub_features:
                _apply(f)

        # feature-specific transl_table precedes the global setting
        if feature.type == "CDS":
            if "transl_table" in feature.qualifiers:
                genetic_code = int(feature.qualifiers["transl_table"][0])
            else:
                genetic_code = transl_table
            feature.qualifiers["transl_table"] = [genetic_code]

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


def _fix_locations(record: SeqRecord, faidx: Optional[io.Faidx]=None) -> None:
    """
    Fix locations of features in records when start/stop codons are absent.
    See DDBJ FAQ for detail.
    https://www.ddbj.nig.ac.jp/faq/en/how-to-fix-error-msg-codon-start-e.html

    Args:
        `record` may contain CDSs as its features to be fixed.
        `faidx` is for FASTA info. If None, get sequence from `record.seq` instead.
    """
    count_fix_codon_start_loc = []

    def _runner(features: List[SeqFeature], seq: Seq) -> None:
        for f in features:
            if f.type == "CDS":
                cs_list = f.qualifiers.get("codon_start", [1])
                codon_start = int(cs_list[0])
                phase = codon_start - 1  # to 0-based phase index
                if f.location is None:
                    logging.error("f.location is None. Something is wrong: {}".format(f))
                    continue
                if f.location.strand is None:
                    logging.error("f.location.strand is None! Something is wrong: {}".format(f))
                    continue

                # CDS must have the qualifier "transl_table" as done in _add_transl_table_to_cds
                genetic_code = int(f.qualifiers["transl_table"][0])

                if not utils.has_start_codon(seq, f.location, genetic_code, phase):
                    f.location = _fix_absent_start_codon(f.location)

                if not utils.has_stop_codon(seq, f.location, genetic_code, phase):
                    f.location = _fix_absent_stop_codon(f.location)

            if hasattr(f, "sub_features"):
                _runner(f.sub_features, seq)

    def _fix_loc(location: FeatureLocation, is_at_start=True) -> FeatureLocation:
        """Credit: EMBLmyGFF3"""

        if location.strand is None:
            logging.error("location.starnd is unavailable!")
            return location

        # Note that elements in location.parts should be already sorted
        # by _join_features in the order of (location.start, location.end)
        # regardless of the strand (+/-).
        # This is why `idx` is set as follows.
        #
        sign = 1 if is_at_start else -1
        if location.strand * sign > 0:
            # left-end (5' terminal):
            #    start of strand (+)   OR   end of strand (-)
            if len(location.parts) > 1:
                idx = 0  # always 5' terminal
                location.parts[idx] = FeatureLocation(
                    BeforePosition(location.parts[idx].start),
                    location.parts[idx].end,
                    strand=location.parts[idx].strand,
                )
            else:
                location = FeatureLocation(
                    BeforePosition(location.start), location.end, strand=location.strand
                )
        else:
            # right-end (3' terminal):
            #    end of strand (+)   OR   start of strand (-)
            if len(location.parts) > 1:
                idx = -1  # always 3' terminal
                location.parts[idx] = FeatureLocation(
                    location.parts[idx].start,
                    AfterPosition(location.parts[idx].end),
                    strand=location.parts[idx].strand,
                )
            else:
                location = FeatureLocation(
                    location.start, AfterPosition(location.end), strand=location.strand
                )
        return location

    def _fix_absent_start_codon(location: FeatureLocation) -> FeatureLocation:
        return _fix_loc(location, is_at_start=True)

    def _fix_absent_stop_codon(location: FeatureLocation) -> FeatureLocation:
        return _fix_loc(location, is_at_start=False)

    s = record.seq if faidx is None else io.get_seq(faidx, record.id)
    if s:
        seq = Seq(s)
        _runner(record.features, seq)

    if count_fix_codon_start_loc:
        msg = "Changed to /codon_start=1    (SeqID: {},\tcount: {})".format(record.id, len(count_fix_codon_start_loc))
        logging.info(msg)


def _merge_rna_and_exons(rec: SeqRecord) -> None:
    """Set .location of joined exons as the location of RNA,
    then rename such exons into __exon to discard.
    """

    def _helper(features: List[SeqFeature]) -> None:
        for f in features:
            if "RNA" in f.type:  # f.type can be mRNA, ncRNA, misc_RNA, tRNA, rRNA, etc.
                joined_exons = [subf for subf in f.sub_features if subf.type == "exon"]
                if len(joined_exons) > 1:
                    logging.warning("Something is wrong with joining exons: {}".format(f))
                    continue
                elif not joined_exons:
                    # when mRNA does not have exons as its children, just leave it as is
                    continue

                joined_exon = joined_exons[0]
                # set the location of joined exons as mRNA's location
                f.location = joined_exon.location
                # rename exon to __exon to supress output
                joined_exon.type = utils.DUMMY_ORIGINALLY_EXON

            if hasattr(f, "sub_features"):
                _helper(f.sub_features)

    # just call the _helper
    _helper(rec.features)


def _handle_source_rec(
    rec: SeqRecord,
    metadata: OrderedDict[str, OrderedDict[str, Any]],
    id_to_seqlen: Optional[OrderedDict[str, int]]=None,
) -> bool:
    """
    Add "source" features in certain cases:
      [source] in metadata will inserts "source" to each entry
      UNLESS either of the following applies.
      [NOTE] GFF3's "region" type corresponds to annotation's "source" feature.
      [NOTE] User-input metadata may contain "[COMMON.source]" items.
    """
    is_inserting = False
    if ("source" in metadata) and ("source" in metadata.get("COMMON", ())):
        msg = "Both [COMMON.source] and [source] exist in the metadata file."
        logging.warning(msg)

    if "source" in metadata:
        if rec.features and rec.features[0].type != "source":
            is_inserting = True
            src_length = id_to_seqlen[rec.id] if id_to_seqlen else len(rec.seq)
            src_qualifiers = metadata["source"]
            src = _get_source(src_length, src_qualifiers)
            rec.features.insert(0, src)
        else:
            source = rec.features[0]
            for qual_key, qual_val in metadata["source"].items():
                if qual_key in source.qualifiers:
                    source.qualifiers[qual_key].extend(qual_val)
                else:
                    source.qualifiers[qual_key] = qual_val

    return is_inserting


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
    """Sort features"""

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
    path_rename_scheme: str,
    metadata: OrderedDict[str, OrderedDict[str, Any]],
    locus_tag_prefix: str,
    transl_table: int,
    joinables: Tuple[str, ...],
) -> Generator[SeqRecord, None, None]:
    """Create a list of `SeqRecord`s and apply various transformations"""

    # Load FASTA as pysam class
    faidx = io.load_fasta_as_faidx(path_fasta)

    # Get SeqIDs and Sequence lengths
    id_to_seqlen = io.get_seqlens(faidx)
    fasta_ids = list(id_to_seqlen.keys())

    # Create record from GFF3 (or dummy if unavailable)
    if path_gff3 is not None:
        records = io.load_gff3_as_seqrecords(path_gff3)

    else:
        # Create dummy SeqRecords with IDs from FASTA
        records = [SeqRecord("", id=seq_id) for seq_id in fasta_ids]

    # Rename feature keys and/or qualifier keys/values
    f = RenameHandler(path_rename_scheme, locus_tag_prefix).run

    for rec in records:
        rec = f(rec)

        # Convert codon_start value to 1-based indexing
        _convert_codon_start_to_1_based(rec)

        # Add "assembly_gap" features
        _set_assembly_gap(rec, metadata, faidx=faidx)

        # Add the transl_table qualifier to CDS feature each
        _add_transl_table_to_cds(rec, transl_table)

        # Insert "source" features from metadata if necessary
        _handle_source_rec(rec, metadata, id_to_seqlen)

        # Check characters in qualifier values
        _regularize_qualifier_value_letters(rec)

        # Join features if the key is in `joinables`
        if joinables:
            rec = _join_features(rec, joinables)

        # Merge exons with their parent mRNA
        _merge_rna_and_exons(rec)

        # Check start and stop codons in CDSs
        _fix_locations(rec, faidx=faidx)

        # Assign single value to /product and put the rest to /inference
        _assign_single_product(rec)

        # Remove duplicates within a qualifier
        _remove_duplicates_in_qualifiers(rec)

        # Sort features
        _sort_features(rec.features)

        # Check if SeqIDs have valid characters
        msg = (
            "\n\n"
            "Found invalid letter(s) in the 1st column of the GFF3: {}\n"
            "Consider running this script to correct entry names in the DDBJ annotation:\n\n"
            "   $ normalize-entry-names <output.ann>\n"
        )
        if utils.is_invalid_as_seqid(rec.id):
            logging.warning(msg.format(rec.id))

        yield rec

    faidx.close()


def run_with_genbank(
    path_genbank: str,
    metadata: OrderedDict[str, OrderedDict[str, Any]],
    transl_table: int,
) -> Generator[SeqRecord, None, None]:
    """Create a list of `SeqRecord`s and apply various transformations"""

    # Load Genbank as a generator of seqrecord
    records = io.load_genbank_as_seqrecords(path_genbank)

    for rec in records:

        # Add "assembly_gap" features
        _set_assembly_gap(rec, metadata)

        # Add the transl_table qualifier to CDS feature each
        _add_transl_table_to_cds(rec, transl_table)

        # Remove duplicates within a qualifier
        _remove_duplicates_in_qualifiers(rec)

        # Check if SeqIDs have valid characters
        msg = (
            "\n\n"
            "Found invalid letter(s) in the 1st column of the GFF3: {}\n"
            "Consider running this script to correct entry names in the DDBJ annotation:\n\n"
            "   $ normalize-entry-names <output.ann>\n"
        )
        if utils.is_invalid_as_seqid(rec.id):
            logging.warning(msg.format(rec.id))

        yield rec
