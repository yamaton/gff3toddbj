from typing import Any, Dict, FrozenSet, Generator, List, Optional, Tuple, Iterable, OrderedDict, DefaultDict
import collections
import re
import logging

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import (
    SimpleLocation,
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
        SimpleLocation(start, end, strand=1)
        for (start, end) in _get_assembly_gap_locations(seq)
    ]

    features = []
    for loc in locs:
        f = SeqFeature(loc, type="assembly_gap", qualifiers=qualifiers)
        length = loc.end - loc.start
        if f.qualifiers.get("estimated_length", None) == ["<COMPUTE>"]:
            f.qualifiers["estimated_length"] = [length]
        if length >= 10:
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
    loc = SimpleLocation(0, length, strand=1)  # 0-based and [0, length)
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
        self.d = utils.to_lowercase_keys(io.load_toml_tables(filepath))

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
                feature.qualifiers = self._run_on_qualifiers(feature.qualifiers, self.d["__ANY__"])

            type_ = feature.type
            type_lower = type_.lower()
            if type_lower in self.d:
                below_type = self.d[type_lower]
                attribute_keys = utils.get_attribute_keys(below_type)
                if not attribute_keys:
                    # ------- level-1: type name matching only -------
                    if "feature_key" not in below_type:
                        logging.error("Renaming config requires feature_key at [{}]".format(type))
                        continue

                    new_type = RenameHandler._DUMMY_PREFIX + below_type["feature_key"]
                    if "qualifier_key" not in below_type.keys():
                        # Rename Type / Feature
                        feature.type = new_type
                    else:
                        # (Case 3) Type --> (Feature key, Qualifier key-value)
                        qual_key = below_type["qualifier_key"]
                        value = below_type.get("qualifier_value", None)
                        if value is None:
                            msg = "[{}] /{} ... qualifier_value is missing in the metadata (case 3)".format(type_, qual_key)
                            logging.error(msg)
                            continue
                        feature.type = new_type
                        if qual_key in feature.qualifiers:
                            if not isinstance(feature.qualifiers[qual_key], list):
                                feature.qualifiers[qual_key] = [feature.qualifiers[qual_key]]
                            if isinstance(value, list):
                                feature.qualifiers[qual_key].extend(value)
                            else:
                                feature.qualifiers[qual_key].append(value)
                        else:
                            feature.qualifiers[qual_key] = [value]
                else:
                    # Going deeper the renaming config like [CDS.*]
                    for (qual_key, qual_vals) in feature.qualifiers.items():
                        qual_vals_lowered = [s.lower() for s in qual_vals]
                        qual_key_lower = qual_key.lower()
                        if qual_key_lower in attribute_keys:
                            below_attribute = below_type[qual_key_lower]
                            attribute_values = utils.get_attribute_values(below_attribute)
                            new_type = below_attribute.get("feature_key", "")
                            if new_type and (not attribute_values):
                                # ------- level-2: (feature_key, qualifier_key) matching
                                feature.type = RenameHandler._DUMMY_PREFIX + new_type
                                if "qualifier_key" in below_attribute:
                                    new_qual_key = below_attribute["qualifier_key"]
                                    if "qualifier_value" in below_attribute:
                                        # set qualifier value
                                        feature.qualifiers[new_qual_key] = below_attribute["qualifier_value"]
                                    else:
                                        # rename qualifier
                                        feature.qualifiers[new_qual_key] = feature.qualifiers[qual_key]
                                        del feature.qualifiers[qual_key]

                            elif set(qual_vals_lowered) & set(attribute_values):
                                # ------- level-3: (feature key, qualifier key, qualifier value) matching
                                intersect = set(qual_vals_lowered) & set(attribute_values)
                                for attribute_value in intersect:
                                    bottom = below_attribute[attribute_value]
                                    if "feature_key" in bottom:
                                        feature.type = bottom["feature_key"]
                                    elif all(x in bottom for x in ("qualifier_key", "qualifier_value")):
                                        new_key = bottom["qualifier_key"]
                                        new_val = bottom["qualifier_value"]
                                        feature.qualifiers[new_key] = [new_val]
                                        matched_val = next(v for v in qual_vals if v.lower() == attribute_value)
                                        feature.qualifiers[qual_key].remove(matched_val)
                                        if not feature.qualifiers[qual_key]:
                                            del feature.qualifiers[qual_key]

        return features


    def _run_on_qualifiers(self, qualifiers: OrderedDict, subtree: OrderedDict) -> OrderedDict:
        """Handle attribute-qualifier renamining
        """
        res = collections.OrderedDict()
        for name, vals in qualifiers.items():
            attr_key = name.lower()
            if attr_key in subtree:
                attr_vals = utils.get_attribute_values(subtree[attr_key])
                vals_lower = [v.lower() for v in vals]
                if not attr_vals:
                    # ---- level-2: matching to the 2nd depth ----
                    new_name = subtree[attr_key].get("qualifier_key", "")
                    if new_name:
                        name = new_name
                        if "qualifier_value" in subtree[attr_key]:
                            # overwrite qualifier value with the setting
                            vals = [subtree[attr_key]["qualifier_value"]]
                        else:
                            prefix = subtree[attr_key].get("qualifier_value_prefix", "")
                            vals = [prefix + v for v in vals]
                    else:
                        # Remove the item from qualifiers if "qualifier_key" is not set, or emtpy
                        continue

                elif set(attr_vals) & set(vals_lower):
                    intersect = set(attr_vals) & set(vals_lower)
                    for attr_val in intersect:
                        # ---- level-3: matching to the 3rd depth ----
                        subsubtree = subtree[attr_key]
                        d = subsubtree[attr_val]
                        if not all(x in d for x in ("qualifier_key", "qualifier_value")):
                            logging.error("Bad format in the renaming config: {} - {}".format(attr_key, attr_val))
                        # following will drop the rest of values under attr_key. ok???
                        name = d["qualifier_key"]
                        vals = [d["qualifier_value"]]
                        # add to res
                        if name not in res:
                            res[name] = []
                        res[name] += vals
                        continue

            # add to res
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
        return (int(f.location.start), int(f.location.end))

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

        if parent is None:
            for f in features:
                if hasattr(f, "sub_features") and f.sub_features:
                    f.sub_features = _runner(f.sub_features, f.type)
            return features

        # Use a dictionary where the key is a Hashable identifier.
        # For non-joinable features, use the object's memory ID.
        # For joinable features, use the (type, product) tuple.
        groups: DefaultDict[Any, List[SeqFeature]] = collections.defaultdict(list)

        # We need to maintain order, so keep track of the keys as they appear
        order = []

        for f in features:
            if f.type not in joinables_:
                key = id(f)  # Use memory address as a unique, hashable key
            else:
                prod = tuple(f.qualifiers.get("product", []))
                key = (f.type, prod)

            if key not in groups:
                order.append(key)
            groups[key].append(f)

        res = []
        for key in order:
            fs = groups[key]

            # If the key is an integer, it's a non-joinable feature (the id)
            if isinstance(key, int):
                res.extend(fs)
            elif len(fs) == 1:
                res.extend(fs)
            elif parent == "gene":
                xs = [f for f in fs if "ribosomal_slippage" in f.qualifiers]
                rest = [f for f in fs if "ribosomal_slippage" not in f.qualifiers]
                if xs:
                    res.append(_join(xs))
                res.extend(rest)
            else:
                res.append(_join(fs))

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
            # skip if CDS has /pseudo or /pseudogene=...
            if f.type == "CDS" and all(x not in f.qualifiers for x in ("pseudo", "pseudogene")):
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

    def _fix_loc(location: SimpleLocation, is_at_start=True) -> SimpleLocation:
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
                location.parts[idx] = SimpleLocation(
                    BeforePosition(location.parts[idx].start),
                    location.parts[idx].end,
                    strand=location.parts[idx].strand,
                )
            else:
                location = SimpleLocation(
                    BeforePosition(location.start), location.end, strand=location.strand
                )
        else:
            # right-end (3' terminal):
            #    end of strand (+)   OR   start of strand (-)
            if len(location.parts) > 1:
                idx = -1  # always 3' terminal
                location.parts[idx] = SimpleLocation(
                    location.parts[idx].start,
                    AfterPosition(location.parts[idx].end),
                    strand=location.parts[idx].strand,
                )
            else:
                location = SimpleLocation(
                    location.start, AfterPosition(location.end), strand=location.strand
                )
        return location

    def _fix_absent_start_codon(location: SimpleLocation) -> SimpleLocation:
        return _fix_loc(location, is_at_start=True)

    def _fix_absent_stop_codon(location: SimpleLocation) -> SimpleLocation:
        return _fix_loc(location, is_at_start=False)

    s = record.seq if faidx is None else io.get_seq(faidx, record.id)
    if s:
        seq = Seq(s)
        _runner(record.features, seq)

    if count_fix_codon_start_loc:
        msg = "Changed to /codon_start=1    (SeqID: {},\tcount: {})".format(record.id, len(count_fix_codon_start_loc))
        logging.info(msg)


def _merge_exons_with_parent(rec: SeqRecord) -> None:
    """Set .location of joined exons as the location of RNA, or *_segment, or *_region
    then rename such exons into __exon to be discarded.
    """

    def _helper(features: List[SeqFeature]) -> None:
        for f in features:
            # f.type can be mRNA, ncRNA, misc_RNA, tRNA, rRNA, V_segment, C_region, etc.
            if any(kwd in f.type for kwd in ["RNA", "_segment", "_region"]):
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
        if (not rec.features) or (rec.features and rec.features[0].type != "source"):
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


def _handle_topology(rec: SeqRecord) -> None:
    """Add TOPOLOGY feature with /circular qualifier if "source" has
    Is_circular=true value. Also handle origin-spanning features.

    https://https.ncbi.nlm.nih.gov/datasets/docs/v1/reference-docs/file-formats/about-ncbi-gff3/#origin-spanning-features
    """
    def _scan_origin_spannig(features: List[SeqFeature], record_end_pos: int) -> None:
        for f in features:
            if f.location is None:
                continue
            try:
                start = int(f.location.start)
                end = int(f.location.end)
            except (AttributeError, TypeError, ValueError):
                continue
            if end > record_end_pos:
                before_origin = SimpleLocation(start, record_end_pos, strand=f.location.strand)
                end_true = end - record_end_pos
                after_origin = SimpleLocation(0, end_true, strand=f.location.strand)
                f.location = CompoundLocation([before_origin, after_origin])

                # scan only the children of matched features
                if hasattr(f, "sub_features"):
                    _scan_origin_spannig(f.sub_features, record_end_pos)


    # Assume "source" always comes to the top if exists
    first_feature = rec.features[0]
    if first_feature.type == "source":
        k = next((k for k in first_feature.qualifiers if k.lower() == "is_circular"), None)
        if k is not None:
            del first_feature.qualifiers[k]
            new_feature = SeqFeature(None, type="TOPOLOGY", qualifiers={"circular": [""]})
            rec.features.insert(1, new_feature)

            record_end_pos = int(first_feature.location.end)
            _scan_origin_spannig(rec.features, record_end_pos)


def _copy_qualifiers_to_children(rec: SeqRecord, qualifier_key_data: Dict[str, FrozenSet[str]]) -> None:
    """Copy qualifiers in data.values() to children if a feature in data.keys() has them.
    """
    def _assign_qual(features: List[SeqFeature], qkey: str, qvals: List[str]):
        for f in features:
            if qkey in f.qualifiers:
                f.qualifiers[qkey].extend(qvals)
            else:
                f.qualifiers[qkey] = qvals
            if hasattr(f, "sub_features"):
                _assign_qual(f.sub_features, qkey, qvals)

    def _runner(features: List[SeqFeature]):
        for f in features:
            if f.type in qualifier_key_data:
                qualifier_keys = qualifier_key_data[f.type]
                for qkey in qualifier_keys:
                    if qkey in f.qualifiers and hasattr(f, "sub_features"):
                        qvals = f.qualifiers[qkey]
                        _assign_qual(f.sub_features, qkey, qvals)

            if hasattr(f, "sub_features"):
                _runner(f.sub_features)

    _runner(rec.features)


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


def _make_gene_have_single_value(rec: SeqRecord) -> None:
    """Take the first item in /gene values as the value of /gene
    and put the rest as /gene_synonym values

    Ref: https://www.ddbj.nig.ac.jp/ddbj/qualifiers-e.html#gene
    """
    def _helper(features: List[SeqFeature]) -> None:
        for f in features:
            if ("gene" in f.qualifiers) and len(f.qualifiers["gene"]) > 1:
                values = f.qualifiers["gene"]
                head = values[0]
                f.qualifiers["gene"] = [head]

                rest = values[1:]
                if "gene_synonym" not in f.qualifiers:
                    f.qualifiers["gene_synonym"] = []
                f.qualifiers["gene_synonym"].extend(rest)

            if hasattr(f, "sub_features"):
                _helper(f.sub_features)

    # just call the _helper
    _helper(rec.features)


def _sort_features(features: List[SeqFeature]) -> None:
    """Sort features"""

    def type_priority(type_name: str) -> float:
        d = {
            "source": -1,
            "TOPOLOGY": -0.9,
            "gene": 0,
            "misc_feature": 0.5,
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
        if f.location:
            res = (f.location.start, type_priority(f.type), f.location.end, f.id)
        else:
            res = (0, type_priority(f.type), 100000000, f.id)
        return res

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
    data_qualifiers_to_children: Dict[str, FrozenSet[str]],
) -> Generator[SeqRecord, None, None]:
    """Create a list of `SeqRecord`s and apply various transformations"""

    # Load FASTA as pysam class
    faidx = io.load_fasta_as_faidx(path_fasta)

    # Get SeqIDs and Sequence lengths
    id_to_seqlen = io.get_seqlens(faidx)
    fasta_ids = list(id_to_seqlen.keys())

    # Create records from GFF3 (or dummy if unavailable)
    if path_gff3 is not None:
        records = io.load_gff3_as_seqrecords(path_gff3)

        # Add extra if FASTA contains more entries
        gff_seqids = {rec.id for rec in records}
        records_extra = [SeqRecord(Seq(None, 0), id=fasta_id) for fasta_id in fasta_ids if fasta_id not in gff_seqids]
        records += records_extra
    else:
        # Create dummy SeqRecords with IDs from FASTA
        records = [SeqRecord(Seq(None, 0), id=seq_id) for seq_id in fasta_ids]

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

        # Add TOPOLOGY feature /circular qualifier if "source" has "Is_circular"
        _handle_topology(rec)

        # Check characters in qualifier values
        _regularize_qualifier_value_letters(rec)

        # Join features if the key is in `joinables`
        if joinables:
            rec = _join_features(rec, joinables)

        # Merge exons with their parent such as mRNA, misc_RNA, V_segment, C_region
        _merge_exons_with_parent(rec)

        # Check start and stop codons in CDSs
        _fix_locations(rec, faidx=faidx)

        # Assign single value to /product and put the rest to /inference
        _assign_single_product(rec)

        # Make /gene have a single value; put the rest to /gene_synonym
        _make_gene_have_single_value(rec)

        # Copy qualifiers to children
        _copy_qualifiers_to_children(rec, data_qualifiers_to_children)

        # Remove duplicates within a qualifier
        _remove_duplicates_in_qualifiers(rec)

        # Sort features
        _sort_features(rec.features)

        # Check if SeqIDs have valid characters
        msg = (
            "\n\n"
            "Found invalid letter(s) in the 1st column of the GFF3 (or FASTA header if fasta-only): {}\n"
            "Consider running this script to correct entry names in the DDBJ annotation:\n\n"
            "   $ normalize-entry-names <output.ann>\n"
        )
        if utils.is_invalid_as_seqid(rec.id):
            logging.warning(msg.format(rec.id))

        yield rec

    faidx.close()
