from typing import Any, Dict, Generator, List, Tuple, Iterable
import collections
import yaml
import re
import pathlib
import gzip
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


def load_fasta_as_seq(filepath) -> Generator[Seq]:
    """
    Load FASTA file as Seq
    """
    p = pathlib.Path(filepath)
    if p.suffix == ".gz":
        with gzip.open(filepath, "rt") as f:
            return Bio.SeqIO.parse(f, "fasta")
    else:
        recs = Bio.SeqIO.parse(filepath, "fasta")
    return recs


def load_yaml_as_dict(filepath) -> Dict[str, Any]:
    """
    Load YAML as python dictionary
    """
    with open(filepath) as fp:
        d = yaml.safe_load(fp)
    return d


def merge_dicts(dict_args: Iterable[Dict]):
    """
    Given any number of dictionaries, shallow copy and merge into a new dict,
    precedence goes to key-value pairs in latter dictionaries.

    >>> merge_dicts([{'us', 1, 'canada, 1}, {'egypt': 20}, {'greece': 30, 'netherlands', 31}])
    {'us', 1, 'canada, 1, 'egypt': 20, 'greece': 30, 'netherlands', 31}
    """
    res = {}
    for d in dict_args:
        res.update(d)
    return res


def get_assembly_gap(seq: Seq) -> List[SeqFeature]:
    """
    Get assembly_gap features from seq.

    [NOTE] A location is of format (begin, end)
    with 1-based indexing, and inclusive in both sides.

    >>> s = Seq("atatnnngattacanccc")
    >>> get_assembly_gap(s)
    [(5, 7), (15, 15)]
    """
    locs = [FeatureLocation(start, end, strand=1) for (start, end) in _get_assembly_gap_locations(seq)]
    features = [SeqFeature(loc, type="assembly_gap") for loc in locs]
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


class TranslateQualifiers(object):
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

    def __init__(self, filepaths: Iterable[str]):
        self.paths = filepaths
        ds = [load_yaml_as_dict(p) for p in filepaths]
        self.trans_table = merge_dicts(ds)

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


class TranslateFeatures(object):
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

    def __init__(self, filepaths: Iterable[str]):
        self.paths = filepaths
        ds = [load_yaml_as_dict(p) for p in filepaths]
        self.d = merge_dicts(ds)

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


def join_features(record: SeqRecord, joinable=["CDS"]) -> SeqRecord:
    """
    Join features
    """

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

    def _helper(features: List[SeqFeature]) -> List[SeqFeature]:
        triples_or_features = collections.defaultdict(list)
        for f in features:
            if f.type not in joinable:
                triples_or_features[f] = [True]
            else:
                prod = f.qualifiers.get("product", None)
                triple = (f.type, f.id, prod)
                triples_or_features[triple].append(f)

        res = []
        seen = set()
        for triple_or_f, fs in triples_or_features.items():
            if isinstance(triple_or_f, SeqFeature):
                res.append(triple_or_f)
            else:
                if triple_or_f not in seen:
                    seen.add(triple_or_f)
                    joined_feature = _join(fs)
                    res.append(joined_feature)

        # join after upper levels
        for f in res:
            if hasattr(f, "sub_features") and f.sub_features:
                f.sub_features = _helper(f.sub_features)

        return res

    record.features = _helper(record.features)
    return record


def fix_codon_start_values(rec: SeqRecord):
    """Convert `codon_start` qualifier value
    from 0-based (in GFF3 'phase' column)
    to   1-based (in INSDC table definition)
    """
    def _fix_feature(feature: SeqFeature):
        trans = {'0': '1', '1': '2', '2': '3'}

        if hasattr(feature, "sub_features"):
            for f in feature.sub_features:
                _fix_feature(f)

        if "codon_start" in feature.qualifiers:
            x = feature.qualifiers["codon_start"][0]
            y = trans[x]
            feature.qualifiers["codon_start"] = [y]

    for f in rec.features:
        _fix_feature(f)


def run(path_to_gff3, path_to_fasta, trans_features, trans_qualifiers, is_joining=False) -> List[SeqRecord]:
    """
    Create SeqRecord and run all translations
    """
    records = load_gff3_as_seqrecords(path_to_gff3)
    f = TranslateFeatures(trans_features).run
    g = TranslateQualifiers(trans_qualifiers).run

    # translate features and qualifiers
    records = [g(f(rec)) for rec in records]

    # fix codon_start
    for rec in records:
        fix_codon_start_values(rec)

    # add assembly_gap
    seqs = load_fasta_as_seq(path_to_fasta)
    gaps = {seq.id: get_assembly_gap(seq) for seq in seqs}
    for rec in records:
        if rec.id in gaps:
            rec.features.extend(gaps[rec.id])

    # join features (such as CDS)
    if is_joining:
        records = [join_features(rec, joinable=["CDS"]) for rec in records]

    return records
