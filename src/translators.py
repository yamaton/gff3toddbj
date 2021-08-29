from typing import Any, Dict, List, Tuple, Iterable
import yaml
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

def load_as_dict(filepath) -> Dict[str, Any]:
    """
    Load YAML as python dictionary
    """
    with open(filepath) as fp:
        d = yaml.load(fp)
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


def get_assembly_gap(seq: Seq) -> List[Tuple[int, int]]:
    """
    Get assembly_gap locations from seq.

    [NOTE] A location is of format (begin, end)
    with 1-based indexing, and inclusive in both sides.

    >>> s = Seq("atatnnngattacanccc")
    [(5, 7), (15, 15)]
    """
    s = str(seq)
    patt = re.compile("n+")
    matches = patt.finditer(s)

    segments = []
    for m in matches:
        a, b = m.span()  # 0-based, left-inclusive, right-exclusive
        tup = (a + 1, b) # 1-based, both-inclusive
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
        ds = [load_as_dict(p) for p in filepaths]
        self.d = merge_dicts(ds)

    def run(self, record: SeqRecord) -> SeqRecord:
        """Modifies record according to the GFF3-attributes-to-qualifiers translation JSON data.
        """
        record.features = self._run_on_features(record.features)
        return record

    def _run_on_features(self, features: List[SeqFeature]) -> List[SeqFeature]:
        for feature in features:
            if hasattr(feature, 'sub_features'):
                feature.sub_features = self._run_on_features(feature.sub_features)
            if hasattr(feature, 'qualifiers'):
                feature.qualifiers = self._run_on_qualifiers(feature.qualifiers)
        return features

    def _run_on_qualifiers(self, qualifiers: Dict) -> Dict:
        res = dict()
        for name, vals in qualifiers.items():
            if name in self.d:
                # Replace the item name
                new_name = self.d[name]["target"]
                if new_name:
                    prefix = self.d[name].get("prefix", "")
                    res[new_name] = [prefix + v for v in vals]
                else:
                    # Remove the item from qualifiers if "target" is not set, or emtpy
                    pass
            else:
                res[name] = vals
        return res


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
        ds = [load_as_dict(p) for p in filepaths]
        self.d = merge_dicts(ds)

    def run(self, record: SeqRecord) -> SeqRecord:
        """Modifies record according to the GFF3-types-to-features translation JSON data.
        """
        record.features = self._run_on_features(record.features)
        return record

    def _run_on_features(self, features: List[SeqFeature]) -> List[SeqFeature]:
        for feature in features:
            if hasattr(feature, 'sub_features'):
                feature.sub_features = self._run_on_features(feature.sub_features)
            name = feature.type
            if name in self.d:
                new_name = self.d[name]["target"]
                feature.type = new_name
        return features

