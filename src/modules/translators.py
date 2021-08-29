from typing import Any, Dict, List, Iterable
import json
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

def load_as_dict(filepath) -> Dict[str, Any]:
    """
    Load JSON file as python dictionary
    """
    with open(filepath) as fp:
        d = json.load(fp)
    return d

def merge_dicts(dict_args: Iterable[Dict]):
    """
    Given any number of dictionaries, shallow copy and merge into a new dict,
    precedence goes to key-value pairs in latter dictionaries.
    """
    res = {}
    for d in dict_args:
        res.update(d)
    return res

class GFF3AttributesToQualifiers(object):
    """Convert GFF3 attributes into qualifiers according to definition in JSON files.

    Paths to JSON files are required at the instantiation.
    JSON must have the schema

    ```
    <feature-name-as-in-type>: {
        "target": <corrected-feature-name>,
        "prefix": <prefix>,  (optional)
    }
    ```

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


class GFF3TypesToFeatures(object):
    """Convert GFF3 types into features according to definition in JSON files

    ```
    <name-as-in-type>: {
        "target": <corrected-feature-name>,
    }
    ```
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
