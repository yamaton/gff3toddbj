"""
Check if parsing (f) DDBJ annotation, followed by formatting (g) into DDBJ equals idnetity.
i.e. (g).(f) = identity

[NOTE] The DDBJ parser ignores COMMON hence the equality does NOT hold if COMMON exists.

"""
import collections
import pathlib
import argparse
import filecmp
import sys
from typing import Optional

from gff3toddbj.parser import load_ddbj
from gff3toddbj.formatter import DDBJFormatter

WORK_DIR = pathlib.Path(__file__).parent


def _test_identity(path_anno_in, path_anno_out=None):

    records = load_ddbj(path_anno_in)
    fmt = DDBJFormatter(collections.OrderedDict(), collections.OrderedDict())
    gen = fmt.run(records, ignore_rules=True)

    print("output: {}".format(str(path_anno_out)))
    with open(path_anno_out, "w") as f:
        for line in gen:
            print(line, file=f)


def check_sanity(path_input) -> Optional[bool]:
    with open(path_input, "r") as f:
        content = f.read()
    if "COMMON" in content:
        return None

    p = pathlib.Path(path_input)
    path_output = p.parent / (p.stem + "__out" + p.suffix)
    _test_identity(path_input, path_output)

    res = filecmp.cmp(path_input, path_output)
    if res:
        path_output.unlink()
    return res


def test_sanity():
    annotation_files = ["augustus.ann", "maker.ann", "prokka.ann", "fasta_only.ann"]
    paths = [WORK_DIR / "golden" / name for name in annotation_files]
    for path in paths:
        assert check_sanity(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to DDBJ")
    args = parser.parse_args()
    path_input = args.input
    res = check_sanity(path_input)

    print()
    print("----------------")
    if res is None:
        print("    SKIPPED   ")
    elif res:
        print("    PASS      ")
    else:
        print("    FAIL      ")
    print("----------------")


if __name__ == "__main__":
    main()
