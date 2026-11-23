"""COLUMNS in fusion_events.py must match the column list its load step reads.

Since 0.4.0 (contract v1) that is the prototype's load-events step in fusions-proto-pg; the F4
template duckdb-fusion-events in duckdb.yaml loads the 0.3.x column set and is not checked.

Run from the repository: python images/gennext-bpview/tests/test_columns_match.py
(skips when the workflow file is not there, e.g. in the image build)."""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
LOADER = os.path.join(HERE, "..", "..", "..", "workflows", "prototypes", "fusions-testing-ground",
                      "fusions-proto-pg.yaml")


def main():
    if not os.path.exists(LOADER):
        print("test_columns_match: skipped (no fusions-proto-pg.yaml)")
        return
    src = open(os.path.join(HERE, "..", "fusion_events.py")).read()
    ns = {}
    exec(re.search(r"(TIER_FLAGS = .*?\nCOLUMNS = \[.*?\n\])", src, re.S).group(1), ns)
    want = [(c, "VARCHAR" if t.endswith("[]") else t) for c, t in ns["COLUMNS"]]
    lists = sorted(c for c, t in ns["COLUMNS"] if t.endswith("[]"))
    text = open(LOADER).read()
    block = text[text.index("- name: load-events"):]
    block = block[:block.index("fusion_events_meta AS")]
    got = re.findall(r"'([a-z0-9_]+)': '([A-Z]+)'", block)
    assert got == want, f"load-events columns differ from COLUMNS:\n{got}\n{want}"
    split = sorted(re.findall(r"string_split\(([a-z0-9_]+), '\|'\)", block))
    assert split == lists, f"list columns split on load {split} != {lists}"
    print(f"test_columns_match: ok ({len(got)} columns, {len(split)} lists)")


if __name__ == "__main__":
    main()
