"""COLUMNS in fusion_events.py must match the column list duckdb-fusion-events loads.

Run from the repository: python images/gennext-bpview/tests/test_columns_match.py
(skips when tasks/base_tasks/duckdb.yaml is not there, e.g. in the image build)."""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
DUCKDB_YAML = os.path.join(HERE, "..", "..", "..", "tasks", "base_tasks", "duckdb.yaml")


def main():
    if not os.path.exists(DUCKDB_YAML):
        print("test_columns_match: skipped (no duckdb.yaml)")
        return
    src = open(os.path.join(HERE, "..", "fusion_events.py")).read()
    ns = {}
    exec(re.search(r"(TIER_FLAGS = .*?\nCOLUMNS = \[.*?\n\])", src, re.S).group(1), ns)
    want = [(c, "VARCHAR" if t.endswith("[]") else t) for c, t in ns["COLUMNS"]]
    lists = sorted(c for c, t in ns["COLUMNS"] if t.endswith("[]"))
    text = open(DUCKDB_YAML).read()
    block = text[text.index("- name: duckdb-fusion-events"):]
    block = block[:block.index("volumeMounts:")]
    got = re.findall(r"'([a-z0-9_]+)': '([A-Z]+)'", block)
    assert got == want, f"duckdb.yaml columns differ from COLUMNS:\n{got}\n{want}"
    split = sorted(re.findall(r"string_split\(([a-z0-9_]+), '\|'\)", block))
    assert split == lists, f"list columns split on load {split} != {lists}"
    print(f"test_columns_match: ok ({len(got)} columns, {len(split)} lists)")


if __name__ == "__main__":
    main()
