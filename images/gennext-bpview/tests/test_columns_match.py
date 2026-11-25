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
LOADERS = [  # (file, the template name of its load step)
    (os.path.join(HERE, "..", "..", "..", "workflows", "prototypes", "fusions-testing-ground", "fusions-proto-pg.yaml"),
     "load-events"),
    (os.path.join(HERE, "..", "..", "..", "tasks", "base_tasks", "fusion-events-v1.yaml"), "fusion-events-load"),
]


def main():
    src = open(os.path.join(HERE, "..", "fusion_events.py")).read()
    ns = {}
    exec(re.search(r"(TIER_FLAGS = .*?\nCOLUMNS = \[.*?\n\])", src, re.S).group(1), ns)
    want = [(c, "VARCHAR" if t.endswith("[]") else t) for c, t in ns["COLUMNS"]]
    lists = sorted(c for c, t in ns["COLUMNS"] if t.endswith("[]"))
    checked = 0
    for path, name in LOADERS:
        if not os.path.exists(path):
            continue
        text = open(path).read()
        block = text[text.index(f"- name: {name}"):]
        block = block[:block.index("fusion_events_meta AS")]
        got = re.findall(r"'([a-z0-9_]+)': '([A-Z]+)'", block)
        assert got == want, f"{os.path.basename(path)} {name}: columns differ from COLUMNS:\n{got}\n{want}"
        split = sorted(re.findall(r"string_split\(([a-z0-9_]+), '\|'\)", block))
        assert split == lists, f"{os.path.basename(path)} {name}: list columns split on load {split} != {lists}"
        checked += 1
    if not checked:
        print("test_columns_match: skipped (no load step found)")
        return
    print(f"test_columns_match: ok ({len(want)} columns, {len(lists)} lists, {checked} load step(s))")


if __name__ == "__main__":
    main()
