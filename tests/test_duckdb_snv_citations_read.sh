#!/usr/bin/env bash
# duckdb-snv reads the OncoKB MAF with read_csv_auto(sample_size=20000). When
# the sampled rows of a citations column hold only single PMIDs, DuckDB types
# it BIGINT, and the read fails at the first "1;2" value further down. (DuckDB
# 0.8.1 samples a small file whole, so it takes a file a bit over 20k rows.)
# This builds results.duckdb with the real duckdb-snv and duckdb-snv-vcf-only
# SQL, in the image those tasks pin, from generated MAFs of 25,000 rows, with
# single PMIDs everywhere except row 24,978, which has several:
#
#   late-multi  every citations column. Must build, with the citations columns
#               VARCHAR and the multi-PMID value intact.
#   tx-only     only TX_CITATIONS, as in an older MAF (types= may only name
#               columns the file has). Must build.
#
# TASKS selects the task file (default: this checkout's). On origin/main's
# task file duckdb-snv fails on late-multi:
#   git show origin/main:tasks/base_tasks/duckdb.yaml > /tmp/main.yaml
#   TASKS=/tmp/main.yaml tests/test_duckdb_snv_citations_read.sh
#
# Needs docker and python3 with PyYAML.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

TASKS=${TASKS:-tasks/base_tasks/duckdb.yaml}
BASE=tests/fixtures/duckdb-snv-oncokb/full
ROWS=25000
MULTI_AT=24978
WORK=$(mktemp -d)
trap '[ -n "${KEEP_WORK:-}" ] || rm -rf "$WORK"' EXIT

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

# Writes the task's script with the output name filled in; prints its image.
extract_task() {
  python3 - "$TASKS" "$1" "$2" <<'PY'
import sys, yaml
tasks, name, out = sys.argv[1:]
doc = yaml.safe_load(open(tasks))
task = next(t for t in doc["spec"]["templates"] if t["name"] == name)
source = task["script"]["source"].replace(
    "{{=sprig.osBase(inputs.parameters['output-duckdb'])}}", "results.duckdb")
if "{{" in source:
    sys.exit(f"{name}: unhandled template expression in script")
open(out, "w").write(source)
print(task["script"]["image"])
PY
}

python3 - "$BASE" "$WORK" "$ROWS" "$MULTI_AT" <<'PY'
import csv, os, sys
base, work, rows, multi_at = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
def read(name):
    with open(os.path.join(base, name)) as f:
        r = list(csv.reader(f, delimiter="\t"))
    return r[0], r[1]
onc_h, onc = read("oncokb_snv.tsv")
cv_h, cv = read("cancervar_snv.tsv")
cites = [i for i, h in enumerate(onc_h) if h.endswith("_CITATIONS")]
o = {h: i for i, h in enumerate(onc_h)}
c = {h: i for i, h in enumerate(cv_h)}
multi = "40332046;40997296;23456;34567;45678"
for layout in ("late-multi", "tx-only"):
    keep = [i for i in range(len(onc_h))
            if layout == "late-multi" or i not in cites or onc_h[i] == "TX_CITATIONS"]
    d = os.path.join(work, layout)
    os.makedirs(d)
    with open(os.path.join(d, "oncokb_snv.tsv"), "w") as fo, open(os.path.join(d, "cancervar_snv.tsv"), "w") as fc:
        fo.write("\t".join(onc_h[i] for i in keep) + "\n")
        fc.write("\t".join(cv_h) + "\n")
        for n in range(1, rows + 1):
            r, k = list(onc), list(cv)
            pos = str(1000000 + n)
            r[o["Start_Position"]] = r[o["End_Position"]] = pos
            k[c["Start"]] = k[c["End"]] = pos
            for i in cites:
                r[i] = multi if n == multi_at else str(10000000 + n)
            fo.write("\t".join(r[i] for i in keep) + "\n")
            fc.write("\t".join(k) + "\n")
PY

query() {
  docker run --rm -v "$1:/mnt/output" "$2" duckdb -noheader -csv /mnt/output/results.duckdb "$3"
}

for task in duckdb-snv duckdb-snv-vcf-only; do
  for layout in late-multi tx-only; do
    dir="$WORK/$task-$layout"
    mkdir -p "$dir/out"
    image=$(extract_task "$task" "$dir/script.sh")
    docker run --rm \
      -v "$WORK/$layout/oncokb_snv.tsv:/mnt/input/oncokb_snv.tsv:ro" \
      -v "$WORK/$layout/cancervar_snv.tsv:/mnt/input/cancervar_snv.tsv:ro" \
      -v "$dir/out:/mnt/output" -v "$dir/script.sh:/script.sh:ro" \
      "$image" sh /script.sh >"$dir/log" 2>&1 \
      || { grep -m3 -i error "$dir/log" >&2; fail "$task/$layout: build failed"; }

    for view in combined_snv advanced_combined_snv; do
      rows=$(query "$dir/out" "$image" "SELECT count(*) FROM $view")
      [ "$rows" = "$ROWS" ] || fail "$task/$layout: $view has $rows rows, expected $ROWS"
    done
    types=$(query "$dir/out" "$image" "SELECT data_type FROM information_schema.columns
      WHERE table_name = 'oncokb_snv' AND column_name = 'TX_CITATIONS'")
    [ "$types" = VARCHAR ] || fail "$task/$layout: oncokb_snv.TX_CITATIONS is $types, expected VARCHAR"
    if [ "$layout" = late-multi ]; then
      types=$(query "$dir/out" "$image" "SELECT DISTINCT data_type FROM information_schema.columns
        WHERE table_name = 'oncokb_snv' AND column_name LIKE '%\_CITATIONS' ESCAPE '\\'")
      [ "$types" = VARCHAR ] || fail "$task: oncokb_snv citations columns are $types, expected VARCHAR"
      got=$(query "$dir/out" "$image" "SELECT DISTINCT concat_ws('|', TX_CITATIONS, DX_CITATIONS,
        PX_CITATIONS, MUTATION_EFFECT_CITATIONS) FROM combined_snv WHERE TX_CITATIONS LIKE '%;%'")
      want="40332046;40997296;23456;34567;45678|40332046;40997296;23456;34567;45678|40332046;40997296;23456;34567;45678|40332046;40997296;23456;34567;45678"
      [ "$got" = "$want" ] || fail "$task: multi-PMID row came back as '$got'"
    fi
    echo "ok  $task/$layout"
  done
done

echo "PASS: duckdb-snv reads citations as text"
