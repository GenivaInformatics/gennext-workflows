#!/usr/bin/env bash
# The SNV views number variants by cancervar_snv's stored order (row_id). With
# preserve_insertion_order off, a parallel read of a large CancerVar file
# stores rows in any order, so the same input got other ids on every build,
# and the CancerVar evidence list (an unordered aggregate) other item orders
# (#53).
#
# This builds results.duckdb twice with the real duckdb-snv and
# duckdb-snv-vcf-only SQL, in the image those tasks pin, from a generated
# sample of 200,000 rows (more than one DuckDB row group, 122,880). Each row
# has a different set of CancerVar evidence items. Checks, per task:
#   - both builds give identical row_id and evidence on every row;
#   - row_id is the CancerVar file's line order;
#   - evidence lists its items in CBP order (CBP1, CBP2, ...).
#
# TASKS selects the task file (default: this checkout's). On origin/main,
# duckdb-snv fails the first check:
#   git show origin/main:tasks/base_tasks/duckdb.yaml > /tmp/main.yaml
#   TASKS=/tmp/main.yaml tests/test_duckdb_snv_row_order.sh
#
# Needs docker and python3 with PyYAML.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

TASKS=${TASKS:-tasks/base_tasks/duckdb.yaml}
BASE=tests/fixtures/duckdb-snv-oncokb/full
ROWS=200000
WORK=$(mktemp -d)
trap '[ -n "${KEEP_WORK:-}" ] || rm -rf "$WORK"' EXIT

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

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

# Rows sorted by chromosome then position, as ANNOVAR writes them; row n has
# the evidence items of n's bits, so most rows have several.
python3 - "$BASE" "$WORK" "$ROWS" <<'PY'
import csv, os, sys
base, work, rows = sys.argv[1], sys.argv[2], int(sys.argv[3])
def read(name):
    r = list(csv.reader(open(os.path.join(base, name)), delimiter="\t"))
    return r[0], r[1]
oh, o = read("oncokb_snv.tsv")
ch, c = read("cancervar_snv.tsv")
oi = {h: i for i, h in enumerate(oh)}
ci = {h: i for i, h in enumerate(ch)}
ev = ci["CancerVar: CancerVar and Evidence"]
with open(os.path.join(work, "oncokb_snv.tsv"), "w") as fo, \
        open(os.path.join(work, "cancervar_snv.tsv"), "w") as fc:
    fo.write("\t".join(oh) + "\n")
    fc.write("\t".join(ch) + "\n")
    for n in range(1, rows + 1):
        chrom = str((n * 22) // (rows + 1) + 1)
        pos = str(1_000_000 + n)
        r, k = list(o), list(c)
        r[oi["Chromosome"]] = "chr" + chrom
        k[ci["#Chr"]] = chrom
        r[oi["Start_Position"]] = r[oi["End_Position"]] = pos
        k[ci["Start"]] = k[ci["End"]] = pos
        evs = [str((n >> b) & 1) for b in range(12)]
        k[ev] = f"CancerVar: 2#Tier_II_potential EVS=[{', '.join(evs)}]"
        fo.write("\t".join(r) + "\n")
        fc.write("\t".join(k) + "\n")
PY

for task in duckdb-snv duckdb-snv-vcf-only; do
  for build in 1 2; do
    dir="$WORK/$task-$build"
    mkdir -p "$dir/out"
    image=$(extract_task "$task" "$dir/script.sh")
    docker run --rm \
      -v "$WORK/oncokb_snv.tsv:/mnt/input/oncokb_snv.tsv:ro" \
      -v "$WORK/cancervar_snv.tsv:/mnt/input/cancervar_snv.tsv:ro" \
      -v "$dir/out:/mnt/output" -v "$dir/script.sh:/script.sh:ro" \
      "$image" sh /script.sh >"$dir/log" 2>&1 \
      || { grep -m3 -i error "$dir/log" >&2; fail "$task: build $build failed"; }
  done
  q() {
    docker run --rm -v "$WORK/$task-$1/out:/mnt/output" "$image" \
      duckdb -readonly -noheader -list /mnt/output/results.duckdb "$2"
  }
  for view in combined_snv advanced_combined_snv; do
    sig="SELECT md5(string_agg(row_id::VARCHAR || ':' || Start::VARCHAR || ':' || coalesce(\"CancerVar: Evidence\", ''),
      ',' ORDER BY Start)) FROM $view"
    [ "$(q 1 "$sig")" = "$(q 2 "$sig")" ] || fail "$task: $view differs between two builds of the same input"
    n=$(q 1 "SELECT count(*) FROM $view")
    [ "$n" = "$ROWS" ] || fail "$task: $view has $n rows, expected $ROWS"
    off=$(q 1 "SELECT count(*) FROM (SELECT row_id, row_number() OVER (ORDER BY Start) AS line FROM $view) WHERE row_id <> line")
    [ "$off" = 0 ] || fail "$task: $view: $off rows whose row_id isn't their line in the CancerVar file"
    # Row n's items are the set bits of n, in CBP order.
    bad=$(q 1 "SELECT count(*) FROM $view WHERE coalesce(\"CancerVar: Evidence\", '') <> coalesce((
      SELECT string_agg('CBP' || (b + 1), ', ' ORDER BY b) FROM range(12) t(b)
      WHERE (CAST(Start AS BIGINT) - 1000000) >> b & 1 = 1), '')")
    [ "$bad" = 0 ] || fail "$task: $view: $bad rows whose evidence isn't in CBP order"
  done
  echo "ok  $task: same row_id and evidence on two builds, row_id in file order, evidence in CBP order"
done

echo "PASS: SNV row ids and evidence don't depend on the build"
