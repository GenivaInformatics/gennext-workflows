#!/usr/bin/env bash
# The OncoKB annotator v4 (image 0.3.0) writes three more MAF columns after
# ONCOGENIC: PATHOGENIC (Germline), PENETRANCE (Germline) and
# GENOMIC_INDICATOR (Germline). They're empty for our somatic MAFs.
#
# This builds results.duckdb with the real duckdb-snv and duckdb-snv-vcf-only
# SQL, in the image those tasks pin, from the full synthetic MAF (v3 layout)
# and from the same MAF in the v4 layout, and checks that combined_snv and
# advanced_combined_snv come out identical: the same columns and types and
# the same rows. The raw oncokb_snv table may differ only by the three empty
# columns.
#
# Needs docker and python3 with PyYAML.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

TASKS=${TASKS:-tasks/base_tasks/duckdb.yaml}
BASE=tests/fixtures/duckdb-snv-oncokb/full
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

# v3 = the fixture as it is; v4 = the germline columns inserted after ONCOGENIC.
python3 - "$BASE" "$WORK" <<'PY'
import csv, os, shutil, sys
base, work = sys.argv[1:]
germline = ["PATHOGENIC (Germline)", "PENETRANCE (Germline)", "GENOMIC_INDICATOR (Germline)"]
for layout in ("v3", "v4"):
    d = os.path.join(work, layout)
    os.makedirs(d)
    shutil.copy(os.path.join(base, "cancervar_snv.tsv"), d)
    with open(os.path.join(base, "oncokb_snv.tsv")) as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if layout == "v4":
        at = rows[0].index("ONCOGENIC") + 1
        rows = [r[:at] + (germline if i == 0 else ["", "", ""]) + r[at:]
                for i, r in enumerate(rows)]
    with open(os.path.join(d, "oncokb_snv.tsv"), "w", newline="") as f:
        csv.writer(f, delimiter="\t", lineterminator="\n").writerows(rows)
PY

for task in duckdb-snv duckdb-snv-vcf-only; do
  for layout in v3 v4; do
    dir="$WORK/$task-$layout"
    mkdir -p "$dir/out"
    image=$(extract_task "$task" "$dir/script.sh")
    docker run --rm \
      -v "$WORK/$layout/oncokb_snv.tsv:/mnt/input/oncokb_snv.tsv:ro" \
      -v "$WORK/$layout/cancervar_snv.tsv:/mnt/input/cancervar_snv.tsv:ro" \
      -v "$dir/out:/mnt/output" -v "$dir/script.sh:/script.sh:ro" \
      "$image" sh /script.sh >"$dir/log" 2>&1 \
      || { grep -m3 -i error "$dir/log" >&2; fail "$task/$layout: build failed"; }
    grep -qi "error" "$dir/log" && { grep -m3 -i error "$dir/log" >&2; fail "$task/$layout: SQL error"; }
  done

  # Views name their tables unqualified, so each file is opened on its own:
  # its column list is read there and each view is exported to Parquet.
  cmp="$WORK/$task-cmp"
  mkdir -p "$cmp"
  for layout in v3 v4; do
    cp "$WORK/$task-$layout/out/results.duckdb" "$cmp/$layout.duckdb"
    docker run --rm -v "$cmp:/mnt/cmp" "$image" duckdb -readonly -noheader -list "/mnt/cmp/$layout.duckdb" -c "
      COPY (SELECT table_name, column_name, data_type, ordinal_position FROM information_schema.columns
            WHERE table_name IN ('combined_snv', 'advanced_combined_snv', 'oncokb_snv'))
        TO '/mnt/cmp/$layout-columns.parquet' (FORMAT PARQUET);
      COPY (SELECT * FROM combined_snv) TO '/mnt/cmp/$layout-combined_snv.parquet' (FORMAT PARQUET);
      COPY (SELECT * FROM advanced_combined_snv) TO '/mnt/cmp/$layout-advanced_combined_snv.parquet' (FORMAT PARQUET);" \
      >"$cmp/$layout.log" 2>&1
    grep -qi error "$cmp/$layout.log" && { cat "$cmp/$layout.log" >&2; fail "$task/$layout: export failed"; }
  done
  q() {
    docker run --rm -v "$cmp:/mnt/cmp" "$image" duckdb -noheader -list -c "$1"
  }
  cols() {
    q "SELECT string_agg(column_name || ':' || data_type, ',' ORDER BY ordinal_position)
       FROM '/mnt/cmp/$1-columns.parquet' WHERE table_name = '$2';"
  }
  for view in combined_snv advanced_combined_snv; do
    [ -n "$(cols v3 $view)" ] || fail "$task: no columns read for $view"
    [ "$(cols v3 $view)" = "$(cols v4 $view)" ] || fail "$task: $view columns or types differ between layouts"
    n=$(q "SELECT count(*) FROM '/mnt/cmp/v3-$view.parquet'")
    [ "$n" -gt 0 ] || fail "$task: $view is empty"
    [ "$n" = "$(q "SELECT count(*) FROM '/mnt/cmp/v4-$view.parquet'")" ] || fail "$task: $view row counts differ"
    for pair in "v3 v4" "v4 v3"; do
      set -- $pair
      d=$(q "SELECT count(*) FROM (SELECT * FROM '/mnt/cmp/$1-$view.parquet' EXCEPT ALL SELECT * FROM '/mnt/cmp/$2-$view.parquet')")
      [ "$d" = 0 ] || fail "$task: $view has $d rows in $1 that $2 lacks"
    done
  done
  extra=$(q "SELECT string_agg(column_name, '|' ORDER BY column_name) FROM (
    SELECT column_name FROM '/mnt/cmp/v4-columns.parquet' WHERE table_name = 'oncokb_snv'
    EXCEPT SELECT column_name FROM '/mnt/cmp/v3-columns.parquet' WHERE table_name = 'oncokb_snv')")
  [ "$extra" = "GENOMIC_INDICATOR (Germline)|PATHOGENIC (Germline)|PENETRANCE (Germline)" ] \
    || fail "$task: oncokb_snv extra columns in v4 are '$extra'"
  echo "ok  $task: both views identical for the v3 and v4 MAF layouts"
done

echo "PASS: duckdb-snv reads the v3 and v4 annotator layouts"
