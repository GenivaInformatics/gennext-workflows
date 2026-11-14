#!/usr/bin/env bash
# When every OncoKB answer in a sample has a single PMID, or none,
# read_csv_auto types MUTATION_EFFECT_CITATIONS as BIGINT. api-edge reads the
# column as text, so the whole grid fails to load for that sample (#41).
#
# This builds results.duckdb with the real duckdb-sv and duckdb-fusions SQL,
# in the image those tasks use, from synthetic inputs that hold every column
# the SQL reads. Citations are single PMIDs. It checks that the raw table
# reads the column as BIGINT, so the case is exercised, and that every view
# returns it as VARCHAR.
#
# TASKS selects the task file (default: this checkout's). On origin/main the
# views return BIGINT:
#   git show origin/main:tasks/base_tasks/duckdb.yaml > /tmp/main.yaml
#   TASKS=/tmp/main.yaml tests/test_duckdb_sv_fusions_citations.sh
#
# Needs docker and python3 with PyYAML.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

TASKS=${TASKS:-tasks/base_tasks/duckdb.yaml}
WORK=$(mktemp -d)
trap '[ -n "${KEEP_WORK:-}" ] || rm -rf "$WORK"' EXIT

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

# Writes each task's script (output name filled in) and synthetic inputs with
# every column its SQL reads; prints "<task> <image>" per task.
python3 - "$TASKS" "$WORK" > "$WORK/tasks.txt" <<'PY'
import csv, os, re, sys, yaml
tasks, work = sys.argv[1:]
doc = yaml.safe_load(open(tasks))
by_name = {t["name"]: t for t in doc["spec"]["templates"]}

def columns(sql, table, unqualified_from=None):
    cols = set(re.findall(rf"\b{table}\.([A-Za-z0-9_]+)", sql))
    if unqualified_from is not None:
        # advanced_combined_fusions reads the table unqualified: CAST(Name AS ...)
        view = sql[sql.index(unqualified_from):]
        cols |= set(re.findall(r"CAST\(([A-Za-z0-9_]+) AS", view))
    return sorted(cols)

def value(col):
    if col.endswith("_CITATIONS"):
        return "12345678"          # a single PMID: read_csv_auto picks BIGINT
    if col in ("AnnotSV_ID",):
        return "sv1"
    if col in ("Gene_name", "GeneA", "GeneB", "LeftGene", "RightGene"):
        return "BRAF"
    if col in ("Fusion",):
        return "BRAF--KIAA1549"
    return "1"

def write(path, cols):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for _ in range(3):
            w.writerow([value(c) for c in cols])

for name, inputs in (("duckdb-sv", {"oncokb_sv": None, "annotsv_sv": None}),
                     ("duckdb-fusions", {"oncokb_fusions": "CREATE VIEW advanced_combined_fusions"})):
    task = by_name[name]
    src = task["script"]["source"].replace(
        "{{=sprig.osBase(inputs.parameters['output-duckdb'])}}", "results.duckdb")
    if "{{" in src:
        sys.exit(f"{name}: unhandled template expression in script")
    d = os.path.join(work, name)
    os.makedirs(os.path.join(d, "out"))
    open(os.path.join(d, "script.sql"), "w").write(src)
    for table, unq in inputs.items():
        write(os.path.join(d, f"{table}.tsv"), columns(src, table, unq))
    print(name, task["script"]["image"])
PY

query() {
  docker run --rm -v "$1:/mnt/output" "$2" duckdb -readonly -noheader -csv /mnt/output/results.duckdb "$3"
}
type_of() {
  query "$1" "$2" "SELECT data_type FROM information_schema.columns
    WHERE table_name = '$3' AND column_name = 'MUTATION_EFFECT_CITATIONS'"
}

while read -r task image; do
  dir="$WORK/$task"
  case "$task" in
    duckdb-sv)
      mounts=(-v "$dir/oncokb_sv.tsv:/mnt/input/oncokb/oncokb_sv.tsv:ro"
              -v "$dir/annotsv_sv.tsv:/mnt/input/annotsv/annotsv_sv.tsv:ro")
      raw=oncokb_sv; views="combined_sv advanced_combined_sv" ;;
    duckdb-fusions)
      mounts=(-v "$dir/oncokb_fusions.tsv:/mnt/input/oncokb/oncokb_fusions.tsv:ro")
      raw=oncokb_fusions; views="combined_fusions advanced_combined_fusions" ;;
  esac
  docker run --rm "${mounts[@]}" -v "$dir/out:/mnt/output" -v "$dir/script.sql:/script.sql:ro" \
    "$image" duckdb -init /script.sql -c "SELECT 1" >"$dir/log" 2>&1 \
    || { cat "$dir/log" >&2; fail "$task: build failed"; }
  grep -qi "error" "$dir/log" && { grep -i -m3 error "$dir/log" >&2; fail "$task: SQL error"; }

  t=$(type_of "$dir/out" "$image" "$raw")
  [ "$t" = BIGINT ] || fail "$task: fixture doesn't exercise the case ($raw.MUTATION_EFFECT_CITATIONS is $t)"
  for view in $views; do
    t=$(type_of "$dir/out" "$image" "$view")
    [ "$t" = VARCHAR ] || fail "$task: $view.MUTATION_EFFECT_CITATIONS is $t, expected VARCHAR"
    v=$(query "$dir/out" "$image" "SELECT DISTINCT MUTATION_EFFECT_CITATIONS FROM $view")
    [ "$v" = 12345678 ] || fail "$task: $view.MUTATION_EFFECT_CITATIONS came back as '$v'"
  done
  echo "ok  $task: MUTATION_EFFECT_CITATIONS is VARCHAR in $views"
done < "$WORK/tasks.txt"

echo "PASS: SV and fusion views return citations as text"
