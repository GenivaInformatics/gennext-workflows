#!/usr/bin/env bash
# Builds results.duckdb with the real duckdb-snv and duckdb-snv-vcf-only SQL
# from tasks/base_tasks/duckdb.yaml, in the image those tasks pin, from two
# synthetic MAFs: one with every OncoKB column, and an older layout without
# Mutation Effect, the Dx/Px levels or VARIANT_IN_ONCOKB. Checks that both
# views build, that the new OncoKB columns have stable types, and that the
# older layout gives NULL instead of failing.
#
# Needs docker and python3 with PyYAML.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

TASKS=${TASKS:-tasks/base_tasks/duckdb.yaml}
FIXTURES=tests/fixtures/duckdb-snv-oncokb
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

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

query() {
  docker run --rm -v "$1:/mnt/output" "$2" \
    duckdb -noheader -csv /mnt/output/results.duckdb "$3"
}

NEW_COLUMNS="'VARIANT_IN_ONCOKB','MUTATION_EFFECT','MUTATION_EFFECT_CITATIONS',\
'LEVEL_Dx1','LEVEL_Dx2','LEVEL_Dx3','HIGHEST_DX_LEVEL','DX_CITATIONS',\
'LEVEL_Px1','LEVEL_Px2','LEVEL_Px3','HIGHEST_PX_LEVEL','PX_CITATIONS'"

for task in duckdb-snv duckdb-snv-vcf-only; do
  for layout in full trimmed; do
    case_dir="$WORK/$task-$layout"
    mkdir -p "$case_dir/out"
    image=$(extract_task "$task" "$case_dir/script.sh")
    case "$image" in
      *@sha256:*) ;;
      *) fail "$task: image is not pinned by digest: $image" ;;
    esac

    docker run --rm \
      -v "$REPO_ROOT/$FIXTURES/$layout/oncokb_snv.tsv:/mnt/input/oncokb_snv.tsv:ro" \
      -v "$REPO_ROOT/$FIXTURES/$layout/cancervar_snv.tsv:/mnt/input/cancervar_snv.tsv:ro" \
      -v "$case_dir/out:/mnt/output" \
      -v "$case_dir/script.sh:/script.sh:ro" \
      "$image" sh /script.sh >"$case_dir/log" 2>&1 \
      || { cat "$case_dir/log" >&2; fail "$task/$layout: build failed"; }

    for view in combined_snv advanced_combined_snv; do
      rows=$(query "$case_dir/out" "$image" "SELECT count(*) FROM $view")
      [ "$rows" = 1 ] || fail "$task/$layout: $view has $rows rows, expected 1"
    done

    types=$(query "$case_dir/out" "$image" \
      "SELECT column_name || '=' || data_type FROM information_schema.columns
       WHERE table_name = 'combined_snv' AND column_name IN ($NEW_COLUMNS)
       ORDER BY column_name")
    [ "$(wc -l <<<"$types")" = 13 ] \
      || fail "$task/$layout: combined_snv lacks some new columns: $types"
    grep -qx 'VARIANT_IN_ONCOKB=BOOLEAN' <<<"$types" \
      || fail "$task/$layout: VARIANT_IN_ONCOKB is not BOOLEAN: $types"
    if grep -v '^VARIANT_IN_ONCOKB=' <<<"$types" | grep -vq '=VARCHAR$'; then
      fail "$task/$layout: expected VARCHAR for the other new columns: $types"
    fi

    got=$(query "$case_dir/out" "$image" \
      "SELECT concat_ws('|', coalesce(CAST(VARIANT_IN_ONCOKB AS VARCHAR), 'NULL'),
              coalesce(MUTATION_EFFECT, 'NULL'), coalesce(HIGHEST_DX_LEVEL, 'NULL'),
              coalesce(DX_CITATIONS, 'NULL'), coalesce(MUTATION_EFFECT_CITATIONS, 'NULL'))
       FROM combined_snv")
    if [ "$layout" = full ]; then
      want='true|Gain-of-function|LEVEL_Dx2|23456|34567;45678'
    else
      want='NULL|NULL|NULL|NULL|NULL'
    fi
    [ "$got" = "$want" ] || fail "$task/$layout: got '$got', expected '$want'"
    echo "ok  $task/$layout"
  done
done

echo "PASS: duckdb-snv OncoKB columns"
