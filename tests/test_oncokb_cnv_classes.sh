#!/usr/bin/env bash
# Checks the CNV path to OncoKB's copy-number classes (#31) end to end, with
# the real task scripts from this repo, in the images those tasks use, on a
# synthetic AnnotSV TSV shaped like CNVkit calls:
#
#   1. annotsv-to-oncokb-cna-format writes the same filtered TSV as
#      annotsv-to-oncokb-format, the same first four OncoKB columns, and the
#      expected Copy_Number_Alteration per row; amplification-min-copies moves
#      the Amplification cut-off.
#   2. oncokb-cnv (CnaAnnotator) against demo.oncokb.org, next to what
#      oncokb-sv runs today (StructuralVariantAnnotator): same OncoKB columns,
#      rows in input order, rows without a class left empty; the Genome kit
#      skip still writes empty OncoKB columns without calling the API.
#   3. duckdb-cnv builds combined_cnv and advanced_combined_cnv from both
#      outputs with the same columns and types, and the same rows apart from
#      the OncoKB values.
#
# The demo API needs no token and knows only BRAF, TP53 and ROS1.
# Needs docker, network access and python3 with PyYAML.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

FIXTURES=$REPO_ROOT/tests/fixtures/oncokb-cnv-classes
DEMO_URL=${DEMO_URL:-https://demo.oncokb.org}
WORK=$(mktemp -d)
trap '[ -n "${KEEP_WORK:-}" ] || rm -rf "$WORK"' EXIT

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

# extract_task FILE TASK OUT [FROM=TO ...]: writes the task's script with the
# given template expressions filled in; prints its image.
extract_task() {
  python3 - "$@" <<'PY'
import sys, yaml
path, name, out, *subs = sys.argv[1:]
doc = yaml.safe_load(open(path))
task = next(t for t in doc["spec"]["templates"] if t["name"] == name)
source = task["script"]["source"]
for sub in subs:
    old, new = sub.rsplit("=", 1)  # values never contain "="
    if old not in source:
        sys.exit(f"{name}: {old!r} not in script")
    source = source.replace(old, new)
if "{{" in source:
    sys.exit(f"{name}: unhandled template expression in script")
open(out, "w").write(source)
print(task["script"]["image"])
PY
}

run_format() { # run_format TASK DIR [AMP]
  local task=$1 dir=$2 amp=${3:-4} subs image
  mkdir -p "$dir/oncokb" "$dir/filtered"
  subs=(
    "{{inputs.parameters.sample-id}}=S1"
    "{{=sprig.osBase(inputs.parameters['annotsv-tsv'])}}=annotsv_cnv.tsv"
    "{{=sprig.osBase(inputs.parameters['annotsv-filtered-out-tsv'])}}=annotsv_cnv_filtered.tsv"
    "{{=sprig.osBase(inputs.parameters['annotsv-oncokb-out-tsv'])}}=cnv_oncokb.tsv"
  )
  if [ "$task" = annotsv-to-oncokb-cna-format ]; then
    subs+=("{{inputs.parameters.amplification-min-copies}}=$amp")
  fi
  image=$(extract_task tasks/base_tasks/filter-annotsv.yaml "$task" "$dir/script.py" "${subs[@]}")
  docker run --rm \
    -v "$FIXTURES/annotsv_cnv.tsv:/mnt/input/annotsv/annotsv_cnv.tsv:ro" \
    -v "$dir/oncokb:/mnt/output/oncokb" -v "$dir/filtered:/mnt/output/filtered" \
    -v "$dir/script.py:/script.py:ro" -w /mnt/input "$image" python /script.py >"$dir/log" 2>&1 \
    || { cat "$dir/log" >&2; fail "$task failed"; }
}

# classes DIR: "AnnotSV_ID<TAB>gene<TAB>class" per OncoKB row.
classes() {
  python3 - "$1/oncokb/cnv_oncokb.tsv" <<'PY'
import csv, sys
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    print(f"{r['AnnotSV_ID']}\t{r['GeneA']}\t{r['Copy_Number_Alteration']}")
PY
}

# --- 1. Mapping --------------------------------------------------------------
run_format annotsv-to-oncokb-format "$WORK/old"
run_format annotsv-to-oncokb-cna-format "$WORK/new"
cmp -s "$WORK/old/filtered/annotsv_cnv_filtered.tsv" "$WORK/new/filtered/annotsv_cnv_filtered.tsv" \
  || fail "filtered TSV differs from annotsv-to-oncokb-format"
cmp -s "$WORK/old/oncokb/cnv_oncokb.tsv" <(cut -f1-4 "$WORK/new/oncokb/cnv_oncokb.tsv") \
  || fail "AnnotSV_ID/GeneA/GeneB/Sv_Type differ from annotsv-to-oncokb-format"
diff <(tail -n +2 "$FIXTURES/expected_classes.tsv") <(classes "$WORK/new") \
  || fail "Copy_Number_Alteration differs from expected_classes.tsv"
echo "ok  mapping (filtered TSV and SV columns unchanged, classes as expected)"

# A stricter cut-off: CN 4 (and log2 0.9, which CNVkit calls CN 4) become Gain.
run_format annotsv-to-oncokb-cna-format "$WORK/new5" 5
diff <(tail -n +2 "$FIXTURES/expected_classes.tsv" | python3 -c '
import sys
for line in sys.stdin:
    i, g, c = line.rstrip("\n").split("\t")
    if c == "Amplification" and i != "7_140719327_140924929_DUP_1":
        c = "Gain"
    print(f"{i}\t{g}\t{c}")') <(classes "$WORK/new5") \
  || fail "amplification-min-copies=5 did not move the cut-off"
echo "ok  amplification-min-copies=5 moves CN 4 to Gain"

mkdir -p "$WORK/bad/oncokb" "$WORK/bad/filtered"
FORMAT_IMAGE=$(extract_task tasks/base_tasks/filter-annotsv.yaml annotsv-to-oncokb-cna-format "$WORK/bad/script.py" \
  "{{inputs.parameters.sample-id}}=S1" \
  "{{=sprig.osBase(inputs.parameters['annotsv-tsv'])}}=annotsv_cnv.tsv" \
  "{{=sprig.osBase(inputs.parameters['annotsv-filtered-out-tsv'])}}=f.tsv" \
  "{{=sprig.osBase(inputs.parameters['annotsv-oncokb-out-tsv'])}}=o.tsv" \
  "{{inputs.parameters.amplification-min-copies}}=2")
if docker run --rm -v "$FIXTURES/annotsv_cnv.tsv:/mnt/input/annotsv/annotsv_cnv.tsv:ro" \
  -v "$WORK/bad/oncokb:/mnt/output/oncokb" -v "$WORK/bad/filtered:/mnt/output/filtered" \
  -v "$WORK/bad/script.py:/script.py:ro" -w /mnt/input "$FORMAT_IMAGE" python /script.py >/dev/null 2>&1; then
  fail "amplification-min-copies=2 was accepted"
fi
echo "ok  amplification-min-copies below 3 is refused"

# --- 2. OncoKB ---------------------------------------------------------------
ANNOTATOR_IMAGE=$(extract_task tasks/base_tasks/oncokb.yaml oncokb-cnv "$WORK/oncokb-cnv.py" \
  "{{=sprig.osBase(inputs.parameters['oncokb-sv-in-f'])}}=cnv_oncokb.tsv" \
  "{{=sprig.osBase(inputs.parameters['oncokb-sv-out-f'])}}=cnv_oncokb_annot.tsv" \
  "{{inputs.parameters.regions-file}}=/ref/kit.bed" \
  "{{inputs.parameters.tumor-type}}=." \
  "{{inputs.parameters.max-attempts}}=1" \
  "{{inputs.parameters.retry-delay-seconds}}=0" \
  "\"-f\", \"individual\", \"-z\"]=\"-f\", \"individual\", \"-z\", \"-u\", \"$DEMO_URL\"]")
case "$ANNOTATOR_IMAGE" in
  *@sha256:*) ;;
  *) fail "oncokb-cnv image is not pinned by digest: $ANNOTATOR_IMAGE" ;;
esac

# Today: oncokb-sv runs filtered_sv_annotator.py, which runs this.
mkdir -p "$WORK/old/annot"
docker run --rm -e ONCOKB_API_TOKEN=demo -v "$WORK/old/oncokb:/mnt/input:ro" -v "$WORK/old/annot:/mnt/output" \
  --entrypoint python3 "$ANNOTATOR_IMAGE" /app/oncokb/StructuralVariantAnnotator.py \
  -i /mnt/input/cnv_oncokb.tsv -o /mnt/output/cnv_oncokb_annot.tsv -b demo -u "$DEMO_URL" \
  >"$WORK/old/annot.log" 2>&1 || { cat "$WORK/old/annot.log" >&2; fail "StructuralVariantAnnotator failed"; }

mkdir -p "$WORK/new/annot"
docker run --rm -e ONCOKB_API_TOKEN=demo -v "$WORK/new/oncokb:/mnt/input:ro" -v "$WORK/new/annot:/mnt/output" \
  -v "$WORK/oncokb-cnv.py:/task.py:ro" --entrypoint python3 "$ANNOTATOR_IMAGE" /task.py \
  >"$WORK/new/annot.log" 2>&1 || { cat "$WORK/new/annot.log" >&2; fail "oncokb-cnv failed"; }

python3 - "$WORK/old/annot/cnv_oncokb_annot.tsv" "$WORK/new/oncokb/cnv_oncokb.tsv" \
  "$WORK/new/annot/cnv_oncokb_annot.tsv" <<'PY' || fail "oncokb-cnv output"
import csv, sys
def read(p):
    rows = list(csv.reader(open(p), delimiter="\t"))
    return rows[0], rows[1:]
old_h, old_rows = read(sys.argv[1])
in_h, in_rows = read(sys.argv[2])
new_h, new_rows = read(sys.argv[3])
assert old_h[:4] == ["AnnotSV_ID", "GeneA", "GeneB", "Sv_Type"], old_h
suffix = old_h[4:]
assert new_h == in_h + suffix, f"header: {new_h}"
assert all(len(r) == len(new_h) for r in new_rows), "rows are not padded to the header"
assert [r[:len(in_h)] for r in new_rows] == in_rows, "rows are not in input order"
new = {(r[0], r[1]): dict(zip(new_h, r)) for r in new_rows}
old = {(r[0], r[1]): dict(zip(old_h, r)) for r in old_rows}
want = {
    ("7_140719327_140924929_DUP_1", "BRAF"): ("Oncogenic", "Gain-of-function"),
    ("6_117288300_117425855_DUP_1", "ROS1"): ("Likely Oncogenic", "Likely Gain-of-function"),
    ("17_7661779_7687538_DEL_1", "TP53"): ("Likely Oncogenic", "Likely Loss-of-function"),
    ("17_7668402_7675493_DEL_1", "TP53"): ("Likely Oncogenic", "Likely Loss-of-function"),
    ("7_140730000_140800000_DUP_1", "BRAF"): ("Oncogenic", "Gain-of-function"),
    ("6_117300000_117400000_CNV_1", "ROS1"): ("", ""),
}
for key, (oncogenic, effect) in want.items():
    got = (new[key]["ONCOGENIC"], new[key]["MUTATION_EFFECT"])
    assert got == (oncogenic, effect), f"{key}: {got}, expected {(oncogenic, effect)}"
unclassified = new[("6_117300000_117400000_CNV_1", "ROS1")]
assert all(unclassified[h] == "" for h in suffix), "row without a class has OncoKB values"
print("    gene   class           today (SV query)                   now (CNA query)")
for key in new:
    o, n = old[key], new[key]
    print(f"    {key[1]:<6} {n['Copy_Number_Alteration'] or '-':<15} "
          f"{(o['ONCOGENIC'] + ' / ' + o['MUTATION_EFFECT']):<34} "
          f"{(n['ONCOGENIC'] + ' / ' + n['MUTATION_EFFECT']) if n['ANNOTATED'] else 'not annotated'}")
PY
echo "ok  oncokb-cnv: same OncoKB columns, input order, classes annotated"

GK=$(extract_task tasks/base_tasks/oncokb.yaml oncokb-cnv "$WORK/oncokb-cnv-gk.py" \
  "{{=sprig.osBase(inputs.parameters['oncokb-sv-in-f'])}}=cnv_oncokb.tsv" \
  "{{=sprig.osBase(inputs.parameters['oncokb-sv-out-f'])}}=cnv_oncokb_annot.tsv" \
  "{{inputs.parameters.regions-file}}=/ref/hg38_canonical_chrom_contig_regions.bed" \
  "{{inputs.parameters.tumor-type}}=." \
  "{{inputs.parameters.max-attempts}}=1" \
  "{{inputs.parameters.retry-delay-seconds}}=0")
mkdir -p "$WORK/gk"
docker run --rm --network none -e ONCOKB_API_TOKEN=unused -v "$WORK/new/oncokb:/mnt/input:ro" -v "$WORK/gk:/mnt/output" \
  -v "$WORK/oncokb-cnv-gk.py:/task.py:ro" --entrypoint python3 "$GK" /task.py >"$WORK/gk.log" 2>&1 \
  || { cat "$WORK/gk.log" >&2; fail "Genome kit skip failed"; }
python3 - "$WORK/new/annot/cnv_oncokb_annot.tsv" "$WORK/gk/cnv_oncokb_annot.tsv" "$WORK/new/oncokb/cnv_oncokb.tsv" <<'PY' \
  || fail "Genome kit skip output"
import csv, sys
a = list(csv.reader(open(sys.argv[1]), delimiter="\t"))
g = list(csv.reader(open(sys.argv[2]), delimiter="\t"))
i = list(csv.reader(open(sys.argv[3]), delimiter="\t"))
assert g[0] == a[0], "Genome kit header differs"
n = len(i[0])
assert [r[:n] for r in g[1:]] == i[1:] and all(set(r[n:]) == {""} for r in g[1:])
PY
echo "ok  Genome kit: empty OncoKB columns, no network"

# --- 3. DuckDB ---------------------------------------------------------------
build_db() { # build_db DIR ONCOKB_TSV
  local dir=$1 image
  mkdir -p "$dir/db"
  image=$(extract_task tasks/base_tasks/duckdb.yaml duckdb-cnv "$dir/db.sql" \
    "{{=sprig.osBase(inputs.parameters['output-duckdb'])}}=results.duckdb" \
    "{{inputs.parameters.table-suffix}}=")
  docker run --rm -v "$2:/mnt/input/oncokb/oncokb_cnv.tsv:ro" \
    -v "$WORK/old/filtered/annotsv_cnv_filtered.tsv:/mnt/input/annotsv/annotsv_cnv.tsv:ro" \
    -v "$dir/db:/mnt/output" -v "$dir/db.sql:/db.sql:ro" "$image" duckdb -init /db.sql \
    >"$dir/db.log" 2>&1 || { cat "$dir/db.log" >&2; fail "duckdb-cnv failed"; }
  # duckdb -init keeps going after a failed statement, so check the log too.
  if grep -q "Error" "$dir/db.log"; then cat "$dir/db.log" >&2; fail "duckdb-cnv failed"; fi
  echo "$image"
}
DUCKDB_IMAGE=$(build_db "$WORK/old" "$WORK/old/annot/cnv_oncokb_annot.tsv")
build_db "$WORK/new" "$WORK/new/annot/cnv_oncokb_annot.tsv" >/dev/null

q() { docker run --rm -v "$1/db:/mnt/output" "$DUCKDB_IMAGE" duckdb -csv /mnt/output/results.duckdb "$2"; }
ONCOKB_VALUES="ONCOGENIC|MUTATION_EFFECT|MUTATION_EFFECT_CITATIONS|GENE_IN_ONCOKB|VARIANT_IN_ONCOKB|ANNOTATED|LEVEL_.*|HIGHEST_.*|TX_CITATIONS|DX_CITATIONS|PX_CITATIONS"
for view in combined_cnv advanced_combined_cnv; do
  schema="SELECT column_name, data_type FROM information_schema.columns WHERE table_name = '$view' ORDER BY ordinal_position"
  diff <(q "$WORK/old" "$schema") <(q "$WORK/new" "$schema") || fail "$view columns or types differ"
  cols=$(q "$WORK/old" "$schema" | tail -n +2 | cut -d, -f1 | grep -Evx "$ONCOKB_VALUES" | sed 's/.*/"&"/' | paste -sd,)
  rows="SELECT $cols FROM $view ORDER BY ALL"
  diff <(q "$WORK/old" "$rows") <(q "$WORK/new" "$rows") || fail "$view rows differ outside the OncoKB values"
  count="SELECT count(*) FROM $view"
  [ "$(q "$WORK/new" "$count")" = "$(q "$WORK/old" "$count")" ] || fail "$view row count differs"
done
echo "ok  duckdb-cnv: same columns, types and rows apart from the OncoKB values"

echo "PASS: OncoKB CNV copy-number classes"
