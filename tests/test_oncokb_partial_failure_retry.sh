#!/usr/bin/env bash
# When a batch request fails, upstream oncokb-annotator asks OncoKB once per
# query. Each query that still fails prints "Error on annotating the url ..."
# and leaves its row unannotated, and the annotator still exits 0 (#42).
#
# This runs the real task scripts from tasks/base_tasks/oncokb.yaml and
# oncokb-fusion-annotation.yaml, in the image those tasks pin, with a fake
# annotator mounted over the upstream script. The fake writes a valid output
# and prints that line on the attempts listed in FAIL_ATTEMPTS. Checks:
#   - a failed call on attempt 1 makes the task retry, and attempt 2 succeeds;
#   - failed calls on every attempt make the task exit non-zero;
#   - a clean run passes on the first attempt;
#   - the token is not printed.
#
# Needs docker and python3 with PyYAML.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

WORK=$(mktemp -d)
trap '[ -n "${KEEP_WORK:-}" ] || rm -rf "$WORK"' EXIT
TOKEN=secret-token-for-test

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

# The fake annotator: takes upstream's -i/-o or our wrappers' --input/--output,
# copies the input with the OncoKB columns appended, and on the attempts in
# FAIL_ATTEMPTS prints upstream's line for a query it gave up on.
cat >"$WORK/fake_annotator.py" <<'PY'
import os, sys
sys.path.insert(0, "/app/oncokb")
from AnnotatorCore import get_oncokb_annotation_column_headers
argv = sys.argv[1:]
def opt(*names):
    for n in names:
        if n in argv:
            return argv[argv.index(n) + 1]
src, dst = opt("-i", "--input"), opt("-o", "--output")
counter = "/tmp/fake-annotator-attempt"
attempt = int(open(counter).read()) + 1 if os.path.exists(counter) else 1
open(counter, "w").write(str(attempt))
fail_on = {int(a) for a in os.environ.get("FAIL_ATTEMPTS", "").split(",") if a}
lines = open(src).read().splitlines()
suffix = get_oncokb_annotation_column_headers(False, False)
with open(dst, "w") as f:
    f.write("\t".join(lines[0].split("\t") + suffix) + "\n")
    for line in lines[1:]:
        f.write("\t".join(line.split("\t") + ["True"] + [""] * (len(suffix) - 1)) + "\n")
if attempt in fail_on:
    print("Error on annotating the url https://www.oncokb.org/api/v1/annotate/x?hugoSymbol=BRAF")
sys.exit(int(os.environ.get("EXIT_CODE", "0")))
PY

# Writes a task's script with its template expressions filled in.
extract_task() {
  python3 - "$1" "$2" "$3" <<'PY'
import re, sys, yaml
tasks, name, out = sys.argv[1:]
doc = yaml.safe_load(open(tasks))
task = next(t for t in doc["spec"]["templates"] if t["name"] == name)
values = {"max-attempts": "2", "retry-delay-seconds": "0", "tumor-type": "MEL",
          "regions-file": "/sample/panel.bed", "ref-version": "hg38"}
src = task["script"]["source"]
src = re.sub(r"\{\{=sprig\.osBase\(inputs\.parameters\['([^']+)'\]\)\}\}",
             lambda m: {"annovar-maf-f": "in.maf", "annovar-oncokb-maf-f": "out.maf",
                        "oncokb-sv-in-f": "in.tsv", "oncokb-sv-out-f": "out.tsv",
                        "oncokb-fusions-in-f": "in.tsv", "oncokb-fusions-out-f": "out.tsv"}[m.group(1)], src)
src = re.sub(r"\{\{inputs\.parameters\.([a-z-]+)\}\}", lambda m: values[m.group(1)], src)
src = src.replace("sleep 10", "sleep 0")
if "{{" in src:
    sys.exit(f"{name}: unhandled template expression: {re.findall(r'{{[^}]*}}', src)}")
open(out, "w").write(src)
print(task["script"]["image"], task["script"]["command"][0])
PY
}

printf 'Hugo_Symbol\tChromosome\tStart_Position\nBRAF\t7\t140753336\n' >"$WORK/in.maf"
printf 'GeneA\tGeneB\tSv_Type\tCopy_Number_Alteration\tHugo_Symbol\tTumor_Sample_Barcode\nBRAF\tBRAF\tDUP\tAmplification\tBRAF\tS1\n' >"$WORK/in.tsv"

# task, upstream script the fake replaces
cases="oncokb-snv:/app/oncokb/filtered_maf_annotator.py
oncokb-snv-on-node:/app/oncokb/filtered_maf_annotator.py
oncokb-sv:/app/oncokb/filtered_sv_annotator.py
oncokb-sv-on-node:/app/oncokb/filtered_sv_annotator.py
oncokb-fusions:/app/oncokb/FusionAnnotator.py
oncokb-cnv:/app/oncokb/CnaAnnotator.py"

run_case() {  # task fake-path fail-attempts exit-code -> prints the status
  local task=$1 fake=$2 fails=$3 code=$4 dir="$WORK/run-$1-$3-$4" image cmd
  mkdir -p "$dir/in" "$dir/out"
  read -r image cmd < <(extract_task tasks/base_tasks/oncokb.yaml "$task" "$dir/script")
  cp "$WORK/in.maf" "$WORK/in.tsv" "$dir/in/"
  set +e
  docker run --rm -e ONCOKB_API_TOKEN="$TOKEN" -e FAIL_ATTEMPTS="$fails" -e EXIT_CODE="$code" \
    -v "$dir/in:/mnt/input:ro" -v "$dir/out:/mnt/output" -v "$dir/script:/script:ro" \
    -v "$WORK/fake_annotator.py:$fake:ro" -w /app/oncokb "$image" "$cmd" /script >"$dir/log" 2>&1
  echo $?
  set -e
}

while IFS=: read -r task fake; do
  status=$(run_case "$task" "$fake" "" 0)
  [ "$status" = 0 ] || fail "$task: a clean run exited $status"
  status=$(run_case "$task" "$fake" 1 0)
  log="$WORK/run-$task-1-0/log"
  [ "$status" = 0 ] || { cat "$log" >&2; fail "$task: a failed call on attempt 1 should be retried and pass on attempt 2 (exit $status)"; }
  grep -q "Attempt 2 of 2" "$log" || fail "$task: no second attempt after a failed call"
  status=$(run_case "$task" "$fake" 1,2 0)
  log="$WORK/run-$task-1,2-0/log"
  [ "$status" != 0 ] || { cat "$log" >&2; fail "$task: failed calls on every attempt should fail the task"; }
  status=$(run_case "$task" "$fake" "" 3)
  [ "$status" != 0 ] || fail "$task: an annotator exit status 3 should fail the task"
  if grep -rq "$TOKEN" "$WORK"/run-"$task"-*/log; then fail "$task: the token appears in the log"; fi
  echo "ok  $task: retries on failed OncoKB calls, fails when they persist"
done <<<"$cases"

# oncokb-fusion-annotation runs FusionAnnotator from its own script.
task=oncokb-fusion-annotation
dir="$WORK/fusion-annotation"
mkdir -p "$dir/in" "$dir/out"
python3 - tasks/base_tasks/oncokb-fusion-annotation.yaml "$dir/script" <<'PY' >"$dir/image"
import re, sys, yaml
doc = yaml.safe_load(open(sys.argv[1]))
task = doc["spec"]["templates"][0]
src = task["script"]["source"]
src = re.sub(r"\{\{[^}]*\}\}", lambda m: {"{{inputs.parameters.tumor-type}}": "MEL"}.get(m.group(0), "X"), src)
src = src.replace("time.sleep(5)", "time.sleep(0)")
open(sys.argv[2], "w").write(src)
print(task["script"]["image"])
PY
# Only the annotation step is exercised: cut the script after it and feed it
# the input it builds in step 1.
python3 - "$dir/script" <<'PY'
import sys
p = sys.argv[1]; s = open(p).read()
start = s.index("# Step 2: Annotate via OncoKB FusionAnnotator")
end = s.index("# Step 3:")
head = "import os, subprocess, sys, time\nTUMOR_TYPE = 'MEL'\n"
open(p, "w").write(head + s[start:end])
PY
printf 'Fusion\tHugo_Symbol\nBRAF-KIAA1549\tBRAF\n' >"$dir/in/fusions_input.tsv"
for fails in "" 1 1,2,3; do
  set +e
  docker run --rm -e ONCOKB_API_TOKEN="$TOKEN" -e FAIL_ATTEMPTS="$fails" \
    -v "$dir/in/fusions_input.tsv:/tmp/fusions_input.tsv:ro" -v "$dir/script:/script.py:ro" \
    -v "$WORK/fake_annotator.py:/app/oncokb/FusionAnnotator.py:ro" -w /app/oncokb \
    "$(cat "$dir/image")" python3 /script.py >"$dir/log-$fails" 2>&1
  status=$?
  set -e
  case "$fails" in
    "") [ $status = 0 ] || fail "$task: a clean run exited $status" ;;
    1) [ $status = 0 ] && grep -q "Attempt 2 of 3" "$dir/log-$fails" \
         || fail "$task: a failed call on attempt 1 should be retried and pass" ;;
    *) [ $status != 0 ] || fail "$task: failed calls on every attempt should fail the task" ;;
  esac
done
if grep -q "$TOKEN" "$dir"/log-*; then fail "$task: the token appears in the log"; fi
echo "ok  $task: retries on failed OncoKB calls, fails when they persist"

echo "PASS: OncoKB tasks retry when OncoKB calls fail partway"
