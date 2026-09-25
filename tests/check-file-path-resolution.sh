#!/bin/sh
set -eu

REPO_DIR=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
CASE_DIR=$(mktemp -d)
trap 'rm -rf -- "$CASE_DIR"' EXIT

RUN_DIR="$CASE_DIR/run"
SAMPLE_UUID="sample-uuid-1234"
OUTPUT_DIR="$CASE_DIR/outputs"
mkdir -p "$RUN_DIR/$SAMPLE_UUID" "$OUTPUT_DIR"

touch "$RUN_DIR/$SAMPLE_UUID/nested_R1.fastq.gz"
touch "$RUN_DIR/$SAMPLE_UUID/nested_R2.fastq.gz"
touch "$RUN_DIR/fallback_regions.bed"

SCRIPT_SOURCE=$(python3 - "$REPO_DIR/tasks/check-file-path.yaml" <<'PYEXTRACT'
import sys
import yaml
with open(sys.argv[1]) as handle:
    data = yaml.safe_load(handle)
print(next(template['script']['source'] for template in data['spec']['templates'] if template['name'] == 'check-file-path'))
PYEXTRACT
)

export TEST_RUN_DIR="$RUN_DIR"
export TEST_SAMPLE_UUID="$SAMPLE_UUID"
export TEST_OUTPUT_DIR="$OUTPUT_DIR"
export TEST_FILE_1="nested_R1.fastq.gz"
export TEST_FILE_2="nested_R2.fastq.gz"
export TEST_FILE_3="fallback_regions.bed"

RENDERED_SOURCE=$(printf '%s\n' "$SCRIPT_SOURCE" | perl -0777 -pe '
  s/\Q{{inputs.parameters.input-rel-dir}}\E/$ENV{TEST_RUN_DIR}/g;
  s/\Q{{inputs.parameters.sample-uuid}}\E/$ENV{TEST_SAMPLE_UUID}/g;
  s/\Q{{inputs.parameters.input-file-1}}\E/$ENV{TEST_FILE_1}/g;
  s/\Q{{inputs.parameters.input-file-2}}\E/$ENV{TEST_FILE_2}/g;
  s/\Q{{inputs.parameters.input-file-3}}\E/$ENV{TEST_FILE_3}/g;
  s#\Q/tmp/output-file-1.txt\E#$ENV{TEST_OUTPUT_DIR}/output-file-1.txt#g;
  s#\Q/tmp/output-file-2.txt\E#$ENV{TEST_OUTPUT_DIR}/output-file-2.txt#g;
  s#\Q/tmp/output-file-3.txt\E#$ENV{TEST_OUTPUT_DIR}/output-file-3.txt#g;
  s#\Q/tmp/input-rel-dir-1.txt\E#$ENV{TEST_OUTPUT_DIR}/input-rel-dir-1.txt#g;
  s#\Q/tmp/input-rel-dir-2.txt\E#$ENV{TEST_OUTPUT_DIR}/input-rel-dir-2.txt#g;
  s#\Q/tmp/input-rel-dir-3.txt\E#$ENV{TEST_OUTPUT_DIR}/input-rel-dir-3.txt#g;
')

printf '%s\n' "$RENDERED_SOURCE" | python3

assert_file() {
  OUTPUT_FILE="$1"
  EXPECTED="$2"
  ACTUAL=$(cat "$OUTPUT_DIR/$OUTPUT_FILE")
  if [ "$ACTUAL" != "$EXPECTED" ]; then
    echo "$OUTPUT_FILE: expected '$EXPECTED', got '$ACTUAL'" >&2
    exit 1
  fi
}

assert_file output-file-1.txt nested_R1.fastq.gz
assert_file output-file-2.txt nested_R2.fastq.gz
assert_file output-file-3.txt fallback_regions.bed
assert_file input-rel-dir-1.txt "$RUN_DIR/$SAMPLE_UUID/"
assert_file input-rel-dir-2.txt "$RUN_DIR/$SAMPLE_UUID/"
assert_file input-rel-dir-3.txt "$RUN_DIR/"

echo "check-file-path sample-directory preference and run-directory fallback: PASS"
