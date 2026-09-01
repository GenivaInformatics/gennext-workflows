#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

BED_VALIDATOR=tasks/base_tasks/bed_validator.yaml
MUTECT2_TEMPLATE=tasks/base_tasks/mutect2.yaml
SHARDED_TEMPLATE=tasks/base_tasks/mutect2-sharded.yaml
WORKFLOWS=(
  workflows/somatic/v2/somatic-dna-bam.yaml
  workflows/somatic/v2/somatic-dna-fastq.yaml
)

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

task_block() {
  local task_name=$1
  local workflow=$2
  awk -v needle="          - name: $task_name" '
    $0 == needle { printing = 1; first = 1 }
    printing && !first && /^          - name: / { exit }
    printing { print; first = 0 }
  ' "$workflow"
}

echo "Checking manifest contract"

[[ $(grep -c '^          # BEGIN PADDED_BED_PYTHON$' "$BED_VALIDATOR") == 1 ]] ||
  fail "production Python start marker is missing or duplicated"
[[ $(grep -c '^          # END PADDED_BED_PYTHON$' "$BED_VALIDATOR") == 1 ]] ||
  fail "production Python end marker is missing or duplicated"

parameter_block=$(grep -A1 'name: generate-padded-artifacts' "$BED_VALIDATOR")
grep -q 'value: "false"' <<< "$parameter_block" ||
  fail "padded artifact generation must default to false"
grep -q 'if ! python3 -' "$BED_VALIDATOR" ||
  fail "production generator is not guarded by a hard-fail condition"
grep -q 'Refusing to use unpadded intervals' "$BED_VALIDATOR" ||
  fail "padding failure does not explain that unpadded fallback is forbidden"
if grep -q 'cp .*validAnalysisBedName' "$BED_VALIDATOR"; then
  fail "analysis BED silently falls back to the unpadded target BED"
fi

yaml_opt_ins=$(grep -R --include='*.yaml' -h 'name: generate-padded-artifacts' tasks workflows | wc -l)
[[ $yaml_opt_ins == 4 ]] ||
  fail "expected one parameter definition and three opt-ins (two v2 + repair), found $yaml_opt_ins"

for workflow in "${WORKFLOWS[@]}"; do
  validator=$(task_block validate-bed-file "$workflow")
  grep -A1 -q 'name: generate-padded-artifacts' <<< "$validator" ||
    fail "$workflow primary validator does not opt in"
  grep -A1 'name: generate-padded-artifacts' <<< "$validator" | grep -q 'value: "true"' ||
    fail "$workflow primary validator opt-in is not true"
  grep -A1 'name: padding' <<< "$validator" | grep -q 'value: "100"' ||
    fail "$workflow primary validator padding is not 100bp"
  grep -A1 'name: reference-fai-file' <<< "$validator" | grep -q '\.fai"' ||
    fail "$workflow primary validator does not provide the reference FASTA index"

  msi_validator=$(task_block validate-msi-bed-file "$workflow")
  if grep -q 'generate-padded-artifacts' <<< "$msi_validator"; then
    fail "$workflow MSI validator must retain the default, non-generating behavior"
  fi

  mutect2=$(task_block mutect2 "$workflow")
  grep -A1 'name: regions-file' <<< "$mutect2" | grep -q 'validated-bed-analysis-file' ||
    fail "$workflow Mutect2 does not consume the analysis BED"

  bcftools=$(task_block bcftools-filter "$workflow")
  grep -q 'mutect2.outputs.parameters.output-vcf' <<< "$bcftools" ||
    fail "$workflow bcftools-filter is not the post-Mutect2 step"
  grep -A1 'name: regions-file' <<< "$bcftools" | grep -q 'validated-bed-analysis-file' ||
    fail "$workflow post-Mutect2 bcftools does not consume the analysis BED"

  tmb=$(task_block collect-tmb-data "$workflow")
  grep -A1 'name: target-bed-file' <<< "$tmb" | grep -q 'validated-bed-file' ||
    fail "$workflow TMB does not retain the original target BED"
  if grep -qE 'validated-bed-(analysis|igv)-file' <<< "$tmb"; then
    fail "$workflow TMB incorrectly consumes a padded or presentation BED"
  fi
done

analysis_consumers=$(grep -R --include='*.yaml' -h 'validated-bed-analysis-file' workflows | wc -l)
[[ $analysis_consumers == 8 ]] ||
  fail "expected exactly eight analysis BED consumers (bam:2 + fastq:2 + repair:4), found $analysis_consumers"
# The IGV BED may only appear in validation contexts (e.g. padding-repair's
# validate-candidate), never as input to an analysis tool like Mutect2/bcftools/TMB.
igv_non_validation=$(grep -R --include='*.yaml' -h 'validated-bed-igv-file' workflows |
  { grep -v 'validate-candidate\|validate-bed\|igv-bed' || true; } | wc -l)
[[ $igv_non_validation == 0 ]] ||
  fail "the presentation-only IGV BED is passed to non-validation analysis tasks ($igv_non_validation occurrences)"

direct_padding=$(grep -c -- '--interval-padding 100' "$MUTECT2_TEMPLATE")
[[ $direct_padding == 2 ]] ||
  fail "direct Mutect2 paths must retain two 100bp padding flags"
if grep -q -- '--interval-padding' "$SHARDED_TEMPLATE"; then
  fail "v2 sharded Mutect2 would double-pad the already padded analysis BED"
fi

echo "Executing the exact embedded production generator"

TEST_DIR=$(mktemp -d)
cleanup() {
  rm -rf -- "$TEST_DIR"
}
trap cleanup EXIT

awk '
  /^          # BEGIN PADDED_BED_PYTHON$/ { capture = 1; next }
  /^          # END PADDED_BED_PYTHON$/ { capture = 0; found_end = 1; next }
  capture { sub(/^          /, ""); print }
  END { if (!found_end) exit 1 }
' "$BED_VALIDATOR" > "$TEST_DIR/padding.py"

cat > "$TEST_DIR/reference.fai" <<'EOF'
chr1	20000	0	0	0
chr2	5000	0	0	0
EOF

cat > "$TEST_DIR/test.bed" <<'EOF'
chr2	4900	5000
chr1	10000	10100
chr1	10	60
chr1	1200	1250
chr1	10300	10400
chr1	1000	1050
EOF

python3 "$TEST_DIR/padding.py" \
  "$TEST_DIR/test.bed" \
  "$TEST_DIR/reference.fai" \
  "$TEST_DIR/analysis.bed" \
  "$TEST_DIR/igv.bed" \
  100

cat > "$TEST_DIR/expected-analysis.bed" <<'EOF'
chr1	0	160
chr1	900	1350
chr1	9900	10500
chr2	4800	5000
EOF

diff -u "$TEST_DIR/expected-analysis.bed" "$TEST_DIR/analysis.bed" ||
  fail "analysis BED does not clip, reference-sort, and merge overlaps/book-ended intervals"
awk 'NF != 3 { exit 1 }' "$TEST_DIR/analysis.bed" ||
  fail "analysis BED is not strict BED3"
if grep -qE '^(track|#)' "$TEST_DIR/analysis.bed"; then
  fail "analysis BED contains a presentation header"
fi

python3 - "$TEST_DIR/igv.bed" <<'PYTHON'
import sys

path = sys.argv[1]
expected_colors = {
    "targeted": "128,128,128",
    "padding": "211,211,211",
}
features = {name: [] for name in expected_colors}

with open(path, encoding="utf-8") as handle:
    header = handle.readline().rstrip("\n")
    if header != 'track name="Targeted and Padding Regions" itemRgb="On"':
        raise SystemExit(f"unexpected IGV track header: {header!r}")

    for line_number, line in enumerate(handle, start=2):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise SystemExit(
                f"IGV line {line_number} has {len(fields)} columns instead of BED9"
            )

        chrom, start_text, end_text, feature_type = fields[:4]
        if feature_type not in expected_colors:
            raise SystemExit(f"unknown IGV feature type: {feature_type!r}")
        if fields[8] != expected_colors[feature_type]:
            raise SystemExit(
                f"{feature_type} uses {fields[8]!r}, expected "
                f"{expected_colors[feature_type]!r}"
            )

        start, end = int(start_text), int(end_text)
        if start >= end:
            raise SystemExit(f"invalid IGV interval: {chrom}:{start}-{end}")
        features[feature_type].append((chrom, start, end))

for feature_type, intervals in features.items():
    if not intervals:
        raise SystemExit(f"IGV BED contains no {feature_type} features")

for target_chrom, target_start, target_end in features["targeted"]:
    for padding_chrom, padding_start, padding_end in features["padding"]:
        overlaps = (
            target_chrom == padding_chrom
            and target_start < padding_end
            and padding_start < target_end
        )
        if overlaps:
            raise SystemExit(
                "targeted and padding-only features overlap: "
                f"{target_chrom}:{target_start}-{target_end} and "
                f"{padding_chrom}:{padding_start}-{padding_end}"
            )
PYTHON

cat > "$TEST_DIR/bad-contig.bed" <<'EOF'
chrMissing	10	20
EOF
if python3 "$TEST_DIR/padding.py" \
  "$TEST_DIR/bad-contig.bed" \
  "$TEST_DIR/reference.fai" \
  "$TEST_DIR/bad-analysis.bed" \
  "$TEST_DIR/bad-igv.bed" \
  100 >/dev/null 2>&1; then
  fail "generator accepted a BED contig absent from authoritative GTF bounds"
fi

touch "$TEST_DIR/empty.fai"
if python3 "$TEST_DIR/padding.py" \
  "$TEST_DIR/test.bed" \
  "$TEST_DIR/empty.fai" \
  "$TEST_DIR/no-bounds-analysis.bed" \
  "$TEST_DIR/no-bounds-igv.bed" \
  100 >/dev/null 2>&1; then
  fail "generator accepted a FASTA index with no authoritative contig bounds"
fi

echo "BED padding regression test passed"
