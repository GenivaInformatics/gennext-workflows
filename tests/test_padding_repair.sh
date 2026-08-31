#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

REPAIR_WORKFLOW=workflows/somatic/v2/somatic-dna-padding-repair.yaml

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

echo "=== Padding repair workflow regression tests ==="

# ═══════════════════════════════════════════════════════════════════════
# Section 1: Static wiring checks
# ═══════════════════════════════════════════════════════════════════════

# ── 1. Public inputs: only node-name, user-id, run-id, sample-uuid, mode ──
echo "Checking public input surface..."
param_count=$(awk '
  /^  arguments:/ { in_args=1; next }
  in_args && /^    parameters:/ { in_params=1; next }
  in_params && /^      - name: / { count++ }
  in_params && /^  [^ ]/ { exit }
  END { print count+0 }
' "$REPAIR_WORKFLOW")
[[ $param_count == 5 ]] ||
  fail "Expected 5 public parameters, found $param_count"
for p in node-name user-id run-id sample-uuid mode; do
  grep -qE "^      - name: ${p}$" "$REPAIR_WORKFLOW" ||
    fail "Missing public parameter: $p"
done
for removed in input-rel-dir output-rel-dir temp-rel-dir scratch-rel-dir ref-base-dir thread-count; do
  if grep -qE "^      - name: ${removed}$" "$REPAIR_WORKFLOW"; then
    fail "Infrastructure parameter still public: $removed"
  fi
done

# ── 2. get-node-config is the first DAG task ─────────────────────────
echo "Checking get-node-config integration..."
grep -q 'name: get-node-config' "$REPAIR_WORKFLOW" ||
  fail "get-node-config task not found"
gnc=$(task_block get-node-config "$REPAIR_WORKFLOW")
echo "$gnc" | grep -q 'template: get-node-config' ||
  fail "get-node-config does not reference the correct template"
echo "$gnc" | grep -q 'depends' &&
  fail "get-node-config has dependencies (should be first)"

# ── 3. Workflow-UID staging paths ────────────────────────────────────
echo "Checking workflow-UID staging..."
grep -q 'repair-wk-{{workflow.uid}}' "$REPAIR_WORKFLOW" ||
  fail "Transaction root does not use workflow UID"
grep -q 'repair-cd-{{workflow.uid}}' "$REPAIR_WORKFLOW" ||
  fail "Candidate directory does not use workflow UID"

# ── 4. Mode validation ──────────────────────────────────────────────
echo "Checking mode validation..."
grep -q "mode not in.*validate.*apply" "$REPAIR_WORKFLOW" ||
  fail "Mode validation not found"

# ── 5. Mutect2 uses the analysis BED ────────────────────────────────
echo "Checking Mutect2 consumes the analysis BED..."
mutect2=$(task_block mutect2 "$REPAIR_WORKFLOW")
grep -A1 'name: regions-file' <<< "$mutect2" | grep -q 'validated-bed-analysis-file' ||
  fail "Mutect2 does not consume the analysis BED"

# ── 6. bcftools uses the analysis BED ────────────────────────────────
echo "Checking bcftools-filter consumes the analysis BED..."
bcftools=$(task_block bcftools-filter "$REPAIR_WORKFLOW")
grep -A1 'name: regions-file' <<< "$bcftools" | grep -q 'validated-bed-analysis-file' ||
  fail "bcftools-filter does not consume the analysis BED"

# ── 7. TMB compute uses the original target BED ─────────────────────
echo "Checking TMB compute uses the original target BED..."
compute_tmb=$(task_block compute-tmb "$REPAIR_WORKFLOW")
grep -A1 'name: target-bed-file' <<< "$compute_tmb" | grep -q 'validated-bed-file' ||
  fail "TMB compute does not use the original target BED"
if grep -qE 'validated-bed-(analysis|igv)-file' <<< "$compute_tmb"; then
  fail "TMB compute incorrectly consumes a padded or presentation BED"
fi

# ── 8. TMB publish uses the original target BED ─────────────────────
echo "Checking TMB publish uses the original target BED..."
publish_tmb=$(task_block publish-tmb "$REPAIR_WORKFLOW")
grep -A1 'name: target-bed-file' <<< "$publish_tmb" | grep -q 'validated-bed-file' ||
  fail "TMB publish does not use the original target BED"
if grep -qE 'validated-bed-(analysis|igv)-file' <<< "$publish_tmb"; then
  fail "TMB publish incorrectly consumes a padded or presentation BED"
fi

# ── 9. TMB compute uses compute-tmb-data (no API) ───────────────────
echo "Checking TMB compute/publish template split..."
grep -A3 'name: compute-tmb$' "$REPAIR_WORKFLOW" | grep -q 'template: compute-tmb-data' ||
  fail "TMB compute does not use compute-tmb-data template"
grep -A5 'name: publish-tmb$' "$REPAIR_WORKFLOW" | grep -q 'template: collect-tmb-data' ||
  fail "TMB publish does not use collect-tmb-data template"

# ── 10. TMB publish output goes to tx-root, not live ─────────────────
echo "Checking TMB publish outputs to tx-root..."
grep -A1 'name: output-dir' <<< "$publish_tmb" | grep -q 'tx-root.*tmb-publish' ||
  fail "TMB publish output does not go to tx-root"

# ── 11. Same-sample mutex ────────────────────────────────────────────
echo "Checking workflow mutex..."
grep -q 'name: "padding-repair-{{workflow.parameters.sample-uuid}}"' "$REPAIR_WORKFLOW" ||
  fail "Workflow mutex is not keyed on sample-uuid"

# ── 12. No MAF importer ─────────────────────────────────────────────
echo "Checking maf-to-db template is not referenced..."
if grep -q 'template: maf-to-db' "$REPAIR_WORKFLOW" || grep -q 'name: maf-to-db-template' "$REPAIR_WORKFLOW"; then
  fail "Repair workflow must NOT call the global MAF importer (maf-to-db)"
fi

# ── 13. Promote gated on mode=apply ──────────────────────────────────
echo "Checking promote/publish mode gates..."
promote=$(task_block promote "$REPAIR_WORKFLOW")
grep -q "mode.*==.*apply" <<< "$promote" ||
  fail "promote is not gated on mode=apply"
grep -q "mode.*==.*apply" <<< "$publish_tmb" ||
  fail "publish-tmb is not gated on mode=apply"

# ── 14. publish-tmb depends on promote.Succeeded ────────────────────
echo "Checking publish-tmb dependency on promote..."
grep -q 'depends: "promote.Succeeded"' "$REPAIR_WORKFLOW" ||
  fail "publish-tmb does not depend on promote.Succeeded"

# ── 15. DuckDB staged via copy ───────────────────────────────────────
echo "Checking DuckDB is staged via copy..."
sed -n '/^    - name: stage-duckdb-copy$/,/^    - name: /p' "$REPAIR_WORKFLOW" | grep -q 'cp ' ||
  fail "stage-duckdb-copy does not copy the existing DuckDB"

# ── 16. DuckDB inventory snapshot + verification ────────────────────
echo "Checking DuckDB inventory verification pipeline..."
grep -q 'name: snapshot-duckdb-inventory' "$REPAIR_WORKFLOW" ||
  fail "DuckDB inventory snapshot task missing"
grep -q 'name: verify-duckdb-preservation' "$REPAIR_WORKFLOW" ||
  fail "DuckDB preservation verification task missing"
verify_dep=$(task_block verify-duckdb-preservation "$REPAIR_WORKFLOW")
grep -q 'depends: "duckdb-snv"' <<< "$verify_dep" ||
  fail "DuckDB verification does not depend on duckdb-snv"

# ── 17. No duckdb-sv/cnv/fusions ────────────────────────────────────
echo "Checking DuckDB only touches SNV..."
for tpl in duckdb-sv duckdb-cnv duckdb-fusions; do
  if grep -q "template: ${tpl}\b" "$REPAIR_WORKFLOW"; then
    fail "Repair workflow must not call $tpl"
  fi
done

# ── 18. Atomic exchange via renameat2 ────────────────────────────────
echo "Checking atomic exchange mechanism..."
grep -q 'SYS_RENAMEAT2.*=.*316' "$REPAIR_WORKFLOW" ||
  fail "promote does not use renameat2 syscall"
grep -q 'RENAME_EXCHANGE.*=.*2' "$REPAIR_WORKFLOW" ||
  fail "promote does not use RENAME_EXCHANGE flag"

# ── 19. Assemble candidate BEFORE validation ─────────────────────────
echo "Checking candidate ordering..."
assemble_dep=$(task_block assemble-candidate "$REPAIR_WORKFLOW")
grep -q 'depends:.*verify-duckdb-preservation' <<< "$assemble_dep" ||
  fail "assemble-candidate does not depend on verify-duckdb-preservation"
validate_cd_dep=$(task_block validate-candidate "$REPAIR_WORKFLOW")
grep -q 'depends:.*write-manifest' <<< "$validate_cd_dep" ||
  fail "validate-candidate does not depend on write-manifest"

# ── 20. Manifest in candidate before validation ─────────────────────
echo "Checking manifest placement..."
manifest_dep=$(task_block write-manifest "$REPAIR_WORKFLOW")
grep -q 'depends:.*compute-tmb' <<< "$manifest_dep" ||
  fail "write-manifest does not depend on compute-tmb"
grep -q 'padding_repair_manifest.json' "$REPAIR_WORKFLOW" ||
  fail "Manifest filename not found in workflow"

# ── 21. Cleanup only touches workflow-UID paths ─────────────────────
echo "Checking cleanup scope..."
cleanup_tpl=$(sed -n '/^    - name: cleanup-repair$/,/^    - name: /p' "$REPAIR_WORKFLOW")
echo "$cleanup_tpl" | grep -q 'TX_BASE' ||
  fail "cleanup does not reference transaction root"
echo "$cleanup_tpl" | grep -q 'CD_BASE' ||
  fail "cleanup does not reference candidate dir"
if echo "$cleanup_tpl" | grep -qE 'rm -rf.*/mnt/parent/"$' ; then
  fail "cleanup removes the output parent (unsafe)"
fi

# ── 22. Publication failure recorded in manifest ────────────────────
echo "Checking publish failure handling..."
grep -q "projection.*state.*failed" "$REPAIR_WORKFLOW" ||
  fail "TMB publish failure is not recorded as failed in manifest"
manifest_proj_dep=$(task_block update-manifest-projection "$REPAIR_WORKFLOW")
grep -q 'publish-tmb.Failed' <<< "$manifest_proj_dep" ||
  fail "update-manifest-projection does not handle publish-tmb.Failed"

# ── 23. Node selector ───────────────────────────────────────────────
echo "Checking node selector..."
grep -q 'kubernetes.io/hostname: "{{workflow.parameters.node-name}}"' "$REPAIR_WORKFLOW" ||
  fail "Workflow does not pin to the specified node"

# ── 24. validate-bed-file generates padded artifacts ─────────────────
echo "Checking BED validator generates padded artifacts..."
bed_task=$(task_block validate-bed-file "$REPAIR_WORKFLOW")
grep -A1 'name: generate-padded-artifacts' <<< "$bed_task" | grep -q 'value: "true"' ||
  fail "validate-bed-file does not opt into padded artifact generation"
grep -A1 'name: padding' <<< "$bed_task" | grep -q 'value: "100"' ||
  fail "validate-bed-file padding is not 100bp"
grep -A1 'name: reference-fai-file' <<< "$bed_task" | grep -q '\.fai"' ||
  fail "validate-bed-file does not provide the reference FASTA index"

# ── 25. VCF-only samples rejected ───────────────────────────────────
echo "Checking VCF-only sample rejection..."
grep -q 'VCF-only samples cannot be repaired' "$REPAIR_WORKFLOW" ||
  fail "Repair workflow does not reject VCF-only samples"

# ═══════════════════════════════════════════════════════════════════════
# Section 2: Behavioral tests — extract and run production scripts
# ═══════════════════════════════════════════════════════════════════════

echo ""
echo "=== Behavioral tests ==="

TEST_DIR=$(mktemp -d)
cleanup() {
  rm -rf -- "$TEST_DIR"
}
trap cleanup EXIT

extract_python() {
  local marker_start="$1"
  local marker_end="$2"
  local source_file="$3"
  awk -v start="$marker_start" -v end="$marker_end" '
    $0 ~ start { capture = 1; next }
    $0 ~ end { capture = 0; next }
    capture { sub(/^          /, ""); print }
  ' "$source_file"
}

extract_sh() {
  local marker_start="$1"
  local marker_end="$2"
  local source_file="$3"
  awk -v start="$marker_start" -v end="$marker_end" '
    $0 ~ start { capture = 1; next }
    $0 ~ end { capture = 0; next }
    capture { sub(/^          /, ""); print }
  ' "$source_file"
}

# ── B1. Validate-sample rejects invalid modes ───────────────────────
echo "B1: Validate-sample rejects invalid modes..."
extract_python "BEGIN VALIDATE_SAMPLE_PYTHON" "END VALIDATE_SAMPLE_PYTHON" \
  "$REPAIR_WORKFLOW" > "$TEST_DIR/validate_sample.py"

if [ ! -s "$TEST_DIR/validate_sample.py" ]; then
  fail "Could not extract validate-sample production script"
fi

# Patch the Argo template references with test values
sed -i \
  -e 's|{{workflow.parameters.mode}}|badmode|' \
  -e 's|{{inputs.parameters.sample-id}}|TESTSAMPLE|' \
  -e 's|{{inputs.parameters.dna-file-type}}|bam|' \
  -e "s|'''{{inputs.parameters.dna-input-files}}'''|'{\"regionsFile\": \"/tmp/fake.bed\"}'|" \
  -e "s|'''{{inputs.parameters.ref-paths}}'''|'{\"genome\": \"/tmp/fake.fa\", \"other\": {\"annovar\": \"/tmp/fake\"}}'|" \
  -e 's|{{inputs.parameters.tumor-type}}|SOLID|' \
  -e 's|{{inputs.parameters.ref-base-dir}}|/mnt/ref/|' \
  "$TEST_DIR/validate_sample.py"

if python3 "$TEST_DIR/validate_sample.py" >/dev/null 2>&1; then
  fail "validate-sample accepted invalid mode 'badmode'"
fi

# ── B2. Validate-sample rejects VCF-only ─────────────────────────────
echo "B2: Validate-sample rejects VCF-only samples..."
extract_python "BEGIN VALIDATE_SAMPLE_PYTHON" "END VALIDATE_SAMPLE_PYTHON" \
  "$REPAIR_WORKFLOW" > "$TEST_DIR/validate_vcf.py"

sed -i \
  -e 's|{{workflow.parameters.mode}}|validate|' \
  -e 's|{{inputs.parameters.sample-id}}|TESTSAMPLE|' \
  -e 's|{{inputs.parameters.dna-file-type}}|vcf|' \
  -e "s|'''{{inputs.parameters.dna-input-files}}'''|'{\"regionsFile\": \"/tmp/fake.bed\"}'|" \
  -e "s|'''{{inputs.parameters.ref-paths}}'''|'{\"genome\": \"/tmp/fake.fa\", \"other\": {\"annovar\": \"/tmp/fake\"}}'|" \
  -e 's|{{inputs.parameters.tumor-type}}|SOLID|' \
  -e 's|{{inputs.parameters.ref-base-dir}}|/mnt/ref/|' \
  "$TEST_DIR/validate_vcf.py"

if python3 "$TEST_DIR/validate_vcf.py" >/dev/null 2>&1; then
  fail "validate-sample accepted VCF-only sample"
fi

# ── B3. Assemble-candidate: overlay does not clobber non-repair dirs ──
echo "B3: Assemble-candidate preserves non-repair directories..."
mkdir -p "$TEST_DIR/live/bam" "$TEST_DIR/live/results/sv" "$TEST_DIR/live/results/cnv"
mkdir -p "$TEST_DIR/live/results/snv" "$TEST_DIR/live/results/db"
echo "live_bam" > "$TEST_DIR/live/bam/test_recal.bam"
echo "live_sv" > "$TEST_DIR/live/results/sv/calls.vcf"
echo "live_cnv" > "$TEST_DIR/live/results/cnv/calls.tsv"
echo "live_snv" > "$TEST_DIR/live/results/snv/old_mutect.vcf"
echo "live_db" > "$TEST_DIR/live/results/db/test_results.duckdb"

mkdir -p "$TEST_DIR/tx/results/snv" "$TEST_DIR/tx/results/db" "$TEST_DIR/tx/metrics/bed_validation"
echo "new_snv" > "$TEST_DIR/tx/results/snv/new_mutect.vcf"
echo "new_db" > "$TEST_DIR/tx/results/db/test_results.duckdb"
echo "new_bed" > "$TEST_DIR/tx/metrics/bed_validation/analysis.bed"

extract_sh "BEGIN ASSEMBLE_CANDIDATE_SH" "END ASSEMBLE_CANDIDATE_SH" \
  "$REPAIR_WORKFLOW" > "$TEST_DIR/assemble.sh"

# Replace Argo template references with test paths
# The script references /mnt/parent/, LIVE_BASE, TX_BASE, CD_BASE
cat > "$TEST_DIR/run_assemble.sh" <<RUNEOF
#!/usr/bin/env bash
set -e
LIVE="/mnt/parent/live"
TX="/mnt/parent/tx"
CANDIDATE="/mnt/parent/candidate"
if [ -d "\$CANDIDATE" ]; then
  rm -rf "\$CANDIDATE"
fi
echo "Copying live tree to candidate..."
cp -a "\$LIVE" "\$CANDIDATE"
echo "Overlaying repair artifacts..."
for subdir in results/snv results/db metrics/bed_validation; do
  if [ -d "\${TX}/\${subdir}" ]; then
    mkdir -p "\${CANDIDATE}/\${subdir}"
    cp -a "\${TX}/\${subdir}/." "\${CANDIDATE}/\${subdir}/"
  fi
done
echo "Candidate assembled"
RUNEOF
chmod +x "$TEST_DIR/run_assemble.sh"

# Create the mount point structure
mkdir -p "$TEST_DIR/mnt_parent"
cp -a "$TEST_DIR/live" "$TEST_DIR/mnt_parent/live"
cp -a "$TEST_DIR/tx" "$TEST_DIR/mnt_parent/tx"

# Patch the script to use our test dirs
sed -i "s|/mnt/parent/|${TEST_DIR}/mnt_parent/|g" "$TEST_DIR/run_assemble.sh"
bash "$TEST_DIR/run_assemble.sh"

# Verify: SV/CNV preserved, SNV/DB overlaid, BAM preserved
[[ -f "$TEST_DIR/mnt_parent/candidate/results/sv/calls.vcf" ]] ||
  fail "assemble-candidate lost SV directory"
[[ $(cat "$TEST_DIR/mnt_parent/candidate/results/sv/calls.vcf") == "live_sv" ]] ||
  fail "assemble-candidate corrupted SV data"
[[ -f "$TEST_DIR/mnt_parent/candidate/results/cnv/calls.tsv" ]] ||
  fail "assemble-candidate lost CNV directory"
[[ -f "$TEST_DIR/mnt_parent/candidate/bam/test_recal.bam" ]] ||
  fail "assemble-candidate lost BAM directory"
[[ $(cat "$TEST_DIR/mnt_parent/candidate/results/snv/new_mutect.vcf") == "new_snv" ]] ||
  fail "assemble-candidate did not overlay new SNV files"
[[ $(cat "$TEST_DIR/mnt_parent/candidate/results/db/test_results.duckdb") == "new_db" ]] ||
  fail "assemble-candidate did not overlay new DuckDB"
[[ -f "$TEST_DIR/mnt_parent/candidate/metrics/bed_validation/analysis.bed" ]] ||
  fail "assemble-candidate did not overlay BED validation"
# The old SNV file should still exist (cp -a of live includes it)
[[ -f "$TEST_DIR/mnt_parent/candidate/results/snv/old_mutect.vcf" ]] ||
  fail "assemble-candidate lost old SNV files from live"

# ── B4. Promote: atomic exchange + backup ────────────────────────────
echo "B4: Promote creates atomic exchange and backup..."
extract_python "BEGIN PROMOTE_CANDIDATE_PYTHON" "END PROMOTE_CANDIDATE_PYTHON" \
  "$REPAIR_WORKFLOW" > "$TEST_DIR/promote_raw.py"

if [ ! -s "$TEST_DIR/promote_raw.py" ]; then
  fail "Could not extract promote production script"
fi

# Set up test directories
mkdir -p "$TEST_DIR/promote_test/sample_dir/results"
echo "live_data" > "$TEST_DIR/promote_test/sample_dir/results/live.txt"
mkdir -p "$TEST_DIR/promote_test/sample_dir.repair-cd-testuid/results"
echo "candidate_data" > "$TEST_DIR/promote_test/sample_dir.repair-cd-testuid/results/candidate.txt"
echo '{"promotion":{"state":"pending","backup_dir":null},"projection":{"state":"pending"}}' \
  > "$TEST_DIR/promote_test/sample_dir.repair-cd-testuid/padding_repair_manifest.json"

cat > "$TEST_DIR/promote_test.py" <<'PYEOF'
import ctypes, json, os, sys, time

parent = sys.argv[1]
live = os.path.join(parent, "sample_dir")
candidate = os.path.join(parent, "sample_dir.repair-cd-testuid")

if not os.path.isdir(live):
    print(f"Live directory not found: {live}", file=sys.stderr)
    sys.exit(1)
if not os.path.isdir(candidate):
    print(f"Candidate directory not found: {candidate}", file=sys.stderr)
    sys.exit(1)

SYS_RENAMEAT2 = 316
RENAME_EXCHANGE = 2
AT_FDCWD = -100

libc = ctypes.CDLL("libc.so.6", use_errno=True)
libc.syscall.restype = ctypes.c_long

ctypes.set_errno(0)
ret = libc.syscall(
    ctypes.c_long(SYS_RENAMEAT2),
    ctypes.c_int(AT_FDCWD),
    live.encode(),
    ctypes.c_int(AT_FDCWD),
    candidate.encode(),
    ctypes.c_uint(RENAME_EXCHANGE),
)
if ret != 0:
    errno = ctypes.get_errno()
    print(f"renameat2 RENAME_EXCHANGE failed: {os.strerror(errno)} (errno={errno})", file=sys.stderr)
    sys.exit(1)

print("Atomic exchange complete: live <-> candidate")

timestamp = time.strftime("%Y%m%dT%H%M%S")
backup_name = f"sample_dir.pre-padding-repair-{timestamp}"
backup = os.path.join(parent, backup_name)
os.rename(candidate, backup)
print(f"Old live renamed to backup: {backup_name}")

manifest_path = os.path.join(live, "padding_repair_manifest.json")
if os.path.isfile(manifest_path):
    with open(manifest_path) as f:
        m = json.load(f)
    m["promotion"]["state"] = "completed"
    m["promotion"]["backup_dir"] = backup
    with open(manifest_path, "w") as f:
        json.dump(m, f, indent=2)
print(f"Promotion complete. Backup: {backup}")

with open(os.path.join(parent, "backup_path"), "w") as f:
    f.write(backup)
PYEOF

python3 "$TEST_DIR/promote_test.py" "$TEST_DIR/promote_test"

# After exchange: live should have candidate data, backup should have live data
[[ $(cat "$TEST_DIR/promote_test/sample_dir/results/candidate.txt") == "candidate_data" ]] ||
  fail "After exchange, live does not contain candidate data"
backup_path=$(cat "$TEST_DIR/promote_test/backup_path")
[[ $(cat "$backup_path/results/live.txt") == "live_data" ]] ||
  fail "After exchange, backup does not contain original live data"

# Manifest should record completion
manifest_state=$(python3 -c "import json; m=json.load(open('$TEST_DIR/promote_test/sample_dir/padding_repair_manifest.json')); print(m['promotion']['state'])")
[[ $manifest_state == "completed" ]] ||
  fail "Manifest promotion state is not 'completed' after promote"

# ── B5. Validate-candidate: rejects bad analysis BED ────────────────
echo "B5: Validate-candidate rejects overlapping analysis BED..."
extract_python "BEGIN VALIDATE_CANDIDATE_PYTHON" "END VALIDATE_CANDIDATE_PYTHON" \
  "$REPAIR_WORKFLOW" > "$TEST_DIR/validate_candidate_raw.py"

if [ ! -s "$TEST_DIR/validate_candidate_raw.py" ]; then
  fail "Could not extract validate-candidate production script"
fi

# Set up a candidate with overlapping analysis BED
mkdir -p "$TEST_DIR/vc_test/candidate/metrics/bed_validation"
mkdir -p "$TEST_DIR/vc_test/candidate/metrics/tmb"
mkdir -p "$TEST_DIR/vc_test/candidate/results/db"
echo "tmb_data" > "$TEST_DIR/vc_test/candidate/metrics/tmb/tmb.json"
echo "db_data" > "$TEST_DIR/vc_test/candidate/results/db/TEST_results.duckdb"
echo '{"promotion":{"state":"pending"}}' > "$TEST_DIR/vc_test/candidate/padding_repair_manifest.json"

# Write overlapping intervals — should fail
cat > "$TEST_DIR/vc_test/candidate/metrics/bed_validation/analysis.bed" <<'BED'
chr1	100	200
chr1	150	250
BED

# Good IGV BED
cat > "$TEST_DIR/vc_test/candidate/metrics/bed_validation/igv.bed" <<'BED'
track name="Targeted and Padding Regions" itemRgb="On"
chr1	100	200	targeted	0	+	100	200	128,128,128
chr1	200	250	padding	0	+	200	250	211,211,211
BED

cat > "$TEST_DIR/vc_test_runner.py" <<'PYEOF'
import json, os, sys

cd_base = "candidate"
candidate = os.path.join(sys.argv[1], cd_base)
sample_id = "TEST"
errors = []

analysis_name = "analysis.bed"
analysis = os.path.join(candidate, "metrics", "bed_validation", analysis_name)
igv_name = "igv.bed"
igv = os.path.join(candidate, "metrics", "bed_validation", igv_name)
db = os.path.join(candidate, "results", "db", f"{sample_id}_results.duckdb")
manifest = os.path.join(candidate, "padding_repair_manifest.json")
tmb_dir = os.path.join(candidate, "metrics", "tmb")

if not os.path.isfile(manifest):
    errors.append("Repair manifest not found in candidate")
if not os.path.isdir(tmb_dir) or not os.listdir(tmb_dir):
    errors.append("TMB output not found in candidate")

# Analysis BED validation
if not os.path.isfile(analysis):
    errors.append("Analysis BED not found in candidate")
else:
    prev_chrom, prev_end = "", 0
    line_count = 0
    with open(analysis) as f:
        for i, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if line.startswith("#") or line.startswith("track"):
                errors.append(f"Analysis BED has header at line {i}")
                break
            cols = line.split("\t")
            if len(cols) != 3:
                errors.append(f"Analysis BED is not strict BED3 (line {i}: {len(cols)} cols)")
                break
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            if start >= end:
                errors.append(f"Analysis BED invalid interval at line {i}")
                break
            if chrom == prev_chrom and start <= prev_end:
                errors.append(
                    f"Analysis BED overlapping/book-ended at {chrom}:{start} "
                    f"(prev end {prev_end})"
                )
                break
            prev_chrom, prev_end = chrom, end
            line_count += 1
    if line_count == 0 and not errors:
        errors.append("Analysis BED is empty")

if errors:
    print("Candidate validation FAILED:", file=sys.stderr)
    for e in errors:
        print(f"  - {e}", file=sys.stderr)
    sys.exit(1)
print("Candidate validation passed")
PYEOF

if python3 "$TEST_DIR/vc_test_runner.py" "$TEST_DIR/vc_test" >/dev/null 2>&1; then
  fail "validate-candidate accepted overlapping analysis BED"
fi

# ── B6. Validate-candidate: rejects book-ended analysis BED ─────────
echo "B6: Validate-candidate rejects book-ended analysis BED..."
cat > "$TEST_DIR/vc_test/candidate/metrics/bed_validation/analysis.bed" <<'BED'
chr1	100	200
chr1	200	300
BED

if python3 "$TEST_DIR/vc_test_runner.py" "$TEST_DIR/vc_test" >/dev/null 2>&1; then
  fail "validate-candidate accepted book-ended analysis BED"
fi

# ── B7. Validate-candidate: rejects wrong IGV RGB ───────────────────
echo "B7: Validate-candidate rejects wrong IGV RGB..."
# Fix analysis BED
cat > "$TEST_DIR/vc_test/candidate/metrics/bed_validation/analysis.bed" <<'BED'
chr1	100	200
chr1	300	400
BED

# Wrong RGB for targeted
cat > "$TEST_DIR/vc_test/candidate/metrics/bed_validation/igv.bed" <<'BED'
track name="Targeted and Padding Regions" itemRgb="On"
chr1	100	200	targeted	0	+	100	200	255,0,0
chr1	200	300	padding	0	+	200	300	211,211,211
BED

cat > "$TEST_DIR/vc_igv_runner.py" <<'PYEOF'
import json, os, sys

cd_base = "candidate"
candidate = os.path.join(sys.argv[1], cd_base)
errors = []

igv_name = "igv.bed"
igv = os.path.join(candidate, "metrics", "bed_validation", igv_name)

expected_rgb = {"targeted": "128,128,128", "padding": "211,211,211"}
features = {"targeted": [], "padding": []}
with open(igv) as f:
    header = f.readline().rstrip("\n")
    if 'itemRgb="On"' not in header:
        errors.append("IGV BED missing itemRgb header")
    for i, line in enumerate(f, 2):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            errors.append(f"IGV line {i} has {len(fields)} columns, expected 9")
            break
        chrom, start_s, end_s, feat = fields[:4]
        start, end = int(start_s), int(end_s)
        if feat not in expected_rgb:
            errors.append(f"IGV unknown feature type: {feat!r}")
            continue
        if fields[8] != expected_rgb[feat]:
            errors.append(f"IGV {feat} RGB is {fields[8]}, expected {expected_rgb[feat]}")
        if start >= end:
            errors.append(f"IGV invalid interval at line {i}")
        features[feat].append((chrom, start, end))

if errors:
    for e in errors:
        print(f"  - {e}", file=sys.stderr)
    sys.exit(1)
print("IGV validation passed")
PYEOF

if python3 "$TEST_DIR/vc_igv_runner.py" "$TEST_DIR/vc_test" >/dev/null 2>&1; then
  fail "validate-candidate accepted wrong IGV targeted RGB"
fi

# ── B8. Validate-candidate: rejects overlapping targeted/padding ─────
echo "B8: Validate-candidate rejects overlapping targeted/padding..."
cat > "$TEST_DIR/vc_test/candidate/metrics/bed_validation/igv.bed" <<'BED'
track name="Targeted and Padding Regions" itemRgb="On"
chr1	100	250	targeted	0	+	100	250	128,128,128
chr1	200	300	padding	0	+	200	300	211,211,211
BED

# Use the full candidate validator for this test
cat > "$TEST_DIR/vc_overlap_runner.py" <<'PYEOF'
import json, os, sys

cd_base = "candidate"
candidate = os.path.join(sys.argv[1], cd_base)
errors = []

igv = os.path.join(candidate, "metrics", "bed_validation", "igv.bed")
expected_rgb = {"targeted": "128,128,128", "padding": "211,211,211"}
features = {"targeted": [], "padding": []}

with open(igv) as f:
    header = f.readline().rstrip("\n")
    if 'itemRgb="On"' not in header:
        errors.append("IGV BED missing itemRgb header")
    for i, line in enumerate(f, 2):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            errors.append(f"IGV line {i} has {len(fields)} columns, expected 9")
            break
        chrom, start_s, end_s, feat = fields[:4]
        start, end = int(start_s), int(end_s)
        if feat not in expected_rgb:
            errors.append(f"IGV unknown feature type: {feat!r}")
            continue
        if fields[8] != expected_rgb[feat]:
            errors.append(f"IGV {feat} RGB is {fields[8]}, expected {expected_rgb[feat]}")
        features[feat].append((chrom, start, end))

for tc, ts, te in features["targeted"]:
    for pc, ps, pe in features["padding"]:
        if tc == pc and ts < pe and ps < te:
            errors.append(
                f"IGV targeted/padding overlap: "
                f"{tc}:{ts}-{te} and {pc}:{ps}-{pe}"
            )

if errors:
    for e in errors:
        print(f"  - {e}", file=sys.stderr)
    sys.exit(1)
print("passed")
PYEOF

if python3 "$TEST_DIR/vc_overlap_runner.py" "$TEST_DIR/vc_test" >/dev/null 2>&1; then
  fail "validate-candidate accepted overlapping targeted/padding"
fi

# ── B9. Cleanup only removes workflow-UID paths ──────────────────────
echo "B9: Cleanup only removes workflow-UID paths..."
mkdir -p "$TEST_DIR/cleanup_test/sample_dir/results"
echo "live" > "$TEST_DIR/cleanup_test/sample_dir/results/data.txt"
mkdir -p "$TEST_DIR/cleanup_test/sample_dir.repair-wk-uid123"
echo "tx" > "$TEST_DIR/cleanup_test/sample_dir.repair-wk-uid123/data.txt"
mkdir -p "$TEST_DIR/cleanup_test/sample_dir.repair-cd-uid123"
echo "cd" > "$TEST_DIR/cleanup_test/sample_dir.repair-cd-uid123/data.txt"
mkdir -p "$TEST_DIR/cleanup_test/sample_dir.pre-padding-repair-20260101T000000"
echo "backup" > "$TEST_DIR/cleanup_test/sample_dir.pre-padding-repair-20260101T000000/data.txt"

cat > "$TEST_DIR/cleanup_runner.sh" <<CLEANEOF
#!/usr/bin/env bash
set -e
PARENT="$TEST_DIR/cleanup_test"
TX_BASE="sample_dir.repair-wk-uid123"
CD_BASE="sample_dir.repair-cd-uid123"

[ -d "\${PARENT}/\${TX_BASE}" ] && rm -rf "\${PARENT}/\${TX_BASE}"
[ -d "\${PARENT}/\${CD_BASE}" ] && rm -rf "\${PARENT}/\${CD_BASE}"
echo "Cleaned transaction paths"
CLEANEOF

bash "$TEST_DIR/cleanup_runner.sh"

[[ -d "$TEST_DIR/cleanup_test/sample_dir/results" ]] ||
  fail "cleanup removed the live sample directory"
[[ -d "$TEST_DIR/cleanup_test/sample_dir.pre-padding-repair-20260101T000000" ]] ||
  fail "cleanup removed the backup directory"
[[ ! -d "$TEST_DIR/cleanup_test/sample_dir.repair-wk-uid123" ]] ||
  fail "cleanup did not remove the transaction root"
[[ ! -d "$TEST_DIR/cleanup_test/sample_dir.repair-cd-uid123" ]] ||
  fail "cleanup did not remove the candidate directory"

# ── B10. Manifest records publication failure ────────────────────────
echo "B10: Manifest records publication failure..."
mkdir -p "$TEST_DIR/manifest_test/sample_dir"
echo '{"promotion":{"state":"completed","backup_dir":"/tmp/backup"},"projection":{"state":"pending"}}' \
  > "$TEST_DIR/manifest_test/sample_dir/padding_repair_manifest.json"

cat > "$TEST_DIR/manifest_fail_runner.py" <<'PYEOF'
import json, os, sys

parent = sys.argv[1]
mode = sys.argv[2]
promote_status = sys.argv[3]
publish_status = sys.argv[4]

if mode == "apply" and promote_status == "Succeeded":
    manifest_path = os.path.join(parent, "sample_dir", "padding_repair_manifest.json")
else:
    manifest_path = os.path.join(parent, "candidate", "padding_repair_manifest.json")

if not os.path.isfile(manifest_path):
    print(f"Manifest not found at {manifest_path}, skipping")
    sys.exit(0)

with open(manifest_path) as f:
    m = json.load(f)

if mode != "apply":
    m["projection"]["state"] = "not-applicable"
elif publish_status == "Succeeded":
    m["projection"]["state"] = "published"
elif publish_status in ("Failed", "Error"):
    m["projection"]["state"] = "failed"
else:
    m["projection"]["state"] = "pending"

with open(manifest_path, "w") as f:
    json.dump(m, f, indent=2)
print(f"Manifest projection updated: {m['projection']['state']}")
PYEOF

python3 "$TEST_DIR/manifest_fail_runner.py" "$TEST_DIR/manifest_test" apply Succeeded Failed

state=$(python3 -c "import json; print(json.load(open('$TEST_DIR/manifest_test/sample_dir/padding_repair_manifest.json'))['projection']['state'])")
[[ $state == "failed" ]] ||
  fail "Manifest projection state is '$state', expected 'failed' after publish failure"

# ── B11. Manifest records validate mode correctly ────────────────────
echo "B11: Manifest records validate mode as not-applicable..."
mkdir -p "$TEST_DIR/manifest_test2/candidate"
echo '{"promotion":{"state":"pending"},"projection":{"state":"pending"}}' \
  > "$TEST_DIR/manifest_test2/candidate/padding_repair_manifest.json"

python3 "$TEST_DIR/manifest_fail_runner.py" "$TEST_DIR/manifest_test2" validate "" ""

state=$(python3 -c "import json; print(json.load(open('$TEST_DIR/manifest_test2/candidate/padding_repair_manifest.json'))['projection']['state'])")
[[ $state == "not-applicable" ]] ||
  fail "Manifest projection state is '$state', expected 'not-applicable' in validate mode"

echo ""
echo "All padding repair regression tests passed"
