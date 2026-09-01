#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"

WORKFLOW=workflows/somatic/v2/somatic-dna-padding-repair.yaml

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

task_block() {
  local task_name=$1
  awk -v needle="          - name: $task_name" '
    $0 == needle { printing = 1; first = 1 }
    printing && !first && /^          - name: / { exit }
    printing { print; first = 0 }
  ' "$WORKFLOW"
}

extract_block() {
  local begin=$1
  local end=$2
  local output=$3
  awk -v begin="$begin" -v end="$end" '
    index($0, begin) { printing = 1; next }
    index($0, end) { printing = 0; exit }
    printing {
      sub(/^          /, "")
      print
    }
  ' "$WORKFLOW" > "$output"
  [[ -s $output ]] || fail "Production block $begin was not extracted"
}

echo "=== Static workflow contract ==="

param_count=$(awk '
  /^  arguments:/ { in_args=1; next }
  in_args && /^    parameters:/ { in_params=1; next }
  in_params && /^      - name: / { count++ }
  in_params && /^  [^ ]/ { exit }
  END { print count+0 }
' "$WORKFLOW")
[[ $param_count == 5 ]] || fail "Expected five public parameters, found $param_count"
for parameter in node-name user-id run-id sample-uuid mode; do
  grep -qE "^      - name: ${parameter}$" "$WORKFLOW" ||
    fail "Missing public parameter $parameter"
done

grep -q 'name: "padding-repair-{{workflow.parameters.sample-uuid}}"' "$WORKFLOW" ||
  fail "Sample mutex is missing"
grep -q 'repair-wk-{{workflow.uid}}' "$WORKFLOW" ||
  fail "Transaction path is not workflow-UID scoped"
grep -q 'repair-cd-{{workflow.uid}}' "$WORKFLOW" ||
  fail "Candidate path is not workflow-UID scoped"

mutect=$(task_block mutect2)
bcftools=$(task_block bcftools-filter)
compute_tmb=$(task_block compute-tmb)
publish_tmb=$(task_block publish-tmb)
finalize=$(task_block finalize)

grep -A1 'name: regions-file' <<< "$mutect" | grep -q validated-bed-analysis-file ||
  fail "Mutect2 does not consume the padded analysis BED"
grep -A1 'name: regions-file' <<< "$bcftools" | grep -q validated-bed-analysis-file ||
  fail "bcftools does not consume the padded analysis BED"
grep -A1 'name: target-bed-file' <<< "$compute_tmb" | grep -q validated-bed-file ||
  fail "TMB compute does not consume the original target BED"
if grep -qE 'validated-bed-(analysis|igv)-file' <<< "$compute_tmb"; then
  fail "TMB compute consumes a padded or presentation BED"
fi
grep -q 'template: publish-tmb-result' <<< "$publish_tmb" ||
  fail "TMB publication does not use the publish-only template"
grep -q 'tmb_api_payload.json' <<< "$publish_tmb" ||
  fail "TMB publication does not consume the staged payload"
if grep -qE 'duckdb-file|target-bed-file|cosmic-vcf|gnomad-vcf' <<< "$publish_tmb"; then
  fail "TMB publication still has recalculation inputs"
fi
grep -q 'test:///mnt/output/tmb_api_payload.json' "$WORKFLOW" ||
  fail "TMB compute does not stage the exact API payload"
grep -q -- '--api-token STAGED' "$WORKFLOW" ||
  fail "TMB compute does not use the non-secret staged token"

for terminal in \
  promote.Failed promote.Errored promote.Omitted \
  publish-tmb.Failed publish-tmb.Errored publish-tmb.Omitted \
  validate-candidate.Failed validate-candidate.Errored validate-candidate.Omitted; do
  grep -q "$terminal" <<< "$finalize" ||
    fail "Finalizer does not cover $terminal"
done

if grep -qE 'template: maf-to-db|name: maf-to-db-template' "$WORKFLOW"; then
  fail "Repair workflow must not call the global MAF importer"
fi
for excluded in duckdb-sv duckdb-cnv duckdb-fusions; do
  if grep -q "template: ${excluded}" "$WORKFLOW"; then
    fail "Repair workflow must not call $excluded"
  fi
done
if grep -q 'pip install.*duckdb' "$WORKFLOW"; then
  fail "DuckDB validation installs a runtime dependency from the network"
fi
[[ $(grep -c 'image: bergun/gennext-tmb-calc:0.5.0' "$WORKFLOW") -ge 4 ]] ||
  fail "DuckDB/TMB validation runtimes are not image-pinned"
grep -q 'ctypes.CDLL(None' "$WORKFLOW" ||
  fail "Promotion is not libc-portable"
grep -q 'pre-padding-repair-{workflow_uid}' "$WORKFLOW" ||
  fail "Backup name is not workflow-UID scoped"
grep -qE 'gennext.bio/source-revision: "[0-9a-f]{40}"' "$WORKFLOW" ||
  fail "Workflow source revision is not pinned to a substantive commit"

echo "=== Production-script behavior ==="

TEST_DIR=$(mktemp -d)
cleanup() {
  if [[ -n ${SERVER_PID:-} ]]; then
    kill "$SERVER_PID" 2>/dev/null || true
  fi
  rm -rf -- "$TEST_DIR"
}
trap cleanup EXIT

mkdir -p "$TEST_DIR/bin"
cat > "$TEST_DIR/bin/python3" <<'PYTHON_WRAPPER'
#!/usr/bin/env bash
set -euo pipefail

env_args=()
for name in \
  DB_PATH INVENTORY_PATH INVENTORY_JSON \
  MOUNT_PARENT CANDIDATE_BASE SAMPLE_ID ANALYSIS_BED_NAME IGV_BED_NAME \
  TARGET_BED_NAME INPUT_BAM_NAME WORKFLOW_UID WORKFLOW_NAME SOURCE_REVISION \
  SAMPLE_UUID RUN_ID USER_ID NODE_NAME REPAIR_MODE REFERENCE_FAI \
  LIVE_BASE OUTPUT_PARENT_HOST BACKUP_OUTPUT_PATH TX_BASE \
  VALIDATION_STATUS PROMOTE_STATUS PUBLISH_STATUS \
  REQUEST_JSON API_ENDPOINT GENNEXT_CLUSTER_TOKEN; do
  if [[ -v $name ]]; then
    env_args+=(-e "$name")
  fi
done

exec docker run --rm -i \
  --add-host=host.docker.internal:host-gateway \
  -v "$TEST_ROOT:$TEST_ROOT" \
  "${env_args[@]}" \
  --entrypoint python \
  bergun/gennext-tmb-calc:0.5.0 \
  "$@"
PYTHON_WRAPPER
chmod +x "$TEST_DIR/bin/python3"
export TEST_ROOT="$TEST_DIR"

extract_block "# BEGIN ASSEMBLE_CANDIDATE_SH" "# END ASSEMBLE_CANDIDATE_SH" "$TEST_DIR/assemble.sh"
extract_block "# BEGIN SNAPSHOT_DUCKDB_PYTHON" "# END SNAPSHOT_DUCKDB_PYTHON" "$TEST_DIR/snapshot.sh"
extract_block "# BEGIN VERIFY_DUCKDB_PYTHON" "# END VERIFY_DUCKDB_PYTHON" "$TEST_DIR/verify.sh"
extract_block "# BEGIN WRITE_MANIFEST_PYTHON" "# END WRITE_MANIFEST_PYTHON" "$TEST_DIR/write_manifest.py"
extract_block "# BEGIN VALIDATE_CANDIDATE_PYTHON" "# END VALIDATE_CANDIDATE_PYTHON" "$TEST_DIR/validate_candidate.py"
extract_block "# BEGIN PROMOTE_CANDIDATE_PYTHON" "# END PROMOTE_CANDIDATE_PYTHON" "$TEST_DIR/promote.py"
extract_block "# BEGIN FINALIZE_REPAIR_PYTHON" "# END FINALIZE_REPAIR_PYTHON" "$TEST_DIR/finalize.py"
extract_block "# BEGIN PUBLISH_TMB_PYTHON" "# END PUBLISH_TMB_PYTHON" "$TEST_DIR/publish.py"

# Candidate assembly must replace all affected trees and preserve unrelated output.
PARENT="$TEST_DIR/assembly"
LIVE_BASE=sample
TX_BASE=sample.repair-wk-test-uid
CANDIDATE_BASE=sample.repair-cd-test-uid
mkdir -p \
  "$PARENT/$LIVE_BASE/results/snv" \
  "$PARENT/$LIVE_BASE/results/db" \
  "$PARENT/$LIVE_BASE/results/cnv" \
  "$PARENT/$LIVE_BASE/metrics/tmb" \
  "$PARENT/$LIVE_BASE/metrics/bed_validation" \
  "$PARENT/$TX_BASE/results/snv" \
  "$PARENT/$TX_BASE/results/db" \
  "$PARENT/$TX_BASE/metrics/bed_validation"
printf old > "$PARENT/$LIVE_BASE/results/snv/stale.vcf"
printf old-tmb > "$PARENT/$LIVE_BASE/metrics/tmb/stale.json"
printf old-bed > "$PARENT/$LIVE_BASE/metrics/bed_validation/stale.bed"
printf keep-cnv > "$PARENT/$LIVE_BASE/results/cnv/cnv.txt"
printf new-snv > "$PARENT/$TX_BASE/results/snv/new.vcf"
printf new-db > "$PARENT/$TX_BASE/results/db/S_results.duckdb"
printf new-bed > "$PARENT/$TX_BASE/metrics/bed_validation/analysis.bed"

MOUNT_PARENT="$PARENT" LIVE_BASE="$LIVE_BASE" TX_BASE="$TX_BASE" \
  CANDIDATE_BASE="$CANDIDATE_BASE" sh "$TEST_DIR/assemble.sh"

CANDIDATE="$PARENT/$CANDIDATE_BASE"
[[ -f $CANDIDATE/results/snv/new.vcf ]] || fail "Fresh SNV output was not assembled"
[[ ! -e $CANDIDATE/results/snv/stale.vcf ]] || fail "Stale SNV output survived assembly"
[[ ! -e $CANDIDATE/metrics/tmb ]] || fail "Stale TMB output survived assembly"
[[ ! -e $CANDIDATE/metrics/bed_validation/stale.bed ]] ||
  fail "Stale BED validation output survived assembly"
[[ $(<"$CANDIDATE/results/cnv/cnv.txt") == keep-cnv ]] ||
  fail "Unrelated CNV output changed during assembly"

# Snapshot/verify must compare content and definitions, not only row counts.
DB_DIR="$TEST_DIR/duckdb"
mkdir -p "$DB_DIR"
DB_PATH="$DB_DIR/results.duckdb"
PATH="$TEST_DIR/bin:$PATH" python3 - "$DB_PATH" <<'PY'
import duckdb
import sys

db = duckdb.connect(sys.argv[1])
db.execute("CREATE TABLE cnv(id INTEGER, value VARCHAR)")
db.execute("INSERT INTO cnv VALUES (1, 'alpha'), (2, 'beta')")
db.execute("CREATE VIEW cnv_view AS SELECT * FROM cnv")
db.execute("CREATE TABLE oncokb_snv(id INTEGER)")
db.execute("CREATE TABLE cancervar_snv(id INTEGER)")
db.execute("CREATE VIEW combined_snv AS SELECT * FROM oncokb_snv")
db.execute("CREATE VIEW advanced_combined_snv AS SELECT * FROM oncokb_snv")
db.close()
PY

INVENTORY_PATH="$TEST_DIR/inventory.json" DB_PATH="$DB_PATH" \
  PATH="$TEST_DIR/bin:$PATH" sh "$TEST_DIR/snapshot.sh"
INVENTORY_JSON=$(<"$TEST_DIR/inventory.json")

PATH="$TEST_DIR/bin:$PATH" python3 - "$DB_PATH" <<'PY'
import duckdb
import sys

db = duckdb.connect(sys.argv[1])
db.execute("INSERT INTO oncokb_snv VALUES (1)")
db.close()
PY
DB_PATH="$DB_PATH" INVENTORY_JSON="$INVENTORY_JSON" \
  PATH="$TEST_DIR/bin:$PATH" sh "$TEST_DIR/verify.sh"

PATH="$TEST_DIR/bin:$PATH" python3 - "$DB_PATH" <<'PY'
import duckdb
import sys

db = duckdb.connect(sys.argv[1])
db.execute("UPDATE cnv SET value = 'changed' WHERE id = 1")
db.close()
PY
if DB_PATH="$DB_PATH" INVENTORY_JSON="$INVENTORY_JSON" \
  PATH="$TEST_DIR/bin:$PATH" sh "$TEST_DIR/verify.sh" >/dev/null 2>&1; then
  fail "DuckDB verification accepted same-count non-SNV content drift"
fi

PATH="$TEST_DIR/bin:$PATH" python3 - "$DB_PATH" <<'PY'
import duckdb
import sys

db = duckdb.connect(sys.argv[1])
db.execute("UPDATE cnv SET value = 'alpha' WHERE id = 1")
db.execute("CREATE OR REPLACE VIEW cnv_view AS SELECT id, value FROM cnv WHERE id > 0")
db.close()
PY
if DB_PATH="$DB_PATH" INVENTORY_JSON="$INVENTORY_JSON" \
  PATH="$TEST_DIR/bin:$PATH" sh "$TEST_DIR/verify.sh" >/dev/null 2>&1; then
  fail "DuckDB verification accepted same-result non-SNV view definition drift"
fi

# Build a complete candidate, write its production manifest, and validate it.
SAMPLE_ID=S
mkdir -p \
  "$CANDIDATE/bam" \
  "$CANDIDATE/results/snv/vcf" \
  "$CANDIDATE/results/db" \
  "$CANDIDATE/metrics/bed_validation" \
  "$CANDIDATE/metrics/tmb"
printf bam > "$CANDIDATE/bam/${SAMPLE_ID}_recal.bam"
printf bai > "$CANDIDATE/bam/${SAMPLE_ID}_recal.bam.bai"
printf db > "$CANDIDATE/results/db/${SAMPLE_ID}_results.duckdb"
for artifact in \
  "vcf/${SAMPLE_ID}_mutect2_marked.vcf.gz" \
  "vcf/${SAMPLE_ID}_mutect2_filtered.vcf.gz" \
  "${SAMPLE_ID}_mutect2_vlod.vcf.gz" \
  "${SAMPLE_ID}.hg38_multianno.txt" \
  "${SAMPLE_ID}.hg38_multianno.txt.cancervar" \
  "${SAMPLE_ID}.maf" \
  "${SAMPLE_ID}_oncokb.maf"; do
  printf artifact > "$CANDIDATE/results/snv/$artifact"
done

TARGET_BED_NAME=target_validated.bed
ANALYSIS_BED_NAME=target_analysis.bed
IGV_BED_NAME=target_igv.bed
cat > "$CANDIDATE/metrics/bed_validation/$TARGET_BED_NAME" <<'EOF_TARGET'
chr1	10	20
EOF_TARGET
cat > "$CANDIDATE/metrics/bed_validation/$ANALYSIS_BED_NAME" <<'EOF_ANALYSIS'
chr1	0	30
chr2	5	10
EOF_ANALYSIS
cat > "$CANDIDATE/metrics/bed_validation/$IGV_BED_NAME" <<'EOF_IGV'
track name="Targeted and Padding Regions" itemRgb="On"
chr1	0	10	padding	0	.	0	10	211,211,211
chr1	10	20	targeted	0	.	10	20	128,128,128
chr1	20	30	padding	0	.	20	30	211,211,211
chr2	5	10	targeted	0	.	5	10	128,128,128
EOF_IGV
cat > "$TEST_DIR/reference.fai" <<'EOF_FAI'
chr1	100	0	0	0
chr2	100	0	0	0
EOF_FAI
cat > "$CANDIDATE/metrics/tmb/tmb_report_${SAMPLE_ID}.json" <<'EOF_TMB'
{"results":{"tmb_per_mb":2.5,"tmb_variants":5,"total_target_bases":2000000,"cds_tmb_per_mb":3.0,"cds_tmb_variants":3,"targeted_cds_bases":1000000}}
EOF_TMB
printf report > "$CANDIDATE/metrics/tmb/tmb_report_${SAMPLE_ID}.txt"
cat > "$CANDIDATE/metrics/tmb/tmb_api_payload.json" <<'EOF_PAYLOAD'
{"runId":"run-1","sampleId":"S","somaticData":{"tmbValue":"2.50","tmbVariantCount":5,"tmbTargetRegionSize":2000000,"cdsTmbValue":"3.00","cdsTmbVariantCount":3,"cdsTmbTargetRegionSize":1000000},"token":"STAGED"}
EOF_PAYLOAD

COMMON_MANIFEST_ENV=(
  MOUNT_PARENT="$PARENT"
  CANDIDATE_BASE="$CANDIDATE_BASE"
  SAMPLE_ID="$SAMPLE_ID"
  ANALYSIS_BED_NAME="$ANALYSIS_BED_NAME"
  IGV_BED_NAME="$IGV_BED_NAME"
  TARGET_BED_NAME="$TARGET_BED_NAME"
  INPUT_BAM_NAME="${SAMPLE_ID}_recal.bam"
  WORKFLOW_UID=test-uid
  WORKFLOW_NAME=test-repair
  SOURCE_REVISION=test-sha
  SAMPLE_UUID=sample-uuid
  RUN_ID=run-1
  USER_ID=user-1
  NODE_NAME=node-1
  REPAIR_MODE=validate
)
env "${COMMON_MANIFEST_ENV[@]}" PATH="$TEST_DIR/bin:$PATH" \
  python3 "$TEST_DIR/write_manifest.py"

env \
  MOUNT_PARENT="$PARENT" \
  REFERENCE_FAI="$TEST_DIR/reference.fai" \
  CANDIDATE_BASE="$CANDIDATE_BASE" \
  SAMPLE_ID="$SAMPLE_ID" \
  ANALYSIS_BED_NAME="$ANALYSIS_BED_NAME" \
  IGV_BED_NAME="$IGV_BED_NAME" \
  TARGET_BED_NAME="$TARGET_BED_NAME" \
  WORKFLOW_UID=test-uid \
  SOURCE_REVISION=test-sha \
  PATH="$TEST_DIR/bin:$PATH" \
  python3 "$TEST_DIR/validate_candidate.py"

cp "$CANDIDATE/metrics/bed_validation/$ANALYSIS_BED_NAME" "$TEST_DIR/analysis.good"
cat > "$CANDIDATE/metrics/bed_validation/$ANALYSIS_BED_NAME" <<'EOF_BAD_ORDER'
chr2	5	10
chr1	0	30
EOF_BAD_ORDER
if env \
  MOUNT_PARENT="$PARENT" REFERENCE_FAI="$TEST_DIR/reference.fai" \
  CANDIDATE_BASE="$CANDIDATE_BASE" SAMPLE_ID="$SAMPLE_ID" \
  ANALYSIS_BED_NAME="$ANALYSIS_BED_NAME" IGV_BED_NAME="$IGV_BED_NAME" \
  TARGET_BED_NAME="$TARGET_BED_NAME" WORKFLOW_UID=test-uid \
  SOURCE_REVISION=test-sha PATH="$TEST_DIR/bin:$PATH" \
  python3 "$TEST_DIR/validate_candidate.py" >/dev/null 2>&1; then
  fail "Candidate validation accepted BED intervals outside FASTA contig order"
fi
mv "$TEST_DIR/analysis.good" "$CANDIDATE/metrics/bed_validation/$ANALYSIS_BED_NAME"
env "${COMMON_MANIFEST_ENV[@]}" PATH="$TEST_DIR/bin:$PATH" \
  python3 "$TEST_DIR/write_manifest.py"

# Atomic promotion must retain a UID backup and be idempotently recoverable.
PROMOTE_PARENT="$TEST_DIR/promote"
mkdir -p "$PROMOTE_PARENT/live" "$PROMOTE_PARENT/candidate"
printf old > "$PROMOTE_PARENT/live/old.txt"
printf new > "$PROMOTE_PARENT/candidate/new.txt"
cat > "$PROMOTE_PARENT/candidate/padding_repair_manifest.json" <<'EOF_PROMOTE_MANIFEST'
{"workflow_uid":"promote-uid","validation":{"state":"validated"},"promotion":{"state":"pending","backup_dir":null},"projection":{"state":"pending"}}
EOF_PROMOTE_MANIFEST
PROMOTE_ENV=(
  MOUNT_PARENT="$PROMOTE_PARENT"
  LIVE_BASE=live
  CANDIDATE_BASE=candidate
  OUTPUT_PARENT_HOST=/host/output/
  WORKFLOW_UID=promote-uid
  BACKUP_OUTPUT_PATH="$TEST_DIR/backup_output"
)
env "${PROMOTE_ENV[@]}" PATH="$TEST_DIR/bin:$PATH" python3 "$TEST_DIR/promote.py"
[[ -f $PROMOTE_PARENT/live/new.txt ]] || fail "New candidate was not promoted"
[[ -f $PROMOTE_PARENT/live.pre-padding-repair-promote-uid/old.txt ]] ||
  fail "Old live tree was not retained as the UID backup"
env "${PROMOTE_ENV[@]}" PATH="$TEST_DIR/bin:$PATH" python3 "$TEST_DIR/promote.py"

# Finalization must record Errored publication and clean only successful UID work.
mkdir -p "$PROMOTE_PARENT/live.repair-wk-promote-uid"
env \
  MOUNT_PARENT="$PROMOTE_PARENT" \
  LIVE_BASE=live \
  CANDIDATE_BASE=candidate \
  TX_BASE=live.repair-wk-promote-uid \
  OUTPUT_PARENT_HOST=/host/output/ \
  WORKFLOW_UID=promote-uid \
  REPAIR_MODE=apply \
  VALIDATION_STATUS=Succeeded \
  PROMOTE_STATUS=Succeeded \
  PUBLISH_STATUS=Errored \
  PATH="$TEST_DIR/bin:$PATH" \
  python3 "$TEST_DIR/finalize.py"
[[ ! -e $PROMOTE_PARENT/live.repair-wk-promote-uid ]] ||
  fail "Successful UID transaction work was not cleaned"
PATH="$TEST_DIR/bin:$PATH" python3 - "$PROMOTE_PARENT/live/padding_repair_manifest.json" <<'PY'
import json
import sys

with open(sys.argv[1]) as handle:
    manifest = json.load(handle)
assert manifest["validation"]["state"] == "validated"
assert manifest["promotion"]["state"] == "completed"
assert manifest["projection"]["state"] == "failed"
assert manifest["terminal_tasks"]["tmb_publication"] == "Errored"
PY

# Publish-only code must post the staged payload with the runtime token unchanged.
cat > "$TEST_DIR/server.py" <<'PY_SERVER'
import pathlib
import sys
from http.server import BaseHTTPRequestHandler, HTTPServer

port_file = pathlib.Path(sys.argv[1])
body_file = pathlib.Path(sys.argv[2])

class Handler(BaseHTTPRequestHandler):
    def do_POST(self):
        length = int(self.headers["Content-Length"])
        body_file.write_bytes(self.rfile.read(length))
        self.send_response(200)
        self.end_headers()
        self.wfile.write(b"ok")

    def log_message(self, *_):
        pass

server = HTTPServer(("0.0.0.0", 0), Handler)
port_file.write_text(str(server.server_port))
server.handle_request()
PY_SERVER
python3 "$TEST_DIR/server.py" "$TEST_DIR/server.port" "$TEST_DIR/posted.json" &
SERVER_PID=$!
for _ in $(seq 1 50); do
  [[ -s $TEST_DIR/server.port ]] && break
  sleep 0.1
done
[[ -s $TEST_DIR/server.port ]] || fail "Local publication test server did not start"
PORT=$(<"$TEST_DIR/server.port")
REQUEST_JSON="$CANDIDATE/metrics/tmb/tmb_api_payload.json" \
  API_ENDPOINT="http://host.docker.internal:${PORT}/somatic" \
  GENNEXT_CLUSTER_TOKEN=runtime-secret \
  PATH="$TEST_DIR/bin:$PATH" \
  python3 "$TEST_DIR/publish.py"
wait "$SERVER_PID"
unset SERVER_PID
PATH="$TEST_DIR/bin:$PATH" python3 - "$TEST_DIR/posted.json" <<'PY'
import json
import sys

with open(sys.argv[1]) as handle:
    payload = json.load(handle)
assert payload["token"] == "runtime-secret"
assert payload["somaticData"]["tmbValue"] == "2.50"
assert payload["somaticData"]["tmbVariantCount"] == 5
PY

echo "All padding repair workflow tests passed."
