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

# ── 1. Mutect2 uses the analysis BED (not the original target BED) ───
echo "Checking Mutect2 consumes the analysis BED..."
mutect2=$(task_block mutect2 "$REPAIR_WORKFLOW")
grep -A1 'name: regions-file' <<< "$mutect2" | grep -q 'validated-bed-analysis-file' ||
  fail "Mutect2 does not consume the analysis BED"

# ── 2. bcftools uses the analysis BED ────────────────────────────────
echo "Checking bcftools-filter consumes the analysis BED..."
bcftools=$(task_block bcftools-filter "$REPAIR_WORKFLOW")
grep -A1 'name: regions-file' <<< "$bcftools" | grep -q 'validated-bed-analysis-file' ||
  fail "bcftools-filter does not consume the analysis BED"

# ── 3. TMB publish uses the original unpadded target BED ─────────────
echo "Checking TMB uses the original target BED..."
tmb=$(task_block publish-tmb "$REPAIR_WORKFLOW")
grep -A1 'name: target-bed-file' <<< "$tmb" | grep -q 'validated-bed-file' ||
  fail "TMB does not use the original target BED"
if grep -qE 'validated-bed-(analysis|igv)-file' <<< "$tmb"; then
  fail "TMB incorrectly consumes a padded or presentation BED"
fi

# ── 4. TMB pre-compute also uses the original target BED ─────────────
echo "Checking TMB pre-compute uses the original target BED..."
compute_tmb=$(task_block compute-tmb "$REPAIR_WORKFLOW")
grep -A1 'name: target-bed-file' <<< "$compute_tmb" | grep -q 'validated-bed-file' ||
  fail "TMB pre-compute does not use the original target BED"
if grep -qE 'validated-bed-(analysis|igv)-file' <<< "$compute_tmb"; then
  fail "TMB pre-compute incorrectly consumes a padded or presentation BED"
fi

# ── 5. Same-sample mutex is configured ───────────────────────────────
echo "Checking workflow mutex..."
grep -q 'name: "padding-repair-{{workflow.parameters.sample-uuid}}"' "$REPAIR_WORKFLOW" ||
  fail "Workflow mutex is not keyed on sample-uuid"

# ── 6. Global MAF importer (maf-to-db) is NOT called ────────────────
echo "Checking maf-to-db template is not referenced..."
if grep -q 'template: maf-to-db' "$REPAIR_WORKFLOW" || grep -q 'name: maf-to-db-template' "$REPAIR_WORKFLOW"; then
  fail "Repair workflow must NOT call the global MAF importer (maf-to-db)"
fi

# ── 7. validate mode does NOT call promote or publish ────────────────
echo "Checking validate mode gates..."
promote=$(task_block promote-and-backup "$REPAIR_WORKFLOW")
grep -q "mode.*==.*apply" <<< "$promote" ||
  fail "promote-and-backup is not gated on mode=apply"
publish=$(task_block publish-tmb "$REPAIR_WORKFLOW")
grep -q "mode.*==.*apply" <<< "$publish" ||
  fail "publish-tmb is not gated on mode=apply"

# ── 8. DuckDB refresh only touches SNV objects ───────────────────────
echo "Checking DuckDB SNV refresh does not drop non-SNV objects..."
# The repair uses the existing duckdb-snv template which only DROPs
# oncokb_snv, cancervar_snv, combined_snv, advanced_combined_snv.
# Verify the repair workflow does NOT call duckdb-sv, duckdb-cnv, or
# duckdb-fusions.
if grep -q 'template: duckdb-sv\b' "$REPAIR_WORKFLOW"; then
  fail "Repair workflow must not call duckdb-sv"
fi
if grep -q 'template: duckdb-cnv\b' "$REPAIR_WORKFLOW"; then
  fail "Repair workflow must not call duckdb-cnv"
fi
if grep -q 'template: duckdb-fusions\b' "$REPAIR_WORKFLOW"; then
  fail "Repair workflow must not call duckdb-fusions"
fi

# ── 9. validate-bed-file generates padded artifacts ──────────────────
echo "Checking BED validator generates padded artifacts..."
bed_task=$(task_block validate-bed-file "$REPAIR_WORKFLOW")
grep -A1 'name: generate-padded-artifacts' <<< "$bed_task" | grep -q 'value: "true"' ||
  fail "validate-bed-file does not opt into padded artifact generation"
grep -A1 'name: padding' <<< "$bed_task" | grep -q 'value: "100"' ||
  fail "validate-bed-file padding is not 100bp"
grep -A1 'name: reference-fai-file' <<< "$bed_task" | grep -q '\.fai"' ||
  fail "validate-bed-file does not provide the reference FASTA index"

# ── 10. VCF-only samples are rejected ────────────────────────────────
echo "Checking VCF-only sample rejection..."
grep -q 'VCF-only samples cannot be repaired' "$REPAIR_WORKFLOW" ||
  fail "Repair workflow does not reject VCF-only samples"

# ── 11. Promotion creates a timestamped backup ──────────────────────
echo "Checking promotion creates timestamped backup..."
grep -q 'pre-padding-repair' "$REPAIR_WORKFLOW" ||
  fail "Promotion does not create a timestamped backup"

# ── 12. Promotion uses atomic rename (not rsync) ────────────────────
echo "Checking promotion uses rename, not rsync..."
if grep -q 'rsync' "$REPAIR_WORKFLOW"; then
  fail "Repair workflow uses rsync for promotion (must use rename)"
fi

# ── 13. DuckDB staged copy preserves SV/CNV (copy, not create new) ──
echo "Checking DuckDB is staged via copy..."
# The stage-duckdb-copy template definition should use cp, not create a fresh db
sed -n '/^    - name: stage-duckdb-copy$/,/^    - name: /p' "$REPAIR_WORKFLOW" | grep -q 'cp ' ||
  fail "stage-duckdb-copy does not copy the existing DuckDB (SV/CNV would be lost)"

# ── 14. Node selector is configured ─────────────────────────────────
echo "Checking node selector..."
grep -q 'kubernetes.io/hostname: "{{workflow.parameters.node-name}}"' "$REPAIR_WORKFLOW" ||
  fail "Workflow does not pin to the specified node"

echo ""
echo "All padding repair regression tests passed"
