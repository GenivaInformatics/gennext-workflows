#!/bin/bash
# Test script for the batch orchestrator prototype.
# Deploys the templates and submits a test workflow to gennext-dev namespace.
#
# Usage: ./test-batch-prototype.sh [max-concurrent]
#   max-concurrent: number of parallel samples (default: 2)

set -euo pipefail

NAMESPACE="gennext-dev"
MAX_CONCURRENT="${1:-2}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== Deploying WorkflowTemplates to ${NAMESPACE} ==="
kubectl apply -f "${SCRIPT_DIR}/sample-worker-prototype.yaml" -n "${NAMESPACE}"
kubectl apply -f "${SCRIPT_DIR}/batch-orchestrator-prototype.yaml" -n "${NAMESPACE}"

echo ""
echo "=== Submitting test batch workflow ==="
echo "Samples: 3 (TUMOR-A, TUMOR-B, TUMOR-C)"
echo "Max concurrent: ${MAX_CONCURRENT}"
echo ""

SAMPLES='[{"id":"aaaa-1111-uuid","sampleId":"TUMOR-A"},{"id":"bbbb-2222-uuid","sampleId":"TUMOR-B"},{"id":"cccc-3333-uuid","sampleId":"TUMOR-C"}]'

argo submit --from workflowtemplate/batch-orchestrator-prototype \
  -n "${NAMESPACE}" \
  -p "samples=${SAMPLES}" \
  -p "max-concurrent=${MAX_CONCURRENT}" \
  --watch

echo ""
echo "=== Child workflows spawned by orchestrator ==="
argo list -n "${NAMESPACE}" -l "gennext.bio/orchestrator" --no-headers
