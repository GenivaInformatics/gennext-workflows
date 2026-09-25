#!/usr/bin/env bash
set -euo pipefail
repo_dir=$(cd "$(dirname "$0")/../.." && pwd)
kubectl --context default -n argo create configmap wgs-cnv-validation-code \
  --from-file="$repo_dir/images/cnv-cohort-builder/cnv_cohort_builder.py" \
  --from-file="$repo_dir/images/cnv-cohort-builder/test_wgs.py" \
  --from-file="$repo_dir/images/cnv-cohort-builder/smoke_wgs.py" \
  --from-file="$repo_dir/tests/wgs-cnv/validate_node.py" \
  --dry-run=client -o yaml | kubectl --context default apply -f -
