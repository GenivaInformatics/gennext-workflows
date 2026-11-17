#!/usr/bin/env bash
# Steps that share a node must declare the host memory they use, or the
# scheduler can place two of them on a node that can't hold both. This checks
# the steps sized so far:
#   - each listed template requests memory;
#   - every fq2bam variant pins --memory-limit (Parabricks otherwise sorts in
#     half of the host's RAM, whatever the pod requests) and requests at least
#     that limit plus 8 GiB for the BWA index, work queues and the process.
#
# Needs python3 with PyYAML.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
cd "$SCRIPT_DIR/.."

python3 - <<'PY'
import re
import sys

import yaml

STEPS = {
    "tasks/base_tasks/fastqc.yaml": ["fastqc"],
    "tasks/base_tasks/fastp.yaml": ["fastp"],
    "tasks/base_tasks/star-align.yaml": ["star-align-arriba"],
    "tasks/base_tasks/gpu_tasks/parabricks-fq2bam.yaml": [
        "parabricks-fq2bamfast-no-bqsr",
        "parabricks-fq2bam-no-bqsr",
        "parabricks-fq2bam-no-bqsr-on-node",
        "parabricks-fq2bamfast-no-bqsr-on-node",
    ],
    "tasks/base_tasks/gpu_tasks/parabricks-rnafq2bam.yaml": ["parabricks-rnafq2bam"],
}
UNITS = {"Mi": 1 / 1024, "Gi": 1}


def gib(quantity):
    m = re.fullmatch(r"(\d+)(Mi|Gi)", str(quantity))
    if not m:
        sys.exit(f"unexpected memory quantity {quantity!r} (use Mi or Gi)")
    return int(m.group(1)) * UNITS[m.group(2)]


failures = []
for path, names in STEPS.items():
    templates = {t["name"]: t for t in yaml.safe_load(open(path))["spec"]["templates"]}
    for name in names:
        if name not in templates:
            failures.append(f"{path}: no template {name}")
            continue
        body = templates[name].get("script") or templates[name].get("container")
        memory = ((body.get("resources") or {}).get("requests") or {}).get("memory")
        if not memory:
            failures.append(f"{name}: no memory request")
            continue
        if "fq2bam" in name and "rnafq2bam" not in name:
            limits = re.findall(r"--memory-limit (\d+)", body["source"])
            if len(limits) != 1:
                failures.append(f"{name}: expected one --memory-limit, found {limits}")
            elif gib(memory) < int(limits[0]) + 8:
                failures.append(f"{name}: requests {memory}, less than --memory-limit {limits[0]} GiB + 8 GiB")
        print(f"ok  {name}: requests {memory}")

if failures:
    print("\n".join("FAIL: " + f for f in failures), file=sys.stderr)
    sys.exit(1)
print("PASS: sized steps request their host memory")
PY
