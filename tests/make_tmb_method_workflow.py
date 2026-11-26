#!/usr/bin/env python3
"""Writes a one-off Argo Workflow that runs one TMB method row on one real sample.

For the checks before a method becomes a default (deploy plan V3, V5b, V6):
what the method computes on a real sample, the label it would get, the exact
arguments that ran, and the pod's time and peak memory. Nothing is posted.

The Workflow is self-contained; it references no WorkflowTemplate:

  v2       the caller task of somatic-dna-<input-type>-v2 (its manifest
           expressions verbatim) -> collect-tmb-data-v2 without its publish
           step -> compute-tmb with its script verbatim. The manifest is a
           workflow argument built from the method row, with the effective
           arguments for the input type (the backend's input_overrides rule).
           One change: compute-tmb's command runs the script under a small
           Python wrapper that records wall time, CPU time, the calculator's
           peak RSS (getrusage of waited-for children) and the container
           cgroup's memory.peak where the kernel has it, in v6.json.
  report   reads the staged payload and v6.json and prints tmbValue,
           tmbTargetRegionSize, the CDS figures, the label by the method's
           thresholds (as TmbStatus::derive does) and the effective args.

It writes only under --scratch-dir. The DuckDB, the BED and the reference
directory are mounted read-only. No step has the cluster token or a backend URL.

The method file is the admin API's create body: calculatorImage, cliArgs,
referenceFiles, targetSource, denominator, inputOverrides, and optionally
label, highThreshold (default 10.0) and minRegionSizeBp (default 1300000).

Usage:
  tests/make_tmb_method_workflow.py --method platform-v2.json --node NODE \\
      --run-id RUN --sample-id SAMPLE --duckdb /path/<sample>_results.duckdb \\
      --target-bed /path/<bed>.valid.bed --scratch-dir /runspace/.../scratch/x/ \\
      [--input-type fastq|bam|vcf] [--ref-base-dir /runspace/gennext/data/ref/] > wf.yaml
"""

import argparse
import copy
import json
import pathlib
import re
import sys

import yaml

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from make_tmb_v0_workflow import COMPARE_IMAGE, caller_task, rebase, templates  # noqa: E402

REMOVABLE = ("--min-depth", "--depth-on-missing")  # the backend's REMOVABLE_OPTIONS

MEASURE = r"""
import json, os, resource, subprocess, sys, time
start = time.monotonic()
rc = subprocess.call([sys.executable, sys.argv[1]])
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
peak = None
for path in ("/sys/fs/cgroup/memory.peak", "/sys/fs/cgroup/memory/memory.max_usage_in_bytes"):
    try:
        with open(path) as handle:
            peak = int(handle.read())
        break
    except (OSError, ValueError):
        pass
v6 = {"exit": rc, "wall_s": round(time.monotonic() - start, 1),
      "user_s": round(usage.ru_utime, 1), "sys_s": round(usage.ru_stime, 1),
      "calculator_max_rss_mib": round(usage.ru_maxrss / 1024, 1),
      "cgroup_memory_peak_mib": None if peak is None else round(peak / 2**20, 1)}
print("V6 " + json.dumps(v6), flush=True)
try:
    with open("/mnt/output/v6.json", "w") as handle:
        json.dump(v6, handle)
except OSError as err:
    print(f"could not write v6.json: {err}", file=sys.stderr)
sys.exit(rc)
"""

REPORT = r"""
import json, os, sys
method = json.loads(os.environ["METHOD"])
effective = json.loads(os.environ["EFFECTIVE_ARGS"])
with open("/mnt/out/tmb_api_payload.json") as handle:
    sd = json.load(handle)["somaticData"]
v6 = {}
try:
    with open("/mnt/out/v6.json") as handle:
        v6 = json.load(handle)
except OSError:
    pass
value, region = float(sd["tmbValue"]), int(sd["tmbTargetRegionSize"])
high, floor = float(method["highThreshold"]), int(method["minRegionSizeBp"])
if region < floor:
    kind = "coding" if method["denominator"] == "cds" else "target"
    label = f"Undeterminable ({kind} region {region} bp < {floor} bp)"
elif value >= high:
    label = f"TMB-H (>= {high:g} mut/Mb)"
else:
    label = f"TMB-L (< {high:g} mut/Mb)"
out = {
    "sample": os.environ["SAMPLE_ID"], "inputType": os.environ["INPUT_TYPE"],
    "method": method.get("label"), "targetSource": method["targetSource"],
    "denominator": method["denominator"],
    "tmbValue": sd["tmbValue"], "tmbVariantCount": sd["tmbVariantCount"],
    "tmbTargetRegionSize": region,
    "cdsTmbValue": sd.get("cdsTmbValue"), "cdsTmbTargetRegionSize": sd.get("cdsTmbTargetRegionSize"),
    "label": label, "effectiveArgs": effective, "v6": v6,
}
print("RESULT " + json.dumps(out))
print(json.dumps(out, indent=1))
if "--min-depth" in effective and "--min-depth" in method.get("inputOverrides", {}).get(os.environ["INPUT_TYPE"], {}).get("remove", []):
    print("V3 FAIL: --min-depth is still in the effective args", file=sys.stderr)
    sys.exit(1)
"""


def effective_args(method, input_type):
    """The backend's effective_args: cli_args with that input type's removals."""
    args = list(method["cliArgs"])
    removals = method.get("inputOverrides", {}).get(input_type, {}).get("remove", [])
    for option in removals:
        if option not in REMOVABLE:
            raise SystemExit(f"{option} cannot be removed by an input override")
        kept, i = [], 0
        while i < len(args):
            if args[i] == option:
                i += 1
                if i < len(args) and not args[i].startswith("--"):
                    i += 1
                continue
            kept.append(args[i])
            i += 1
        args = kept
    return args


def build(args, method):
    for key in ("calculatorImage", "cliArgs", "referenceFiles", "targetSource", "denominator"):
        if key not in method:
            raise SystemExit(f"method file: {key} is missing")
    method = dict(method)
    method.setdefault("inputOverrides", {})
    method.setdefault("highThreshold", 10.0)
    method.setdefault("minRegionSizeBp", 1300000)
    if not re.fullmatch(r"[a-z0-9./_-]+@sha256:[0-9a-f]{64}", method["calculatorImage"]):
        raise SystemExit("calculatorImage must be <repository>@sha256:<64 hex>")
    effective = effective_args(method, args.input_type)
    sample = {
        "tmbMethodId": args.method_id,
        "tmbMethodScope": "kit",
        "tmbMethodVersion": 2,
        "tmbCalculatorImage": method["calculatorImage"],
        "tmbInputType": args.input_type,
        "tmbMethodArgs": json.dumps(effective),
        "tmbReferenceFiles": json.dumps(method["referenceFiles"]),
        "tmbTargetSource": method["targetSource"],
    }

    v2 = templates("tasks/collect_tmb_data_v2.yaml")
    wrapper = copy.deepcopy(v2["collect-tmb-data-v2"])
    wrapper["steps"] = [[s for s in g if s["name"] != "publish-tmb"] for g in wrapper["steps"]]
    wrapper["steps"] = [g for g in wrapper["steps"] if g]
    if [s["name"] for g in wrapper["steps"] for s in g] != ["compute-tmb"]:
        raise SystemExit("collect-tmb-data-v2: unexpected steps")
    wrapper["inputs"]["parameters"] = [p for p in wrapper["inputs"]["parameters"]
                                       if p["name"] != "api-base-url"]
    compute = copy.deepcopy(v2["compute-tmb"])
    if compute["script"]["command"] != ["python"]:
        raise SystemExit("compute-tmb: unexpected command")
    compute["script"]["command"] = ["python", "-c", MEASURE]
    task = caller_task(args.input_type)
    if task["templateRef"] != {"name": "collect-tmb-data-v2", "template": "collect-tmb-data-v2"}:
        raise SystemExit("the caller does not call collect-tmb-data-v2")

    scratch = args.scratch_dir.rstrip("/") + "/"
    out_dir = scratch + args.input_type
    report = {
        "name": "report",
        "volumes": [{"name": "out", "hostPath": {"path": out_dir + "/metrics/tmb/",
                                                 "type": "Directory"}}],
        "script": {
            "image": COMPARE_IMAGE,
            "command": ["python"],
            "env": [{"name": "METHOD", "value": json.dumps(method)},
                    {"name": "EFFECTIVE_ARGS", "value": json.dumps(effective)},
                    {"name": "SAMPLE_ID", "value": args.sample_id},
                    {"name": "INPUT_TYPE", "value": args.input_type}],
            "volumeMounts": [{"name": "out", "mountPath": "/mnt/out", "readOnly": True}],
            "source": REPORT,
        },
    }
    inputs = ["run-id", "sample-id", "sample-run-data", "ref-base-dir",
              "duckdb-file", "target-bed-file"]
    main = {
        "name": "main",
        "inputs": {"parameters": [{"name": n} for n in inputs]},
        "dag": {"tasks": [
            {"name": "v2", "template": "collect-tmb-data-v2",
             "arguments": {"parameters": [
                 {"name": p["name"], "value": rebase(p["value"], out_dir)}
                 for p in task["arguments"]["parameters"]]}},
            {"name": "report", "template": "report", "depends": "v2"},
        ]},
    }
    values = {"run-id": args.run_id, "sample-id": args.sample_id,
              "sample-run-data": json.dumps({"samples": [sample]}),
              "ref-base-dir": args.ref_base_dir,
              "duckdb-file": args.duckdb, "target-bed-file": args.target_bed}
    return {
        "apiVersion": "argoproj.io/v1alpha1",
        "kind": "Workflow",
        "metadata": {"generateName": args.generate_name, "namespace": "argo",
                     "labels": {"gennext.bio/purpose": "tmb-method-check"}},
        "spec": {
            "serviceAccountName": "gennext-workflows",
            "entrypoint": "main",
            "nodeSelector": {"kubernetes.io/hostname": args.node},
            "arguments": {"parameters": [{"name": k, "value": v} for k, v in values.items()]},
            "templates": [main, wrapper, compute, report],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--method", required=True, type=pathlib.Path)
    parser.add_argument("--node", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--duckdb", required=True)
    parser.add_argument("--target-bed", required=True)
    parser.add_argument("--scratch-dir", required=True)
    parser.add_argument("--input-type", default="fastq", choices=["fastq", "bam", "vcf"])
    parser.add_argument("--ref-base-dir", default="/runspace/gennext/data/ref/")
    parser.add_argument("--method-id", default="00000000-0000-4000-8000-00000000000b",
                        help="any non-empty id: compute only requires one, nothing posts")
    parser.add_argument("--generate-name", default="tmb-method-")
    args = parser.parse_args()
    method = json.loads(args.method.read_text())
    yaml.safe_dump(build(args, method), sys.stdout, sort_keys=False, width=1000)


if __name__ == "__main__":
    main()
