#!/usr/bin/env python3
"""Writes a one-off Argo Workflow for V0 on one real sample (nothing posts).

V0 in the deploy plan: on a few existing samples, collect-tmb-data-v2 with the
kit's v1 must stage a payload byte-identical to what today's command writes.
tests/test_tmb_v2_byte_identity.py proves this on synthetic data, with Argo's
substitution emulated. This runs it in Argo itself, so the quoted JSON in
tmbMethodArgs and tmbReferenceFiles goes through real Argo, from the
expression to the parameter to the env value.

The Workflow is self-contained; it references no WorkflowTemplate:

  today    collect-tmb-data from tasks/collect_tmb_data.yaml, verbatim, with
           the callers' arguments as before this change. Only two things
           differ: --api-endpoint is the test:// path, and the token is the
           literal STAGED (no secret is mounted).
  v2       the branch's caller task (somatic-dna-<input-type>-v2) ->
           collect-tmb-data-v2 -> compute-tmb, all verbatim from this
           checkout, except that the wrapper's publish-tmb step is removed.
           The manifest is a workflow argument carrying the v1 method row
           from tests/fixtures/tmb-v2/methods.json (real checksums).
  compare  byte comparison of the two staged payloads; fails if they differ.

Both write only under --scratch-dir (today/ and v2/). The DuckDB, the BED and
the reference directory are mounted read-only. No step has the cluster token
or a backend URL.

Usage:
  tests/make_tmb_v0_workflow.py --node NODE --run-id RUN --sample-id SAMPLE \\
      --duckdb /path/<sample>_results.duckdb --target-bed /path/<bed>.valid.bed \\
      --scratch-dir /runspace/.../scratch/tmb-v0-<sample>/ \\
      [--input-type fastq|bam|vcf] [--target-source kit_bed|cds_bed] \\
      [--ref-base-dir /runspace/gennext/data/ref/] [--method-id UUID] > wf.yaml
  argo submit -n argo wf.yaml --watch
"""

import argparse
import copy
import json
import pathlib
import sys

import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent
FIX = REPO / "tests/fixtures/tmb-v2"
PROD_SOMATIC = ("http://gennext-backend.gennext-prod.svc.cluster.local:3030"
                "/api/v0/runs/{{inputs.parameters.run-id}}/somatic")
STAGED_PATH = "test:///mnt/output/tmb_api_payload.json"
COMPARE_IMAGE = ("python:3.11-alpine3.18@sha256:"
                 "4799cdc47a1b7ecc2829e71e74d06ce7c0af54dec4f5ff139f869cc2c4283b1d")


def templates(path):
    out = {}
    for doc in yaml.safe_load_all((REPO / path).read_text()):
        for tmpl in (doc or {}).get("spec", {}).get("templates", []):
            out[tmpl["name"]] = tmpl
    return out


def caller_task(input_type):
    path = f"workflows/somatic/v2/somatic-dna-{input_type}.yaml"
    for tmpl in templates(path).values():
        for task in (tmpl.get("dag") or {}).get("tasks", []):
            if task["name"] == "collect-tmb-data":
                return task
    raise SystemExit(f"{path}: no collect-tmb-data task")


def rebase(value, out_dir):
    """A caller argument, with its upstream task outputs replaced by inputs."""
    return (value
            .replace("{{tasks.duckdb-snv.outputs.parameters.output-duckdb}}",
                     "{{inputs.parameters.duckdb-file}}")
            .replace("{{tasks.validate-bed-file.outputs.parameters.validated-bed-file}}",
                     "{{inputs.parameters.target-bed-file}}")
            .replace("{{inputs.parameters.sample-temp-dir}}", out_dir))


def build(args):
    method = json.loads((FIX / "methods.json").read_text())["v1"]
    sample = {
        "tmbMethodId": args.method_id,
        "tmbMethodScope": "kit",
        "tmbMethodVersion": 1,
        "tmbCalculatorImage": method["calculator_image"],
        "tmbInputType": args.input_type,
        # v1 has no input overrides: its effective args are its cli_args.
        "tmbMethodArgs": json.dumps(method["cli_args"]),
        "tmbReferenceFiles": json.dumps(method["reference_files"]),
        "tmbTargetSource": args.target_source,
    }
    manifest = json.dumps({"samples": [sample]})

    # today: collect-tmb-data verbatim, endpoint and token staged.
    today = copy.deepcopy(templates("tasks/collect_tmb_data.yaml")["collect-tmb-data"])
    source = today["script"]["source"]
    if source.count(PROD_SOMATIC) != 1:
        raise SystemExit("collect-tmb-data: the prod endpoint is not where expected")
    today["script"]["source"] = source.replace(PROD_SOMATIC, STAGED_PATH)
    today["script"]["env"] = [{"name": "GENNEXT_CLUSTER_TOKEN", "value": "STAGED"}]
    today_args = json.loads((FIX / "today-callers.json").read_text())["callers"][args.input_type]

    # v2: wrapper without publish, compute-tmb verbatim.
    v2 = templates("tasks/collect_tmb_data_v2.yaml")
    wrapper = copy.deepcopy(v2["collect-tmb-data-v2"])
    wrapper["steps"] = [[step for step in group if step["name"] != "publish-tmb"]
                        for group in wrapper["steps"]]
    wrapper["steps"] = [group for group in wrapper["steps"] if group]
    if [s["name"] for g in wrapper["steps"] for s in g] != ["compute-tmb"]:
        raise SystemExit("collect-tmb-data-v2: unexpected steps")
    # Only publish used the backend URL.
    wrapper["inputs"]["parameters"] = [p for p in wrapper["inputs"]["parameters"]
                                       if p["name"] != "api-base-url"]
    compute = copy.deepcopy(v2["compute-tmb"])
    task = caller_task(args.input_type)
    if task["templateRef"] != {"name": "collect-tmb-data-v2", "template": "collect-tmb-data-v2"}:
        raise SystemExit("the caller does not call collect-tmb-data-v2")

    scratch = args.scratch_dir.rstrip("/") + "/"
    compare = {
        "name": "compare",
        "volumes": [{"name": "scratch", "hostPath": {"path": scratch, "type": "Directory"}}],
        "script": {
            "image": COMPARE_IMAGE,
            "command": ["python"],
            "volumeMounts": [{"name": "scratch", "mountPath": "/mnt/scratch", "readOnly": True}],
            "source": (
                "import hashlib, json, sys\n"
                "a = open('/mnt/scratch/today/metrics/tmb/tmb_api_payload.json', 'rb').read()\n"
                "b = open('/mnt/scratch/v2/metrics/tmb/tmb_api_payload.json', 'rb').read()\n"
                "for name, data in (('today', a), ('v2', b)):\n"
                "    print(name, hashlib.sha256(data).hexdigest(), len(data), 'bytes')\n"
                "print(json.dumps(json.loads(a)['somaticData'], indent=1)[:2000])\n"
                "if a != b:\n"
                "    print('V0 FAIL: the staged payloads differ', file=sys.stderr)\n"
                "    sys.exit(1)\n"
                "print('V0 PASS: byte-identical')\n"
            ),
        },
    }

    inputs = ["run-id", "sample-id", "sample-run-data", "ref-base-dir",
              "duckdb-file", "target-bed-file"]
    main = {
        "name": "main",
        "inputs": {"parameters": [{"name": n} for n in inputs]},
        "dag": {"tasks": [
            {"name": "today", "template": "collect-tmb-data",
             "arguments": {"parameters": [
                 {"name": p["name"], "value": rebase(p["value"], scratch + "today")}
                 for p in today_args["parameters"]]}},
            {"name": "v2", "template": "collect-tmb-data-v2",
             "arguments": {"parameters": [
                 {"name": p["name"], "value": rebase(p["value"], scratch + "v2")}
                 for p in task["arguments"]["parameters"]]}},
            {"name": "compare", "template": "compare", "depends": "today && v2"},
        ]},
    }
    values = {"run-id": args.run_id, "sample-id": args.sample_id,
              "sample-run-data": manifest, "ref-base-dir": args.ref_base_dir,
              "duckdb-file": args.duckdb, "target-bed-file": args.target_bed}
    return {
        "apiVersion": "argoproj.io/v1alpha1",
        "kind": "Workflow",
        "metadata": {"generateName": "tmb-v0-", "namespace": "argo",
                     "labels": {"gennext.bio/purpose": "tmb-v0-check"}},
        "spec": {
            "serviceAccountName": "gennext-workflows",
            "entrypoint": "main",
            "nodeSelector": {"kubernetes.io/hostname": args.node},
            "arguments": {"parameters": [{"name": k, "value": v} for k, v in values.items()]},
            "templates": [main, today, wrapper, compute, compare],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--node", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--duckdb", required=True)
    parser.add_argument("--target-bed", required=True)
    parser.add_argument("--scratch-dir", required=True)
    parser.add_argument("--input-type", default="fastq", choices=["fastq", "bam", "vcf"])
    parser.add_argument("--target-source", default="kit_bed", choices=["kit_bed", "cds_bed"])
    parser.add_argument("--ref-base-dir", default="/runspace/gennext/data/ref/")
    parser.add_argument("--method-id", default="00000000-0000-4000-8000-00000000000a",
                        help="any non-empty id: compute only requires one, nothing posts")
    args = parser.parse_args()
    if args.target_source == "cds_bed" and "hg38_canonical_chrom_contig_regions" not in args.target_bed:
        print("warning: cds_bed is the Genome kit's v1 only", file=sys.stderr)
    yaml.safe_dump(build(args), sys.stdout, sort_keys=False, width=1000)


if __name__ == "__main__":
    main()
