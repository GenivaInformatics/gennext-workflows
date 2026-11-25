#!/usr/bin/env python3
"""collect-tmb-data-v2 with v1 methods stages byte-identical TMB payloads (V0).

Today, somatic-dna-{fastq,bam,vcf}-v2 call collect-tmb-data with fixed
arguments. With versioning they call collect-tmb-data-v2, which runs whatever
method the run manifest names. Before the cutover every kit is on its v1, so
for v1 the new path must stage exactly the bytes the old one would have
posted. This runs both from the templates, in the real images, on synthetic
fixtures (tests/fixtures/tmb-v2):

  today   collect-tmb-data's script (tasks/collect_tmb_data.yaml) with the
          callers' arguments as they were before this change
          (fixtures/tmb-v2/today-callers.json), --api-endpoint swapped for the
          same test:// path and GENNEXT_CLUSTER_TOKEN=STAGED (both stage the
          token as STAGED; publish replaces it);
  v2      the callers in this branch -> collect-tmb-data-v2 -> compute-tmb,
          with a manifest that carries the v1 method (methods.json).

Checks:
  V0   for fastq, bam and vcf callers: the two payloads are byte-identical;
       also for the Genome kit (today: BED filename -> CDS target; v2:
       tmbTargetSource cds_bed);
  V0b  somatic-dna-padding-repair-v1's compute-tmb-data stages the same bytes
       as v2 with v1 (its posts are recorded as v1);
  v2   the 0.6.0 method runs through the same template; vcf input with the
       effective args (no --min-depth) keeps the low-depth variant;
  the compute stage fails on a missing or wrong reference file, an empty
  method id, an unpinned image or a method that sets a fixed option, and
  never leaves a stale payload; arguments are never shell-interpreted;
  publish posts the staged payload plus tmbMethodId, tmbInputType and
  tmbEffectiveArgs with the cluster token, and fails on a mismatch or non-200.

hostPath volumes are emulated under a temporary "host" root, so every mount
is the template's own. Needs python3 with PyYAML, and Docker (pulls the
pinned images on first run).
"""

import hashlib
import http.server
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time

import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent
FIX = REPO / "tests/fixtures/tmb-v2"
METHODS = json.loads((FIX / "methods.json").read_text())
V1_IMAGE = METHODS["v1"]["calculator_image"]
V2_IMAGE = METHODS["v2"]["calculator_image"]
PROD_SOMATIC = "http://gennext-backend.gennext-prod.svc.cluster.local:3030/api/v0/runs/{run}/somatic"
STAGED_PATH = "test:///mnt/output/tmb_api_payload.json"
RUN_ID = "11111111-2222-4333-8444-555555555555"
SAMPLE_ID = "TMBV2TEST01"
METHOD_ID = "0b7f2c9e-1d3a-4c55-9e1f-6a2b3c4d5e6f"
KIT_BP = sum(int(f[2]) - int(f[1]) for f in
             (line.split("\t") for line in (FIX / "kit.bed").read_text().splitlines()))
FAILURES = []


def check(ok, what):
    print(("ok   " if ok else "FAIL ") + what, flush=True)
    if not ok:
        FAILURES.append(what)


def templates(path):
    out = {}
    for doc in yaml.safe_load_all((REPO / path).read_text()):
        for tmpl in (doc or {}).get("spec", {}).get("templates", []):
            out[tmpl["name"]] = tmpl
    return out


def dag_task(path, name):
    for tmpl in templates(path).values():
        for task in (tmpl.get("dag") or {}).get("tasks", []):
            if task["name"] == name:
                return task
    raise SystemExit(f"{path}: no task {name}")


# ── Argo substitution, for the forms these templates use ────────────────

MANIFEST_EXPR = re.compile(
    r"\{\{= sprig\.default\('', jsonpath\(inputs\.parameters\['sample-run-data'\], "
    r"'\$\.samples\[0\]'\)\['(\w+)'\]\) \}\}"
)


def go_clean(path):
    return os.path.normpath(path) if path else "."


def render(text, params, manifest=None):
    """Argo's {{...}} substitution for the forms used here; anything else fails."""
    def manifest_field(match):
        if manifest is None:
            raise SystemExit(f"manifest expression without a manifest: {match.group(0)}")
        value = manifest["samples"][0].get(match.group(1))
        return "" if value is None else str(value)

    text = MANIFEST_EXPR.sub(manifest_field, text)
    text = re.sub(
        r"\{\{=sprig\.default\('([^']*)', inputs\.parameters\['([\w-]+)'\]\)\}\}",
        lambda m: params[m.group(2)] or m.group(1), text)
    funcs = {"osBase": lambda p: os.path.basename(p.rstrip("/")),
             "osDir": lambda p: go_clean(os.path.dirname(p)),
             "osClean": go_clean}
    text = re.sub(
        r"\{\{=sprig\.(osBase|osDir|osClean)\(inputs\.parameters\['([\w-]+)'\]\)\}\}",
        lambda m: funcs[m.group(1)](params[m.group(2)]), text)
    text = re.sub(r"\{\{\s*([\w.-]+)\s*\}\}",
                  lambda m: params[m.group(1)] if m.group(1) in params else m.group(0),
                  text)
    left = re.findall(r"\{\{.*?\}\}", text)
    if left:
        raise SystemExit(f"unresolved: {left}")
    return text


def inputs_for(tmpl, args):
    """A template's inputs: its defaults, overridden by the caller's args."""
    params = {}
    for p in tmpl.get("inputs", {}).get("parameters", []):
        if "value" in p:
            params[p["name"]] = str(p["value"])
    params.update(args)
    missing = [p["name"] for p in tmpl.get("inputs", {}).get("parameters", [])
               if p["name"] not in params]
    if missing:
        raise SystemExit(f"{tmpl['name']}: inputs not supplied: {missing}")
    return {f"inputs.parameters.{k}": v for k, v in params.items()} | params


# ── A container run with the template's own volumes ────────────────────

class Host:
    """A temporary root that stands in for the node's filesystem."""

    def __init__(self, root):
        self.root = pathlib.Path(root)

    def path(self, abs_path):
        return self.root / go_clean(abs_path).lstrip("/")

    def put(self, abs_path, src):
        dest = self.path(abs_path)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if src.is_dir():
            shutil.copytree(src, dest, dirs_exist_ok=True)
        else:
            shutil.copy(src, dest)


def run_script(host, tmpl, params, secrets, image=None, network="none", extra=()):
    """Runs a script template the way Argo would; returns the CompletedProcess."""
    script = tmpl["script"]
    image = image or render(script["image"], params)
    staging = host.root / "staging" / f"{tmpl['name']}-{os.urandom(4).hex()}"
    staging.parent.mkdir(parents=True, exist_ok=True)
    staging.write_text(render(script["source"], params))
    cmd = ["docker", "run", "--rm", f"--network={network}",
           "--user", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp",
           "-v", f"{staging}:/argo/staging/script:ro"]
    if script.get("workingDir"):
        cmd += ["-w", script["workingDir"]]
    volumes = {}
    for vol in tmpl.get("volumes", []):
        hp = vol["hostPath"]
        local = host.path(render(hp["path"], params))
        kind = hp.get("type", "")
        if kind == "DirectoryOrCreate":
            local.mkdir(parents=True, exist_ok=True)
        elif kind == "Directory" and not local.is_dir():
            return subprocess.CompletedProcess(cmd, 99, "", f"hostPath {local} is not a directory")
        elif kind == "File" and not local.is_file():
            return subprocess.CompletedProcess(cmd, 99, "", f"hostPath {local} is not a file")
        volumes[vol["name"]] = local
    for mount in script.get("volumeMounts", []):
        ro = ":ro" if mount.get("readOnly") else ""
        cmd += ["-v", f"{volumes[mount['name']]}:{render(mount['mountPath'], params)}{ro}"]
    for src, dest in extra:
        cmd += ["-v", f"{src}:{dest}:ro"]
    for env in script.get("env", []):
        if "valueFrom" in env:
            value = secrets[env["valueFrom"]["secretKeyRef"]["key"]]
        else:
            value = render(env["value"], params)
        cmd += ["-e", f"{env['name']}={value}"]
    # Argo's command replaces the image's ENTRYPOINT.
    cmd += ["--entrypoint", script["command"][0], image]
    cmd += script["command"][1:] + ["/argo/staging/script"]
    started = time.monotonic()
    run = subprocess.run(cmd, capture_output=True, text=True, check=False)
    print(f"     ran {tmpl['name']} in {image.split('@')[0]}@{image.split('@sha256:')[-1][:8]} "
          f"({time.monotonic() - started:.0f} s, exit {run.returncode})", flush=True)
    return run


# ── Fixtures ──────────────────────────────────────────────────────────

def build_duckdb(host, abs_path):
    """combined_snv from variants.json, written by the calculator's own duckdb."""
    rows = json.loads((FIX / "variants.json").read_text())["rows"]
    out = host.path(abs_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    (out.parent / "rows.json").write_text(json.dumps(rows))
    code = r"""
import json, duckdb
rows = json.load(open('/w/rows.json'))
con = duckdb.connect('/w/DBNAME')
con.execute('''CREATE TABLE combined_snv ("Chrom" VARCHAR, "Start" BIGINT, "End" BIGINT,
  "Ref" VARCHAR, "Alt" VARCHAR, "Func.refGene" VARCHAR, "Gene.refGene" VARCHAR,
  "ExonicFunc.refGene" VARCHAR, "AF" DOUBLE, "Filter" VARCHAR, "CancerVar: Tier" VARCHAR,
  "DETS" DOUBLE, "DP" VARCHAR, "aaChange" VARCHAR, "txChange" VARCHAR)''')
genes = {'chr17': 'TP53', 'chr12': 'KRAS', 'chr22': 'CHEK2', 'chr3': 'PIK3CA', 'chr13': 'BRCA2'}
for i, r in enumerate(rows):
    con.execute('INSERT INTO combined_snv VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', [
        r['Chrom'], r['Start'], r['Start'] + max(len(r['Ref']), 1) - 1, r['Ref'], r['Alt'],
        r.get('Func.refGene', 'exonic'), genes[r['Chrom']], r['ExonicFunc.refGene'],
        r['AF'], r.get('Filter', 'PASS'), 'Tier III Unknown', r['DETS'], r['DP'],
        'p.X%dY' % i, 'c.%d' % i])
con.close()
""".replace("DBNAME", out.name)
    subprocess.run(["docker", "run", "--rm", "--network=none",
                    "--user", f"{os.getuid()}:{os.getgid()}",
                    "-v", f"{out.parent}:/w", "--entrypoint", "python", V1_IMAGE,
                    "-c", code], check=True)
    (out.parent / "rows.json").unlink()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def effective_args(method, input_type):
    """The backend's effective_args (design 1.1), to build test manifests."""
    args = list(method["cli_args"])
    for option in method["input_overrides"].get(input_type, {}).get("remove", []):
        while option in args:
            i = args.index(option)
            drop = 2 if i + 1 < len(args) and not args[i + 1].startswith("--") else 1
            del args[i:i + drop]
    return args


def manifest(method_key, input_type, target_source="kit_bed", **overrides):
    method = METHODS[method_key]
    refs = [{"path": r["path"], "sha256": sha256(FIX / "ref" / r["path"])}
            for r in method["reference_files"]]
    sample = {
        "tmbMethodId": METHOD_ID,
        "tmbMethodScope": "kit",
        "tmbMethodVersion": 1 if method_key == "v1" else 2,
        "tmbCalculatorImage": method["calculator_image"],
        "tmbInputType": input_type,
        "tmbMethodArgs": json.dumps(effective_args(method, input_type)),
        "tmbReferenceFiles": json.dumps(refs),
        "tmbTargetSource": target_source,
        "msiMinCoverage": "500",
    }
    sample.update(overrides)
    return {"samples": [{k: v for k, v in sample.items() if v is not None}]}


REF_BASE = "/runspace/data/ref/"           # trailing slash, as on the nodes
TEMP = f"/runspace/data/tmp/{SAMPLE_ID}/"


def caller_context(bed_name):
    return {
        "inputs.parameters.run-id": RUN_ID,
        "inputs.parameters.sample-id": SAMPLE_ID,
        "inputs.parameters.ref-base-dir": REF_BASE,
        "inputs.parameters.sample-temp-dir": TEMP,
        "tasks.duckdb-snv.outputs.parameters.output-duckdb":
            f"{TEMP}results/db/{SAMPLE_ID}_results.duckdb",
        "tasks.validate-bed-file.outputs.parameters.validated-bed-file":
            f"{TEMP}metrics/bed_validation/{bed_name}",
    }


def new_host(tmp, bed_name):
    host = Host(tmp)
    ctx = caller_context(bed_name)
    host.put(REF_BASE + "hg38/other", FIX / "ref")
    host.put(ctx["tasks.validate-bed-file.outputs.parameters.validated-bed-file"],
             FIX / bed_name)
    build_duckdb(host, ctx["tasks.duckdb-snv.outputs.parameters.output-duckdb"])
    return host, ctx


def tmb_dir(host):
    return host.path(f"{TEMP}/metrics/tmb/")


# ── The two paths ─────────────────────────────────────────────────────

# The Genome kit's target is the image's whole CDS BED: one tabix query per
# interval (195k) per VCF takes minutes. Both Genome-kit runs read this
# five-gene cut of it instead, mounted over the image's file.
CDS_SUBSET = [(FIX / "cds_subset.bed", "/app/data/gencode_v48_canonical_cds.bed")]


def run_today(host, ctx, caller, extra=()):
    today = json.loads((FIX / "today-callers.json").read_text())["callers"][caller]
    tmpl = templates(f"tasks/collect_tmb_data.yaml")[today["templateRef"]["template"]]
    args = {p["name"]: render(p["value"], ctx) for p in today["parameters"]}
    params = inputs_for(tmpl, args)
    endpoint = PROD_SOMATIC.format(run=RUN_ID)
    source = render(tmpl["script"]["source"], params)
    check(source.count(endpoint) == 1, f"today/{caller}: the script posts to the prod endpoint once")
    tmpl = dict(tmpl, script=dict(tmpl["script"],
                                  source=tmpl["script"]["source"].replace(
                                      PROD_SOMATIC.format(run="{{inputs.parameters.run-id}}"),
                                      STAGED_PATH)))
    run = run_script(host, tmpl, params, {"gennext_cluster_token": "STAGED"}, extra=extra)
    return run, tmb_dir(host) / "tmb_api_payload.json"


def v2_steps(host, ctx, caller, mf):
    """The branch's caller -> collect-tmb-data-v2 -> its steps' inputs."""
    task = dag_task(f"workflows/somatic/v2/somatic-dna-{caller}.yaml", "collect-tmb-data")
    ref = task["templateRef"]
    v2 = templates("tasks/collect_tmb_data_v2.yaml")
    args = {p["name"]: render(p["value"], ctx, mf) for p in task["arguments"]["parameters"]}
    outer = inputs_for(v2[ref["template"]], args)
    steps = {}
    for group in v2[ref["template"]]["steps"]:
        for step in group:
            step_args = {p["name"]: render(p["value"], outer)
                         for p in step["arguments"]["parameters"]}
            steps[step["name"]] = (v2[step["template"]],
                                   inputs_for(v2[step["template"]], step_args))
    return ref, steps


def run_v2_compute(host, ctx, caller, mf, image=None, extra=()):
    _, steps = v2_steps(host, ctx, caller, mf)
    tmpl, params = steps["compute-tmb"]
    run = run_script(host, tmpl, params, {}, image=image, extra=extra)
    return run, tmb_dir(host) / "tmb_api_payload.json"


def run_padding_repair(host, ctx):
    tmpl = templates("workflows/somatic/v2/somatic-dna-padding-repair.yaml")["compute-tmb-data"]
    other = REF_BASE + "/hg38/other/"
    args = {
        "run-id": RUN_ID, "sample-id": SAMPLE_ID,
        "duckdb-file": ctx["tasks.duckdb-snv.outputs.parameters.output-duckdb"],
        "target-bed-file": ctx["tasks.validate-bed-file.outputs.parameters.validated-bed-file"],
        "output-dir": f"{TEMP}/metrics/tmb/",
        "cosmic-vcf": other + "cosmic/Cosmic_CompleteTargetedScreensMutant_v103_GRCh38.vcf.gz",
        "gnomad-vcf": other + "common_variants/gnomad4/gnomad4_1pct.vcf.gz",
        "turkish-variome-vcf": other + "common_variants/turkish_variome/turkish_variome_1pct.vcf.gz",
    }
    # The live template names the tag, which resolves to v1's digest.
    check(tmpl["script"]["image"] == "bergun/gennext-tmb-calc:0.5.0",
          "padding-repair: compute-tmb-data is still on tmb-calc 0.5.0")
    run = run_script(host, tmpl, inputs_for(tmpl, args), {}, image=V1_IMAGE)
    return run, tmb_dir(host) / "tmb_api_payload.json"


def show(run):
    return f"exit {run.returncode}\n{run.stdout[-2000:]}\n{run.stderr[-2000:]}"


def compare(label, run_a, pa, run_b, pb):
    ok = run_a.returncode == 0 and run_b.returncode == 0 and pa.is_file() and pb.is_file()
    check(ok, f"{label}: both runs stage a payload")
    if not ok:
        print(show(run_a), "\n", show(run_b))
        return None
    same = pa.read_bytes() == pb.read_bytes()
    check(same, f"{label}: payloads are byte-identical ({sha256(pa)[:16]})")
    if not same:
        print(pa.read_text(), "\n---\n", pb.read_text())
    return json.loads(pa.read_text())


# ── Static checks on the templates ────────────────────────────────────

def static_checks():
    fields = {"tmb-method-id": "tmbMethodId", "tmb-input-type": "tmbInputType",
              "tmb-calculator-image": "tmbCalculatorImage", "tmb-method-args": "tmbMethodArgs",
              "tmb-reference-files": "tmbReferenceFiles", "tmb-target-source": "tmbTargetSource"}
    for caller in ("fastq", "bam", "vcf"):
        path = f"workflows/somatic/v2/somatic-dna-{caller}.yaml"
        task = dag_task(path, "collect-tmb-data")
        check(task["templateRef"] == {"name": "collect-tmb-data-v2",
                                      "template": "collect-tmb-data-v2"},
              f"{caller}: collect-tmb-data calls collect-tmb-data-v2")
        args = {p["name"]: p["value"] for p in task["arguments"]["parameters"]}
        for param, field in fields.items():
            m = MANIFEST_EXPR.fullmatch(args.get(param, ""))
            check(bool(m) and m.group(1) == field, f"{caller}: {param} is the manifest's {field}")
        check(not {"cosmic-vcf", "gnomad-vcf", "turkish-variome-vcf"} & args.keys(),
              f"{caller}: no cosmic, gnomad or turkish arguments")
        inputs = [p["name"] for p in templates(path)[f"somatic-testing-dna-gpu-{caller}-v2"]
                  ["inputs"]["parameters"]]
        check("sample-run-data" in inputs, f"{caller}: takes sample-run-data")
    for path in ("workflows/per-sample-somatic-dna.yaml", "workflows/per-sample-somatic-dna-rna.yaml"):
        task = dag_task(path, "run-dna-vcf")
        names = [p["name"] for p in task["arguments"]["parameters"]]
        check("sample-run-data" in names, f"{path}: passes sample-run-data to the vcf caller")
    v2 = templates("tasks/collect_tmb_data_v2.yaml")
    for name in ("compute-tmb", "publish-tmb"):
        vols = v2[name].get("volumes", [])
        check(all("osDir" not in v["hostPath"]["path"] for v in vols)
              and not any(n in v["name"] for v in vols for n in ("cosmic", "gnomad", "turkish")),
              f"collect-tmb-data-v2/{name}: no per-database volumes")
        check("@sha256:" in v2[name]["script"]["image"] or "{{" in v2[name]["script"]["image"],
              f"collect-tmb-data-v2/{name}: image pinned by digest or from the manifest")


# ── Publish, against a stub backend ───────────────────────────────────

class Backend(http.server.BaseHTTPRequestHandler):
    posts = []
    status = 200

    def do_POST(self):  # noqa: N802 (http.server's name)
        body = self.rfile.read(int(self.headers["Content-Length"]))
        Backend.posts.append((self.path, json.loads(body)))
        self.send_response(Backend.status)
        self.end_headers()
        self.wfile.write(b"stub")

    def log_message(self, *args):
        pass


def publish(host, ctx, caller, mf, status=200):
    server = http.server.HTTPServer(("127.0.0.1", 0), Backend)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    Backend.posts, Backend.status = [], status
    try:
        _, steps = v2_steps(host, ctx, caller, mf)
        tmpl, params = steps["publish-tmb"]
        params = dict(params, **{"inputs.parameters.api-base-url":
                                 f"http://127.0.0.1:{server.server_port}"})
        run = run_script(host, tmpl, params, {"gennext_cluster_token": "cluster-token"},
                         network="host")
        return run, list(Backend.posts)
    finally:
        server.shutdown()


def v0_workflow_checks():
    """The gennext-dev V0 Workflow can't post and runs this branch's code."""
    sys.path.insert(0, str(REPO / "tests"))
    import make_tmb_v0_workflow as gen  # noqa: E402 (path set just above)

    args = type("Args", (), dict(
        node="n1", run_id=RUN_ID, sample_id=SAMPLE_ID, duckdb="/d/s_results.duckdb",
        target_bed="/d/kit.valid.bed", scratch_dir="/s/tmb-v0/", input_type="vcf",
        target_source="kit_bed", ref_base_dir=REF_BASE, method_id=METHOD_ID))
    wf = gen.build(args)
    text = json.dumps(wf)
    check("secretKeyRef" not in text and "gennext-backend" not in text
          and "publish-tmb" not in text and "http://" not in text,
          "V0 workflow: no secret, no backend URL, no publish step")
    by_name = {t["name"]: t for t in wf["spec"]["templates"]}
    check(by_name["compute-tmb"] == templates("tasks/collect_tmb_data_v2.yaml")["compute-tmb"],
          "V0 workflow: compute-tmb is this branch's, verbatim")
    caller = dag_task("workflows/somatic/v2/somatic-dna-vcf.yaml", "collect-tmb-data")
    v2_task = [t for t in by_name["main"]["dag"]["tasks"] if t["name"] == "v2"][0]
    check([p["value"] for p in v2_task["arguments"]["parameters"] if p["name"].startswith("tmb-")]
          == [p["value"] for p in caller["arguments"]["parameters"] if p["name"].startswith("tmb-")],
          "V0 workflow: the vcf caller's manifest expressions, verbatim")
    mounts = [m for t in wf["spec"]["templates"] for m in (t.get("script") or {}).get("volumeMounts", [])
              if not m.get("readOnly")]
    check(all(m["name"] == "output-dir" for m in mounts),
          "V0 workflow: only the scratch output directories are writable")
    sample = json.loads(wf["spec"]["arguments"]["parameters"][2]["value"])["samples"][0]
    check(json.loads(sample["tmbMethodArgs"]) == METHODS["v1"]["cli_args"]
          and sample["tmbCalculatorImage"] == V1_IMAGE, "V0 workflow: carries the v1 method row")


def main():
    static_checks()
    v0_workflow_checks()
    with tempfile.TemporaryDirectory() as tmp:
        # V0, per caller, kit BED.
        v0 = {}
        for caller in ("fastq", "bam", "vcf"):
            host, ctx = new_host(f"{tmp}/{caller}-today", "kit.bed")
            a = run_today(host, ctx, caller)
            host2, ctx2 = new_host(f"{tmp}/{caller}-v2", "kit.bed")
            b = run_v2_compute(host2, ctx2, caller, manifest("v1", caller))
            v0[caller] = compare(f"V0 {caller}", *a, *b)
            if caller == "fastq":
                v0_host, v0_ctx = host2, ctx2
        kit = v0["fastq"]
        if kit:
            sd = kit["somaticData"]
            check(kit["token"] == "STAGED" and kit["runId"] == RUN_ID
                  and kit["sampleId"] == SAMPLE_ID, "V0: payload carries run, sample, STAGED")
            # 2 TP53 + the TP53 deletion + KRAS at AF 0.001 + 2 PIK3CA, over the kit BED
            check(sd["tmbVariantCount"] == 6 and sd["tmbTargetRegionSize"] == KIT_BP,
                  f"V0: every filter of today's method acted on the fixture "
                  f"({sd['tmbVariantCount']} variants over {sd['tmbTargetRegionSize']} bp)")
            steps = [s["stepName"] for s in sd.get("processingSteps", [])]
            print("     steps:", " > ".join(f"{s['stepName']} {s['variantCount']}"
                                           for s in sd.get("processingSteps", [])))
            check(any("COSMIC" in s for s in steps) and any("lacklist" in s for s in steps),
                  "V0: the COSMIC whitelist and the gene blacklist ran")

        # V0, Genome kit: today switches on the BED filename, v2 on tmbTargetSource.
        genome_bed = "hg38_canonical_chrom_contig_regions.bed"
        host, ctx = new_host(f"{tmp}/genome-today", genome_bed)
        a = run_today(host, ctx, "fastq", extra=CDS_SUBSET)
        host, ctx = new_host(f"{tmp}/genome-v2", genome_bed)
        b = run_v2_compute(host, ctx, "fastq", manifest("v1", "fastq", "cds_bed"), extra=CDS_SUBSET)
        genome = compare("V0 Genome kit (CDS target)", *a, *b)
        if genome and kit:
            check(genome["somaticData"]["tmbVariantCount"] == 7,
                  "V0 Genome kit: the CDS target keeps the BRCA2 variant outside the kit BED")

        # V0b: padding repair stages the same bytes as v2 with v1.
        host, ctx = new_host(f"{tmp}/padding", "kit.bed")
        a = run_padding_repair(host, ctx)
        compare("V0b padding repair = v2 with v1", *a,
                run_v2_compute(v0_host, v0_ctx, "fastq", manifest("v1", "fastq"))[0],
                tmb_dir(v0_host) / "tmb_api_payload.json")

        # v2 through the same template: fastq vs vcf effective args.
        counts = {}
        for caller in ("fastq", "vcf"):
            host, ctx = new_host(f"{tmp}/v2-{caller}", "kit.bed")
            run, path = run_v2_compute(host, ctx, caller, manifest("v2", caller))
            ok = run.returncode == 0 and path.is_file()
            check(ok, f"v2 {caller}: the 0.6.0 method runs through collect-tmb-data-v2")
            if not ok:
                print(show(run))
                continue
            counts[caller] = json.loads(path.read_text())["somaticData"]
            check("--min-depth" in run.stdout if caller == "fastq" else "--min-depth" not in run.stdout,
                  f"v2 {caller}: the command line {'has' if caller == 'fastq' else 'lacks'} --min-depth")
        if len(counts) == 2:
            check(counts["vcf"]["tmbVariantCount"] == counts["fastq"]["tmbVariantCount"] + 1,
                  "v2: vcf input keeps the DP 20 variant that --min-depth 50 removes "
                  f"({counts['fastq']['tmbVariantCount']} vs {counts['vcf']['tmbVariantCount']})")
            check(counts["fastq"]["tmbTargetRegionSize"] < KIT_BP,
                  "v2: the denominator is the CDS part of the kit BED")

        # The compute stage refuses a bad method and leaves no payload behind.
        host, ctx = new_host(f"{tmp}/neg", "kit.bed")
        payload = tmb_dir(host) / "tmb_api_payload.json"
        good = manifest("v1", "fastq")
        bad_ref = json.loads(good["samples"][0]["tmbReferenceFiles"])
        bad_ref[2]["sha256"] = "0" * 64
        listed = json.loads(good["samples"][0]["tmbReferenceFiles"])
        missing_ref = listed + [{"path": "common_variants/gnomad4/not_installed_here.vcf.gz",
                                 "sha256": "0" * 64}]
        args = METHODS["v1"]["cli_args"]
        cases = [
            ("a reference file with the wrong checksum", {"tmbReferenceFiles": json.dumps(bad_ref)}, "checksum mismatch"),
            ("a reference file missing on the node", {"tmbReferenceFiles": json.dumps(missing_ref)}, "missing on this node"),
            ("a reference path outside /mnt/ref", {"tmbReferenceFiles": json.dumps(listed + [{"path": "../x", "sha256": "0" * 64}])}, "relative"),
            ("an empty tmbMethodId", {"tmbMethodId": None}, "no tmbMethodId"),
            ("an image not pinned by digest", {"tmbCalculatorImage": "bergun/gennext-tmb-calc:0.5.0"}, "not pinned"),
            ("an image with a tag (the backend refuses it)", {"tmbCalculatorImage": "bergun/gennext-tmb-calc:0.5.0@" + V1_IMAGE.split("@")[1]}, "not pinned"),
            ("args naming a /mnt/ref file that is not listed", {"tmbMethodArgs": json.dumps(args + ["--cosmic-vcf", "/mnt/ref/cosmic/other.vcf.gz"])}, "not a listed reference file"),
            ("args naming /mnt/ref outside the mount", {"tmbMethodArgs": json.dumps(args + ["--cosmic-vcf", "/mnt/refx/a.vcf.gz"])}, "must name a file under"),
            ("a short option other than -v", {"tmbMethodArgs": json.dumps(args + ["-o", "/tmp"])}, "short option"),
            ("a method that sets --api-endpoint", {"tmbMethodArgs": json.dumps(args + ["--api-endpoint", "http://x"])}, "may not set"),
            ("a method that sets a prefix of --output-dir", {"tmbMethodArgs": json.dumps(args + ["--output", "/tmp"])}, "may not set"),
            ("an unknown target source", {"tmbTargetSource": "bed"}, "target source"),
            ("an unknown input type", {"tmbInputType": "cram"}, "input type"),
            ("arguments that are not a JSON array", {"tmbMethodArgs": "--verbose"}, "not JSON"),
        ]
        for what, override, message in cases:
            payload.parent.mkdir(parents=True, exist_ok=True)
            payload.write_text("stale")
            run, _ = run_v2_compute(host, ctx, "fastq", manifest("v1", "fastq", **override),
                                    image=V1_IMAGE)
            check(run.returncode != 0 and message in run.stderr,
                  f"compute fails on {what}")
        # A backend that predates TMB methods: none of the six fields. The
        # image falls back to the pinned python image, so the pod starts and
        # the script says why it fails (an empty image would be refused by
        # Kubernetes before any script ran).
        old = manifest("v1", "fastq", tmbMethodId=None, tmbInputType=None,
                       tmbCalculatorImage=None, tmbMethodArgs=None,
                       tmbReferenceFiles=None, tmbTargetSource=None)
        _, steps = v2_steps(host, ctx, "fastq", old)
        tmpl, params = steps["compute-tmb"]
        fallback = render(tmpl["script"]["image"], params)
        check(fallback.startswith("python:3.11-alpine3.18@sha256:"),
              "compute: a manifest without tmb* fields still gets a pinned image")
        run, _ = run_v2_compute(host, ctx, "fastq", old)
        check(run.returncode != 0 and "the run manifest has no tmbMethodId" in run.stderr,
              "compute: a manifest without tmb* fields fails with its own message")
        check(render(tmpl["script"]["image"], dict(params, **{"tmb-calculator-image": V1_IMAGE}))
              == V1_IMAGE, "compute: with the field set, the manifest's image is used")
        payload.write_text("stale")
        canary = "$(touch /mnt/output/shell-ran)"
        run, _ = run_v2_compute(host, ctx, "fastq",
                                manifest("v1", "fastq", tmbMethodArgs=json.dumps(args + [canary])))
        check(run.returncode != 0 and not payload.exists(),
              "compute: a failed calculator run leaves no payload (the stale one is removed)")
        check(not (tmb_dir(host) / "shell-ran").exists(),
              "compute: method arguments are never shell-interpreted")

        # Publish the V0 fastq payload.
        run_v2_compute(v0_host, v0_ctx, "fastq", manifest("v1", "fastq"))
        staged = json.loads((tmb_dir(v0_host) / "tmb_api_payload.json").read_text())
        run, posts = publish(v0_host, v0_ctx, "fastq", manifest("v1", "fastq"))
        check(run.returncode == 0 and len(posts) == 1, "publish: one POST, exit 0")
        if not (run.returncode == 0 and len(posts) == 1):
            print(show(run))
        else:
            path, body = posts[0]
            expect = json.loads(json.dumps(staged))
            expect["token"] = "cluster-token"
            expect["somaticData"].update(tmbMethodId=METHOD_ID, tmbInputType="fastq",
                                         tmbEffectiveArgs=METHODS["v1"]["cli_args"])
            check(path == f"/api/v0/runs/{RUN_ID}/somatic", "publish: posts to the run's somatic route")
            check(body == expect, "publish: staged payload + tmbMethodId, tmbInputType, "
                                  "tmbEffectiveArgs, cluster token")
        run, _ = publish(v0_host, v0_ctx, "fastq", manifest("v1", "fastq"), status=422)
        check(run.returncode != 0 and "422" in run.stderr, "publish: a non-200 answer fails the step")
        p = tmb_dir(v0_host) / "tmb_api_payload.json"
        p.write_text(json.dumps(dict(staged, sampleId="OTHER")))
        run, posts = publish(v0_host, v0_ctx, "fastq", manifest("v1", "fastq"))
        check(run.returncode != 0 and not posts, "publish: a payload for another sample is refused")
        p.write_text(json.dumps(dict(staged, token="real")))
        run, posts = publish(v0_host, v0_ctx, "fastq", manifest("v1", "fastq"))
        check(run.returncode != 0 and not posts, "publish: a payload without the STAGED token is refused")

    if FAILURES:
        print(f"\n{len(FAILURES)} check(s) failed")
        sys.exit(1)
    print("\nall checks passed")


if __name__ == "__main__":
    main()
