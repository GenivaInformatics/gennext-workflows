#!/usr/bin/env python3
"""HRD that could not be computed is posted as null, never as 0.

When CNVkit's call.cns has no allele-specific copy numbers (cn1/cn2), the
scarhrd step cannot run scarHRD. It writes a fallback _hrd.json so the sample
still completes. That fallback used to hold zeros, which the backend stored
and the frontend showed as a real HRD score of 0. It now holds nulls and a
reason.

This runs the real scripts from the templates:
  - the scarhrd step's fallback branch, from both per-sample CNV workflows,
    on a call.cns without cn1/cn2;
  - collect-hrd-data, from tasks/collect_somatic_data.yaml, against a stub
    backend.
Checks:
  - both fallbacks write null scores and a reason, and exit 0;
  - collect-hrd-data posts the nulls as JSON null and exits 0;
  - real scores are still posted as integers (floats cast to int);
  - a missing key still fails the step.

Needs python3 with PyYAML and requests, and bash.
"""

import http.server
import json
import os
import pathlib
import subprocess
import sys
import tempfile
import threading

import yaml

REPO = pathlib.Path(__file__).resolve().parent.parent
WORKFLOWS = [
    "workflows/per-sample-cnv-analysis.yaml",
    "workflows/per-sample-cnv-cohort-analysis.yaml",
]
FAILURES = []


def check(ok, what):
    print(("ok   " if ok else "FAIL ") + what)
    if not ok:
        FAILURES.append(what)


def template(path, name):
    for doc in yaml.safe_load_all((REPO / path).read_text()):
        for tmpl in (doc or {}).get("spec", {}).get("templates", []):
            if tmpl.get("name") == name:
                return tmpl
    raise SystemExit(f"{path}: no template {name}")


def scarhrd_fallback(path):
    script = template(path, "scarhrd")["script"]["args"][0]
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp)
        (root / "output").mkdir()
        # A plain call.cns: CNVkit's columns without cn1/cn2.
        (root / "output" / "S1.call.cns").write_text(
            "chromosome\tstart\tend\tgene\tlog2\tcn\tdepth\tprobes\tweight\n"
            "chr1\t1000\t2000\tA\t0.0\t2\t100\t5\t1\n"
        )
        rendered = script.replace("{{inputs.parameters.sample-id}}", "S1").replace(
            "/mnt/", f"{tmp}/"
        )
        run = subprocess.run(
            ["bash", "-euo", "pipefail", "-c", rendered],
            capture_output=True,
            text=True,
            check=False,
        )
        check(run.returncode == 0, f"{path}: the fallback exits 0")
        check(
            "Writing null HRD fallback" in run.stderr,
            f"{path}: the fallback logs its WARN line",
        )
        data = json.loads((root / "hrd-output" / "S1_hrd.json").read_text())
        check(
            all(data[k] is None for k in ("HRD", "Telomeric_AI", "LST", "HRD_sum")),
            f"{path}: the fallback writes null scores, not 0 ({data})",
        )
        check(bool(data.get("reason")), f"{path}: the fallback gives a reason")


class Backend(http.server.BaseHTTPRequestHandler):
    posts = []

    def do_POST(self):  # noqa: N802 (http.server's name)
        body = self.rfile.read(int(self.headers["Content-Length"]))
        Backend.posts.append((self.path, json.loads(body)))
        self.send_response(200)
        self.end_headers()

    def log_message(self, *args):
        pass


def collect(hrd_json):
    """Runs collect-hrd-data on hrd_json; returns (exit code, posted bodies)."""
    tmpl = template("tasks/collect_somatic_data.yaml", "collect-hrd-data")
    source = tmpl["script"]["source"]
    server = http.server.HTTPServer(("127.0.0.1", 0), Backend)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    Backend.posts = []
    try:
        with tempfile.TemporaryDirectory() as tmp:
            (pathlib.Path(tmp) / "S1_hrd.json").write_text(json.dumps(hrd_json))
            rendered = (
                source.replace("{{inputs.parameters.run-id}}", "run-1")
                .replace("{{inputs.parameters.sample-id}}", "S1")
                .replace("{{=sprig.osBase(inputs.parameters['hrd-json-f'])}}", "S1_hrd.json")
                .replace("/mnt/input/", f"{tmp}/")
                .replace(
                    "gennext-backend.gennext-prod.svc.cluster.local:3030",
                    f"127.0.0.1:{server.server_port}",
                )
            )
            check("{{" not in rendered, "collect-hrd-data: every parameter rendered")
            run = subprocess.run(
                [sys.executable, "-c", rendered],
                capture_output=True,
                text=True,
                env={**os.environ, "GENNEXT_CLUSTER_TOKEN": "test-token"},
                check=False,
            )
        return run.returncode, [body for _, body in Backend.posts], run
    finally:
        server.shutdown()


def main():
    for path in WORKFLOWS:
        scarhrd_fallback(path)

    nulls = {"HRD": None, "Telomeric_AI": None, "LST": None, "HRD_sum": None}
    code, posts, run = collect({**nulls, "reason": "no cn1/cn2"})
    check(code == 0, f"collect-hrd-data: null scores exit 0 ({run.stderr.strip()[-200:]})")
    check(len(posts) == 1, "collect-hrd-data: null scores are posted once")
    if posts:
        sent = posts[0]["somaticData"]
        check(
            sent == {"hrdValue": None, "hrdTelomericAi": None, "hrdLst": None, "hrdSum": None},
            f"collect-hrd-data: posted as JSON null, not 0 ({sent})",
        )
    check("HRD not computed: no cn1/cn2" in run.stdout, "collect-hrd-data: logs the reason")

    code, posts, _ = collect({"HRD": 13.0, "Telomeric_AI": 13.0, "LST": 7.0, "HRD_sum": 33.0})
    check(code == 0 and len(posts) == 1, "collect-hrd-data: real scores exit 0")
    if posts:
        sent = posts[0]["somaticData"]
        check(
            sent == {"hrdValue": 13, "hrdTelomericAi": 13, "hrdLst": 7, "hrdSum": 33},
            f"collect-hrd-data: real scores are posted as integers ({sent})",
        )

    code, posts, _ = collect({"HRD": 1.0, "Telomeric_AI": 2.0, "LST": 3.0})
    check(code != 0 and not posts, "collect-hrd-data: a missing key still fails")

    if FAILURES:
        print(f"{len(FAILURES)} failed")
        sys.exit(1)
    print("all passed")


if __name__ == "__main__":
    main()
