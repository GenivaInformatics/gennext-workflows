"""End-to-end check of fusions-proto-pg's load-events step: fusion-events build output -> the
template's own load script, run in the pinned DuckDB 0.8.1 image.

Run from the repository (needs Docker and PyYAML, and the gennext-bpview:ci image that
`mise run test:gennext-bpview` builds): python3 images/gennext-bpview/tests/test_load_events.py
"""
import os
import re
import subprocess
import sys
import tempfile

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..", "..", "..")
TEMPLATE = os.path.join(REPO, "workflows", "prototypes", "fusions-testing-ground", "fusions-proto-pg.yaml")
BPVIEW_IMAGE = os.environ.get("BPVIEW_IMAGE", "gennext-bpview:ci")

ARRIBA = ("#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\t"
          "type\tsplit_reads1\tsplit_reads2\tdiscordant_mates\tconfidence\treading_frame\ttags\tdirection1\t"
          "direction2\tarriba_status\tfilters\n"
          "A1\tB1\t+/+\t+/+\tchr1:1000\tchr2:500\tsplice-site\tsplice-site\ttranslocation\t3\t2\t1\thigh\t"
          "in-frame\tMitelman\tdownstream\tupstream\tkept\tduplicates(4)\n"
          "A2\tB2\t+/+\t+/+\tchr3:1000\tchr4:500\tsplice-site\tsplice-site\ttranslocation\t1\t0\t0\tlow\t"
          "in-frame\tMitelman\tdownstream\tupstream\trescued\tin_vitro(1)\n")


def sh(*cmd, check=True):
    return subprocess.run(cmd, check=check, capture_output=True, text=True)


def main():
    doc = yaml.safe_load(open(TEMPLATE))
    t = [x for x in doc["spec"]["templates"] if x["name"] == "load-events"][0]
    image = t["script"]["image"]
    script = t["script"]["source"].replace("{{inputs.parameters.dev-sample-id}}", "T")
    assert "{{" not in script, "unsubstituted parameter in load-events"
    uid = f"{os.getuid()}:{os.getgid()}"
    with tempfile.TemporaryDirectory() as d:
        os.makedirs(f"{d}/results/db")
        os.makedirs(f"{d}/results/fusions/fusion-events")
        open(f"{d}/arriba.tsv", "w").write(ARRIBA)
        open(f"{d}/load.sh", "w").write(script)
        sh("docker", "run", "--rm", "--network", "none", "-u", uid, "-v", f"{d}:/p", "-w", "/p", BPVIEW_IMAGE,
           "fusion-events", "build", "--sample-id", "T", "--arriba-tsv", "arriba.tsv", "--protein-domains", ".",
           "--prototype", "--arriba-settings", "filters_restored",
           "--out-tsv", "results/fusions/fusion-events/T_fusion_events.tsv",
           "--out-meta-tsv", "results/fusions/fusion-events/T_fusion_events_meta.tsv")
        # an existing table in the copied DuckDB must survive the load
        sh("docker", "run", "--rm", "--network", "none", "-u", uid, "-v", f"{d}/results:/mnt/results", image,
           "duckdb", "/mnt/results/db/T_fusions.duckdb", "CREATE TABLE fusions_with_oncokb AS SELECT 42 AS x")
        run = ["docker", "run", "--rm", "--network", "none", "-u", uid, "-v", f"{d}/results:/mnt/results",
               "-v", f"{d}/load.sh:/load.sh:ro", image, "sh", "/load.sh"]
        for attempt in (1, 2):  # a retried load must give the same result
            r = sh(*run, check=False)
            assert r.returncode == 0, (attempt, r.stdout, r.stderr)
            assert "loaded into" in r.stdout and ": 2 events" in r.stdout, r.stdout
        q = sh("docker", "run", "--rm", "--network", "none", "-u", uid, "-v", f"{d}/results:/mnt/results", image,
               "duckdb", "-csv", "-noheader", "/mnt/results/db/T_fusions.duckdb",
               "SELECT (SELECT x FROM fusions_with_oncokb), (SELECT typeof(tier_flags) FROM fusion_events LIMIT 1), "
               "(SELECT count(*) FROM fusion_events WHERE list_contains(tier_flags, 'rescued')), "
               "(SELECT value FROM fusion_events_meta WHERE key = 'arriba_settings'), "
               "(SELECT min(event_rank) FROM fusion_events)").stdout.strip()
        assert q == "42,VARCHAR[],1,filters_restored,1", q
        # a meta count that doesn't match the table fails the step
        meta = f"{d}/results/fusions/fusion-events/T_fusion_events_meta.tsv"
        text = open(meta).read()
        assert "events\t2" in text, text
        open(meta, "w").write(re.sub(r"events\t2", "events\t3", text))
        r = sh(*run, check=False)
        assert r.returncode != 0 and "differs from the meta count" in r.stdout, (r.stdout, r.stderr)
    print("test_load_events: ok")


if __name__ == "__main__":
    sys.exit(main())
