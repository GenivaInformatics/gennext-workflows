"""Synthetic test of fusion_events.py: merge, annotation join, support, PoN recurrence, tiers.

Run: python tests/test_fusion_events.py   (fusion_events.py and bpview.py one directory up)
"""
import csv
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
FE = os.path.join(HERE, "..", "fusion_events.py")
sys.path.insert(0, os.path.join(HERE, ".."))
import fusion_events  # noqa: E402

ARRIBA_HDR = ["#gene1", "gene2", "strand1(gene/fusion)", "strand2(gene/fusion)", "breakpoint1", "breakpoint2",
              "site1", "site2", "type", "split_reads1", "split_reads2", "discordant_mates", "coverage1",
              "coverage2", "confidence", "reading_frame", "tags", "direction1", "direction2"]
SF_HDR = ["#FusionName", "JunctionReadCount", "SpanningFragCount", "SpliceType", "LeftGene", "LeftBreakpoint",
          "RightGene", "RightBreakpoint", "LargeAnchorSupport", "FFPM", "annots"]


def tsv(path, hdr, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(h.lstrip("#"), ".")) for h in hdr) + "\n")


def arriba_row(g1, g2, bp1, bp2, sr1, sr2, conf, frame, tags="."):
    return {"gene1": g1, "gene2": g2, "strand1(gene/fusion)": "+/+", "strand2(gene/fusion)": "+/+",
            "breakpoint1": bp1, "breakpoint2": bp2, "site1": "splice-site", "site2": "splice-site",
            "type": "translocation", "split_reads1": sr1, "split_reads2": sr2, "discordant_mates": 0,
            "coverage1": 10, "coverage2": 10, "confidence": conf, "reading_frame": frame, "tags": tags,
            "direction1": "downstream", "direction2": "upstream"}


def sf_row(g1, g2, bp1, bp2, j, s):
    return {"FusionName": f"{g1}--{g2}", "JunctionReadCount": j, "SpanningFragCount": s,
            "SpliceType": "ONLY_REF_SPLICE", "LeftGene": f"{g1}^ENSG1", "LeftBreakpoint": bp1 + ":+",
            "RightGene": f"{g2}^ENSG2", "RightBreakpoint": bp2 + ":+", "LargeAnchorSupport": "YES_LDAS",
            "FFPM": 0.5, "annots": '["Mitelman"]'}


def oncokb(level="", onc="Unknown", dbs=()):
    return {"Databases": list(dbs), "Fusion Indication Index (FII)": 0.5, "clinical_tier": "Tier 3",
            "OncoKB": {"ONCOGENIC": onc, "HIGHEST_LEVEL": level, "HIGHEST_DX_LEVEL": "",
                       "MUTATION_EFFECT": "Gain-of-function"}}


def run_build(d, sample, arriba, sf, pairs, extra=()):
    a, s, j, out = (os.path.join(d, f"{sample}.{x}") for x in ("arriba.tsv", "sf.tsv", "oncokb.json", "ev.tsv"))
    tsv(a, ARRIBA_HDR, arriba)
    tsv(s, SF_HDR, sf)
    with open(j, "w") as fh:
        json.dump([dict(v, Fusion=k) for k, v in pairs.items()], fh)
    subprocess.run([sys.executable, FE, "build", "--sample-id", sample, "--arriba-tsv", a, "--starfusion-tsv", s,
                    "--oncokb-json", j, "--out-tsv", out, *extra], check=True)
    with open(out) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(out) as fh:
        assert fh.readline().rstrip("\n").split("\t") == [c for c, _ in fusion_events.COLUMNS]
    return out, {r["fusion"] + "@" + r["breakpoint5"]: r for r in rows}


def main():
    with tempfile.TemporaryDirectory() as d:
        pairs = {"ALK1--PART1": oncokb("LEVEL_1", "Oncogenic"),         # actionable
                 "NTRKX--JUNK": oncokb("LEVEL_1", "Likely Oncogenic"),  # partner-agnostic OncoKB hit
                 "KNOWN--PAIR": oncokb(dbs=["Mitelman"]),               # curated, weak support
                 "RECUR--PAIR": oncokb(dbs=["Mitelman"])}               # curated but recurrent
        arriba = [arriba_row("ALK1", "PART1", "chr1:1000", "chr2:500", 3, 4, "high", "in-frame"),
                  arriba_row("RECUR", "PAIR", "chr5:1000", "chr6:500", 2, 2, "high", "in-frame")]
        sf = [sf_row("ALK1", "PART1", "chr1:1003", "chr2:498", 6, 1),   # same junction within 10 bp
              sf_row("NTRKX", "JUNK", "chr3:1000", "chr4:500", 1, 0),   # 1 molecule, one caller
              sf_row("KNOWN", "PAIR", "chr7:1000", "chr8:500", 1, 0)]
        # PoN: RECUR--PAIR in 3 other samples, ALK1--PART1 in 1 other sample.
        pon_rows = [arriba_row("RECUR", "PAIR", "chr5:1000", "chr6:500", 2, 2, "low", ".")]
        evs = []
        for i, sample in enumerate(("P1", "P2", "P3")):
            evs.append(run_build(d, sample, pon_rows + ([arriba[0]] if i == 0 else []), [], {})[0])
        pon = os.path.join(d, "pon.tsv")
        subprocess.run([sys.executable, FE, "pon", "--events-tsv", *evs, "--out-tsv", pon], check=True)

        _, rows = run_build(d, "S1", arriba, sf, pairs, ["--pon", pon])
        assert len(rows) == 4, sorted(rows)
        alk = rows["ALK1--PART1@chr1:1000:+"]
        assert alk["callers"] == "arriba,starfusion" and alk["n_callers"] == "2", alk
        assert alk["breakpoint3"] == "chr2:500:+" and alk["reading_frame"] == "in-frame"
        assert alk["arriba_split_reads"] == "7" and alk["starfusion_junction_reads"] == "6"
        assert alk["unique_molecules"] == "7" and alk["support_source"] == "caller" and alk["pcr_copies"] == ""
        assert alk["tier"] == "1" and alk["oncokb_highest_level"] == "LEVEL_1", alk["tier_reason"]
        assert alk["pon_gene_pair_count"] == "1" and alk["pon_junction_count"] == "1" and alk["pon_size"] == "3"
        assert alk["event_id"] == "S1:E0001", alk["event_id"]  # tier 1 sorts first

        junk = rows["NTRKX--JUNK@chr3:1000:+"]
        assert junk["tier"] == "3" and "but one caller" in junk["tier_reason"], junk["tier_reason"]
        known = rows["KNOWN--PAIR@chr7:1000:+"]
        assert known["tier"] == "2" and known["tier_reason"] == "known fusion: Mitelman", known["tier_reason"]
        recur = rows["RECUR--PAIR@chr5:1000:+"]
        assert recur["pon_gene_pair_count"] == "3" and recur["tier"] == "3", recur
        assert "recurrent in PoN (3/3 samples): capped at tier 3" in recur["tier_reason"], recur["tier_reason"]

        # bpview counts replace the caller counts when bpview drew the event.
        bp = os.path.join(d, "bpview.events.json")
        with open(bp, "w") as fh:
            json.dump({"events": [{"breakpoint1": "chr1:1000", "breakpoint2": "chr2:500", "keep1": "le",
                                   "keep2": "ge", "contig": "ALK1--PART1__chr1_1000__chr2_500",
                                   "fragments_nondup": {"junction": 2, "spanning": 1},
                                   "fragments_dup": {"junction": 40, "spanning": 3}}]}, fh)
        _, rows = run_build(d, "S1", arriba, sf, pairs, ["--bpview-events", bp, "--tier1-min-unique-molecules", "4"])
        alk = rows["ALK1--PART1@chr1:1000:+"]
        assert alk["unique_molecules"] == "3" and alk["pcr_copies"] == "43" and alk["support_source"] == "bpview"
        assert alk["bpview_contig"] == "ALK1--PART1__chr1_1000__chr2_500"
        assert alk["tier"] == "2" and "(not tier 1: 3 unique molecule(s) < 4)" in alk["tier_reason"], alk["tier_reason"]

        # Empty caller outputs give an empty table with a header.
        _, rows = run_build(d, "S2", [], [], {})
        assert rows == {}
    print("test_fusion_events: ok")


if __name__ == "__main__":
    main()
