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


def arriba_row(g1, g2, bp1, bp2, sr1, sr2, conf, frame, tags=".", strands="+/+"):
    return {"gene1": g1, "gene2": g2, "strand1(gene/fusion)": strands, "strand2(gene/fusion)": strands,
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
    if "--protein-domains" not in extra:
        extra = ["--protein-domains", ".", *extra]
    subprocess.run([sys.executable, FE, "build", "--sample-id", sample, "--arriba-tsv", a, "--starfusion-tsv", s,
                    "--oncokb-json", j, "--out-tsv", out, *extra], check=True)
    with open(out) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(out) as fh:
        assert fh.readline().rstrip("\n").split("\t") == [c for c, _ in fusion_events.COLUMNS]
    return out, {r["fusion"] + "@" + r["breakpoint5"]: r for r in rows}


def test_contract_v1(d):
    """0.4.0 additions: rescue, low support, exons, kinase, start sites, meta, the Arriba-first fallback."""
    # Gene K5 (+ strand, chr1) exons 1..3; gene K3 (- strand, chr2) exons 1..3 (exon 1 rightmost).
    gtf = os.path.join(d, "mane.gtf")
    with open(gtf, "w") as fh:
        def ex(chrom, a, b, strand, gene, n, tag='tag "MANE_Select";'):
            fh.write(f'{chrom}\tT\texon\t{a}\t{b}\t.\t{strand}\t.\tgene_name "{gene}"; '
                     f'transcript_id "{gene}-201"; exon_number {n}; {tag}\n')
        for n, (a, b) in enumerate([(100, 200), (400, 500), (800, 900)], 1):
            ex("chr1", a, b, "+", "K5", n)
        for n, (a, b) in enumerate([(800, 900), (400, 500), (100, 200)], 1):
            ex("chr2", a, b, "-", "K3", n)
        ex("chr1", 100, 900, "+", "K5", 1, tag='tag "basic";')  # a non-MANE transcript is ignored
    dom = os.path.join(d, "domains.gff3")
    with open(dom, "w") as fh:
        # K3's kinase domain at 120-480 (genomic, - strand gene); K5 has none.
        for a, b in ((120, 200), (400, 480)):
            fh.write(f"2\tpfam\tprotein_domain\t{a}\t{b}\t0\t-\t.\tName=Protein kinase domain;"
                     f"gene_name=K3;protein_domain_id=PF00069\n")
    # 5' K5 breaks at the end of its exon 2 (500), 3' K3 (- strand: keeps <= bp) at the start... of
    # exon 2 read 3'->5': K3 exon 2 spans 400-500, its first kept base in transcript order is 500.
    a = arriba_row("K5", "K3", "chr1:500", "chr2:500", 2, 1, "high", "in-frame")
    a["direction2"] = "downstream"  # the - strand 3' partner keeps coordinates <= 500
    _, rows = run_build(d, "K1", [a], [], {}, ["--gtf", gtf, "--protein-domains", dom, "--out-meta-tsv",
                                                os.path.join(d, "K1.meta.tsv"), "--arriba-settings",
                                                "filters_restored", "--prototype"])
    r = list(rows.values())[0]
    assert (r["transcript5"], r["exon5"], r["site5"]) == ("K5-201", "2", "exon_boundary"), r
    assert (r["transcript3"], r["exon3"], r["site3"]) == ("K3-201", "2", "exon_boundary"), r
    assert r["kinase_gene"] == "K3" and r["kinase_retained"] == "true" and r["kinase_domain_fraction"] == "1.0", r
    assert r["support_source"] == "arriba" and r["unique_molecules"] == "3", r
    meta = dict(l.rstrip("\n").split("\t") for l in open(os.path.join(d, "K1.meta.tsv")))
    assert meta["contract_version"] == "1" and meta["arriba_settings"] == "filters_restored", meta
    assert meta["prototype"] == "true" and meta["pon"] == "none" and meta["events"] == "1", meta
    # Intronic 5' breakpoint (600: between exons 2 and 3); 3' at 300 cuts K3's domain (keeps <= 300).
    a = arriba_row("K5", "K3", "chr1:600", "chr2:300", 2, 1, "high", "in-frame")
    a["direction2"] = "downstream"
    _, rows = run_build(d, "K2", [a], [], {}, ["--gtf", gtf, "--protein-domains", dom])
    r = list(rows.values())[0]
    assert (r["exon5"], r["site5"]) == ("2", "intron") and (r["exon3"], r["site3"]) == ("3", "intron"), r
    assert r["kinase_retained"] == "false" and r["kinase_domain_fraction"] == "0.5", r
    # Antisense kept part: no exon numbers; unknown gene: nothing.
    a = arriba_row("K5", "NOPE", "chr1:500", "chr3:500", 2, 1, "high", "in-frame")
    a["direction1"] = "upstream"
    _, rows = run_build(d, "K4", [a], [], {}, ["--gtf", gtf])
    r = list(rows.values())[0]
    assert r["exon5"] == "" and r["site5"] == "" and r["exon3"] == "" and r["kinase_gene"] == "", r

    # A minus-strand 5' partner (K3 keeps >= 400: its exons 1-2) and a plus-strand 3' partner (K5 from 400).
    a = arriba_row("K3", "K5", "chr2:400", "chr1:400", 2, 1, "high", "in-frame")
    a["direction1"] = "upstream"
    _, rows = run_build(d, "K5M", [a], [], {}, ["--gtf", gtf])
    r = list(rows.values())[0]
    assert (r["exon5"], r["site5"], r["exon3"], r["site3"]) == ("2", "exon_boundary", "2", "exon_boundary"), r
    # PF07714 counts as a kinase domain too.
    dom2 = os.path.join(d, "domains2.gff3")
    with open(dom2, "w") as fh:
        fh.write("1\tpfam\tprotein_domain\t150\t190\t0\t+\t.\tName=Protein tyrosine and serine/threonine kinase;"
                 "gene_name=K5;protein_domain_id=PF07714\n")
    a = arriba_row("K5", "NOPE", "chr1:500", "chr3:500", 2, 1, "high", "in-frame")
    _, rows = run_build(d, "K6", [a], [], {}, ["--protein-domains", dom2])
    r = list(rows.values())[0]
    assert (r["kinase_gene"], r["kinase_retained"], r["kinase_domain_fraction"]) == ("K5", "true", "1.0"), r
    # Two Arriba rows at one junction stay two events; the second ID gets #2.
    a = arriba_row("D1", "D2", "chr1:50", "chr2:60", 3, 0, "high", "in-frame")
    out, _ = run_build(d, "DUP", [a, dict(a, split_reads1=1)], [], {})
    ids = [r["event_id"] for r in csv.DictReader(open(out), delimiter="\t")]
    assert ids == ["chr1:50:+>chr2:60:+", "chr1:50:+>chr2:60:+#2"], ids

    # bpview start sites and molecules.
    bp = os.path.join(d, "K.events.json")
    with open(bp, "w") as fh:
        json.dump({"events": [{"breakpoint1": "chr1:500", "breakpoint2": "chr2:500", "keep1": "le", "keep2": "le",
                               "contig": "c1", "fragments_nondup": {"junction": 4, "spanning": 2},
                               "fragments_dup": {"junction": 9}, "junction_start_sites": 3}]}, fh)
    a = arriba_row("K5", "K3", "chr1:500", "chr2:500", 2, 1, "high", "in-frame")
    a["direction2"] = "downstream"
    _, rows = run_build(d, "K1", [a], [], {}, ["--bpview-events", bp])
    r = list(rows.values())[0]
    assert (r["unique_molecules"], r["molecules_junction"], r["molecules_spanning"], r["start_sites"],
            r["pcr_copies"], r["support_source"]) == ("6", "4", "2", "3", "9", "bpview"), r

    # Rescue: discarded calls come back only when support-filtered and known or STAR-Fusion-confirmed.
    kept = [arriba_row("KEEP", "ME", "chr20:1000", "chr21:500", 5, 5, "high", "in-frame")]
    disc_hdr = ARRIBA_HDR + ["filters"]

    def drow(g1, g2, bp1, bp2, sr, filters, tags="."):
        r = arriba_row(g1, g2, bp1, bp2, sr, 0, "low", "in-frame", tags=tags)
        r["filters"] = filters
        return r
    disc = [drow("EML", "ALKK", "chr2:1000", "chr2:5000", 1, "duplicates(300),in_vitro(1)", tags="Mitelman"),  # known
            drow("TRIM", "BRAFF", "chr7:1000", "chr7:9000", 1, "relative_support(1)"),   # STAR-Fusion confirms
            drow("NOVEL", "ONE", "chr8:1000", "chr9:500", 1, "relative_support(1)"),     # novel: stays out
            drow("ART", "IFACT", "chr10:1000", "chr11:500", 3, "blacklist(3)", tags="Mitelman"),  # artefact filter
            drow("ZERO", "READS", "chr12:1000", "chr13:500", 0, "min_support(1)", tags="Mitelman"),  # no read
            drow("KEEP", "ME", "chr20:1001", "chr21:501", 1, "relative_support(1)", tags="Mitelman"),  # copy of kept
            drow("ONC", "HIT", "chr14:1000", "chr15:500", 1, "min_support(1)"),          # OncoKB hit for the pair
            drow("MIX", "ED", "chr16:1000", "chr17:500", 2, "in_vitro(1),blacklist(1)", tags="Mitelman")]  # mixed
    ka, da, sfp, js = (os.path.join(d, x) for x in ("r.arriba.tsv", "r.disc.tsv", "r.sf.tsv", "r.json"))
    tsv(ka, ARRIBA_HDR + ["filters"], [dict(k, filters="duplicates(2)") for k in kept])
    tsv(da, disc_hdr, disc)
    tsv(sfp, SF_HDR, [sf_row("TRIM", "BRAFF", "chr7:1004", "chr7:9000", 1, 0)])
    with open(js, "w") as fh:
        json.dump([dict(oncokb("LEVEL_1", "Oncogenic"), Fusion="ONC--HIT")], fh)
    fl, outa, outl = (os.path.join(d, x) for x in ("r.list.tsv", "r.out.tsv", "r.out.list.tsv"))
    with open(fl, "w") as fh:
        fh.write("KEEP--ME\n")
    subprocess.run([sys.executable, FE, "rescue", "--arriba-tsv", ka, "--arriba-discarded-tsv", da,
                    "--starfusion-tsv", sfp, "--oncokb-json", js, "--fusion-list", fl,
                    "--out-arriba-tsv", outa, "--out-fusion-list", outl], check=True)
    got = list(csv.DictReader(open(outa), delimiter="\t"))
    assert [(g["#gene1"], g["arriba_status"]) for g in got] == [("KEEP", "kept"), ("EML", "rescued"),
                                                               ("TRIM", "rescued"), ("ONC", "rescued")], got
    assert open(outl).read().split() == ["KEEP--ME", "EML--ALKK", "TRIM--BRAFF", "ONC--HIT"]
    out = os.path.join(d, "r.ev.tsv")
    subprocess.run([sys.executable, FE, "build", "--sample-id", "R", "--arriba-tsv", outa, "--starfusion-tsv", sfp,
                    "--oncokb-json", js, "--protein-domains", ".", "--out-tsv", out], check=True)
    ev = {r["gene5"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    assert ev["EML"]["arriba_status"] == "rescued" and ev["EML"]["arriba_filters"] == "in_vitro", ev["EML"]
    assert ev["EML"]["tier_flags"] == "weak_caller_evidence|low_unique_molecules|low_support|rescued", ev["EML"]
    assert ev["EML"]["tier"] == "2" and ev["EML"]["prototype_rules"] == "rescue|low_support", ev["EML"]
    assert ev["TRIM"]["callers"] == "arriba|starfusion" and ev["TRIM"]["breakpoint_offset"] == "4", ev["TRIM"]
    assert ev["KEEP"]["arriba_status"] == "kept" and ev["KEEP"]["arriba_filters"] == "", ev["KEEP"]
    # ONC is a rescued OncoKB LEVEL_1 pair with one molecule: low support keeps it out of tier 1.
    assert ev["ONC"]["oncokb_highest_level"] == "LEVEL_1" and ev["ONC"]["tier"] != "1", ev["ONC"]
    # Order within tier 3: supported (KEEP, 10 molecules) before low support (TRIM).
    assert int(ev["KEEP"]["event_rank"]) < int(ev["TRIM"]["event_rank"]), ev
    # A missing OncoKB file (a directory, or '.') is a warning, not a crash, in rescue and build.
    subprocess.run([sys.executable, FE, "rescue", "--arriba-tsv", ka, "--arriba-discarded-tsv", da,
                    "--starfusion-tsv", sfp, "--oncokb-json", ".", "--out-arriba-tsv", outa], check=True)
    subprocess.run([sys.executable, FE, "build", "--sample-id", "R", "--arriba-tsv", outa, "--starfusion-tsv", sfp,
                    "--oncokb-json", d, "--protein-domains", ".", "--out-tsv", out], check=True)
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    # no OncoKB annotation at all; Arriba's Mitelman tag still makes EML a curated tier 2
    assert rows and all(r["oncokb_highest_level"] == "" for r in rows), rows
    assert {r["gene5"]: r["tier"] for r in rows}["EML"] == "2"
    # An unoriented call's ID has '.' strands.
    a = arriba_row("U1", "U2", "chr1:10", "chr2:20", 2, 2, "high", ".", strands="./.")
    _, rows = run_build(d, "U", [a], [], {})
    assert list(rows.values())[0]["event_id"] == "chr1:10:.>chr2:20:.", rows


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

        # ALK1--PART1 is in 1 of 3 other PoN samples (33%): not recurrent at 50%; RECUR--PAIR is in 3/3.
        _, rows = run_build(d, "S1", arriba, sf, pairs, ["--pon", pon, "--pon-min-fraction", "0.5"])
        assert len(rows) == 4, sorted(rows)
        alk = rows["ALK1--PART1@chr1:1000:+"]
        assert alk["callers"] == "arriba|starfusion" and alk["n_callers"] == "2", alk
        assert alk["caller_arriba"] == "true" and alk["caller_starfusion"] == "true"
        assert alk["breakpoint3"] == "chr2:500:+" and alk["reading_frame"] == "in-frame"
        assert alk["arriba_split_reads"] == "7" and alk["starfusion_junction_reads"] == "6"
        assert alk["unique_molecules"] == "7" and alk["support_source"] == "arriba" and alk["pcr_copies"] == ""
        assert alk["tier"] == "1" and alk["oncokb_highest_level"] == "LEVEL_1", alk["tier_reason"]
        assert alk["tier_basis"] == "oncokb_actionable" and alk["tier_flags"] == "" and alk["pon_recurrent"] == "false"
        assert alk["pon_gene_pair_count"] == "1" and alk["pon_junction_count"] == "1" and alk["pon_size"] == "3"
        assert alk["event_id"] == "chr1:1000:+>chr2:500:+" and alk["event_rank"] == "1", alk  # tier 1 first
        assert alk["breakpoint_offset"] == "3" and alk["low_support"] == "false", alk

        junk = rows["NTRKX--JUNK@chr3:1000:+"]
        assert junk["tier"] == "3" and "but one caller" in junk["tier_reason"], junk["tier_reason"]
        assert junk["tier_basis"] == "none" and junk["tier_flags"] == "weak_caller_evidence|low_support", junk["tier_flags"]
        assert junk["low_support"] == "true" and junk["prototype_rules"] == "low_support"
        assert junk["support_source"] == "starfusion" and junk["breakpoint_offset"] == "", junk
        assert junk["caller_arriba"] == "false" and junk["caller_starfusion"] == "true"
        known = rows["KNOWN--PAIR@chr7:1000:+"]
        assert known["tier"] == "2" and known["tier_reason"].startswith("known fusion: Mitelman; evidence:"), known["tier_reason"]
        assert known["tier_flags"] == "weak_caller_evidence|low_unique_molecules|low_support", known["tier_flags"]
        assert known["tier_basis"] == "curated_known" and known["curated_known"] == "true"
        assert known["known_db"] == "Mitelman", known["known_db"]
        recur = rows["RECUR--PAIR@chr5:1000:+"]
        assert recur["pon_gene_pair_count"] == "3" and recur["tier"] == "3", recur
        assert recur["pon_recurrent"] == "true" and recur["tier_flags"] == "pon_recurrent|pon_capped", recur["tier_flags"]
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
        assert alk["tier_basis"] == "oncokb_actionable" and alk["tier_flags"] == "low_unique_molecules"
        assert alk["pon_recurrent"] == "" and alk["pon_size"] == ""  # no PoN given

        # Annotation comes only from the gene a breakpoint is in. Arriba between genes: ALK1 is
        # 400 kb away, so ALK1--PART1's LEVEL_1 must not apply, unless the distance is allowed.
        far = [arriba_row("ALK1(400000),OTHER(500000)", "PART1", "chr9:1000", "chr10:500", 6, 6, "high", "in-frame")]
        _, rows = run_build(d, "S3", far, [], pairs)
        r = rows["ALK1--PART1@chr9:1000:+"] if "ALK1--PART1@chr9:1000:+" in rows else list(rows.values())[0]
        assert r["tier"] == "3" and r["oncokb_highest_level"] == "" and r["fusion"] == "ALK1--PART1", r
        _, rows = run_build(d, "S3", far, [], pairs, ["--max-intergenic-distance", "400000"])
        assert list(rows.values())[0]["tier"] == "1", list(rows.values())[0]["tier_reason"]

        # An unoriented Arriba call (./.) that matches only the reversed pair: no tier 1, and no
        # curated tier 2 without evidence.
        rev_pairs = {"ALK2--PART2": oncokb("LEVEL_1", "Oncogenic"), "CUR1--CUR2": oncokb(dbs=["Mitelman"])}
        rev = [arriba_row("PART2", "ALK2", "chr11:1000", "chr12:500", 6, 6, "high", "in-frame", strands="./."),
               arriba_row("CUR2", "CUR1", "chr13:1000", "chr14:500", 1, 0, "low", ".", strands="./.")]
        _, rows = run_build(d, "S4", rev, [], rev_pairs)
        act = [r for r in rows.values() if r["fusion"] == "ALK2--PART2"][0]
        assert act["orientation_known"] == "false" and act["tier"] == "2", act
        assert act["tier_flags"] == "orientation_unknown" and act["tier_basis"] == "oncokb_actionable", act["tier_flags"]
        cur = [r for r in rows.values() if r["fusion"] == "CUR1--CUR2"][0]
        assert cur["tier"] == "3" and "orientation_unknown" in cur["tier_flags"], cur

        test_contract_v1(d)

        # Empty caller outputs give an empty table with a header.
        _, rows = run_build(d, "S2", [], [], {})
        assert rows == {}
    print("test_fusion_events: ok")


if __name__ == "__main__":
    main()
