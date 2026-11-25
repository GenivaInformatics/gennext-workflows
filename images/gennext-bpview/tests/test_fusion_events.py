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
    # 0.6.0 rescue rule: a rescued single molecule stays tier III unless OncoKB curates the exact pair
    # (a Mitelman tag alone no longer keeps it in tier II)
    assert ev["EML"]["tier_flags"] == "weak_caller_evidence|low_unique_molecules|rescue_capped|low_support|rescued", ev["EML"]
    assert ev["EML"]["tier"] == "3" and ev["EML"]["prototype_rules"] == "rescue|low_support", ev["EML"]
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
    # no OncoKB annotation at all: no crash, and the rescued single-molecule EML stays tier III (0.6.0)
    assert rows and all(r["oncokb_highest_level"] == "" for r in rows), rows
    assert {r["gene5"]: r["tier"] for r in rows}["EML"] == "3"
    # An unoriented call's ID has '.' strands.
    a = arriba_row("U1", "U2", "chr1:10", "chr2:20", 2, 2, "high", ".", strands="./.")
    _, rows = run_build(d, "U", [a], [], {})
    assert list(rows.values())[0]["event_id"] == "chr1:10:.>chr2:20:.", rows


def test_caller_details(d):
    """0.5.0 detail columns come from each caller's own row and the OncoKB entry, unchanged."""
    hdr = ARRIBA_HDR + ["retained_protein_domains", "closest_genomic_breakpoint1", "closest_genomic_breakpoint2",
                        "gene_id1", "gene_id2", "transcript_id1", "transcript_id2", "filters",
                        "fusion_transcript", "peptide_sequence"]
    a = arriba_row("ALK1", "PART1", "chr1:1000", "chr2:500", 3, 26, "high", "in-frame", tags="Mitelman,CCLE")
    a.update({"site1": "CDS/splice-site", "coverage1": 33, "coverage2": 26,
              "retained_protein_domains": "Ig_domain(100%),Protein_kinase_domain(100%)|TACC_C-terminal(100%)",
              "closest_genomic_breakpoint1": ".", "gene_id1": "ENSG1.1", "transcript_id1": "ENST1.1",
              "filters": "duplicates(448)", "fusion_transcript": "ACGT|TTGA", "peptide_sequence": "Q" * 2500})
    sfh = SF_HDR + ["est_J", "est_S", "LeftBreakDinuc", "LeftBreakEntropy", "RightBreakDinuc", "RightBreakEntropy"]
    sf = dict(sf_row("ALK1", "PART1", "chr1:1000", "chr2:500", 29, 0), est_J="29.00", est_S="0.00",
              LeftBreakDinuc="GT", LeftBreakEntropy="1.95", RightBreakDinuc="AG", RightBreakEntropy="1.8")
    pairs = {"ALK1--PART1": dict(oncokb("LEVEL_3B", "Oncogenic", dbs=["Mitelman", "FusionGDB2"]),
                                 **{"Explained FII": "0.5 * (2 / 2) + 0.5 * (0.5 / 0.5)",
                                    # the OncoKB JSON's own arriba block is from another run: must not leak
                                    "arriba": {"split_reads1": "166", "split_reads2": "300", "coverage1": "999"}})}
    pairs["ALK1--PART1"]["OncoKB"].update({"LEVEL_3B": "Erdafitinib", "LEVEL_4": "Fexagratinib", "LEVEL_Dx2": "",
                                           "HIGHEST_SENSITIVE_LEVEL": "LEVEL_3B", "HIGHEST_RESISTANCE_LEVEL": "",
                                           "HIGHEST_PX_LEVEL": "", "TX_CITATIONS": "37870920;31088831",
                                           "MUTATION_EFFECT_CITATIONS": "22837387", "GENE_IN_ONCOKB": "True",
                                           "VARIANT_IN_ONCOKB": "False"})
    ap, sp, jp, out = (os.path.join(d, x) for x in ("cd.arriba.tsv", "cd.sf.tsv", "cd.json", "cd.ev.tsv"))
    tsv(ap, hdr, [a])
    tsv(sp, sfh, [sf])
    with open(jp, "w") as fh:
        json.dump([dict(v, Fusion=k) for k, v in pairs.items()], fh)
    subprocess.run([sys.executable, FE, "build", "--sample-id", "CD", "--arriba-tsv", ap, "--starfusion-tsv", sp,
                    "--oncokb-json", jp, "--protein-domains", ".", "--out-tsv", out], check=True)
    r = list(csv.DictReader(open(out), delimiter="\t"))[0]
    want = {"arriba_split_reads1": "3", "arriba_split_reads2": "26", "arriba_split_reads": "29",
            "arriba_coverage1": "33", "arriba_coverage2": "26", "arriba_site1": "CDS/splice-site",
            "arriba_filters_raw": "duplicates(448)", "arriba_duplicates": "448", "arriba_tags": "Mitelman|CCLE",
            "arriba_retained_domains1": "Ig_domain(100%)|Protein_kinase_domain(100%)",
            "arriba_retained_domains2": "TACC_C-terminal(100%)", "arriba_closest_genomic_breakpoint1": "",
            "arriba_gene_id1": "ENSG1.1", "arriba_transcript_id1": "ENST1.1", "arriba_gene1": "ALK1",
            "arriba_breakpoint1": "chr1:1000", "arriba_strand1": "+/+",
            "arriba_fusion_transcript": "ACGT|TTGA",
            "arriba_peptide_sequence": "Q" * 2000 + "...[truncated, 2500 characters]",
            "starfusion_left_gene": "ALK1^ENSG1", "starfusion_left_breakpoint": "chr1:1000:+",
            "starfusion_est_junction": "29.0", "starfusion_splice_type": "ONLY_REF_SPLICE",
            "starfusion_large_anchor_support": "true", "starfusion_left_break_dinuc": "GT",
            "starfusion_right_break_entropy": "1.8", "starfusion_annots": "Mitelman",
            "oncokb_highest_sensitive_level": "LEVEL_3B", "oncokb_highest_resistance_level": "",
            "oncokb_treatments": "LEVEL_3B: Erdafitinib|LEVEL_4: Fexagratinib",
            "oncokb_treatment_citations": "37870920|31088831", "oncokb_mutation_effect_citations": "22837387",
            "oncokb_dx_citations": "", "oncokb_gene_in_oncokb": "true", "oncokb_variant_in_oncokb": "false",
            "pair_fii_explained": "0.5 * (2 / 2) + 0.5 * (0.5 / 0.5)",
            "unique_molecules": "29", "support_source": "arriba"}
    bad = {k: (r[k], v) for k, v in want.items() if r[k] != v}
    assert not bad, bad
    # A STAR-Fusion-only event has no Arriba details (null, lists empty).
    tsv(ap, hdr, [])
    subprocess.run([sys.executable, FE, "build", "--sample-id", "CD", "--arriba-tsv", ap, "--starfusion-tsv", sp,
                    "--oncokb-json", jp, "--protein-domains", ".", "--out-tsv", out], check=True)
    r = list(csv.DictReader(open(out), delimiter="\t"))[0]
    assert r["arriba_split_reads1"] == "" and r["arriba_duplicates"] == "" and r["arriba_tags"] == "", r
    assert r["starfusion_left_break_dinuc"] == "GT", r


def test_caller_sides():
    """arriba_split_reads_5p/_3p follow the event's orientation; *_bam_* come from bpview's callers_bam."""
    raw = {"breakpoint1": "chr2:500", "breakpoint2": "chr1:1000", "split_reads1": "7", "split_reads2": "2"}
    ev = dict(chr1="chr1", pos1=1000, chr2="chr2", pos2=500, callers=["arriba", "starfusion"],
              arriba={"raw": raw})  # Arriba's side 1 is the event's 3' (a flipped call)
    b = {"callers_bam": {"arriba": {"fragments_nondup": {"junction": 5, "spanning": 1, "junction_weak": 2},
                                    "fragments_dup": {"junction": 40}, "split_5p": 3, "split_3p": 2,
                                    "cited": 50, "cited_missing": 0},
                         "starfusion": {"fragments_nondup": {"junction": 4}, "fragments_dup": {},
                                        "split_5p": 4, "split_3p": 0, "cited": 4, "cited_missing": 0}}}
    out = fusion_events.caller_sides(ev, b)
    assert (out["arriba_split_reads_5p"], out["arriba_split_reads_3p"]) == (2, 7), out
    assert (out["arriba_bam_molecules"], out["arriba_bam_pcr_copies"]) == (6, 40), out  # junction_weak not counted
    assert (out["arriba_bam_split_reads_5p"], out["arriba_bam_split_reads_3p"]) == (3, 2), out
    assert (out["starfusion_bam_molecules"], out["starfusion_bam_split_reads_5p"], out["starfusion_cited_missing"]) == (4, 4, 0)
    # not flipped; and a caller that didn't call the event gets no BAM columns
    raw.update(breakpoint1="chr1:1000", breakpoint2="chr2:500")
    out = fusion_events.caller_sides(dict(ev, callers=["arriba"]), b)
    assert (out["arriba_split_reads_5p"], out["arriba_split_reads_3p"]) == (7, 2), out
    assert "starfusion_bam_molecules" not in out, out
    # a caller that cites no names (a rescued Arriba call: discarded calls carry no read_identifiers)
    b["callers_bam"]["arriba"]["cited"] = 0
    out = fusion_events.caller_sides(ev, b)
    assert "arriba_bam_molecules" not in out and out["arriba_split_reads_5p"] == 7, out


def test_tier_rules_0_6(d):
    """0.6.0 (Berk, 25 Nov), on the shapes of the real samples: the rescue rule, the order within a tier,
    and the MANE-derived reading frame."""
    gtf = os.path.join(d, "t6.gtf")
    with open(gtf, "w") as fh:
        def gene(chrom, name):
            for n, (a, b) in enumerate([(100, 199), (300, 399), (500, 599)], 1):
                for feat in ("exon", "CDS"):
                    fh.write(f'{chrom}\tT\t{feat}\t{a}\t{b}\t.\t+\t0\tgene_name "{name}"; transcript_id "{name}-201"; '
                             f'exon_number {n}; tag "MANE_Select";\n')
        for i, g in enumerate(("EK", "AK", "NT", "MY", "FG", "CI", "CL", "AR", "FR", "TA")):
            gene(f"chr{i + 1}", g)
    hdr = ARRIBA_HDR + ["filters", "arriba_status"]

    def row(g5, g3, c5, c3, p5, p3, split, conf, frame, status, filters=".", tags="."):
        r = arriba_row(g5, g3, f"{c5}:{p5}", f"{c3}:{p3}", split, 0, conf, frame, tags=tags)
        r.update(filters=filters, arriba_status=status)
        return r
    arr = [
        # EML4-ALK shape: rescued, 1 molecule, exact OncoKB pair, both at exon boundaries, no Arriba frame
        row("EK", "AK", "chr1", "chr2", 199, 300, 1, "low", ".", "rescued", "in_vitro(1)"),
        # partner-agnostic kinase hit: rescued, 2 molecules (at the tier I bar), partner-agnostic LEVEL_1, both callers
        row("NT", "MY", "chr3", "chr4", 150, 350, 2, "low", ".", "rescued", "relative_support(1)"),
        # a partner-agnostic single molecule, rescued, both callers
        row("FG", "CI", "chr5", "chr6", 150, 350, 1, "low", ".", "rescued", "relative_support(1)"),
        # CLDN18-ARHGAP26 shape: kept, 84 molecules, curated (Mitelman), 5' breakpoint mid-exon
        row("CL", "AR", "chr7", "chr8", 350, 300, 84, "high", "in-frame", "kept", tags="Mitelman"),
        # FGFR3-TACC3 shape: kept, 29 molecules, exact OncoKB pair, canonical junction
        row("FR", "TA", "chr9", "chr10", 399, 500, 29, "high", "in-frame", "kept"),
    ]
    sf = [sf_row("NT", "MY", "chr3:150", "chr4:350", 2, 0), sf_row("FG", "CI", "chr5:150", "chr6:350", 1, 0)]
    for r in sf:
        r["annots"] = "."

    def okb(level, onc, exact, dbs=()):
        e = oncokb(level, onc, dbs)
        e["OncoKB"]["VARIANT_IN_ONCOKB"] = "True" if exact else "False"
        return e
    pairs = {"EK--AK": okb("LEVEL_1", "Oncogenic", True, ["Mitelman"]),
             "NT--MY": okb("LEVEL_1", "Likely Oncogenic", False),
             "FG--CI": okb("LEVEL_2", "Likely Oncogenic", False),
             "CL--AR": okb("", "Unknown", False, ["Mitelman"]),
             "FR--TA": okb("LEVEL_3B", "Oncogenic", True)}
    ap, sp, jp, out = (os.path.join(d, x) for x in ("t6.arriba.tsv", "t6.sf.tsv", "t6.json", "t6.ev.tsv"))
    tsv(ap, hdr, arr)
    tsv(sp, SF_HDR, sf)
    with open(jp, "w") as fh:
        json.dump([dict(v, Fusion=k) for k, v in pairs.items()], fh)
    subprocess.run([sys.executable, FE, "build", "--sample-id", "T6", "--arriba-tsv", ap, "--starfusion-tsv", sp,
                    "--oncokb-json", jp, "--gtf", gtf, "--protein-domains", ".", "--out-tsv", out], check=True)
    ev = {r["fusion"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    ek, nt, fg, cl, fr = (ev[k] for k in ("EK--AK", "NT--MY", "FG--CI", "CL--AR", "FR--TA"))
    # EML4-ALK shape: stays tier II (OncoKB curates the exact pair), frame derived from MANE (n5 = 100, n3 = 100)
    assert ek["tier"] == "2" and "rescue_capped" not in ek["tier_flags"], ek
    assert (ek["reading_frame"], ek["reading_frame_source"]) == ("in-frame", "mane"), ek
    assert ek["oncokb_pair_specific"] == "true" and ek["canonical_junction"] == "true", ek
    # partner-agnostic kinase hit with 2 molecules: rescued without an exact OncoKB pair -> tier III
    assert nt["tier"] == "3" and "rescue_capped" in nt["tier_flags"], nt
    assert "without OncoKB curation of this exact gene pair: tier III" in nt["tier_reason"], nt["tier_reason"]
    assert nt["oncokb_pair_specific"] == "false" and nt["reading_frame"] == "", nt  # mid-exon: not computable
    # a partner-agnostic single molecule, rescued: tier III
    assert fg["tier"] == "3" and "rescue_capped" in fg["tier_flags"], fg
    assert "without OncoKB curation of this exact gene pair: tier III" in fg["tier_reason"], fg["tier_reason"]
    # CLDN18-ARHGAP26 shape: Arriba's frame, source arriba; not canonical (mid-exon), so it doesn't lead
    assert cl["tier"] == "2" and (cl["reading_frame"], cl["reading_frame_source"]) == ("in-frame", "arriba"), cl
    assert cl["site5"] == "exon" and cl["canonical_junction"] == "false", cl
    # Order in tier II: canonical exact pairs first (FGFR3-TACC3 shape, then EML4-ALK shape), then the rest
    # by the existing order (CLDN18-ARHGAP26 shape)
    order = sorted((r for r in ev.values() if r["tier"] == "2"), key=lambda r: int(r["event_rank"]))
    assert [r["fusion"] for r in order] == ["FR--TA", "EK--AK", "CL--AR"], [r["fusion"] for r in order]
    # an exact OncoKB pair that is rescued with enough molecules still never reaches tier I
    arr2 = [row("EK", "AK", "chr1", "chr2", 199, 300, 5, "high", ".", "rescued", "relative_support(1)")]
    tsv(ap, hdr, arr2)
    subprocess.run([sys.executable, FE, "build", "--sample-id", "T6", "--arriba-tsv", ap, "--starfusion-tsv", sp,
                    "--oncokb-json", jp, "--gtf", gtf, "--protein-domains", ".", "--out-tsv", out], check=True)
    r = [x for x in csv.DictReader(open(out), delimiter="\t") if x["fusion"] == "EK--AK"][0]
    assert r["tier"] == "2" and "never tier I" in r["tier_reason"], r
    # no counted molecule: tier III through any path, a kept curated call included (Berk, 25 Nov)
    bp = os.path.join(d, "t6.events.json")
    with open(bp, "w") as fh:  # bpview drew the CLDN18-ARHGAP26-shaped call and found no molecule
        json.dump({"events": [{"breakpoint1": "chr7:350", "breakpoint2": "chr8:300", "keep1": "le", "keep2": "ge",
                               "contig": "CL--AR", "fragments_nondup": {"wt5": 6}, "fragments_dup": {"wt5": 141}}]}, fh)
    tsv(ap, hdr, arr)
    subprocess.run([sys.executable, FE, "build", "--sample-id", "T6", "--arriba-tsv", ap, "--starfusion-tsv", sp,
                    "--oncokb-json", jp, "--gtf", gtf, "--protein-domains", ".", "--bpview-events", bp,
                    "--out-tsv", out], check=True)
    r = [x for x in csv.DictReader(open(out), delimiter="\t") if x["fusion"] == "CL--AR"][0]
    assert r["unique_molecules"] == "0" and r["tier"] == "3" and "no_molecules" in r["tier_flags"], r
    assert "no counted molecule: tier III" in r["tier_reason"], r["tier_reason"]
    # mane_frame: out of frame when the phases differ; UTR or intronic breakpoints give None
    t = fe_tx = fusion_events.load_transcripts(gtf)
    e = dict(pos1=199, keep1="le", pos2=100, keep2="ge")
    assert fusion_events.mane_frame(t["EK"], t["AK"], e, "exon_boundary", "exon_boundary") == "out-of-frame"
    assert fusion_events.mane_frame(t["EK"], t["AK"], dict(e, pos1=250), "intron", "exon_boundary") is None
    assert fusion_events.mane_frame(t["EK"], t["AK"], dict(e, pos2=50), "exon_boundary", "exon_boundary") is None


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
        test_caller_details(d)
        test_caller_sides()
        test_tier_rules_0_6(d)

        # Empty caller outputs give an empty table with a header.
        _, rows = run_build(d, "S2", [], [], {})
        assert rows == {}
    print("test_fusion_events: ok")


if __name__ == "__main__":
    main()
