"""Synthetic end-to-end test: two genes on two chromosomes, a fusion G1(<=1000) + G2(>=501).

Run: python tests/test_bpview.py   (needs pysam; bpview.py one directory up)
"""
import json
import os
import random
import subprocess
import sys
import tempfile

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
BPVIEW = os.path.join(HERE, "..", "bpview.py")

random.seed(1)
CHR_A = "".join(random.choice("ACGT") for _ in range(3000))
CHR_B = "".join(random.choice("ACGT") for _ in range(3000))
BP1, BP2 = 1000, 501  # 1-based: G1 keeps 1..1000, G2 keeps 501..


def write_inputs(d):
    fa = os.path.join(d, "genome.fa")
    with open(fa, "w") as fh:
        fh.write(">chrA\n" + CHR_A + "\n>chrB\n" + CHR_B + "\n")
    pysam.faidx(fa)
    gtf = os.path.join(d, "genes.gtf")
    with open(gtf, "w") as fh:
        for chrom, g in (("chrA", "G1"), ("chrB", "G2")):
            attr = f'gene_id "{g}"; gene_name "{g}"; transcript_id "{g}.1";'
            fh.write(f"{chrom}\tt\tgene\t100\t2900\t.\t+\t.\t{attr}\n")
            fh.write(f"{chrom}\tt\texon\t100\t2900\t.\t+\t.\t{attr}\n")
    hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                           "SQ": [{"SN": "chrA", "LN": 3000}, {"SN": "chrB", "LN": 3000}]})
    recs = []

    def rec(name, tid, start, cigar, seq, r1, rev, mate=None, dup=False, sa=None, supp=False):
        r = pysam.AlignedSegment(hdr)
        r.query_name, r.reference_id, r.reference_start = name, tid, start
        r.cigarstring, r.query_sequence = cigar, seq
        r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        r.mapping_quality = 60
        r.is_paired, r.is_read1, r.is_read2 = True, r1, not r1
        r.is_reverse, r.is_duplicate, r.is_supplementary = rev, dup, supp
        if mate:
            r.next_reference_id, r.next_reference_start = mate
        if sa:
            r.set_tag("SA", sa)
        recs.append(r)

    # junction read: 50 bp of chrA 951..1000 then 50 bp of chrB 501..550
    jseq = CHR_A[950:1000] + CHR_B[500:550]
    for k, dup in ((0, False), (1, True)):
        n = f"junc{k}"
        rec(n, 0, 950, "50M50S", jseq, True, False, mate=(1, 700), dup=dup, sa="chrB,501,+,50S50M,60,0;")
        rec(n, 1, 500, "50H50M", CHR_B[500:550], True, False, mate=(1, 700), dup=dup,
            sa="chrA,951,+,50M50S,60,0;", supp=True)
        rec(n, 1, 700, "100M", CHR_B[700:800], False, True, mate=(0, 950), dup=dup)
    # a fragment whose two mates both cross the junction (overlapping mates): one molecule, one start site
    s1 = CHR_A[940:1000] + CHR_B[500:540]
    rec("jj", 0, 940, "60M40S", s1, True, False, mate=(0, 960), sa="chrB,501,+,60S40M,60,0;")
    rec("jj", 1, 500, "60H40M", CHR_B[500:540], True, False, mate=(0, 960), sa="chrA,941,+,60M40S,60,0;", supp=True)
    s2 = CHR_A[960:1000] + CHR_B[500:560]
    rec("jj", 0, 960, "40M60S", s2, False, True, mate=(0, 940), sa="chrB,501,-,40S60M,60,0;")
    rec("jj", 1, 500, "40H60M", CHR_B[500:560], False, True, mate=(0, 940), sa="chrA,961,-,40M60S,60,0;", supp=True)
    # spanning pair: chrA 801..900 / chrB 601..700
    rec("span", 0, 800, "100M", CHR_A[800:900], True, False, mate=(1, 600))
    rec("span", 1, 600, "100M", CHR_B[600:700], False, True, mate=(0, 800))
    # wild-type G2 read crossing bp2: chrB 451..550
    rec("wt", 1, 450, "100M", CHR_B[450:550], True, False, mate=(1, 650))
    rec("wt", 1, 650, "100M", CHR_B[650:750], False, True, mate=(1, 450))
    bam = os.path.join(d, "in.bam")
    with pysam.AlignmentFile(bam + ".u", "wb", header=hdr) as o:
        for r in recs:
            o.write(r)
    pysam.sort("-o", bam, bam + ".u")
    pysam.index(bam)
    arriba = os.path.join(d, "arriba.tsv")
    with open(arriba, "w") as fh:
        fh.write("#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\t"
                 "split_reads1\tsplit_reads2\tdiscordant_mates\tconfidence\tdirection1\tdirection2\ttags\t"
                 "read_identifiers\n")
        # Arriba cites junc0, its copy junc1, the spanning pair, the wild-type read and a name not in the BAM
        fh.write(f"G1\tG2\t+/+\t+/+\tchrA:{BP1}\tchrB:{BP2}\t1\t0\t1\thigh\tdownstream\tupstream\tMitelman\t"
                 "junc0,junc1,span,wt,ghost\n")
    sf = os.path.join(d, "sf.tsv")
    with open(sf, "w") as fh:
        fh.write("#FusionName\tJunctionReadCount\tSpanningFragCount\tSpliceType\tLeftGene\tLeftBreakpoint\t"
                 "RightGene\tRightBreakpoint\tFFPM\tLargeAnchorSupport\tannots\tJunctionReads\tSpanningFrags\n")
        # STAR-Fusion cites junc0 and jj, in its own name format (UMI suffix, '&...@' prefix)
        fh.write(f"G1--G2\t1\t1\tONLY_REF_SPLICE\tG1^X\tchrA:{BP1 + 3}:+\tG2^Y\tchrB:{BP2 - 2}:+\t1.0\tYES\t"
                 '["G1:Oncogene"];["Cosmic","Mitelman"];INTERCHROMOSOMAL[chrA--chrB]\t'
                 "&x/1@junc0:UMI_ACGTA_TTGAC,jj/1\t.\n")
    return fa, gtf, bam, arriba, sf


def test_gene_names():
    sys.path.insert(0, os.path.join(HERE, ".."))
    import bpview
    assert bpview.gene_names("ENSG0001(731),ABL1(120799)") == ["ENSG0001", "ABL1"]
    assert bpview.gene_names("ABL1(120799),ENSG0001(731)") == ["ENSG0001", "ABL1"]
    assert bpview.gene_names("FGFR3") == ["FGFR3"]
    assert bpview.gene_names("TACC3^ENSG00000013810.20") == ["TACC3"]
    assert bpview.gene_names("RP11-1.2") == ["RP11-1.2"]


def test_merge():
    """Nearest-call merge and one event per STAR-Fusion junction."""
    sys.path.insert(0, os.path.join(HERE, ".."))
    import bpview

    def a(p1, p2, oriented=True):
        return dict(gene1="A", gene2="B", names1=["A"], names2=["B"], intergenic1=False, intergenic2=False,
                    chr1="chr1", pos1=p1, keep1="le", chr2="chr2", pos2=p2, keep2="ge", oriented=oriented,
                    arriba=dict(confidence="high", split_reads1=1, split_reads2=1, discordant_mates=0))

    def s(p1, p2, g1="A", j=1):
        return dict(gene1=g1, gene2="B", names1=[g1], names2=["B"], intergenic1=False, intergenic2=False,
                    chr1="chr1", pos1=p1, keep1="le", chr2="chr2", pos2=p2, keep2="ge", oriented=True,
                    starfusion=dict(junction_reads=j, spanning_frags=0, ffpm=0.1,
                                    left_breakpoint=f"chr1:{p1}:+", right_breakpoint=f"chr2:{p2}:+"))

    # Two Arriba calls 8 bp apart; the STAR-Fusion call is 1 bp from the second one. The first
    # match within 10 bp (the old rule) would have joined it to the first Arriba call.
    ev = bpview.merge_events([a(1000, 500), a(1008, 500)], [s(1009, 500)])
    assert len(ev) == 2, ev
    assert ev[0]["callers"] == ["arriba"] and ev[1]["callers"] == ["arriba", "starfusion"], ev
    # STAR-Fusion's two rows for one junction (alternative gene annotation) are one event, labelled
    # by the row with more support, with the other name as an alias; not a second, STAR-Fusion-only event.
    ev = bpview.merge_events([a(1000, 500)], [s(1000, 500, "A", j=1), s(1000, 500, "A-AS1", j=3)])
    assert len(ev) == 1 and ev[0]["callers"] == ["arriba", "starfusion"], ev
    assert ev[0]["starfusion"]["junction_reads"] == 3 and ev[0]["names1"] == ["A-AS1", "A"], ev[0]["names1"]
    ev = bpview.merge_events([], [s(1000, 500, "A"), s(1000, 500, "A2")])
    assert len(ev) == 1 and ev[0]["names1"] == ["A", "A2"], ev
    # One Arriba call, two STAR-Fusion junctions near it: the nearer one joins, the other stays apart.
    ev = bpview.merge_events([a(1000, 500)], [s(1006, 500), s(1002, 500)])
    assert len(ev) == 2 and ev[0]["starfusion"]["left_breakpoint"] == "chr1:1002:+", ev


def test_single_record_junctions():
    """0.4.1: a read-through read (one spliced record, N gap from bp1 to bp2) is a junction read; an ITD's
    wild-type reads, which overlap both segments, are not; an ITD's chimeric split read still is."""
    sys.path.insert(0, os.path.join(HERE, ".."))
    import bpview
    random.seed(7)
    C = "".join(random.choice("ACGT") for _ in range(4000))
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "g.fa")
        with open(fa, "w") as fh:
            fh.write(">chrC\n" + C + "\n")
        pysam.faidx(fa)
        gtf = os.path.join(d, "g.gtf")
        with open(gtf, "w") as fh:
            for g, a, b in (("RT5", 100, 1000), ("RT3", 2001, 3400), ("ITD", 500, 2500)):
                attr = f'gene_id "{g}"; gene_name "{g}"; transcript_id "{g}.1";'
                fh.write(f"chrC\tt\tgene\t{a}\t{b}\t.\t+\t.\t{attr}\n")
                fh.write(f"chrC\tt\texon\t{a}\t{b}\t.\t+\t.\t{attr}\n")
        hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                               "SQ": [{"SN": "chrC", "LN": 4000}]})
        recs = []

        def rec(name, start, cigar, seq, r1=True, rev=False, sa=None, supp=False):
            r = pysam.AlignedSegment(hdr)
            r.query_name, r.reference_id, r.reference_start = name, 0, start
            r.cigarstring, r.query_sequence = cigar, seq
            r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            r.mapping_quality, r.is_paired, r.is_read1, r.is_read2 = 255, False, r1, False
            r.is_reverse, r.is_supplementary = rev, supp
            if sa:
                r.set_tag("SA", sa)
            recs.append(r)
        # read-through RT5(..1000) -> RT3(2001..): one spliced record each, like STAR's 8M23291N106M6143N31M
        rec("rt_a", 950, "50M1000N50M", C[950:1000] + C[2000:2050])   # 50 + 50 aligned: side 5p (tie)
        rec("rt_b", 970, "30M1000N70M", C[970:1000] + C[2000:2070])   # 30 + 70: side 3p
        rec("rt_wt", 920, "100M", C[920:1020])                        # wild-type RT5 read: not a junction
        # ITD inside ITD (500..2500): 5' part up to 2000, 3' part again from 1001 (a duplication)
        rec("itd_wt1", 950, "100M", C[950:1050])                      # overlaps both segments: wild-type
        rec("itd_wt2", 1950, "100M", C[1950:2050])
        s1 = C[1950:2000] + C[1000:1050]                               # chimeric split read across the ITD
        rec("itd_j", 1950, "50M50S", s1, sa="chrC,1001,+,50S50M,60,0;")
        rec("itd_j", 1000, "50H50M", C[1000:1050], sa="chrC,1951,+,50M50S,60,0;", supp=True)
        bam = os.path.join(d, "in.bam")
        with pysam.AlignmentFile(bam + ".u", "wb", header=hdr) as o:
            for r in recs:
                o.write(r)
        pysam.sort("-o", bam, bam + ".u")
        pysam.index(bam)
        arriba = os.path.join(d, "a.tsv")
        with open(arriba, "w") as fh:
            fh.write("#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\t"
                     "split_reads1\tsplit_reads2\tdiscordant_mates\tconfidence\tdirection1\tdirection2\ttype\n")
            fh.write("RT5\tRT3\t+/+\t+/+\tchrC:1000\tchrC:2001\t2\t0\t0\tlow\tdownstream\tupstream\tdeletion/read-through\n")
            fh.write("ITD\tITD\t+/+\t+/+\tchrC:2000\tchrC:1001\t1\t0\t0\tlow\tdownstream\tupstream\tduplication/ITD\n")
        out = os.path.join(d, "out")
        subprocess.run([sys.executable, BPVIEW, "--sample-id", "T", "--bam", bam, "--fasta", fa, "--gtf", gtf,
                        "--arriba-tsv", arriba, "--out-prefix", out], check=True)
        doc = json.load(open(out + ".events.json"))
        ev = {e["gene1"]: e for e in doc["events"]}
        rt, itd = ev["RT5"], ev["ITD"]
        assert rt["fragments_nondup"].get("junction") == 2, rt["fragments_nondup"]
        # (the ITD's reads also cross the read-through breakpoints, so the wild-type counts include them)
        assert itd["fragments_nondup"].get("junction") == 1, itd["fragments_nondup"]
        assert not itd["fragments_nondup"].get("junction_weak"), itd["fragments_nondup"]
        b = pysam.AlignmentFile(out + ".bam")
        by = {}
        for r in b.fetch(rt["contig"]):
            by[r.query_name] = r
        j = rt["junction_1based"]
        # drawn crossing the junction, as one alignment: 50M50M merged to 100M at the junction
        assert by["rt_a"].get_tag("XF") == "junction" and by["rt_a"].get_tag("XG") == "1junction"
        assert by["rt_a"].reference_start == j - 50 and by["rt_a"].cigarstring == "100M", by["rt_a"].cigarstring
        assert by["rt_a"].get_tag("XJ") == "5p" and by["rt_b"].get_tag("XJ") == "3p"
        assert by["rt_wt"].get_tag("XF") == "wt5", by["rt_wt"].get_tag("XF")
        ib = {r.query_name: r.get_tag("XF") for r in b.fetch(itd["contig"])}
        assert ib.get("itd_j") == "junction" and ib.get("itd_wt1") != "junction" and ib.get("itd_wt2") != "junction", ib


def main():
    test_gene_names()
    test_single_record_junctions()
    test_merge()
    with tempfile.TemporaryDirectory() as d:
        fa, gtf, bam, arriba, sf = write_inputs(d)
        out = os.path.join(d, "out")
        subprocess.run([sys.executable, BPVIEW, "--sample-id", "T", "--bam", bam, "--fasta", fa, "--gtf", gtf,
                        "--arriba-tsv", arriba, "--starfusion-tsv", sf, "--out-prefix", out], check=True)
        doc = json.load(open(out + ".events.json"))
        assert doc["n_events"] == 1, doc
        ev = doc["events"][0]
        assert ev["callers"] == ["arriba", "starfusion"], ev["callers"]
        assert ev["junction_1based"] == BP1, ev["junction_1based"]
        nd = ev["fragments_nondup"]
        assert nd.get("junction") == 2, nd  # junc0 and jj
        assert nd.get("spanning") == 1, nd
        assert nd.get("wt3") == 1, nd
        assert ev["fragments_dup"].get("junction") == 1, ev["fragments_dup"]
        # junc0 starts at j-50; jj's mates start at j-60 and j-40 but count once (the leftmost)
        assert ev["junction_start_sites"] == 2, ev["junction_start_sites"]
        # per-caller counts from the cited reads (callers_bam), checkable in the BAM by XC
        ca, cs = ev["callers_bam"]["arriba"], ev["callers_bam"]["starfusion"]
        assert ca["fragments_nondup"] == {"junction": 1, "spanning": 1, "wt3": 1}, ca
        assert ca["fragments_dup"] == {"junction": 1} and ca["cited"] == 5 and ca["cited_missing"] == 1, ca
        assert (ca["split_5p"], ca["split_3p"]) == (1, 0), ca  # junc0's primary record is on G1 (5')
        assert cs["fragments_nondup"] == {"junction": 2} and cs["cited_missing"] == 0, cs
        assert (cs["split_5p"], cs["split_3p"]) == (2, 0), cs
        # contig sequence = chrA[..1000] + chrB[501..]
        contig = pysam.FastaFile(out + ".fa").fetch(ev["contig"])
        j = ev["junction_1based"]
        assert contig[j - 50:j] == CHR_A[950:1000] and contig[j:j + 50] == CHR_B[500:550]
        recs = list(pysam.AlignmentFile(out + ".bam"))
        junc = [r for r in recs if r.get_tag("XF") == "junction" and r.query_name == "junc0" and r.is_read1]
        assert len(junc) == 1 and junc[0].cigarstring == "100M" and junc[0].reference_start == j - 50, junc
        assert junc[0].query_sequence == contig[j - 50:j + 50]
        wt = [r for r in recs if r.get_tag("XF") == "wt3"]
        assert wt and wt[0].cigarstring == "50S50M" and wt[0].reference_start == j, [r.cigarstring for r in wt]
        assert all(r.has_tag("YC") and r.has_tag("XK") and r.has_tag("XG") for r in recs)
        # XG is one class per fragment, on both mates: junc0's mates are a junction read and an
        # ordinary mate, so XF/XK differ per mate but XG is "1junction" on both.
        by_name = {}
        for r in recs:
            by_name.setdefault(r.query_name, []).append(r)
        assert all(len({r.get_tag("XG") for r in rs}) == 1 for rs in by_name.values()), \
            {q: [r.get_tag("XG") for r in rs] for q, rs in by_name.items()}
        j0 = by_name["junc0"]
        assert {r.get_tag("XF") for r in j0} != {"junction"}, [r.get_tag("XF") for r in j0]
        assert {r.get_tag("XG") for r in j0} == {"1junction"}, [r.get_tag("XG") for r in j0]
        assert {r.get_tag("XG") for r in by_name["span"]} == {"2spanning"}
        xc = {q: {r.get_tag("XC") if r.has_tag("XC") else None for r in rs} for q, rs in by_name.items()}
        assert xc["junc0"] == {"arriba,starfusion"} and xc["junc1"] == {"arriba"} and xc["jj"] == {"starfusion"}, xc
        assert xc["span"] == {"arriba"} and xc["wt"] == {"arriba"}, xc
        # XJ: the split read's anchor side, on junction-class records only
        assert {r.get_tag("XJ") for r in by_name["junc0"] if r.get_tag("XF") == "junction"} == {"5p"}
        assert not any(r.has_tag("XJ") for r in recs if r.get_tag("XF") != "junction")
        sys.path.insert(0, os.path.join(HERE, ".."))
        import bpview
        assert bpview.starfusion_qname("&V1/1@V2:UMI_ACGTA_TTGAC") == "V2" and bpview.starfusion_qname("V3/2") == "V3"
        # fields the IGV view used to read from FusionInspector's abridged TSV
        assert ev["fusion_name"] == "G1--G2" and ev["left_gene"] == "G1" and ev["right_gene"] == "G2"
        assert (ev["left_local_breakpoint"], ev["right_local_breakpoint"]) == (j, j + 1)
        assert ev["left_breakpoint"] == f"chrA:{BP1}:+" and ev["right_breakpoint"] == f"chrB:{BP2}:+", ev
        assert ev["junction_read_count"] == 2 and ev["spanning_frag_count"] == 1
        assert ev["junction_read_count_with_dups"] == 3
        assert ev["splice_type"] == "ONLY_REF_SPLICE"
        # both callers' own coordinates are kept, so rows keyed by either caller can be matched
        assert ev["starfusion"]["left_breakpoint"] == f"chrA:{BP1 + 3}:+" and ev["arriba"]["breakpoint1"] == f"chrA:{BP1}"
        assert ev["gene_names1"] == ["G1"] and ev["gene_names2"] == ["G2"]
        assert ev["annots"] == ["G1:Oncogene", "Cosmic", "Mitelman", "INTERCHROMOSOMAL[chrA--chrB]"], ev["annots"]
        assert ev["left_gene_locus"] == "chrA:100-2900:+" and ev["right_gene_locus"] == "chrB:100-2900:+"
        # tabix-indexed tracks, queried per contig like gennext-api-edge does
        for ext in (".bed.gz", ".bed.gz.tbi", ".gtf.gz", ".gtf.gz.tbi"):
            assert os.path.exists(out + ext), ext
        gtf = list(pysam.TabixFile(out + ".gtf.gz").fetch(ev["contig"]))
        assert len(gtf) == 2 and all('orig_coord_info "chr' in l for l in gtf), gtf
        g1 = [l.split("\t") for l in gtf if 'gene_name "G1"' in l][0]
        assert (int(g1[3]), int(g1[4])) == (100, j) and 'orig_coord_info "chrA,100,2900,+"' in g1[8], g1
        assert len(list(pysam.TabixFile(out + ".bed.gz").fetch(ev["contig"]))) == 2
    print("test_bpview: OK")


if __name__ == "__main__":
    main()
