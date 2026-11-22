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
                 "split_reads1\tsplit_reads2\tdiscordant_mates\tconfidence\tdirection1\tdirection2\ttags\n")
        fh.write(f"G1\tG2\t+/+\t+/+\tchrA:{BP1}\tchrB:{BP2}\t1\t0\t1\thigh\tdownstream\tupstream\tMitelman\n")
    sf = os.path.join(d, "sf.tsv")
    with open(sf, "w") as fh:
        fh.write("#FusionName\tJunctionReadCount\tSpanningFragCount\tSpliceType\tLeftGene\tLeftBreakpoint\t"
                 "RightGene\tRightBreakpoint\tFFPM\tLargeAnchorSupport\tannots\n")
        fh.write(f"G1--G2\t1\t1\tONLY_REF_SPLICE\tG1^X\tchrA:{BP1 + 3}:+\tG2^Y\tchrB:{BP2 - 2}:+\t1.0\tYES\t"
                 '["G1:Oncogene"];["Cosmic","Mitelman"];INTERCHROMOSOMAL[chrA--chrB]\n')
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


def main():
    test_gene_names()
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
        assert all(r.has_tag("YC") and r.has_tag("XK") for r in recs)
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
