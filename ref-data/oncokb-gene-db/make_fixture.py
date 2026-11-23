"""Make the small oncokbGene test fixture (tests/fixtures/oncokb-gene-db) from a full build.

usage: make_fixture.py <build-dir> <mane-genepred> <mane-mrna.fa> <hg38.fa> <out-dir>

  build-dir      build.sh's work-dir/build (hg38_oncokbGene.txt, ...Mrna.fa, .isoforms.tsv)
  mane-*         a MANE-only genePred and its mRNA (build_genepred.py with an empty curated
                 list, then retrieve_seq_from_fasta.pl), standing in for refGene/ensGene

The variants are synthetic: well-known hotspots and changes made at chosen codons of the
fixture transcripts (no sample data). Each has an id in the VCF's ID column.
"""
import os
import sys

build, mane_gp, mane_fa, fasta, out = sys.argv[1:6]
GENES = ["BRAF", "EGFR", "KRAS", "TP53", "ROS1", "NF1", "SRC", "PPP2R1A", "DIS3L2", "CDKN2A",
         "CDKN2A_p14", "ACTB", "RET"]
os.makedirs(out, exist_ok=True)


def read_fa(path, keep):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in keep:
                seqs[name] = [line]
            else:
                name = None
        elif name:
            seqs[name].append(line)
    return seqs


def subset(gp_path, fa_path, gp_out, fa_out):
    rows = [l for l in open(gp_path) if l.split("\t")[12] in GENES]
    open(gp_out, "w").writelines(rows)
    names = [l.split("\t")[1] for l in rows]
    fa = read_fa(fa_path, set(names))
    with open(fa_out, "w") as fh:
        for n in names:
            fh.writelines(fa[n])
    return {l.split("\t")[12]: l.rstrip("\n").split("\t") for l in rows}


okb = subset(f"{build}/hg38_oncokbGene.txt", f"{build}/hg38_oncokbGeneMrna.fa",
             f"{out}/hg38_oncokbGene.txt", f"{out}/hg38_oncokbGeneMrna.fa")
subset(mane_gp, mane_fa, f"{out}/hg38_refGene.txt", f"{out}/hg38_refGeneMrna.fa")
subset(mane_gp, mane_fa, f"{out}/hg38_ensGene.txt", f"{out}/hg38_ensGeneMrna.fa")
with open(f"{build}/hg38_oncokbGene.isoforms.tsv") as fh, open(f"{out}/hg38_oncokbGene.isoforms.tsv", "w") as w:
    w.write(fh.readline())
    for line in fh:
        if line.split("\t")[1] in GENES:
            w.write(line)
open(f"{out}/hg38_oncokbGene.version", "w").write("name=oncokbGene\nversion=fixture\n")

# reference bases through the .fai
fai = {l.split("\t")[0]: [int(x) for x in l.split("\t")[1:5]] for l in open(fasta + ".fai")}
fh = open(fasta, "rb")


def ref(chrom, pos, n=1):  # 1-based
    length, offset, lb, lbytes = fai[chrom]
    out = []
    for p in range(pos, pos + n):
        i = p - 1
        fh.seek(offset + (i // lb) * lbytes + i % lb)
        out.append(fh.read(1).decode().upper())
    return "".join(out)


COMP = str.maketrans("ACGT", "TGCA")


def cds_positions(f):
    """Genomic 1-based positions of the CDS, in transcript order."""
    strand, cs, ce = f[3], int(f[6]), int(f[7])
    starts = [int(x) for x in f[9].rstrip(",").split(",")]
    ends = [int(x) for x in f[10].rstrip(",").split(",")]
    pos = [p + 1 for s, e in zip(starts, ends) for p in range(max(s, cs), min(e, ce))]
    return pos if strand == "+" else pos[::-1]


def codon(gene, n):
    f = okb[gene]
    pos = cds_positions(f)[3 * (n - 1):3 * n]
    b = "".join(ref(f[2], p) for p in pos)
    return pos, (b if f[3] == "+" else b.translate(COMP)), f


def snv_at(gene, n, i, new_tx_base):
    """Change base i (0..2) of codon n to new_tx_base (transcript strand)."""
    pos, _, f = codon(gene, n)
    p = pos[i]
    alt = new_tx_base if f[3] == "+" else new_tx_base.translate(COMP)
    return f[2], p, ref(f[2], p), alt


V = []
# hotspots
V += [("braf_v600e", "chr7", 140753336, "A", "T"), ("egfr_l858r", "chr7", 55191822, "T", "G"),
      ("kras_g12d", "chr12", 25245350, "C", "T"), ("tp53_r175h", "chr17", 7675088, "C", "T"),
      ("egfr_e746_a750del", "chr7", 55174771, "AGGAATTAAGAGAAGC", "A")]
# OncoKB isoform differs from MANE: missense at codon 1000 / 300
for gene, n, tag in (("ROS1", 1000, "ros1_missense"), ("NF1", 1000, "nf1_missense"), ("SRC", 300, "src_missense")):
    pos, b, f = codon(gene, n)
    nb = {"A": "C", "C": "A", "G": "T", "T": "G"}[b[0]]
    V.append((tag,) + snv_at(gene, n, 0, nb))
# stop-loss (ext*N) and synonymous stop (*N=)
pp = codon("PPP2R1A", len(cds_positions(okb["PPP2R1A"])) // 3)
V.append(("ppp2r1a_stoploss",) + snv_at("PPP2R1A", len(cds_positions(okb["PPP2R1A"])) // 3, 0, "C"))
nd = len(cds_positions(okb["DIS3L2"])) // 3
_, b, _ = codon("DIS3L2", nd)
i, nb = next((i, x) for i in (2, 1) for x in "ACGT"
             if x != b[i] and b[:i] + x + b[i + 1:] in ("TAA", "TAG", "TGA"))
V.append(("dis3l2_syn_stop",) + snv_at("DIS3L2", nd, i, nb))
# nonsense in RET: first CGA codon after 700
for n in range(700, 1100):
    if codon("RET", n)[1] == "CGA":
        V.append(("ret_nonsense",) + snv_at("RET", n, 0, "T"))
        break
# dup: copy EGFR codons 770-772 after codon 772 (+ strand)
pos, _, f = codon("EGFR", 772)
dup = "".join(ref("chr7", p) for p in cds_positions(okb["EGFR"])[3 * 769:3 * 772])
V.append(("egfr_dup", "chr7", pos[2], ref("chr7", pos[2]), ref("chr7", pos[2]) + dup))
# frameshift: delete one base in KRAS codon 20
pos, _, f = codon("KRAS", 20)
V.append(("kras_fs", "chr12", pos[1] - 1, ref("chr12", pos[1] - 1, 2), ref("chr12", pos[1] - 1)))
# CDKN2A exon 2, read in both frames (p16 and p14ARF)
V.append(("cdkn2a_ex2",) + snv_at("CDKN2A", 80, 0, "T" if codon("CDKN2A", 80)[1][0] != "T" else "A"))
# exon edge: the last coding base of a BRAF internal exon, and +5 into the intron after it
f = okb["BRAF"]
cp = cds_positions(f)
starts = [int(x) for x in f[9].rstrip(",").split(",")]
edge = next(i for i in range(3, len(cp) - 1) if abs(cp[i] - cp[i + 1]) > 1 and i % 3 == 2)
_, b, _ = codon("BRAF", edge // 3 + 1)
syn = {"A": "G", "G": "A", "C": "T", "T": "C"}[b[2]]
V.append(("braf_exon_edge_syn",) + snv_at("BRAF", edge // 3 + 1, 2, syn))
p = cp[edge] - 5  # - strand: intron after the exon is at lower coordinates
V.append(("braf_intron_plus5", "chr7", p, ref("chr7", p), "A" if ref("chr7", p) != "A" else "G"))
# not curated by OncoKB
pos, b, f = codon("ACTB", 100)
V.append(("actb_missense",) + snv_at("ACTB", 100, 0, "A" if b[0] != "A" else "C"))

V.sort(key=lambda v: (v[1], v[2]))
with open(f"{out}/fixture.vcf", "w") as w:
    w.write("##fileformat=VCFv4.2\n")
    for c in sorted({v[1] for v in V}):
        w.write(f"##contig=<ID={c},length={fai[c][0]}>\n")
    w.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH\n")
    for tag, c, p, r, a in V:
        w.write(f"{c}\t{p}\t{tag}\t{r}\t{a}\t.\tPASS\t.\tGT\t0/1\n")
print(f"{len(V)} variants")
