"""Build the oncokbGene sidecar: one row per OncoKB isoform, for annovar-to-maf.

usage: build_isoforms.py <hg38_oncokbGene.txt> <hg38_oncokbGeneMrna.fa> <pep.all.fa.gz>
                         <allCuratedGenes.json> <out.tsv>

Columns (tab-separated, with a header):
  hugo_symbol     OncoKB's hugoSymbol, as OncoKB writes it ("CDKN2A (p14)")
  annovar_gene    the name ANNOVAR reports in Gene.oncokbGene ("CDKN2A_p14")
  transcript      grch38Isoform with its Ensembl version (ENST00000288602.11)
  refseq          grch38RefSeq, as OncoKB lists it (may be empty)
  chrom strand cds_start cds_end exon_starts exon_ends
                  from the genePred (0-based half-open, comma-separated exons)
  protein         the Ensembl protein of the transcript (no trailing "*")
  utr3            the mRNA after the stop codon, for stop-loss extension lengths

Checks (exit 1 if one fails):
  - every OncoKB isoform is in the genePred;
  - every coding isoform's CDS, translated, equals its Ensembl protein, except at
    residue 1 (non-ATG starts) or where the CDS is marked incomplete.
"""
import csv
import gzip
import json
import sys

gp_path, mrna_path, pep_path, curated_path, out_path = sys.argv[1:6]

CODON = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODON[a + b + c] = aas[i]
            i += 1


def translate(s):
    return "".join(CODON.get(s[j:j + 3], "X") for j in range(0, len(s) - 2, 3))


def acc(t):
    return t.split(".")[0]


gp = {}
for line in open(gp_path):
    f = line.rstrip("\n").split("\t")
    gp[acc(f[1])] = f

mrna = {}
name = None
for line in open(mrna_path):
    if line.startswith(">"):
        name = acc(line[1:].split()[0])
        mrna[name] = []
    elif name:
        mrna[name].append(line.strip().upper())
mrna = {k: "".join(v) for k, v in mrna.items()}

pep = {}
name = None
for line in gzip.open(pep_path, "rt"):
    if line.startswith(">"):
        name = next((acc(t.split(":", 1)[1]) for t in line.split() if t.startswith("transcript:")), None)
        pep[name] = []
    elif name:
        pep[name].append(line.strip())
pep = {k: "".join(v) for k, v in pep.items()}

curated = json.load(open(curated_path))
failures, notes, rows = [], [], []
seen = set()
for g in sorted(curated, key=lambda g: g["hugoSymbol"]):
    iso = g.get("grch38Isoform")
    if not iso:
        continue
    f = gp.get(acc(iso))
    if f is None:
        failures.append(f"{g['hugoSymbol']} {iso}: not in the genePred")
        continue
    chrom, strand = f[2], f[3]
    tx_start, tx_end, cds_start, cds_end = (int(x) for x in f[4:8])
    starts = [int(x) for x in f[9].rstrip(",").split(",")]
    ends = [int(x) for x in f[10].rstrip(",").split(",")]
    protein, utr3 = pep.get(acc(iso), "").rstrip("*"), ""
    if cds_start < cds_end:
        seq = mrna.get(acc(iso), "")
        exons = list(zip(starts, ends))
        # mRNA offsets of the CDS, in transcript orientation
        before = sum(max(0, min(e, cds_start) - s) for s, e in exons)
        cds_len = sum(max(0, min(e, cds_end) - max(s, cds_start)) for s, e in exons)
        after = sum(max(0, e - max(s, cds_end)) for s, e in exons)
        off = before if strand == "+" else after
        cds = seq[off:off + cds_len]
        utr3 = seq[off + cds_len:]
        incomplete = "incmpl" in (f[13], f[14])
        tr = translate(cds).rstrip("*")
        if tr == protein:
            pass
        elif tr[1:] == protein[1:] or incomplete:
            notes.append(f"{g['hugoSymbol']} {f[1]}: translation differs at residue 1 or the CDS is incomplete")
        else:
            failures.append(f"{g['hugoSymbol']} {f[1]}: translation differs from the Ensembl protein")
    key = (g["hugoSymbol"], f[1])
    if key in seen:
        continue
    seen.add(key)
    rows.append([g["hugoSymbol"], f[12], f[1], g.get("grch38RefSeq") or "", chrom, strand,
                 str(cds_start), str(cds_end), f[9], f[10], protein, utr3])

with open(out_path, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["hugo_symbol", "annovar_gene", "transcript", "refseq", "chrom", "strand", "cds_start",
                "cds_end", "exon_starts", "exon_ends", "protein", "utr3"])
    w.writerows(rows)

for n in notes:
    print("note:", n)
for x in failures:
    print("FAIL:", x)
print(f"{len(rows)} isoform rows; {len(notes)} notes; {len(failures)} failures")
sys.exit(1 if failures else 0)
