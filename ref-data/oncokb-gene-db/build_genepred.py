"""Build ANNOVAR's oncokbGene database (16-column genePred, UCSC chr names) from an
Ensembl GTF: OncoKB's grch38Isoform for each curated gene, and MANE Select for every
other gene, so the database can stand next to refGene/ensGene as a full gene set.

usage: build_genepred.py <gtf.gz> <allCuratedGenes.json> <out.txt> <report.json> [--oncokb-only]

Gene names are OncoKB's hugoSymbol (MANE genes: the GTF gene_name), with spaces and
parentheses replaced by "_" (ANNOVAR's coding_change.pl cannot parse a space:
"CDKN2A (p14)" becomes "CDKN2A_p14").
"""
import collections
import gzip
import json
import re
import sys

gtf, curated, out, report = sys.argv[1:5]
oncokb_only = "--oncokb-only" in sys.argv
ATTR = re.compile(r'(\S+) "([^"]*)"')


def annovar_name(sym):
    return re.sub(r"[ ()]+", "_", sym).strip("_")

cg = json.load(open(curated))
want = {g["grch38Isoform"]: g["hugoSymbol"] for g in cg if g["grch38Isoform"]}

tx = {}
exons = collections.defaultdict(list)
cds = collections.defaultdict(list)
mane = set()
for line in gzip.open(gtf, "rt"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    feat = f[2]
    if feat not in ("transcript", "exon", "CDS", "stop_codon", "start_codon"):
        continue
    a = {}
    tags = []
    for k, v in ATTR.findall(f[8]):
        if k == "tag":
            tags.append(v)
        else:
            a[k] = v
    t = a["transcript_id"]
    if feat == "transcript":
        tx[t] = {"chrom": f[0], "strand": f[6], "start": int(f[3]) - 1, "end": int(f[4]),
                 "version": a.get("transcript_version"), "gene_id": a["gene_id"],
                 "gene_name": a.get("gene_name") or a["gene_id"], "biotype": a.get("transcript_biotype"),
                 "tags": tags}
        if "MANE_Select" in tags:
            mane.add(t)
    elif feat == "exon":
        exons[t].append((int(f[3]) - 1, int(f[4])))
    elif feat in ("CDS", "stop_codon"):
        cds[t].append((int(f[3]) - 1, int(f[4]), feat, f[7]))

primary = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
chosen = {}
missing = []
for enst, sym in want.items():
    if enst in tx:
        chosen[enst] = sym
    else:
        missing.append((sym, enst))
curated_gene_ids = {tx[t]["gene_id"] for t in chosen}
curated_syms = set(want.values())
n_mane = 0
if not oncokb_only:
    for t in mane:
        d = tx[t]
        if d["gene_id"] in curated_gene_ids or d["gene_name"] in curated_syms:
            continue
        chosen[t] = d["gene_name"]
        n_mane += 1

rows = []
noncoding = []
for t, sym in chosen.items():
    d = tx[t]
    if d["chrom"] not in primary:
        continue
    ex = sorted(exons[t])
    cd = sorted(cds[t])
    if cd:
        cs, ce = min(c[0] for c in cd), max(c[1] for c in cd)
    else:
        cs = ce = d["end"]
        noncoding.append((sym, t))
    frames = []
    for s, e in ex:
        fr = [c for c in cd if c[2] == "CDS" and c[0] < e and c[1] > s]
        frames.append(str((3 - int(fr[0][3])) % 3) if fr else "-1")
    start_stat = "incmpl" if "cds_start_NF" in d["tags"] else ("cmpl" if cd else "none")
    end_stat = "incmpl" if "cds_end_NF" in d["tags"] else ("cmpl" if cd else "none")
    chrom = "chrM" if d["chrom"] == "MT" else "chr" + d["chrom"]
    name = f"{t}.{d['version']}" if d["version"] else t
    rows.append((chrom, d["start"], [
        "0", name, chrom, d["strand"], str(d["start"]), str(d["end"]), str(cs), str(ce), str(len(ex)),
        "".join(f"{s}," for s, _ in ex), "".join(f"{e}," for _, e in ex), "0", annovar_name(sym),
        start_stat, end_stat, "".join(f"{x}," for x in frames)]))

rows.sort(key=lambda r: (r[0], r[1], r[2][1]))
with open(out, "w") as fh:
    for _, _, r in rows:
        fh.write("\t".join(r) + "\n")
json.dump({"oncokb_isoforms_wanted": len(want), "oncokb_isoforms_in_ensembl111": len(want) - len(missing),
           "missing_from_ensembl111": missing, "mane_select_other_genes": n_mane, "rows": len(rows),
           "oncokb_noncoding": [x for x in noncoding if x[0] in curated_syms]},
          open(report, "w"), indent=1)
print(f"rows {len(rows)}; OncoKB isoforms {len(want) - len(missing)}/{len(want)}; MANE others {n_mane}; missing {missing}")
