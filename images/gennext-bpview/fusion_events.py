#!/usr/bin/env python3
"""fusion-events: one row per fusion breakpoint event, for the per-sample fusions DuckDB.

  fusion-events rescue --arriba-tsv A --arriba-discarded-tsv D --starfusion-tsv F --oncokb-json J
                       [--fusion-list L] --out-arriba-tsv A2 [--out-fusion-list L2]
  fusion-events build  --sample-id S --arriba-tsv A --starfusion-tsv F --oncokb-json J
                       [--bpview-events E] [--pon P] [--gtf G] [--protein-domains P]
                       --out-tsv OUT [--out-meta-tsv META]
  fusion-events pon    --events-tsv E1 E2 ... --out-tsv PON

The output follows the event data contract, docs/fusions/event-contract.md (contract_version 1).

`rescue` (prototype rule, Q1 on #30) brings calls back from Arriba's discarded list: at least 1
unique read, removed only by support-level filters (in_vitro, relative_support, min_support), and
a known pair (Arriba's Mitelman tag, a database or OncoKB hit) or a STAR-Fusion call of the same
junction (+-10 bp). It writes Arriba's kept calls plus the rescued ones, with an `arriba_status`
column, for bpview and `build`, and appends the rescued pairs to bpview's fusion list.

`build` merges Arriba and STAR-Fusion calls by breakpoint with bpview's merge (same loaders,
same +-10 bp tolerance, same gene aliases), so every event lines up with its bpview contig.
Gene-pair annotation and PoN pair counts use only the gene each breakpoint is in (STAR-Fusion's
gene, or an Arriba gene without a distance); an intergenic neighbour counts only within
--max-intergenic-distance (default 0).
Gene-pair annotation (OncoKB, fusion-report databases, the existing pair tier) comes from
oncokb-fusion-annotation's fusions_with_oncokb.json: no second OncoKB query. Support comes
from bpview's projected, non-duplicate fragments when bpview drew the event, otherwise from
Arriba's count after its dedup, or STAR-Fusion's (support_source says which).
Exon numbers come from the gene's MANE Select transcript in --gtf, and kinase retention from the
Pfam protein kinase domain (PF00069) in Arriba's protein-domain annotation (--protein-domains).

`pon` turns event tables of many samples into a recurrence panel: one row per sample and
junction. `build --pon` then counts, for each event, how many *other* samples carry the same
gene pair (either orientation) and the same junction (+-tolerance).

Tiers (see tier_event; every threshold is a flag):
  1  actionable: OncoKB LEVEL_1/LEVEL_2/LEVEL_R1 (or Dx1) for the pair in this orientation,
     with fusion evidence: not out-of-frame/stop-codon, at least --tier1-min-unique-molecules,
     and two callers or Arriba confidence high/medium
  2  any OncoKB hit (actionable, Oncogenic / Likely Oncogenic, 3A/3B/4/R2/Dx2/Dx3) with frame
     and caller evidence but not tier 1 (e.g. too few molecules),
     or a curated known fusion (Mitelman, COSMIC) whatever its evidence; tier_flags still records
     the weak evidence
  3  everything else. OncoKB levels are often partner-agnostic (any NTRK fusion), so an OncoKB
     hit without fusion evidence and without a curated pair stays here; tier_reason says so.
An unoriented Arriba call whose annotation matched only the reversed pair gets no tier 1 and no
curated tier 2 without evidence (flag orientation_unknown).
low_support (prototype rule, Q1): fewer than --min-unique-molecules (default 2) unique molecules.
Shown, never dropped; capped at tier 2; ordered after supported events of the same tier.
A gene pair recurrent in the PoN (>= --pon-min-fraction of the other samples, and at least
--pon-min-samples, default 1) is capped at tier 2 (if it would be 1) or 3. tier_reason lists
every rule that applied.
"""
import argparse
import csv
import gzip
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bpview  # noqa: E402  (same image, same directory)

__version__ = "0.6.0"
CONTRACT_VERSION = "1"

ACTIONABLE = {"LEVEL_1", "LEVEL_2", "LEVEL_R1"}
ACTIONABLE_DX = {"LEVEL_Dx1"}
ONCOGENIC_LEVELS = {"LEVEL_3A", "LEVEL_3B", "LEVEL_4", "LEVEL_R2"}
ONCOGENIC_DX = {"LEVEL_Dx2", "LEVEL_Dx3"}
ONCOGENIC = {"Oncogenic", "Likely Oncogenic"}
# The curated databases fusion-report 4.0.1 reports (FusionGDB2 is not curated).
CURATED_DB = {"Mitelman", "COSMIC"}

# Column order and DuckDB types of the output; duckdb-fusion-events reads it with these types.
# Every field a user may filter on is its own typed column. List columns (VARCHAR[]) are written
# '|'-separated in the TSV and split on load; an empty cell is an empty list.
TIER_FLAGS = ("out_of_frame", "weak_caller_evidence", "low_unique_molecules", "low_support", "rescued",
              "rescue_capped", "no_molecules", "readthrough_capped", "orientation_unknown", "pon_recurrent",
              "pon_capped")
COLUMNS = [
    # identity and order
    ("sample_id", "VARCHAR"), ("event_id", "VARCHAR"), ("event_rank", "INTEGER"), ("fusion", "VARCHAR"),
    ("gene5", "VARCHAR"), ("gene3", "VARCHAR"), ("gene5_aliases", "VARCHAR[]"), ("gene3_aliases", "VARCHAR[]"),
    ("breakpoint5", "VARCHAR"), ("breakpoint3", "VARCHAR"), ("orientation_known", "BOOLEAN"),
    # structure
    ("transcript5", "VARCHAR"), ("transcript3", "VARCHAR"), ("exon5", "INTEGER"), ("exon3", "INTEGER"),
    ("site5", "VARCHAR"), ("site3", "VARCHAR"), ("reading_frame", "VARCHAR"), ("splice_type", "VARCHAR"),
    ("kinase_gene", "VARCHAR"), ("kinase_retained", "BOOLEAN"), ("kinase_domain_fraction", "DOUBLE"),
    # callers
    ("callers", "VARCHAR[]"), ("n_callers", "INTEGER"), ("caller_arriba", "BOOLEAN"), ("caller_starfusion", "BOOLEAN"),
    ("arriba_status", "VARCHAR"), ("arriba_filters", "VARCHAR[]"), ("arriba_confidence", "VARCHAR"),
    ("arriba_type", "VARCHAR"), ("arriba_split_reads", "INTEGER"), ("arriba_discordant_mates", "INTEGER"),
    ("starfusion_junction_reads", "INTEGER"), ("starfusion_spanning_frags", "INTEGER"),
    ("starfusion_ffpm", "DOUBLE"), ("breakpoint_offset", "INTEGER"),
    # support
    ("unique_molecules", "INTEGER"), ("molecules_junction", "INTEGER"), ("molecules_spanning", "INTEGER"),
    ("start_sites", "INTEGER"), ("pcr_copies", "INTEGER"), ("support_source", "VARCHAR"),
    ("low_support", "BOOLEAN"),
    # annotation
    ("known_db", "VARCHAR[]"), ("curated_known", "BOOLEAN"),
    ("oncokb_oncogenic", "VARCHAR"), ("oncokb_mutation_effect", "VARCHAR"),
    ("oncokb_highest_level", "VARCHAR"), ("oncokb_highest_dx_level", "VARCHAR"),
    ("pair_clinical_tier", "VARCHAR"), ("pair_fii", "DOUBLE"),
    # recurrence
    ("pon_gene_pair_count", "INTEGER"), ("pon_junction_count", "INTEGER"), ("pon_size", "INTEGER"),
    ("pon_recurrent", "BOOLEAN"),
    # tier_basis: the rule that set the tier (oncokb_actionable | oncokb_oncogenic | curated_known | none);
    # tier_flags: what held it back (TIER_FLAGS); tier_reason: the same in words.
    ("tier", "INTEGER"), ("tier_basis", "VARCHAR"), ("tier_flags", "VARCHAR[]"), ("tier_reason", "VARCHAR"),
    # links and provenance
    ("bpview_contig", "VARCHAR"), ("builder_version", "VARCHAR"), ("prototype_rules", "VARCHAR[]"),
    # --- 0.5.0 additions (contract v1, additive): caller and annotator details for the detail panel.
    # Caller and annotator values, never the event's read count (contract §4 "Support").
    # Arriba, from the arm's own Arriba row; sides 1/2 are Arriba's gene1/gene2 as written (the
    # reverse of 5'/3' when an unoriented call was flipped).
    ("arriba_gene1", "VARCHAR"), ("arriba_gene2", "VARCHAR"),
    ("arriba_breakpoint1", "VARCHAR"), ("arriba_breakpoint2", "VARCHAR"),
    ("arriba_strand1", "VARCHAR"), ("arriba_strand2", "VARCHAR"),
    ("arriba_split_reads1", "INTEGER"), ("arriba_split_reads2", "INTEGER"),
    ("arriba_coverage1", "INTEGER"), ("arriba_coverage2", "INTEGER"),
    ("arriba_site1", "VARCHAR"), ("arriba_site2", "VARCHAR"),
    ("arriba_filters_raw", "VARCHAR"), ("arriba_duplicates", "INTEGER"), ("arriba_tags", "VARCHAR[]"),
    ("arriba_retained_domains1", "VARCHAR[]"), ("arriba_retained_domains2", "VARCHAR[]"),
    ("arriba_closest_genomic_breakpoint1", "VARCHAR"), ("arriba_closest_genomic_breakpoint2", "VARCHAR"),
    ("arriba_gene_id1", "VARCHAR"), ("arriba_gene_id2", "VARCHAR"),
    ("arriba_transcript_id1", "VARCHAR"), ("arriba_transcript_id2", "VARCHAR"),
    ("arriba_fusion_transcript", "VARCHAR"), ("arriba_peptide_sequence", "VARCHAR"),
    # STAR-Fusion, from the event's STAR-Fusion row (left = 5', right = 3')
    ("starfusion_left_gene", "VARCHAR"), ("starfusion_right_gene", "VARCHAR"),
    ("starfusion_left_breakpoint", "VARCHAR"), ("starfusion_right_breakpoint", "VARCHAR"),
    ("starfusion_est_junction", "DOUBLE"), ("starfusion_est_spanning", "DOUBLE"),
    ("starfusion_splice_type", "VARCHAR"), ("starfusion_large_anchor_support", "BOOLEAN"),
    ("starfusion_left_break_dinuc", "VARCHAR"), ("starfusion_right_break_dinuc", "VARCHAR"),
    ("starfusion_left_break_entropy", "DOUBLE"), ("starfusion_right_break_entropy", "DOUBLE"),
    ("starfusion_annots", "VARCHAR[]"),
    # OncoKB (gene pair, fusions_with_oncokb.json)
    ("oncokb_highest_sensitive_level", "VARCHAR"), ("oncokb_highest_resistance_level", "VARCHAR"),
    ("oncokb_highest_px_level", "VARCHAR"), ("oncokb_treatments", "VARCHAR[]"),
    ("oncokb_mutation_effect_citations", "VARCHAR[]"), ("oncokb_treatment_citations", "VARCHAR[]"),
    ("oncokb_dx_citations", "VARCHAR[]"), ("oncokb_px_citations", "VARCHAR[]"),
    ("oncokb_gene_in_oncokb", "BOOLEAN"), ("oncokb_variant_in_oncokb", "BOOLEAN"),
    # OncoKB's variantSummary and matched alteration names, when the annotator stored them (gennext-workflows#75)
    ("oncokb_variant_summary", "VARCHAR"), ("oncokb_matched_alterations", "VARCHAR[]"),
    # fusion-report
    ("pair_fii_explained", "VARCHAR"),
    # Both callers side by side, on the event's 5'/3'. arriba_split_reads_5p/_3p: Arriba's own
    # split_reads1/2 mapped to 5'/3'. *_bam_*: counted in the bpview BAM from the fragments that caller
    # cites (XC tag): split reads by the partner of their primary alignment, molecules and copies of
    # XG 1junction + 2spanning. Null when the caller didn't call the event or bpview didn't draw it.
    ("arriba_split_reads_5p", "INTEGER"), ("arriba_split_reads_3p", "INTEGER"),
    ("arriba_bam_split_reads_5p", "INTEGER"), ("arriba_bam_split_reads_3p", "INTEGER"),
    ("arriba_bam_molecules", "INTEGER"), ("arriba_bam_pcr_copies", "INTEGER"), ("arriba_cited_missing", "INTEGER"),
    ("starfusion_bam_split_reads_5p", "INTEGER"), ("starfusion_bam_split_reads_3p", "INTEGER"),
    ("starfusion_bam_molecules", "INTEGER"), ("starfusion_bam_pcr_copies", "INTEGER"),
    ("starfusion_cited_missing", "INTEGER"),
    # --- 0.6.0 additions (contract v1, additive)
    # reading_frame_source: "arriba" (Arriba's call) or "mane" (derived from the MANE Select CDS phase of
    # both partners when no caller gives a frame); null when reading_frame is null.
    ("reading_frame_source", "VARCHAR"),
    # OncoKB's answer for this specific gene pair (not only one gene's fusions, not as a loss-of-function event),
    # decided from the OncoKB API response stored at annotation time (oncokb_exact_pair).
    ("oncokb_pair_specific", "BOOLEAN"),
    # Both breakpoints at exon boundaries of the MANE Select transcripts.
    ("canonical_junction", "BOOLEAN"),
    # Same chromosome and strand, < READTHROUGH_MAX_DISTANCE apart, no other gene between the breakpoints.
    ("adjacent_readthrough", "BOOLEAN"),
]
READTHROUGH_MAX_DISTANCE = 50_000  # bp: an adjacent-gene read-through is capped at tier 3
LOF_EFFECTS = {"Loss-of-function", "Likely Loss-of-function"}
EXACT_PAIR_RULE = ("OncoKB API answer: a matched alteration is 'GENE5-GENE3 Fusion' (either order), when OncoKB "
                   "names its matches; else VARIANT_IN_ONCOKB=True, two different genes, MUTATION_EFFECT not "
                   "(Likely) Loss-of-function, ONCOGENIC not Unknown")
ONCOKB_LEVELS = ("LEVEL_1", "LEVEL_2", "LEVEL_3A", "LEVEL_3B", "LEVEL_4", "LEVEL_R1", "LEVEL_R2",
                 "LEVEL_Dx1", "LEVEL_Dx2", "LEVEL_Dx3", "LEVEL_Px1", "LEVEL_Px2", "LEVEL_Px3")
SUPPORT_FILTERS = {"in_vitro", "relative_support", "min_support"}
# Pfam "Protein kinase domain" and "Protein tyrosine and serine/threonine kinase". Arriba 2.5.1's file
# lists only PF00069 (ALK, FGFR3, RET, NTRK1-3, ROS1, BRAF, EGFR, MET all carry it); PF07714 is
# accepted too, should a later file list it separately.
KINASE_PFAM = ("PF00069", "PF07714")
PON_COLUMNS = ["sample_id", "gene5", "gene3", "chr5", "pos5", "keep5", "chr3", "pos3", "keep3"]


# ----------------------------------------------------------------------------- inputs

def load_pairs(path):
    """fusions_with_oncokb.json -> {'G1--G2': entry}. Missing or empty file -> {} with a warning."""
    if not path or not os.path.isfile(path) or not os.path.getsize(path):
        print(f"WARNING: gene-pair annotation {path!r} missing or empty: no OncoKB or database "
              "annotation, every event is tier 3", file=sys.stderr)
        return {}
    with open(path) as fh:
        data = json.load(fh)
    return {e["Fusion"]: e for e in data if e.get("Fusion")}


def event_key(chr1, pos1, keep1, chr2, pos2, keep2):
    return (chr1, int(pos1), keep1, chr2, int(pos2), keep2)


def load_bpview(path):
    """bpview events.json -> {junction key: event summary}."""
    if not path or not os.path.isfile(path) or not os.path.getsize(path):
        return {}
    with open(path) as fh:
        doc = json.load(fh)
    out = {}
    for e in doc.get("events", []):
        c1, p1 = e["breakpoint1"].rsplit(":", 1)
        c2, p2 = e["breakpoint2"].rsplit(":", 1)
        out[event_key(c1, p1, e["keep1"], c2, p2, e["keep2"])] = e
    return out


def load_pon(path):
    if not path or not os.path.isfile(path) or not os.path.getsize(path):
        return None
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["pos5"], r["pos3"] = int(r["pos5"]), int(r["pos3"])
            rows.append(r)
    pairs, by_contig = {}, {}
    for r in rows:
        pairs.setdefault(frozenset((r["gene5"], r["gene3"])), set()).add(r["sample_id"])
        by_contig.setdefault((r["chr5"], r["keep5"], r["chr3"], r["keep3"]), []).append(r)
    return dict(pairs=pairs, by_contig=by_contig, samples={r["sample_id"] for r in rows})


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_transcripts(path):
    """GTF -> {gene_name: transcript}, the MANE Select transcript (else Ensembl canonical) of each gene:
    {id, chrom, strand, exons: [(start, end, number)]}. No file -> {}."""
    if not path or path == "." or not os.path.isfile(path):
        if path and path != ".":
            print(f"WARNING: GTF {path!r} not found: exon5/exon3 will be empty", file=sys.stderr)
        return {}
    tx = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            a = f[8]
            g, t, n = (re.search(k + r' "?([^";]+)"?;', a) for k in ("gene_name", "transcript_id", "exon_number"))
            if not (g and t and n):
                continue
            rank = 2 if 'tag "MANE_Select"' in a else 1 if 'tag "Ensembl_canonical"' in a else 0
            if not rank:
                continue
            e = tx.setdefault((g.group(1), t.group(1)),
                              dict(id=t.group(1), chrom=f[0], strand=f[6], rank=rank, exons=[], cds=[]))
            if f[2] == "exon":
                e["exons"].append((int(f[3]), int(f[4]), int(n.group(1))))
            else:
                e["cds"].append((int(f[3]), int(f[4])))
    best = {}
    for (gene, _), t in tx.items():
        if gene not in best or t["rank"] > best[gene]["rank"]:
            best[gene] = t
    return best


def oncokb_exact_pair(entry, fusion):
    """Is OncoKB's answer for this fusion about this specific gene pair? Decided from the OncoKB API response
    stored by the annotator at annotation time (a ruling over the API answer, never a copy of
    OncoKB's data). -> True / False, or None when the pair has no OncoKB annotation.

    1. If the response records the alterations OncoKB matched (ONCOKB_MATCHED_ALTERATIONS, written by the
       annotator from gennext-workflows#75: ';'-joined, e.g. "EML4-ALK Fusion;Fusions"; MATCHED_ALTERATIONS
       is accepted too) and the field is non-empty, the pair is specific when one entry is exactly
       "GENE5-GENE3 Fusion" or "GENE3-GENE5 Fusion".
    2. Otherwise (field absent: outputs before #75; or empty: OncoKB named no alteration): VARIANT_IN_ONCOKB is True, the partners are two different genes, MUTATION_EFFECT is not
       Loss-of-function or Likely Loss-of-function, and ONCOGENIC is not Unknown. OncoKB also sets
       VARIANT_IN_ONCOKB for tumour-suppressor and intragenic fusions it matches as truncating; the effect test
       drops those (rollout check of 9,736 stored rows, 25 Nov)."""
    onc = (entry or {}).get("OncoKB")
    if not onc:
        return None
    g5, _, g3 = (fusion or "").partition("--")
    if not g5 or not g3:
        return False
    matched = onc.get("ONCOKB_MATCHED_ALTERATIONS", onc.get("MATCHED_ALTERATIONS"))
    if matched:
        alts = matched if isinstance(matched, list) else matched.split(";")
        alts = {a.strip() for a in alts if a and a.strip()}
        if alts:
            return f"{g5}-{g3} Fusion" in alts or f"{g3}-{g5} Fusion" in alts
    return (onc.get("VARIANT_IN_ONCOKB") == "True" and g5 != g3
            and onc.get("MUTATION_EFFECT") not in LOF_EFFECTS
            and (onc.get("ONCOGENIC") or "Unknown") != "Unknown")


def gene_spans(transcripts):
    """{chrom: [(start, end, gene)]} from the selected transcripts (MANE Select)."""
    out = {}
    for g, t in transcripts.items():
        if t["exons"]:
            out.setdefault(t["chrom"], []).append((min(a for a, _, _ in t["exons"]), max(b for _, b, _ in t["exons"]), g))
    return out


def adjacent_readthrough(ev, spans):
    """Partner distance (bp) when the event is an adjacent-gene read-through: same chromosome and strand,
    closer than READTHROUGH_MAX_DISTANCE, and no other gene overlapping the interval between the breakpoints;
    else None. Genes are those with a MANE Select transcript in --gtf."""
    if ev["chr1"] != ev["chr2"] or not ev["oriented"] or strand5(ev) != strand3(ev):
        return None
    d = abs(ev["pos1"] - ev["pos2"])
    if d >= READTHROUGH_MAX_DISTANCE:
        return None
    lo, hi = sorted((ev["pos1"], ev["pos2"]))
    partners = set(ev["names1"]) | set(ev["names2"])
    if any(a < hi and b > lo and g not in partners for a, b, g in spans.get(ev["chr1"], [])):
        return None
    return d


def load_kinase_domains(path):
    """Arriba's protein-domain GFF3 -> {gene_name: (chrom, [(start, end)] merged)} for the Pfam protein
    kinase domain. No file -> {}."""
    if not path or path == "." or not os.path.isfile(path):
        if path and path != ".":
            print(f"WARNING: protein domains {path!r} not found: kinase_* will be empty", file=sys.stderr)
        return {}
    segs = {}
    with _open(path) as fh:
        for line in fh:
            if not any(k in line for k in KINASE_PFAM):
                continue
            f = line.rstrip("\n").split("\t")
            m = re.search(r"protein_domain_id=([^;\s]+)", f[8]) if len(f) >= 9 else None
            if not m or m.group(1) not in KINASE_PFAM:
                continue
            g = re.search(r"gene_name=([^;]+)", f[8])
            if not g:
                continue
            chrom = f[0] if f[0].startswith("chr") else "chr" + f[0]
            segs.setdefault(g.group(1), (chrom, []))[1].append((int(f[3]), int(f[4])))
    out = {}
    for gene, (chrom, iv) in segs.items():
        merged = []
        for a, b in sorted(iv):
            if merged and a <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        out[gene] = (chrom, merged)
    return out


# ----------------------------------------------------------------------------- per event

def strand5(ev):
    return "+" if ev["keep1"] == "le" else "-"


def strand3(ev):
    return "+" if ev["keep2"] == "ge" else "-"


def arriba_distances(ev):
    """Distance from the breakpoint to every gene Arriba names (0 = the breakpoint is in it).
    Arriba writes 'GENE' inside a gene and 'GENE1(dist),GENE2(dist)' between genes."""
    out = {}
    a = ev.get("arriba") or {}
    for raw in (a.get("gene1") or "", a.get("gene2") or ""):
        for part in raw.split(","):
            m = re.match(r"\s*([^()]+?)\s*(?:\((\d+)\))?\s*$", part)
            if m and m.group(1):
                d = int(m.group(2)) if m.group(2) else 0
                out[m.group(1)] = min(out.get(m.group(1), d), d)
    return out


def genic_names(ev, side, max_distance):
    """The names of one side that may carry gene-pair annotation and PoN recurrence: the gene the
    breakpoint is in (STAR-Fusion's gene, or an Arriba name without a distance), or an intergenic
    neighbour at most max_distance away. A neighbour further away is a label, not the gene."""
    names, dist = ev["names" + side], arriba_distances(ev)
    sf_gene = names[0] if "starfusion" in ev["callers"] else None
    return [n for n in names if n == sf_gene or dist.get(n, 0) <= max_distance]


def pair_entry(ev, pairs, n1, n2):
    """The fusion-report/OncoKB entry for the gene pair (n1 x n2 genic names) in the event's
    orientation -> (pair, entry, reversed). Unoriented Arriba calls also try the reverse
    orientation; such a match is reversed=True and can't give tier 1 (see tier_event)."""
    for a in n1:
        for b in n2:
            if f"{a}--{b}" in pairs:
                return f"{a}--{b}", pairs[f"{a}--{b}"], False
    if not ev["oriented"]:
        for a in n1:
            for b in n2:
                if f"{b}--{a}" in pairs:
                    return f"{b}--{a}", pairs[f"{b}--{a}"], True
    return None, None, False


def support(ev, bp):
    """-> (unique molecules, junction, spanning, start sites, PCR copies, source). bpview's projected
    non-duplicate fragments when bpview drew the event; otherwise Arriba's count after its dedup, or
    STAR-Fusion's (deduplicated only when STAR-Fusion read a deduplicated junction file)."""
    a, s = ev.get("arriba") or {}, ev.get("starfusion") or {}
    if bp:
        nd, dup = bp.get("fragments_nondup") or {}, bp.get("fragments_dup") or {}
        j, sp = nd.get("junction", 0), nd.get("spanning", 0)
        return (j + sp, j, sp, bp.get("junction_start_sites"),
                dup.get("junction", 0) + dup.get("spanning", 0), "bpview")
    if a:
        return (a.get("split_reads1", 0) + a.get("split_reads2", 0) + a.get("discordant_mates", 0),
                None, None, None, None, "arriba")
    return s.get("junction_reads", 0) + s.get("spanning_frags", 0), None, None, None, None, "starfusion"


def side_sense(side, strand):
    """The keep value under which the kept part of a partner runs with its transcript."""
    return {("5", "+"): "le", ("5", "-"): "ge", ("3", "+"): "ge", ("3", "-"): "le"}.get((side, strand))


def exon_site(tx, side, chrom, pos, keep):
    """-> (transcript id, exon number, site) for one side: the exon of the kept part next to the
    breakpoint (the 5' partner's last, the 3' partner's first) and where the breakpoint lies:
    exon_boundary | exon | intron | outside. (None, None, None) when the kept part runs antisense."""
    if not tx or tx["chrom"] != chrom or keep != side_sense(side, tx["strand"]):
        return None, None, None
    le = keep == "le"
    kept = [x for x in tx["exons"] if (x[0] <= pos if le else x[1] >= pos)]
    if not kept:
        return tx["id"], None, "outside"
    ex = max(kept, key=lambda x: x[1]) if le else min(kept, key=lambda x: x[0])
    if ex[0] <= pos <= ex[1]:
        site = "exon_boundary" if pos == (ex[1] if le else ex[0]) else "exon"
    elif len(kept) < len(tx["exons"]):
        site = "intron"
    else:
        site = "outside"
    return tx["id"], ex[2], site


def _kept_bases(segs, pos, keep):
    """Bases of [start, end] segments on the kept side of a breakpoint (keep 'le': <= pos; 'ge': >= pos)."""
    if keep == "le":
        return sum(max(0, min(b, pos) - a + 1) for a, b in segs)
    return sum(max(0, b - max(a, pos) + 1) for a, b in segs)


def mane_frame(t5, t3, ev, site5, site3):
    """Reading frame from the MANE Select CDS phase of both partners, or None when not computable.
    Both breakpoints must be at exon boundaries inside the partner's CDS: a splice junction is exact, while a
    mid-exon junction can carry inserted or microhomology bases that coordinates can't show. n5 = coding bases
    of the 5' partner kept up to and including its breakpoint; n3 = coding bases of the 3' partner
    before its breakpoint. In frame when n5 and n3 have the same phase (mod 3). Stop codons created at
    the junction are not detected, so this never says stop-codon."""
    if not t5 or not t3 or not t5.get("cds") or not t3.get("cds"):
        return None
    if site5 != "exon_boundary" or site3 != "exon_boundary":
        return None
    for t, pos in ((t5, ev["pos1"]), (t3, ev["pos2"])):
        if not (min(a for a, _ in t["cds"]) <= pos <= max(b for _, b in t["cds"])):
            return None  # in a UTR: the frame depends on which start codon is used
    n5 = _kept_bases(t5["cds"], ev["pos1"], ev["keep1"])
    n3 = sum(b - a + 1 for a, b in t3["cds"]) - _kept_bases(t3["cds"], ev["pos2"], ev["keep2"])
    return "in-frame" if n5 % 3 == n3 % 3 else "out-of-frame"


def gene_tx(names, transcripts, chrom):
    for n in names:
        t = transcripts.get(n)
        if t and t["chrom"] == chrom:
            return t
    return None


def kinase_info(ev, domains, min_fraction):
    """-> (kinase_gene, kinase_retained, kinase_domain_fraction). The partner with a protein kinase
    domain, the 3' partner when both have one; the fraction of its domain on the kept side."""
    hits = []
    for side, names, chrom, pos, keep in (("5", ev["names1"], ev["chr1"], ev["pos1"], ev["keep1"]),
                                          ("3", ev["names2"], ev["chr2"], ev["pos2"], ev["keep2"])):
        # the display gene first, then the other names the callers gave this side
        gene = next((n for n in names if n in domains and domains[n][0] == chrom), None)
        if gene is None:
            continue
        d = domains[gene]
        total = sum(b - a + 1 for a, b in d[1])
        kept = sum(max(0, (min(b, pos) if keep == "le" else b) - (a if keep == "le" else max(a, pos)) + 1)
                   for a, b in d[1])
        hits.append((side, gene, kept / total if total else 0.0))
    if not hits:
        return None, None, None
    side, gene, frac = hits[-1]
    return gene, frac >= min_fraction, round(frac, 3)


def breakpoint_offset(ev):
    s = ev.get("starfusion") or {}
    if "arriba" not in ev.get("callers", []) or not s.get("left_breakpoint"):
        return None
    try:
        p5 = int(s["left_breakpoint"].split(":")[1])
        p3 = int(s["right_breakpoint"].split(":")[1])
    except (IndexError, ValueError, KeyError):
        return None
    return max(abs(p5 - ev["pos1"]), abs(p3 - ev["pos2"]))


def arriba_filters(a):
    """Arriba's filters column -> filter names, without duplicates(n) (it only reports removed copies)."""
    out = []
    for f in (a.get("filters") or "").split(","):
        name = re.sub(r"\(.*\)", "", f).strip()
        if name and name not in (".", "duplicates") and name not in out:
            out.append(name)
    return out


def known_dbs(ev, entry):
    dbs = list((entry or {}).get("Databases") or [])
    for t in ((ev.get("arriba") or {}).get("tags") or "").split(","):
        name = t.split("=")[0].strip()
        if name == "Mitelman" and name not in dbs:
            dbs.append(name)
    return dbs


def pon_counts(ev, pon, sample_id, tol, n1, n2):
    """PoN recurrence: other samples with the same genic gene pair (either orientation; 0 when a
    side has no genic name) and with the same junction (+-tol)."""
    if pon is None:
        return None, None, None
    others = pon["samples"] - {sample_id}
    pair = set()
    for a in n1:
        for b in n2:
            pair |= pon["pairs"].get(frozenset((a, b)), set())
    near = pon["by_contig"].get((ev["chr1"], ev["keep1"], ev["chr2"], ev["keep2"]), [])
    junction = {r["sample_id"] for r in near
                if r["sample_id"] != sample_id
                and abs(r["pos5"] - ev["pos1"]) <= tol and abs(r["pos3"] - ev["pos2"]) <= tol}
    return len(pair - {sample_id}), len(junction), len(others)


def evidence_fails(row, opts):
    """Why the fusion evidence is too weak to carry an OncoKB hit: (basic, support).
    basic: frame and caller evidence, needed for tiers 1 and 2 through OncoKB.
    support: unique molecules, needed for tier 1 only."""
    basic, support = {}, {}
    # only the caller's frame counts as evidence; a MANE-derived frame (reading_frame_source "mane") is shown only
    if (row.get("reading_frame_source") in (None, "arriba") and row["reading_frame"] in ("out-of-frame", "stop-codon")
            and not opts.allow_out_of_frame):
        basic["out_of_frame"] = f"reading frame {row['reading_frame']}"
    if row["n_callers"] < 2 and row["arriba_confidence"] not in ("high", "medium"):
        basic["weak_caller_evidence"] = "one caller and Arriba confidence not high/medium"
    if row["unique_molecules"] < opts.tier1_min_unique_molecules:
        support["low_unique_molecules"] = (f"{row['unique_molecules']} unique molecule(s) < "
                                           f"{opts.tier1_min_unique_molecules}")
    return basic, support


def tier_event(row, opts):
    level, dx = row["oncokb_highest_level"] or "", row["oncokb_highest_dx_level"] or ""
    onc = row["oncokb_oncogenic"] or ""
    curated = sorted(CURATED_DB & set(row["known_db"]))
    basic, low = evidence_fails(row, opts)
    if row["pair_reversed"]:
        # The annotation matched the reversed pair of a call whose orientation is unknown, so
        # it may not describe this event: no tier 1, and no curated tier 2 without evidence.
        low = dict(low, orientation_unknown="orientation unknown; annotation matched the reversed pair")
    # OncoKB levels are often partner-agnostic (any NTRK or ALK fusion), so an OncoKB hit
    # only counts with fusion evidence; without it, only a curated known pair keeps tier 2.
    actionable = level in ACTIONABLE or dx in ACTIONABLE_DX
    oncogenic = actionable or onc in ONCOGENIC or level in ONCOGENIC_LEVELS or dx in ONCOGENIC_DX
    lab = "OncoKB " + ", ".join(x for x in (onc, level, dx) if x and x != "Unknown")
    reasons, flags = [], []
    if actionable and not basic and not low:
        tier, basis = 1, "oncokb_actionable"
        reasons.append(lab)
    elif oncogenic and not basic:
        tier, basis = 2, "oncokb_actionable" if actionable else "oncokb_oncogenic"
        if actionable and low:
            flags += list(low)
        elif "orientation_unknown" in low:
            flags.append("orientation_unknown")
        reasons.append(lab + (" (not tier 1: " + "; ".join(low.values()) + ")" if actionable and low else ""))
    elif curated and not (row["pair_reversed"] and basic):
        tier, basis = 2, "curated_known"
        flags += list(basic) + list(low)
        reasons.append("known fusion: " + ", ".join(curated))
        if basic or low:
            reasons.append("evidence: " + "; ".join(list(basic.values()) + list(low.values())))
        if oncogenic:
            reasons.append(lab + " (not tier 1)")
    else:
        tier, basis = 3, "none"
        if oncogenic or curated:
            flags += list(basic) + ([k for k in low if k == "orientation_unknown"])
            known = (lab if oncogenic else "known fusion: " + ", ".join(curated))
            reasons.append(known + ", but " + "; ".join(list(basic.values())
                                                      + ([low["orientation_unknown"]] if "orientation_unknown" in low else [])))
        else:
            reasons.append("not known or actionable")
    if row["arriba_status"] == "rescued":
        # Rescue rule: a rescued call stays above tier 3 only when OncoKB recognises that exact
        # gene pair, whatever its molecule count (partner-agnostic levels and database tags don't count),
        # and even then it never reaches tier 1.
        if not row["oncokb_pair_specific"] and tier < 3:
            flags.append("rescue_capped")
            reasons.append("rescued call, and OncoKB doesn't recognise this exact gene pair: tier III")
            tier = 3
        elif tier == 1:
            flags.append("rescue_capped")
            reasons.append("rescued call: never tier I")
            tier = 2
    if row["adjacent_readthrough"] and not row["oncokb_pair_specific"] and tier < 3:
        flags.append("readthrough_capped")
        reasons.append(f"adjacent-gene read-through ({row['readthrough_distance']:,} bp apart, under "
                       f"{READTHROUGH_MAX_DISTANCE:,} bp, no gene between), and OncoKB doesn't recognise this exact "
                       "gene pair: tier III")
        tier = 3
    if not row["unique_molecules"] and tier < 3:
        # No counted molecule, no tier above III, whatever the path (a kept curated call too)
        flags.append("no_molecules")
        reasons.append("no counted molecule: tier III")
        tier = 3
    if row["low_support"]:
        flags.append("low_support")
        if tier == 1:
            tier = 2
            reasons.append(f"low support ({row['unique_molecules']} unique molecule(s)): capped at tier 2")
        else:
            reasons.append(f"low support ({row['unique_molecules']} unique molecule(s))")
    if row["arriba_status"] == "rescued":
        flags.append("rescued")
        reasons.append("rescued from Arriba's discarded calls (" + ", ".join(row["arriba_filters"] or []) + ")")
    n, size = row["pon_gene_pair_count"], row["pon_size"]
    recurrent = None if n is None or not size else (n >= opts.pon_min_samples and n / size >= opts.pon_min_fraction)
    if recurrent:
        flags.append("pon_recurrent")
        capped = 2 if tier == 1 else 3
        if capped != tier:
            flags.append("pon_capped")
            reasons.append(f"recurrent in PoN ({n}/{size} samples): capped at tier {capped}")
            tier = capped
        else:
            reasons.append(f"recurrent in PoN ({n}/{size} samples)")
    return dict(tier=tier, tier_basis=basis, tier_flags=flags, tier_reason="; ".join(reasons),
                pon_recurrent=recurrent)


def _has(path):
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def load_arriba_raw(path):
    """bpview's Arriba events, each with its source row (arriba["raw"]) for the detail columns."""
    out = []
    for r in bpview._tsv(path):
        e = bpview.arriba_event(r)
        if e:
            e["arriba"]["raw"] = r
            out.append(e)
    return out


def load_starfusion_raw(path):
    """bpview's STAR-Fusion events, each with its source row (starfusion["raw"]). load_starfusion
    skips rows without three-part breakpoints; the same rows are skipped here, so the two line up."""
    rows = [r for r in bpview._tsv(path)
            if len(r.get("LeftBreakpoint", "").split(":")) == 3 and len(r.get("RightBreakpoint", "").split(":")) == 3]
    evs = bpview.load_starfusion(path)
    assert len(rows) == len(evs), (len(rows), len(evs))
    for r, e in zip(rows, evs):
        e["starfusion"]["raw"] = r
    return evs


def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _str(v):
    return None if v in (None, "", ".") else v


SEQ_CAP = 2000  # characters; one pathological Arriba row must not bloat the events response


def _cap(v):
    v = _str(v)
    if v is None or len(v) <= SEQ_CAP:
        return v
    return v[:SEQ_CAP] + f"...[truncated, {len(v)} characters]"


def _split(v, sep):
    return [x.strip() for x in (v or "").split(sep) if x.strip() and x.strip() != "."]


def caller_details(ev, entry):
    """The 0.5.0 detail columns: raw caller and annotator values (labelled as such, never support)."""
    a = ((ev.get("arriba") or {}).get("raw")) or {}
    s = ((ev.get("starfusion") or {}).get("raw")) or {}
    onc = (entry or {}).get("OncoKB") or {}
    dom = (a.get("retained_protein_domains") or "").split("|") if _str(a.get("retained_protein_domains")) else []
    dup = re.search(r"duplicates\((\d+)\)", a.get("filters") or "")
    tf = {"True": True, "False": False}
    return dict(
        arriba_gene1=_str(a.get("gene1")), arriba_gene2=_str(a.get("gene2")),
        arriba_breakpoint1=_str(a.get("breakpoint1")), arriba_breakpoint2=_str(a.get("breakpoint2")),
        arriba_strand1=_str(a.get("strand1(gene/fusion)")), arriba_strand2=_str(a.get("strand2(gene/fusion)")),
        arriba_split_reads1=_int(a.get("split_reads1")), arriba_split_reads2=_int(a.get("split_reads2")),
        arriba_coverage1=_int(a.get("coverage1")), arriba_coverage2=_int(a.get("coverage2")),
        arriba_site1=_str(a.get("site1")), arriba_site2=_str(a.get("site2")),
        arriba_filters_raw=_str(a.get("filters")), arriba_duplicates=int(dup.group(1)) if dup else (0 if a else None),
        arriba_tags=_split(a.get("tags"), ","),
        arriba_retained_domains1=_split(dom[0], ",") if len(dom) > 0 else [],
        arriba_retained_domains2=_split(dom[1], ",") if len(dom) > 1 else [],
        arriba_closest_genomic_breakpoint1=_str(a.get("closest_genomic_breakpoint1")),
        arriba_closest_genomic_breakpoint2=_str(a.get("closest_genomic_breakpoint2")),
        arriba_gene_id1=_str(a.get("gene_id1")), arriba_gene_id2=_str(a.get("gene_id2")),
        arriba_transcript_id1=_str(a.get("transcript_id1")), arriba_transcript_id2=_str(a.get("transcript_id2")),
        arriba_fusion_transcript=_cap(a.get("fusion_transcript")), arriba_peptide_sequence=_cap(a.get("peptide_sequence")),
        starfusion_left_gene=_str(s.get("LeftGene")), starfusion_right_gene=_str(s.get("RightGene")),
        starfusion_left_breakpoint=_str(s.get("LeftBreakpoint")), starfusion_right_breakpoint=_str(s.get("RightBreakpoint")),
        starfusion_est_junction=_float(s.get("est_J")), starfusion_est_spanning=_float(s.get("est_S")),
        starfusion_splice_type=_str(s.get("SpliceType")),
        # YES_LDAS / NO_LDAS -> true / false; anything else unknown
        starfusion_large_anchor_support={"YES_LDAS": True, "NO_LDAS": False}.get(s.get("LargeAnchorSupport")),
        starfusion_left_break_dinuc=_str(s.get("LeftBreakDinuc")), starfusion_right_break_dinuc=_str(s.get("RightBreakDinuc")),
        starfusion_left_break_entropy=_float(s.get("LeftBreakEntropy")),
        starfusion_right_break_entropy=_float(s.get("RightBreakEntropy")),
        starfusion_annots=bpview.parse_annots(s.get("annots")) if s else [],
        oncokb_highest_sensitive_level=_str(onc.get("HIGHEST_SENSITIVE_LEVEL")),
        oncokb_highest_resistance_level=_str(onc.get("HIGHEST_RESISTANCE_LEVEL")),
        oncokb_highest_px_level=_str(onc.get("HIGHEST_PX_LEVEL")),
        oncokb_treatments=[f"{lvl}: {onc[lvl]}" for lvl in ONCOKB_LEVELS if _str(onc.get(lvl))],
        oncokb_mutation_effect_citations=_split(onc.get("MUTATION_EFFECT_CITATIONS"), ";"),
        oncokb_treatment_citations=_split(onc.get("TX_CITATIONS"), ";"),
        oncokb_dx_citations=_split(onc.get("DX_CITATIONS"), ";"),
        oncokb_px_citations=_split(onc.get("PX_CITATIONS"), ";"),
        oncokb_gene_in_oncokb=tf.get(onc.get("GENE_IN_ONCOKB")),
        oncokb_variant_in_oncokb=tf.get(onc.get("VARIANT_IN_ONCOKB")),
        oncokb_variant_summary=_str(onc.get("ONCOKB_VARIANT_SUMMARY")),
        oncokb_matched_alterations=_split(onc.get("ONCOKB_MATCHED_ALTERATIONS"), ";"),
        pair_fii_explained=_str((entry or {}).get("Explained FII")))


def caller_sides(ev, b):
    """The per-side and per-caller BAM columns (see COLUMNS)."""
    out = {}
    a = (ev.get("arriba") or {}).get("raw") or {}
    if a:
        s1, s2 = _int(a.get("split_reads1")), _int(a.get("split_reads2"))
        try:
            c1, p1 = a["breakpoint1"].rsplit(":", 1)
            flipped = not (c1 == ev["chr1"] and int(p1) == ev["pos1"])
        except (KeyError, ValueError):
            flipped = False
        out.update(arriba_split_reads_5p=s2 if flipped else s1, arriba_split_reads_3p=s1 if flipped else s2)
    cb = (b or {}).get("callers_bam") or {}
    for c in ("arriba", "starfusion"):
        v = cb.get(c) if c in ev["callers"] else None
        # No names to look up (Arriba writes read_identifiers only for kept calls, so a rescued call
        # cites none): unknown, not zero.
        if v is None or not v.get("cited"):
            continue
        nd, du = v.get("fragments_nondup") or {}, v.get("fragments_dup") or {}
        out.update({f"{c}_bam_split_reads_5p": v.get("split_5p"), f"{c}_bam_split_reads_3p": v.get("split_3p"),
                    f"{c}_bam_molecules": nd.get("junction", 0) + nd.get("spanning", 0),
                    f"{c}_bam_pcr_copies": du.get("junction", 0) + du.get("spanning", 0),
                    f"{c}_cited_missing": v.get("cited_missing")})
    return out


def build_rows(opts):
    arriba = load_arriba_raw(opts.arriba_tsv) if _has(opts.arriba_tsv) else []
    starf = load_starfusion_raw(opts.starfusion_tsv) if _has(opts.starfusion_tsv) else []
    events = bpview.merge_events(arriba, starf, opts.merge_tolerance)
    pairs, bp, pon = load_pairs(opts.oncokb_json), load_bpview(opts.bpview_events), load_pon(opts.pon)
    transcripts, domains = load_transcripts(opts.gtf), load_kinase_domains(opts.protein_domains)
    spans = gene_spans(transcripts)
    rows = []
    for ev in events:
        a, s = ev.get("arriba") or {}, ev.get("starfusion") or {}
        key = event_key(ev["chr1"], ev["pos1"], ev["keep1"], ev["chr2"], ev["pos2"], ev["keep2"])
        b = bp.get(key)
        n1, n2 = (genic_names(ev, i, opts.max_intergenic_distance) for i in ("1", "2"))
        fusion, entry, reversed_match = pair_entry(ev, pairs, n1, n2)
        onc = (entry or {}).get("OncoKB") or {}
        u, mj, ms, sites, dup, src = support(ev, b)
        npair, njunc, size = pon_counts(ev, pon, opts.sample_id, opts.merge_tolerance, n1, n2)
        t5, x5, site5 = exon_site(gene_tx(ev["names1"], transcripts, ev["chr1"]), "5", ev["chr1"], ev["pos1"], ev["keep1"])
        t3, x3, site3 = exon_site(gene_tx(ev["names2"], transcripts, ev["chr2"]), "3", ev["chr2"], ev["pos2"], ev["keep2"])
        kg, kr, kf = kinase_info(ev, domains, opts.kinase_min_fraction)
        caller_frame = a.get("reading_frame") if a.get("reading_frame") not in (None, ".", "") else None
        frame, frame_src = (caller_frame, "arriba") if caller_frame else (None, None)
        if frame is None:
            tx5, tx3 = (gene_tx(ev["names" + i], transcripts, ev["chr" + i]) for i in ("1", "2"))
            frame = mane_frame(tx5, tx3, ev, site5, site3)
            frame_src = "mane" if frame else None
        status = a.get("status") if a else None
        bp5 = f"{ev['chr1']}:{ev['pos1']}:{strand5(ev)}"
        bp3 = f"{ev['chr2']}:{ev['pos2']}:{strand3(ev)}"
        o5, o3 = (strand5(ev), strand3(ev)) if ev["oriented"] else (".", ".")
        row = dict(
            sample_id=opts.sample_id,
            event_id=f"{ev['chr1']}:{ev['pos1']}:{o5}>{ev['chr2']}:{ev['pos2']}:{o3}",
            fusion=fusion or f"{ev['gene1']}--{ev['gene2']}",
            gene5=ev["gene1"], gene3=ev["gene2"],
            gene5_aliases=list(ev["names1"]), gene3_aliases=list(ev["names2"]),
            breakpoint5=bp5, breakpoint3=bp3, orientation_known=bool(ev["oriented"]),
            transcript5=t5, transcript3=t3, exon5=x5, exon3=x3, site5=site5, site3=site3,
            # tiers use the caller's frame only; a MANE-derived frame is shown, marked by reading_frame_source
            reading_frame=frame, reading_frame_source=frame_src,
            canonical_junction=(site5 == "exon_boundary" and site3 == "exon_boundary") if site5 and site3 else None,
            splice_type=bpview.splice_type(ev),
            kinase_gene=kg, kinase_retained=kr, kinase_domain_fraction=kf,
            callers=list(ev["callers"]), n_callers=len(ev["callers"]),
            caller_arriba="arriba" in ev["callers"], caller_starfusion="starfusion" in ev["callers"],
            arriba_status=status,
            arriba_filters=arriba_filters(a) if status == "rescued" else None,
            arriba_confidence=a.get("confidence"), arriba_type=a.get("type"),
            arriba_split_reads=(a.get("split_reads1", 0) + a.get("split_reads2", 0)) if a else None,
            arriba_discordant_mates=a.get("discordant_mates") if a else None,
            starfusion_junction_reads=s.get("junction_reads") if s else None,
            starfusion_spanning_frags=s.get("spanning_frags") if s else None,
            starfusion_ffpm=s.get("ffpm") if s else None,
            breakpoint_offset=breakpoint_offset(ev),
            unique_molecules=u, molecules_junction=mj, molecules_spanning=ms, start_sites=sites,
            pcr_copies=dup, support_source=src, low_support=u < opts.min_unique_molecules,
            known_db=known_dbs(ev, entry),
            oncokb_oncogenic=onc.get("ONCOGENIC") or None,
            oncokb_mutation_effect=onc.get("MUTATION_EFFECT") or None,
            oncokb_highest_level=onc.get("HIGHEST_LEVEL") or None,
            oncokb_highest_dx_level=onc.get("HIGHEST_DX_LEVEL") or None,
            pair_clinical_tier=(entry or {}).get("clinical_tier"),
            pair_fii=(entry or {}).get("Fusion Indication Index (FII)"),
            pon_gene_pair_count=npair, pon_junction_count=njunc, pon_size=size,
            bpview_contig=b["contig"] if b else None, builder_version=__version__)
        row.update(caller_details(ev, entry))
        row["oncokb_pair_specific"] = oncokb_exact_pair(entry, fusion)
        rt = adjacent_readthrough(ev, spans) if spans else None
        row["adjacent_readthrough"] = (rt is not None) if spans else None
        row["readthrough_distance"] = rt
        row.update({c: None for c, _ in COLUMNS if c not in row})
        row.update(caller_sides(ev, b))
        row["curated_known"] = bool(CURATED_DB & set(row["known_db"]))
        row["pair_reversed"] = reversed_match
        row.update(tier_event(row, opts))
        row["prototype_rules"] = ((["rescue"] if status == "rescued" else [])
                                  + (["low_support"] if row["low_support"] else []))
        rows.append((row, bpview.priority(ev)))
    # Default order: tier; then exact curated or OncoKB-curated actionable/oncogenic pairs at a canonical
    # junction (both breakpoints at exon boundaries) first; then supported before low support,
    # unique molecules, and caller priority as before.
    def leads(r):
        exact = r["curated_known"] or r["oncokb_pair_specific"]
        # an event with no supporting molecule in the BAM doesn't lead
        return bool(r["canonical_junction"]) and bool(exact) and (r["unique_molecules"] or 0) >= 1

    rows.sort(key=lambda x: (x[0]["tier"], not leads(x[0]), x[0]["low_support"], -x[0]["unique_molecules"],
                             tuple(-v for v in x[1])))
    out, seen = [], {}
    for i, (row, _) in enumerate(rows):
        row["event_rank"] = i + 1
        n = seen.get(row["event_id"], 0)
        seen[row["event_id"]] = n + 1
        if n:  # one junction is one event after the merge; keep IDs unique if a caller repeats one
            row["event_id"] += f"#{n + 1}"
        out.append(row)
    return out


def write_meta(path, opts, rows):
    pon = load_pon(opts.pon)
    meta = [("contract_version", CONTRACT_VERSION), ("builder_version", __version__),
            ("prototype", "true" if opts.prototype else "false"),
            ("arriba_settings", opts.arriba_settings or ""),
            ("pon", f"{os.path.basename(opts.pon)} ({len(pon['samples'])} samples)" if pon else "none"),
            ("oncokb_exact_pair", EXACT_PAIR_RULE),
            ("events", str(len(rows)))]
    with open(path, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in meta:
            fh.write(f"{k}\t{v}\n")


# ----------------------------------------------------------------------------- rescue

def rescue(opts):
    """Arriba's kept calls plus the discarded ones the rescue rule brings back (see the docstring).
    -> (kept, rescued) counts."""
    def read(path):
        if not path or not os.path.isfile(path) or not os.path.getsize(path):
            return None, []
        with open(path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            cols = [h.lstrip("#") for h in hdr]
            return hdr, [dict(zip(cols, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]

    hdr, kept = read(opts.arriba_tsv)
    dhdr = None
    if opts.arriba_discarded_tsv and os.path.isfile(opts.arriba_discarded_tsv) and os.path.getsize(opts.arriba_discarded_tsv):
        with open(opts.arriba_discarded_tsv) as fh:
            dhdr = fh.readline().rstrip("\n").split("\t")
    hdr = hdr or dhdr
    if hdr is None:
        raise SystemExit("rescue: no Arriba header in either input")
    pairs = load_pairs(opts.oncokb_json)
    sf = (bpview.collapse_starfusion(bpview.load_starfusion(opts.starfusion_tsv))
          if opts.starfusion_tsv and os.path.isfile(opts.starfusion_tsv) and os.path.getsize(opts.starfusion_tsv) else [])
    tol = opts.merge_tolerance

    def index(events):
        by = {}
        for o in events:
            by.setdefault((o["chr1"], o["chr2"]), []).append(o)
        return by

    def near(e, by):
        f = bpview._flip(e)
        cands = by.get((e["chr1"], e["chr2"]), []) + (by.get((f["chr1"], f["chr2"]), []) if not e["oriented"] else [])
        return any(bpview._same(e, o, tol) or (not e["oriented"] and bpview._same(f, o, tol)) for o in cands)

    def discarded_rows():
        if dhdr is None:
            return
        cols = [h.lstrip("#") for h in dhdr]
        with open(opts.arriba_discarded_tsv) as fh:
            fh.readline()
            for l in fh:
                if l.strip():
                    yield dict(zip(cols, l.rstrip("\n").split("\t")))

    kept_ev = index([e for e in (bpview.arriba_event(r) for r in kept) if e])
    sf = index(sf)
    rescued, keys = [], set()
    for r in discarded_rows():
        e = bpview.arriba_event(r)
        if not e:
            continue
        a = e["arriba"]
        if a["split_reads1"] + a["split_reads2"] + a["discordant_mates"] < 1:
            continue
        fs = set(arriba_filters(a))
        if not fs or not fs <= SUPPORT_FILTERS:
            continue
        k = (e["chr1"], e["pos1"], e["keep1"], e["chr2"], e["pos2"], e["keep2"])
        if k in keys or near(e, kept_ev):
            continue
        e["callers"] = ["arriba"]
        n1, n2 = genic_names(e, "1", 0), genic_names(e, "2", 0)
        _, entry, _ = pair_entry(dict(e, oriented=False), pairs, n1, n2)
        onc = (entry or {}).get("OncoKB") or {}
        known = ("Mitelman" in (a.get("tags") or "") or bool((entry or {}).get("Databases"))
                 or (onc.get("ONCOGENIC") or "") in ONCOGENIC or bool(onc.get("HIGHEST_LEVEL")))
        if not (known or near(e, sf)):
            continue
        keys.add(k)
        rescued.append((r, e))
    with open(opts.out_arriba_tsv, "w") as fh:
        cols = [h.lstrip("#") for h in hdr]
        fh.write("\t".join(hdr + ["arriba_status"]) + "\n")
        for status, rs in (("kept", kept), ("rescued", [r for r, _ in rescued])):
            for r in rs:
                fh.write("\t".join([r.get(c, ".") for c in cols] + [status]) + "\n")
    if opts.out_fusion_list:
        lines = []
        if opts.fusion_list and os.path.isfile(opts.fusion_list):
            lines = [l.rstrip("\n") for l in open(opts.fusion_list) if l.strip()]
        listed = {l.split("\t")[0] for l in lines}
        with open(opts.out_fusion_list, "w") as fh:
            for l in lines:
                fh.write(l + "\n")
            for _, e in rescued:
                name = f"{e['gene1']}--{e['gene2']}"
                if name not in listed:
                    listed.add(name)
                    fh.write(name + "\n")
    return len(kept), len(rescued)


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (list, tuple)):
        return "|".join(str(x).replace("|", "/") for x in v)
    return str(v).replace("\t", " ").replace("\n", " ")


def write_tsv(path, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(c for c, _ in COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r[c]) for c, _ in COLUMNS) + "\n")


# ----------------------------------------------------------------------------- PoN

def build_pon(paths, out):
    n = 0
    with open(out, "w") as fh:
        fh.write("\t".join(PON_COLUMNS) + "\n")
        for p in paths:
            with open(p) as ih:
                for r in csv.DictReader(ih, delimiter="\t"):
                    c5, p5, s5 = r["breakpoint5"].rsplit(":", 2)
                    c3, p3, s3 = r["breakpoint3"].rsplit(":", 2)
                    fh.write("\t".join([r["sample_id"], r["gene5"], r["gene3"],
                                        c5, p5, "le" if s5 == "+" else "ge",
                                        c3, p3, "ge" if s3 == "+" else "le"]) + "\n")
                    n += 1
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("rescue")
    r.add_argument("--arriba-tsv", required=True, help="Arriba fusions.tsv")
    r.add_argument("--arriba-discarded-tsv", required=True, help="Arriba fusions.discarded.tsv")
    r.add_argument("--starfusion-tsv")
    r.add_argument("--oncokb-json")
    r.add_argument("--fusion-list", help="bpview's fusion list; rescued pairs are appended")
    r.add_argument("--out-arriba-tsv", required=True)
    r.add_argument("--out-fusion-list")
    r.add_argument("--merge-tolerance", type=int, default=10)
    b = sub.add_parser("build")
    b.add_argument("--sample-id", required=True)
    b.add_argument("--arriba-tsv")
    b.add_argument("--starfusion-tsv")
    b.add_argument("--oncokb-json", help="fusions_with_oncokb.json from oncokb-fusion-annotation")
    b.add_argument("--bpview-events", help="bpview events.json (optional)")
    b.add_argument("--pon", help="recurrence panel from `fusion-events pon` (optional)")
    b.add_argument("--gtf", help="GTF with MANE Select transcripts, for exon5/exon3 (optional)")
    b.add_argument("--protein-domains", default=os.path.join(os.path.dirname(os.path.realpath(__file__)), "ref",
                                                            "protein_domains_hg38_GRCh38_v2.5.1.gff3"),
                   help="Arriba's protein-domain GFF3, for kinase_* ('.' for none)")
    b.add_argument("--out-tsv", required=True)
    b.add_argument("--out-meta-tsv", help="fusion_events_meta as key/value TSV (optional)")
    b.add_argument("--arriba-settings", help="recorded in the meta table, e.g. max_sensitivity")
    b.add_argument("--prototype", action="store_true", help="mark the output as built by prototype rules")
    b.add_argument("--merge-tolerance", type=int, default=10)
    b.add_argument("--tier1-min-unique-molecules", type=int, default=2)
    b.add_argument("--min-unique-molecules", type=int, default=2,
                   help="below this, an event is low_support: shown, capped at tier 2")
    b.add_argument("--kinase-min-fraction", type=float, default=0.9,
                   help="kinase_retained when at least this fraction of the kinase domain is kept")
    b.add_argument("--allow-out-of-frame", action="store_true",
                   help="do not hold out-of-frame / stop-codon events back from tiers 1-2")
    b.add_argument("--max-intergenic-distance", type=int, default=0,
                   help="an intergenic breakpoint takes a neighbouring gene's annotation only within this "
                        "many bp (0: only the gene the breakpoint is in)")
    b.add_argument("--pon-min-samples", type=int, default=1)
    b.add_argument("--pon-min-fraction", type=float, default=0.10)
    p = sub.add_parser("pon")
    p.add_argument("--events-tsv", nargs="+", required=True)
    p.add_argument("--out-tsv", required=True)
    opts = ap.parse_args(argv)
    if opts.cmd == "pon":
        n = build_pon(opts.events_tsv, opts.out_tsv)
        print(f"fusion-events {__version__}: PoN with {n} junctions from {len(opts.events_tsv)} samples")
        return
    if opts.cmd == "rescue":
        k, n = rescue(opts)
        print(f"fusion-events {__version__}: {k} kept Arriba calls, {n} rescued")
        return
    rows = build_rows(opts)
    write_tsv(opts.out_tsv, rows)
    if opts.out_meta_tsv:
        write_meta(opts.out_meta_tsv, opts, rows)
    tiers = {t: sum(1 for r in rows if r["tier"] == t) for t in (1, 2, 3)}
    print(f"fusion-events {__version__}: {len(rows)} events; tier 1/2/3 = {tiers[1]}/{tiers[2]}/{tiers[3]}; "
          f"{sum(1 for r in rows if r['low_support'])} low support, "
          f"{sum(1 for r in rows if r['arriba_status'] == 'rescued')} rescued")


if __name__ == "__main__":
    main()
