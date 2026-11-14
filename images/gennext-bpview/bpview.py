#!/usr/bin/env python3
"""bpview: per-breakpoint fusion contigs built by *projecting* existing alignments.

For every fusion breakpoint called by Arriba and/or STAR-Fusion, build a synthetic contig

    [5' partner genomic segment ending at bp1] + [3' partner genomic segment starting at bp2]

in fusion-transcript orientation, and lift the STAR genome alignments of the sample into
contig coordinates. Nothing is realigned, so what is shown is exactly what the callers saw:

  * split reads (primary + SA supplementary) whose parts land on the two segments are merged
    into one alignment that runs straight through the junction           -> XF:Z:junction
    (anchors shorter than --min-anchor on either side, or a gap > 5 bp)    -> XF:Z:junction_weak
  * mates on opposite segments                                           -> XF:Z:spanning
  * reads running past the breakpoint in the wild-type configuration are
    soft-clipped at the breakpoint                                       -> XF:Z:wt5 / wt3
  * chimeric reads ending at the breakpoint with a partner elsewhere     -> XF:Z:chimeric_other
  * everything else in the segments                                      -> XF:Z:other

Other tags: XB:i distance of the nearest block end to the junction, YC:Z RGB colour by class
(igv.js colorBy tag), XK:Z sort key (evidence classes first). QNAME, duplicate flag, MAPQ, RX,
MI and NH are copied from the source alignment. Counts in the events JSON are recomputed from the
projected records, so the number reported equals the number of fragments drawn.

Outputs (one set per sample, one contig per event):
  <prefix>.fa/.fa.fai    contigs
  <prefix>.bam/.bam.bai  projected alignments (coordinate sorted)
  <prefix>.bed[.gz/.tbi] partner exons (all transcripts in the GTF) in contig coordinates
  <prefix>.gtf[.gz/.tbi] the same exons as GTF, with gene_name, transcript_id and
                         orig_coord_info "chrom,start,end,strand" (FusionInspector convention)
  <prefix>.events.json   event list with breakpoints, callers and support per class

Each event also carries the fields the FusionInspector abridged TSV provided to the IGV view,
under FusionInspector-like names: left/right_local_breakpoint, left/right_breakpoint
(chr:pos:strand), junction_read_count, spanning_frag_count, splice_type, annots and
left/right_gene_locus. The counts are unique (non-duplicate) fragments; the *_with_dups
counts include PCR duplicates, which is closer to what FusionInspector reported.
"""
import argparse
import collections
import gzip
import hashlib
import json
import os
import re
import sys

import pysam

__version__ = "0.2.0"

REF_CONSUME = {0, 2, 3, 7, 8}
QRY_CONSUME = {0, 1, 4, 7, 8}
EVIDENCE = {"junction", "junction_weak", "spanning", "chimeric_other"}
SORT_ORDER = {"junction": "1", "spanning": "2", "junction_weak": "3", "wt3": "4", "wt5": "5",
              "chimeric_other": "6", "other": "7"}
COLOURS = {"junction": "21,128,92", "junction_weak": "110,190,160", "spanning": "140,70,170",
           "wt5": "52,88,200", "wt3": "214,92,30", "chimeric_other": "215,160,0", "other": "175,178,185"}
CONFIDENCE_RANK = {"high": 3, "medium": 2, "low": 1}


# ----------------------------------------------------------------------------- utilities

def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def gene_names(name):
    """All gene names in a caller's gene field, nearest first.
    Arriba: 'GENE', or 'GENE1(1234),GENE2(567)' for a breakpoint between genes (distances);
    STAR-Fusion: 'GENE^ENSG...'."""
    out = []
    for part in name.split("^")[0].split(","):
        m = re.match(r"\s*([^()]+?)\s*(?:\((\d+)\))?\s*$", part)
        if m and m.group(1):
            out.append((int(m.group(2)) if m.group(2) else -1, m.group(1)))
    return [n for _, n in sorted(out, key=lambda x: x[0])] or [name]


def clean_gene(name):
    return gene_names(name)[0]


def safe(s):
    return re.sub(r"[^A-Za-z0-9._-]", "_", s)


class Segment:
    """Genomic interval [start, end] (1-based inclusive) mapped onto contig [off, off+len)."""

    def __init__(self, chrom, start, end, rc, off):
        self.chrom, self.start, self.end, self.rc, self.off = chrom, start, end, rc, off
        self.len = end - start + 1

    def to_contig(self, g):  # genomic 1-based -> contig 0-based
        return self.off + ((self.end - g) if self.rc else (g - self.start))

    def as_dict(self):
        return dict(chrom=self.chrom, start=self.start, end=self.end, rc=self.rc, offset=self.off)


def read_gtf(gtf):
    """gene_name -> {chrom, start, end, strand, exons[(start, end, transcript_id)]}."""
    out = {}
    op = gzip.open if gtf.endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            m = re.search(r'gene_name "([^"]+)"', f[8])
            if not m:
                continue
            g = out.setdefault(m.group(1), {"chrom": f[0], "start": int(f[3]), "end": int(f[4]),
                                            "strand": f[6], "exons": []})
            if f[2] == "gene":
                g.update(chrom=f[0], start=int(f[3]), end=int(f[4]), strand=f[6])
            else:
                g["start"], g["end"] = min(g["start"], int(f[3])), max(g["end"], int(f[4]))
                tid = re.search(r'transcript_id "([^"]+)"', f[8])
                g["exons"].append((int(f[3]), int(f[4]), tid.group(1) if tid else ""))
    return out


# ----------------------------------------------------------------------------- events

def _tsv(path):
    with open(path) as fh:
        hdr = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= len(hdr):
                yield dict(zip(hdr, f))


def _cp(s):
    p = s.split(":")
    return p[0], int(p[1])


def load_arriba(path):
    """keep = 'le' when the partner's retained sequence lies at coordinates <= breakpoint.
    Arriba direction 'downstream' means the other partner is joined above the breakpoint."""
    out = []
    for r in _tsv(path):
        try:
            (c1, p1), (c2, p2) = _cp(r["breakpoint1"]), _cp(r["breakpoint2"])
        except (KeyError, ValueError, IndexError):
            continue
        d1, d2 = r.get("direction1", ""), r.get("direction2", "")
        if d1 not in ("upstream", "downstream") or d2 not in ("upstream", "downstream"):
            continue
        s1 = r.get("strand1(gene/fusion)", "./.").split("/")[-1]
        s2 = r.get("strand2(gene/fusion)", "./.").split("/")[-1]
        n1, n2 = gene_names(r["gene1"]), gene_names(r["gene2"])
        out.append(dict(gene1=n1[0], gene2=n2[0], names1=n1, names2=n2,
                        intergenic1="(" in r["gene1"], intergenic2="(" in r["gene2"],
                        chr1=c1, pos1=p1, keep1="le" if d1 == "downstream" else "ge",
                        chr2=c2, pos2=p2, keep2="le" if d2 == "downstream" else "ge",
                        oriented=s1 in "+-" and s2 in "+-" and s1 != "." and s2 != ".",
                        arriba=dict(confidence=r.get("confidence"), type=r.get("type"),
                                    split_reads1=int(r.get("split_reads1") or 0),
                                    split_reads2=int(r.get("split_reads2") or 0),
                                    discordant_mates=int(r.get("discordant_mates") or 0),
                                    reading_frame=r.get("reading_frame"),
                                    site1=r.get("site1"), site2=r.get("site2"),
                                    tags=r.get("tags"),
                                    gene1=r["gene1"], gene2=r["gene2"],
                                    breakpoint1=r["breakpoint1"], breakpoint2=r["breakpoint2"])))
    return out


def load_starfusion(path):
    """LeftBreakpoint is the 5' partner. chr:pos:strand; '+' 5' keeps <= bp, '+' 3' keeps >= bp."""
    out = []
    for r in _tsv(path):
        try:
            c1, p1, s1 = r["LeftBreakpoint"].split(":")
            c2, p2, s2 = r["RightBreakpoint"].split(":")
        except (KeyError, ValueError):
            continue
        g1, g2 = clean_gene(r["LeftGene"]), clean_gene(r["RightGene"])
        out.append(dict(gene1=g1, gene2=g2, names1=[g1], names2=[g2],
                        intergenic1=False, intergenic2=False,
                        chr1=c1, pos1=int(p1), keep1="le" if s1 == "+" else "ge",
                        chr2=c2, pos2=int(p2), keep2="ge" if s2 == "+" else "le",
                        oriented=True,
                        starfusion=dict(junction_reads=int(r.get("JunctionReadCount") or 0),
                                        spanning_frags=int(r.get("SpanningFragCount") or 0),
                                        ffpm=float(r.get("FFPM") or 0), splice_type=r.get("SpliceType"),
                                        large_anchor_support=r.get("LargeAnchorSupport"),
                                        annots=r.get("annots"),
                                        left_breakpoint=r["LeftBreakpoint"],
                                        right_breakpoint=r["RightBreakpoint"])))
    return out


def _same(a, b, tol):
    return (a["chr1"] == b["chr1"] and a["chr2"] == b["chr2"] and a["keep1"] == b["keep1"]
            and a["keep2"] == b["keep2"] and abs(a["pos1"] - b["pos1"]) <= tol
            and abs(a["pos2"] - b["pos2"]) <= tol)


def _flip(e):
    """Same junction described from the other partner's side (unoriented Arriba calls)."""
    inv = {"le": "ge", "ge": "le"}
    f = dict(e)
    f.update(gene1=e["gene2"], gene2=e["gene1"], names1=e["names2"], names2=e["names1"],
             intergenic1=e["intergenic2"], intergenic2=e["intergenic1"],
             chr1=e["chr2"], pos1=e["pos2"], keep1=inv[e["keep2"]],
             chr2=e["chr1"], pos2=e["pos1"], keep2=inv[e["keep1"]])
    return f


def merge_events(arriba, starfusion, tol=10):
    events = [dict(e, callers=["arriba"]) for e in arriba]
    for s in starfusion:
        hit = None
        for e in events:
            if "starfusion" in e:
                continue
            if _same(e, s, tol):
                hit = e
                break
            if not e["oriented"] and _same(_flip(e), s, tol):
                # STAR-Fusion knows the orientation: adopt it
                e.update({k: v for k, v in _flip(e).items() if k not in ("arriba", "callers")})
                e["oriented"] = True
                hit = e
                break
        if hit:
            hit["starfusion"] = s["starfusion"]
            hit["callers"].append("starfusion")
            for i in ("1", "2"):
                sg = s["gene" + i]
                hit["names" + i] = [sg] + [n for n in hit["names" + i] if n != sg]
                if hit["intergenic" + i]:
                    # Arriba only knows the neighbouring genes; STAR-Fusion's annotated gene is the better label
                    hit["gene" + i] = sg
        else:
            events.append(dict(s, callers=["starfusion"]))
    return events


def priority(e):
    a, s = e.get("arriba") or {}, e.get("starfusion") or {}
    support = (a.get("split_reads1", 0) + a.get("split_reads2", 0) + a.get("discordant_mates", 0)
               + s.get("junction_reads", 0) + s.get("spanning_frags", 0))
    return (len(e["callers"]), CONFIDENCE_RANK.get(a.get("confidence"), 0), support)


# ----------------------------------------------------------------------------- projection

def build_segments(ev, genes, cap, pad, off):
    segs = []
    for i, (chrom, pos, keep, gname) in enumerate([(ev["chr1"], ev["pos1"], ev["keep1"], ev["gene1"]),
                                                    (ev["chr2"], ev["pos2"], ev["keep2"], ev["gene2"])]):
        g = genes.get(gname)
        on_gene = g is not None and g["chrom"] == chrom and g["start"] - pad <= pos <= g["end"] + pad
        if keep == "le":
            far = (g["start"] - pad) if on_gene else pos - 10000
            start, end = max(1, far, pos - cap + 1), pos
        else:
            far = (g["end"] + pad) if on_gene else pos + 10000
            start, end = pos, min(far, pos + cap - 1)
        # 5' segment must END at the junction on the contig; 3' segment must START there.
        rc = (keep == "ge") if i == 0 else (keep == "le")
        s = Segment(chrom, start, end, rc, off)
        off += s.len
        segs.append(s)
    return segs


def _expand(rec):
    pos = rec.reference_start + 1
    for op, ln in rec.cigartuples:
        if op in REF_CONSUME:
            yield op, ln, pos
            pos += ln
        else:
            yield op, ln, None


def project(rec, seg):
    """Project one aligned record into a segment; ref blocks outside it become soft clips."""
    if rec.reference_name != seg.chrom:
        return None
    ops = []
    for op, ln, p in _expand(rec):
        if op in (4, 5):
            ops.append([op, ln, None, False])
            continue
        if p is None:  # insertion
            ops.append([op, ln, None, None])
            continue
        a, b = p, p + ln - 1
        cuts = sorted({a, b + 1} | {x for x in (seg.start, seg.end + 1) if a < x <= b})
        for x, y in zip(cuts, cuts[1:]):
            ops.append([op, y - x, x, seg.start <= x <= seg.end])
    for i, o in enumerate(ops):
        if o[3] is None:
            prev = next((ops[j][3] for j in range(i - 1, -1, -1)
                         if ops[j][3] is not None and ops[j][0] not in (4, 5)), False)
            nxt = next((ops[j][3] for j in range(i + 1, len(ops))
                        if ops[j][3] is not None and ops[j][0] not in (4, 5)), False)
            o[3] = prev and nxt
    inside = [i for i, o in enumerate(ops) if o[3] and o[0] in (0, 7, 8)]
    if not inside:
        return None
    first, last = inside[0], inside[-1]
    lead = sum(o[1] for o in ops[:first] if o[0] in QRY_CONSUME or o[0] == 5)
    trail = sum(o[1] for o in ops[last + 1:] if o[0] in QRY_CONSUME or o[0] == 5)
    hlead = sum(o[1] for o in ops[:first] if o[0] == 5)
    htrail = sum(o[1] for o in ops[last + 1:] if o[0] == 5)
    core = [(o[0], o[1]) for o in ops[first:last + 1]]
    clipped_ref = any(o[0] in (0, 2, 3, 7, 8) and not o[3] for o in ops if o[2] is not None)
    gstart = ops[first][2]
    gend = gstart + sum(l for op, l in core if op in REF_CONSUME) - 1
    if seg.rc:
        core = core[::-1]
        lead, trail = trail, lead
        hlead, htrail = htrail, hlead
        cstart = seg.to_contig(gend)
        rev = not rec.is_reverse
    else:
        cstart = seg.to_contig(gstart)
        rev = rec.is_reverse
    return dict(cstart=cstart, core=core, lead=lead, trail=trail, hlead=hlead, htrail=htrail, rev=rev,
                cend=cstart + sum(l for op, l in core if op in REF_CONSUME), clipped=clipped_ref)


def _full_seq(recs):
    for r in recs:
        if r.query_sequence and not any(op == 5 for op, _ in r.cigartuples):
            seq, qual = r.query_sequence, r.query_qualities
            if r.is_reverse:
                seq, qual = revcomp(seq), (qual[::-1] if qual is not None else None)
            return seq, qual
    return None, None


def _merge_ops(cig):
    out = []
    for o, l in cig:
        if o == 5 or l <= 0:
            continue
        if out and out[-1][0] == o:
            out[-1] = (o, out[-1][1] + l)
        else:
            out.append((o, l))
    return out


def classify_mate(rl, segs, junction, min_anchor):
    """Return (chosen, cigar, seq, qual, tag) for one mate's records, or None."""
    sq, qu = _full_seq(rl)
    proj = []
    for r in rl:
        for si, s in enumerate(segs):
            p = project(r, s)
            if p:
                p["seg"], p["rec"] = si, r
                proj.append(p)
    if not proj:
        return None
    partial = sq is None  # every record of this mate is hard-clipped
    tag, chosen = "other", None
    if len({p["seg"] for p in proj}) == 2 and not partial:
        A = max((p for p in proj if p["seg"] == 0), key=lambda p: p["cend"])
        B = min((p for p in proj if p["seg"] == 1), key=lambda p: p["cstart"])
        if A["rev"] == B["rev"] and A["rec"] is not B["rec"]:
            L = len(sq)
            qa0, qa1 = A["lead"], L - A["trail"]
            qb0, qb1 = B["lead"], L - B["trail"]
            if qa0 <= qb0:
                coreB, bstart, ov = list(B["core"]), B["cstart"], qa1 - qb0
                while ov > 0 and coreB:  # trim query overlap (microhomology) from B's left
                    op, ln = coreB[0]
                    if op in (0, 7, 8):
                        k = min(ln, ov)
                        coreB[0] = (op, ln - k)
                        bstart += k
                        ov -= k
                        qb0 += k
                    elif op == 1:
                        k = min(ln, ov)
                        coreB[0] = (op, ln - k)
                        ov -= k
                        qb0 += k
                    elif op in (2, 3):
                        bstart += ln
                        coreB[0] = (op, 0)
                    if coreB and coreB[0][1] == 0:
                        coreB.pop(0)
                gap = bstart - A["cend"]
                ins = qb0 - qa1
                if gap >= 0 and coreB and coreB[0][0] in (0, 7, 8):
                    mid = ([(1, ins)] if ins > 0 else []) + ([(3 if gap > 30 else 2, gap)] if gap > 0 else [])
                    cig = ([(4, qa0)] if qa0 else []) + A["core"] + mid + coreB + ([(4, L - qb1)] if L - qb1 else [])
                    anchor_a = junction - A["cstart"]
                    anchor_b = bstart + sum(l for o, l in coreB if o in REF_CONSUME) - junction
                    tag = "junction" if (anchor_a >= min_anchor and anchor_b >= min_anchor and gap <= 5) else "junction_weak"
                    chosen = dict(cstart=A["cstart"], cig=cig, rev=A["rev"], src=A["rec"], dist=gap,
                                  mapq=min(A["rec"].mapping_quality, B["rec"].mapping_quality))
    if chosen is None:
        prim = [p for p in proj if not p["rec"].is_supplementary] or proj
        p = max(prim, key=lambda p: sum(l for o, l in p["core"] if o in (0, 7, 8)))
        lead, trail = p["lead"], p["trail"]
        if partial:  # only hard-clipped records: use this record's stored (clipped) sequence
            r0 = p["rec"]
            sq = revcomp(r0.query_sequence) if r0.is_reverse else r0.query_sequence
            qu = r0.query_qualities[::-1] if (r0.is_reverse and r0.query_qualities is not None) else r0.query_qualities
            lead, trail = lead - p["hlead"], trail - p["htrail"]
        cig = ([(4, lead)] if lead else []) + p["core"] + ([(4, trail)] if trail else [])
        near = min(abs(p["cstart"] - junction), abs(p["cend"] - junction))
        if p["clipped"] and near <= 3:
            tag = "wt5" if p["seg"] == 0 else "wt3"
        elif p["rec"].has_tag("SA") and near <= 3:
            tag = "chimeric_other"
        chosen = dict(cstart=p["cstart"], cig=cig, rev=p["rev"], src=p["rec"], dist=near,
                      mapq=p["rec"].mapping_quality)
    s_out, q_out = (revcomp(sq), qu[::-1] if qu is not None else None) if chosen["rev"] else (sq, qu)
    cig = _merge_ops(chosen["cig"])
    if sum(l for o, l in cig if o in QRY_CONSUME) != len(s_out):
        return "len_mismatch"
    return chosen, cig, s_out, q_out, tag


def project_event(bam, ev, segs, tid, opts, header):
    """Project all reads overlapping the event's segments. Returns (records, counts)."""
    junction = segs[0].off + segs[0].len  # contig 0-based index of the first 3' base
    by_q = collections.defaultdict(list)
    for s, bp in zip(segs, (ev["pos1"], ev["pos2"])):
        for r in bam.fetch(s.chrom, s.start - 1, s.end):
            if r.is_unmapped or r.is_secondary:
                continue
            # duplicates only matter for junction/chimeric evidence near the breakpoint
            if r.is_duplicate and not (r.reference_start - opts.dup_window <= bp <= r.reference_end + opts.dup_window):
                continue
            by_q[r.query_name].append(r)

    frags = {}
    skipped = 0
    for qn, recs in by_q.items():
        uniq = {}
        for r in recs:
            uniq[(r.is_read1, r.is_supplementary, r.reference_name, r.reference_start, r.cigarstring)] = r
        recs = list(uniq.values())
        mates = {True: [r for r in recs if r.is_read1 or not r.is_paired], False: [r for r in recs if r.is_read2]}
        final = {}
        for is_r1, rl in mates.items():
            if not rl:
                continue
            res = classify_mate(rl, segs, junction, opts.min_anchor)
            if res == "len_mismatch":
                skipped += 1
            elif res:
                final[is_r1] = list(res)
        if not final:
            continue
        tags = {v[4] for v in final.values()}
        if len(final) == 2 and "junction" not in tags:
            sides = {0 if v[0]["cstart"] < junction else 1 for v in final.values()}
            if sides == {0, 1} and all(v[4] in ("other", "chimeric_other") for v in final.values()):
                for v in final.values():
                    v[4] = "spanning"
        frags[qn] = final

    # fragment class = best class among its mates
    def frag_class(final):
        return min((v[4] for v in final.values()), key=lambda t: SORT_ORDER[t])

    counts = dict(nondup=collections.Counter(), dup=collections.Counter())
    keep = set()
    shown = collections.Counter()
    for qn in sorted(frags, key=lambda q: hashlib.md5(q.encode()).hexdigest()):
        final = frags[qn]
        cls = frag_class(final)
        dup = any(v[0]["src"].is_duplicate for v in final.values())
        counts["dup" if dup else "nondup"][cls] += 1
        if cls in EVIDENCE:
            if not dup or opts.keep_dups in ("evidence", "all"):
                keep.add(qn)
        elif not dup or opts.keep_dups == "all":
            shown[cls] += 1
            if not opts.max_context or shown[cls] <= opts.max_context:
                keep.add(qn)

    out = []
    for qn in keep:
        final = frags[qn]
        recs = {}
        for is_r1, (ch, cig, s_out, q_out, tag) in final.items():
            src = ch["src"]
            r = pysam.AlignedSegment(header)
            r.query_name = qn
            r.query_sequence = s_out
            if q_out is not None:
                r.query_qualities = q_out
            r.flag = 0
            r.is_paired = src.is_paired
            r.is_read1, r.is_read2 = src.is_read1, src.is_read2
            r.is_duplicate = src.is_duplicate
            r.is_qcfail = src.is_qcfail
            r.is_reverse = ch["rev"]
            r.reference_id = tid
            r.reference_start = ch["cstart"]
            r.mapping_quality = ch["mapq"]
            r.cigartuples = cig
            r.set_tag("XF", tag)
            r.set_tag("XB", int(ch["dist"]))
            r.set_tag("YC", COLOURS[tag])
            r.set_tag("XK", SORT_ORDER[tag] + tag)
            for t in ("RX", "MI", "NH"):
                if src.has_tag(t):
                    r.set_tag(t, src.get_tag(t))
            recs[is_r1] = r
        for k, r in recs.items():
            if not r.is_paired:
                continue
            m = recs.get(not k)
            r.next_reference_id = tid
            if m is None:
                r.mate_is_unmapped = True
                r.next_reference_start = r.reference_start
            else:
                r.next_reference_start = m.reference_start
                r.mate_is_reverse = m.is_reverse
                r.is_proper_pair = r.is_reverse != m.is_reverse
                lo, hi = min(r.reference_start, m.reference_start), max(r.reference_end, m.reference_end)
                r.template_length = (hi - lo) * (1 if r.reference_start <= m.reference_start else -1)
        out.extend(recs.values())
    return out, counts, skipped


def parse_annots(s):
    """FusionAnnotator string -> list. Same rules as gennext-api-edge parse_annotations:
    '["A:Oncogene","B"];["Cosmic"];INTRACHROMOSOMAL[chr4:0.05Mb]' -> [A:Oncogene, B, Cosmic, INTRA...]."""
    out = []
    for chunk in (s or "").split(";"):
        chunk = chunk.strip()
        if not chunk or chunk == ".":
            continue
        if chunk.startswith("[") and chunk.endswith("]"):
            out += [i.strip().strip("\"'") for i in chunk[1:-1].split(",") if i.strip().strip("\"'")]
        else:
            out.append(chunk)
    return out


def splice_type(ev):
    """STAR-Fusion's SpliceType when available; otherwise derived from Arriba's sites."""
    sf, a = ev.get("starfusion") or {}, ev.get("arriba") or {}
    if sf.get("splice_type"):
        return sf["splice_type"]
    if a.get("site1") or a.get("site2"):
        both = a.get("site1") == "splice-site" and a.get("site2") == "splice-site"
        return "ONLY_REF_SPLICE" if both else "INCL_NON_REF_SPLICE"
    return None


def gene_locus(g, chrom):
    if not g or g["chrom"] != chrom:
        return None
    return f"{g['chrom']}:{g['start']}-{g['end']}:{g['strand']}"


def fi_fields(ev, segs, nd, dup, genes):
    """The fields the IGV fusion view used to read from FusionInspector's abridged TSV."""
    j = segs[0].len
    s5 = "+" if ev["keep1"] == "le" else "-"
    s3 = "+" if ev["keep2"] == "ge" else "-"
    annots = parse_annots((ev.get("starfusion") or {}).get("annots"))
    for t in ((ev.get("arriba") or {}).get("tags") or "").split(","):
        t = t.strip()
        if t and t != "." and t not in annots:
            annots.append(t)
    return dict(
        fusion_name=f"{ev['gene1']}--{ev['gene2']}",
        left_gene=ev["gene1"], right_gene=ev["gene2"],
        left_local_breakpoint=j, right_local_breakpoint=j + 1,
        left_breakpoint=f"{ev['chr1']}:{ev['pos1']}:{s5}",
        right_breakpoint=f"{ev['chr2']}:{ev['pos2']}:{s3}",
        junction_read_count=nd.get("junction", 0),
        spanning_frag_count=nd.get("spanning", 0),
        junction_read_count_with_dups=nd.get("junction", 0) + dup.get("junction", 0),
        spanning_frag_count_with_dups=nd.get("spanning", 0) + dup.get("spanning", 0),
        splice_type=splice_type(ev), annots=annots,
        left_gene_locus=gene_locus(genes.get(ev["gene1"]), ev["chr1"]),
        right_gene_locus=gene_locus(genes.get(ev["gene2"]), ev["chr2"]))


def exon_gtf(name, segs, ev, genes):
    """Exons of both partners in contig coordinates, one record per transcript exon."""
    lines = []
    for s, gname in zip(segs, (ev["gene1"], ev["gene2"])):
        g = genes.get(gname)
        if not g or g["chrom"] != s.chrom:
            continue
        strand = g["strand"] if not s.rc else {"+": "-", "-": "+"}.get(g["strand"], ".")
        for es, ee, tid in g["exons"]:
            if ee < s.start or es > s.end:
                continue
            c1, c2 = sorted((s.to_contig(max(es, s.start)), s.to_contig(min(ee, s.end))))
            attrs = (f'gene_id "{gname}"; transcript_id "{tid or gname}"; gene_name "{gname}"; '
                     f'orig_coord_info "{g["chrom"]},{es},{ee},{g["strand"]}";')
            lines.append((c1, c2, f"{name}\tbpview\texon\t{c1 + 1}\t{c2 + 1}\t.\t{strand}\t.\t{attrs}\n"))
    return [l for *_, l in sorted(lines)]


def bgzip_index(path, preset):
    """Write path.gz + path.gz.tbi next to path (the plain file is kept)."""
    pysam.tabix_index(path, preset=preset, force=True, keep_original=True)


def exon_bed(name, segs, ev, genes):
    lines = []
    for s, gname in zip(segs, (ev["gene1"], ev["gene2"])):
        g = genes.get(gname)
        if not g or g["chrom"] != s.chrom:
            continue
        seen = set()
        for es, ee, _tid in g["exons"]:
            if ee < s.start or es > s.end:
                continue
            c1, c2 = sorted((s.to_contig(max(es, s.start)), s.to_contig(min(ee, s.end))))
            if (c1, c2) in seen:
                continue
            seen.add((c1, c2))
            strand = g["strand"] if not s.rc else {"+": "-", "-": "+"}.get(g["strand"], ".")
            lines.append((c1, f"{name}\t{c1}\t{c2 + 1}\t{gname}\t0\t{strand}\n"))
    return [l for _, l in sorted(lines)]


# ----------------------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--bam", required=True, help="STAR genome BAM with duplicate flags (indexed)")
    ap.add_argument("--fasta", required=True, help="genome FASTA the BAM was aligned to (indexed)")
    ap.add_argument("--gtf", required=True, help="annotation for partner extents and exon track")
    ap.add_argument("--arriba-tsv")
    ap.add_argument("--starfusion-tsv")
    ap.add_argument("--fusion-list", help="optional ranked allow-list of GENE1--GENE2 pairs (first column); "
                                          "list order decides which events survive --max-events")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--max-events", type=int, default=300)
    ap.add_argument("--merge-tolerance", type=int, default=10)
    ap.add_argument("--cap", type=int, default=60000, help="max bp per partner segment")
    ap.add_argument("--pad", type=int, default=500)
    ap.add_argument("--min-anchor", type=int, default=10)
    ap.add_argument("--max-context", type=int, default=40,
                    help="max non-evidence fragments drawn per class per event (0 = all)")
    ap.add_argument("--keep-dups", choices=("none", "evidence", "all"), default="evidence")
    ap.add_argument("--dup-window", type=int, default=500)
    ap.add_argument("--version", action="version", version=__version__)
    opts = ap.parse_args(argv)

    arriba = load_arriba(opts.arriba_tsv) if opts.arriba_tsv and os.path.getsize(opts.arriba_tsv) else []
    starf = load_starfusion(opts.starfusion_tsv) if opts.starfusion_tsv and os.path.getsize(opts.starfusion_tsv) else []
    events = merge_events(arriba, starf, opts.merge_tolerance)
    n_called = len(events)
    if opts.fusion_list:
        # The list is ranked strongest-first upstream (filter-artifact-fusions: database-recorded
        # fusions, then in-frame, then concordant). Its order decides which events survive the cap;
        # bpview's own priority only orders events within one gene pair.
        rank = {}
        for l in open(opts.fusion_list):
            if l.strip() and not l.startswith("#"):
                rank.setdefault(l.strip().split("\t")[0], len(rank))

        def list_rank(e):
            r = [rank[k] for a in e["names1"] for b in e["names2"]
                 for k in (f"{a}--{b}", f"{b}--{a}") if k in rank]
            return min(r) if r else None

        for e in events:
            e["list_rank"] = list_rank(e)
        events = [e for e in events if e["list_rank"] is not None]
        events.sort(key=priority, reverse=True)
        events.sort(key=lambda e: e["list_rank"])  # stable: priority breaks ties within a pair
    else:
        events.sort(key=priority, reverse=True)
    n_allowed = len(events)
    events = events[:opts.max_events]

    genes = read_gtf(opts.gtf)
    fa = pysam.FastaFile(opts.fasta)
    bam = pysam.AlignmentFile(opts.bam)
    fa_contigs = set(fa.references)

    contigs, used = [], set()
    for i, ev in enumerate(events):
        if ev["chr1"] not in fa_contigs or ev["chr2"] not in fa_contigs:
            continue
        segs = build_segments(ev, genes, opts.cap, opts.pad, 0)
        seq = "".join(revcomp(fa.fetch(s.chrom, s.start - 1, s.end).upper()) if s.rc
                      else fa.fetch(s.chrom, s.start - 1, s.end).upper() for s in segs)
        name = safe(f"{ev['gene1']}--{ev['gene2']}__{ev['chr1']}_{ev['pos1']}__{ev['chr2']}_{ev['pos2']}")
        if name in used:
            continue
        used.add(name)
        contigs.append(dict(id=f"E{len(contigs) + 1:03d}", name=name, ev=ev, segs=segs, seq=seq))

    with open(opts.out_prefix + ".fa", "w") as fh:
        for c in contigs:
            fh.write(f">{c['name']}\n")
            for j in range(0, len(c["seq"]), 60):
                fh.write(c["seq"][j:j + 60] + "\n")
    if contigs:
        pysam.faidx(opts.out_prefix + ".fa")
    else:
        open(opts.out_prefix + ".fa.fai", "w").close()

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c["name"], "LN": len(c["seq"])} for c in contigs] or [{"SN": "empty", "LN": 1}],
        "PG": [{"ID": "bpview", "PN": "bpview", "VN": __version__, "CL": " ".join(sys.argv)}]})
    tmp = opts.out_prefix + ".unsorted.bam"
    summaries = []
    with pysam.AlignmentFile(tmp, "wb", header=header) as ob, open(opts.out_prefix + ".bed", "w") as bed, \
            open(opts.out_prefix + ".gtf", "w") as gtf:
        for tid, c in enumerate(contigs):
            recs, counts, skipped = project_event(bam, c["ev"], c["segs"], tid, opts, header)
            for r in recs:
                ob.write(r)
            for line in exon_bed(c["name"], c["segs"], c["ev"], genes):
                bed.write(line)
            for line in exon_gtf(c["name"], c["segs"], c["ev"], genes):
                gtf.write(line)
            nd = counts["nondup"]
            j = nd["junction"] + nd["spanning"]
            ev = c["ev"]
            summaries.append(dict(
                id=c["id"], contig=c["name"], contig_length=len(c["seq"]),
                junction_1based=c["segs"][0].len,
                gene1=ev["gene1"], gene2=ev["gene2"],
                breakpoint1=f"{ev['chr1']}:{ev['pos1']}", breakpoint2=f"{ev['chr2']}:{ev['pos2']}",
                keep1=ev["keep1"], keep2=ev["keep2"], orientation_known=ev["oriented"],
                gene_names1=ev["names1"], gene_names2=ev["names2"], fusion_list_rank=ev.get("list_rank"),
                callers=ev["callers"], arriba=ev.get("arriba"), starfusion=ev.get("starfusion"),
                segments=[s.as_dict() for s in c["segs"]],
                fragments_nondup=dict(nd), fragments_dup=dict(counts["dup"]),
                fusion_fraction_5p=round(j / (j + nd["wt5"]), 3) if j + nd["wt5"] else None,
                fusion_fraction_3p=round(j / (j + nd["wt3"]), 3) if j + nd["wt3"] else None,
                records_drawn=len(recs), skipped_length_mismatch=skipped,
                **fi_fields(ev, c["segs"], nd, counts["dup"], genes)))
    pysam.sort("-o", opts.out_prefix + ".bam", tmp)
    os.remove(tmp)
    pysam.index(opts.out_prefix + ".bam")
    bgzip_index(opts.out_prefix + ".bed", "bed")
    bgzip_index(opts.out_prefix + ".gtf", "gff")

    doc = dict(bpview_version=__version__, sample_id=opts.sample_id,
               inputs=dict(bam=opts.bam, arriba_tsv=opts.arriba_tsv, starfusion_tsv=opts.starfusion_tsv,
                           fusion_list=opts.fusion_list, gtf=opts.gtf),
               settings=dict(max_events=opts.max_events, merge_tolerance=opts.merge_tolerance, cap=opts.cap,
                             pad=opts.pad, min_anchor=opts.min_anchor, max_context=opts.max_context,
                             keep_dups=opts.keep_dups),
               n_events_called=n_called, n_events_allowed=n_allowed, n_events=len(summaries),
               events=summaries)
    with open(opts.out_prefix + ".events.json", "w") as fh:
        json.dump(doc, fh, indent=1)
    print(f"bpview {__version__}: {n_called} merged events, {n_allowed} allowed, {len(summaries)} written")


if __name__ == "__main__":
    main()
