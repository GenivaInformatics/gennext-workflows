#!/usr/bin/env python3
"""fusion-events: one row per fusion breakpoint event, for the per-sample fusions DuckDB.

  fusion-events build  --sample-id S --arriba-tsv A --starfusion-tsv F --oncokb-json J
                       [--bpview-events E] [--pon P] --out-tsv OUT
  fusion-events pon    --events-tsv E1 E2 ... --out-tsv PON

`build` merges Arriba and STAR-Fusion calls by breakpoint with bpview's merge (same loaders,
same +-10 bp tolerance, same gene aliases), so every event lines up with its bpview contig.
Gene-pair annotation (OncoKB, fusion-report databases, the existing pair tier) comes from
oncokb-fusion-annotation's fusions_with_oncokb.json: no second OncoKB query. Support comes
from bpview's projected, non-duplicate fragments when bpview drew the event, otherwise from
the callers' own counts (support_source says which).

`pon` turns event tables of many samples into a recurrence panel: one row per sample and
junction. `build --pon` then counts, for each event, how many *other* samples carry the same
gene pair (either orientation) and the same junction (+-tolerance).

Tiers (see tier_event; every threshold is a flag):
  1  actionable: OncoKB LEVEL_1/LEVEL_2/LEVEL_R1 (or Dx1) for the pair in this orientation,
     with fusion evidence: not out-of-frame/stop-codon, at least --tier1-min-unique-molecules,
     and two callers or Arriba confidence high/medium
  2  any OncoKB hit (actionable, Oncogenic / Likely Oncogenic, 3A/3B/4/R2/Dx2/Dx3) with frame
     and caller evidence but not tier 1 (e.g. too few molecules),
     or a curated known fusion (Mitelman, COSMIC, CIViC) whatever its evidence
  3  everything else. OncoKB levels are often partner-agnostic (any NTRK fusion), so an OncoKB
     hit without fusion evidence and without a curated pair stays here; tier_reason says so.
A gene pair recurrent in the PoN (>= --pon-min-samples and >= --pon-min-fraction of the other
samples) is capped at tier 2 (if it would be 1) or 3. tier_reason lists every rule that applied.
"""
import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bpview  # noqa: E402  (same image, same directory)

__version__ = "0.3.0"

ACTIONABLE = {"LEVEL_1", "LEVEL_2", "LEVEL_R1"}
ACTIONABLE_DX = {"LEVEL_Dx1"}
ONCOGENIC_LEVELS = {"LEVEL_3A", "LEVEL_3B", "LEVEL_4", "LEVEL_R2"}
ONCOGENIC_DX = {"LEVEL_Dx2", "LEVEL_Dx3"}
ONCOGENIC = {"Oncogenic", "Likely Oncogenic"}
CURATED_DB = {"Mitelman", "COSMIC", "CIViC"}

# Column order and DuckDB types of the output; duckdb-fusion-events reads it with these types.
# Every field a user may filter on is its own typed column. List columns (VARCHAR[]) are written
# '|'-separated in the TSV and split on load; an empty cell is an empty list.
TIER_FLAGS = ("out_of_frame", "weak_caller_evidence", "low_unique_molecules", "pon_recurrent", "pon_capped")
COLUMNS = [
    ("sample_id", "VARCHAR"), ("event_id", "VARCHAR"), ("fusion", "VARCHAR"),
    ("gene5", "VARCHAR"), ("gene3", "VARCHAR"), ("gene5_aliases", "VARCHAR[]"), ("gene3_aliases", "VARCHAR[]"),
    ("breakpoint5", "VARCHAR"), ("breakpoint3", "VARCHAR"), ("orientation_known", "BOOLEAN"),
    ("callers", "VARCHAR[]"), ("n_callers", "INTEGER"), ("caller_arriba", "BOOLEAN"), ("caller_starfusion", "BOOLEAN"),
    ("arriba_status", "VARCHAR"), ("arriba_confidence", "VARCHAR"), ("arriba_type", "VARCHAR"),
    ("arriba_split_reads", "INTEGER"), ("arriba_discordant_mates", "INTEGER"),
    ("starfusion_junction_reads", "INTEGER"), ("starfusion_spanning_frags", "INTEGER"),
    ("starfusion_ffpm", "DOUBLE"),
    ("reading_frame", "VARCHAR"), ("splice_type", "VARCHAR"),
    ("unique_molecules", "INTEGER"), ("pcr_copies", "INTEGER"), ("support_source", "VARCHAR"),
    ("known_db", "VARCHAR[]"), ("curated_known", "BOOLEAN"),
    ("oncokb_oncogenic", "VARCHAR"), ("oncokb_mutation_effect", "VARCHAR"),
    ("oncokb_highest_level", "VARCHAR"), ("oncokb_highest_dx_level", "VARCHAR"),
    ("pair_clinical_tier", "VARCHAR"), ("pair_fii", "DOUBLE"),
    ("pon_gene_pair_count", "INTEGER"), ("pon_junction_count", "INTEGER"), ("pon_size", "INTEGER"),
    ("pon_recurrent", "BOOLEAN"),
    # tier_basis: the rule that set the tier (oncokb_actionable | oncokb_oncogenic | curated_known | none);
    # tier_flags: what held it back (TIER_FLAGS); tier_reason: the same in words.
    ("tier", "INTEGER"), ("tier_basis", "VARCHAR"), ("tier_flags", "VARCHAR[]"), ("tier_reason", "VARCHAR"),
    ("bpview_contig", "VARCHAR"), ("builder_version", "VARCHAR"),
]
PON_COLUMNS = ["sample_id", "gene5", "gene3", "chr5", "pos5", "keep5", "chr3", "pos3", "keep3"]


# ----------------------------------------------------------------------------- inputs

def load_pairs(path):
    """fusions_with_oncokb.json -> {'G1--G2': entry}. Missing or empty file -> {}."""
    if not path or not os.path.exists(path) or not os.path.getsize(path):
        return {}
    with open(path) as fh:
        data = json.load(fh)
    return {e["Fusion"]: e for e in data if e.get("Fusion")}


def event_key(chr1, pos1, keep1, chr2, pos2, keep2):
    return (chr1, int(pos1), keep1, chr2, int(pos2), keep2)


def load_bpview(path):
    """bpview events.json -> {junction key: event summary}."""
    if not path or not os.path.exists(path) or not os.path.getsize(path):
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
    if not path or not os.path.exists(path) or not os.path.getsize(path):
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


# ----------------------------------------------------------------------------- per event

def strand5(ev):
    return "+" if ev["keep1"] == "le" else "-"


def strand3(ev):
    return "+" if ev["keep2"] == "ge" else "-"


def pair_entry(ev, pairs):
    """The fusion-report/OncoKB entry for this event's gene pair, in the event's orientation.
    Unoriented Arriba calls also try the reverse orientation."""
    for a in ev["names1"]:
        for b in ev["names2"]:
            if f"{a}--{b}" in pairs:
                return f"{a}--{b}", pairs[f"{a}--{b}"]
    if not ev["oriented"]:
        for a in ev["names1"]:
            for b in ev["names2"]:
                if f"{b}--{a}" in pairs:
                    return f"{b}--{a}", pairs[f"{b}--{a}"]
    return None, None


def support(ev, bp):
    a, s = ev.get("arriba") or {}, ev.get("starfusion") or {}
    if bp:
        nd, dup = bp.get("fragments_nondup") or {}, bp.get("fragments_dup") or {}
        u = nd.get("junction", 0) + nd.get("spanning", 0)
        return u, dup.get("junction", 0) + dup.get("spanning", 0), "bpview"
    caller = max(a.get("split_reads1", 0) + a.get("split_reads2", 0) + a.get("discordant_mates", 0),
                 s.get("junction_reads", 0) + s.get("spanning_frags", 0))
    return caller, None, "caller"


def known_dbs(ev, entry):
    dbs = list((entry or {}).get("Databases") or [])
    for t in ((ev.get("arriba") or {}).get("tags") or "").split(","):
        name = t.split("=")[0].strip()
        if name == "Mitelman" and name not in dbs:
            dbs.append(name)
    return dbs


def pon_counts(ev, pon, sample_id, tol):
    if pon is None:
        return None, None, None
    others = pon["samples"] - {sample_id}
    pair = set()
    for a in ev["names1"]:
        for b in ev["names2"]:
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
    if row["reading_frame"] in ("out-of-frame", "stop-codon") and not opts.allow_out_of_frame:
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
        reasons.append(lab + (" (not tier 1: " + "; ".join(low.values()) + ")" if actionable and low else ""))
    elif curated:
        tier, basis = 2, "curated_known"
        reasons.append("known fusion: " + ", ".join(curated))
        if oncogenic:
            flags += list(basic)
            reasons.append(lab + ", but " + "; ".join(basic.values()))
    else:
        tier, basis = 3, "none"
        if oncogenic:
            flags += list(basic)
            reasons.append(lab + ", but " + "; ".join(basic.values()))
        else:
            reasons.append("not known or actionable")
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


def build_rows(opts):
    arriba = bpview.load_arriba(opts.arriba_tsv) if opts.arriba_tsv and os.path.getsize(opts.arriba_tsv) else []
    starf = (bpview.load_starfusion(opts.starfusion_tsv)
             if opts.starfusion_tsv and os.path.getsize(opts.starfusion_tsv) else [])
    events = bpview.merge_events(arriba, starf, opts.merge_tolerance)
    pairs, bp, pon = load_pairs(opts.oncokb_json), load_bpview(opts.bpview_events), load_pon(opts.pon)
    rows = []
    for ev in events:
        a, s = ev.get("arriba") or {}, ev.get("starfusion") or {}
        key = event_key(ev["chr1"], ev["pos1"], ev["keep1"], ev["chr2"], ev["pos2"], ev["keep2"])
        b = bp.get(key)
        fusion, entry = pair_entry(ev, pairs)
        onc = (entry or {}).get("OncoKB") or {}
        u, dup, src = support(ev, b)
        npair, njunc, size = pon_counts(ev, pon, opts.sample_id, opts.merge_tolerance)
        row = dict(
            sample_id=opts.sample_id, event_id=None,
            fusion=fusion or f"{ev['gene1']}--{ev['gene2']}",
            gene5=ev["gene1"], gene3=ev["gene2"],
            gene5_aliases=list(ev["names1"]), gene3_aliases=list(ev["names2"]),
            breakpoint5=f"{ev['chr1']}:{ev['pos1']}:{strand5(ev)}",
            breakpoint3=f"{ev['chr2']}:{ev['pos2']}:{strand3(ev)}",
            orientation_known=bool(ev["oriented"]),
            callers=list(ev["callers"]), n_callers=len(ev["callers"]),
            caller_arriba="arriba" in ev["callers"], caller_starfusion="starfusion" in ev["callers"],
            arriba_status="kept" if a else None,
            arriba_confidence=a.get("confidence"), arriba_type=a.get("type"),
            arriba_split_reads=(a.get("split_reads1", 0) + a.get("split_reads2", 0)) if a else None,
            arriba_discordant_mates=a.get("discordant_mates") if a else None,
            starfusion_junction_reads=s.get("junction_reads") if s else None,
            starfusion_spanning_frags=s.get("spanning_frags") if s else None,
            starfusion_ffpm=s.get("ffpm") if s else None,
            reading_frame=(a.get("reading_frame") if a.get("reading_frame") not in (None, ".") else None),
            splice_type=bpview.splice_type(ev),
            unique_molecules=u, pcr_copies=dup, support_source=src,
            known_db=known_dbs(ev, entry),
            oncokb_oncogenic=onc.get("ONCOGENIC") or None,
            oncokb_mutation_effect=onc.get("MUTATION_EFFECT") or None,
            oncokb_highest_level=onc.get("HIGHEST_LEVEL") or None,
            oncokb_highest_dx_level=onc.get("HIGHEST_DX_LEVEL") or None,
            pair_clinical_tier=(entry or {}).get("clinical_tier"),
            pair_fii=(entry or {}).get("Fusion Indication Index (FII)"),
            pon_gene_pair_count=npair, pon_junction_count=njunc, pon_size=size,
            bpview_contig=b["contig"] if b else None, builder_version=__version__)
        row["curated_known"] = bool(CURATED_DB & set(row["known_db"]))
        row.update(tier_event(row, opts))
        rows.append((row, bpview.priority(ev)))
    rows.sort(key=lambda x: x[1], reverse=True)
    rows.sort(key=lambda x: x[0]["tier"])  # stable: caller priority orders events within a tier
    out = []
    for i, (row, _) in enumerate(rows):
        row["event_id"] = f"{opts.sample_id}:E{i + 1:04d}"
        out.append(row)
    return out


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
    b = sub.add_parser("build")
    b.add_argument("--sample-id", required=True)
    b.add_argument("--arriba-tsv")
    b.add_argument("--starfusion-tsv")
    b.add_argument("--oncokb-json", help="fusions_with_oncokb.json from oncokb-fusion-annotation")
    b.add_argument("--bpview-events", help="bpview events.json (optional)")
    b.add_argument("--pon", help="recurrence panel from `fusion-events pon` (optional)")
    b.add_argument("--out-tsv", required=True)
    b.add_argument("--merge-tolerance", type=int, default=10)
    b.add_argument("--tier1-min-unique-molecules", type=int, default=2)
    b.add_argument("--allow-out-of-frame", action="store_true",
                   help="do not hold out-of-frame / stop-codon events back from tiers 1-2")
    b.add_argument("--pon-min-samples", type=int, default=3)
    b.add_argument("--pon-min-fraction", type=float, default=0.10)
    p = sub.add_parser("pon")
    p.add_argument("--events-tsv", nargs="+", required=True)
    p.add_argument("--out-tsv", required=True)
    opts = ap.parse_args(argv)
    if opts.cmd == "pon":
        n = build_pon(opts.events_tsv, opts.out_tsv)
        print(f"fusion-events {__version__}: PoN with {n} junctions from {len(opts.events_tsv)} samples")
        return
    rows = build_rows(opts)
    write_tsv(opts.out_tsv, rows)
    tiers = {t: sum(1 for r in rows if r["tier"] == t) for t in (1, 2, 3)}
    print(f"fusion-events {__version__}: {len(rows)} events; tier 1/2/3 = {tiers[1]}/{tiers[2]}/{tiers[3]}")


if __name__ == "__main__":
    main()
