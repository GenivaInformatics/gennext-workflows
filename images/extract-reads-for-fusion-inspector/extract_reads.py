#!/usr/bin/env python3
"""Pre-filter the upstream genome-aligned BAM for FusionInspector visualization.

Drops candidate fusions matching well-known artifact classes (mitochondrial,
IG@/IGH/IGL/IGK, sno/snRNA, pseudogenes, antisense, unannotated clones,
self-fusions) and keeps reads that either:
  - have a supplementary alignment (junction reads / split-reads), or
  - align within a ±N-kb window around any surviving fusion breakpoint
    (spanning fragments) when --breakpoint-window > 0.

Setting --breakpoint-window=0 disables the BED-based spanning-fragment
collection and yields a chimeric-only read set: smaller, faster for FI to
realign, but FI's SpanningFragCount will be ~0. Pipeline-wide counts
should be sourced from fusion-report's fusions.json (caller-native).

Dup-flagged reads ARE retained — they are needed downstream so the post-FI
samtools markdup pass on the realigned consolidated.bam can re-flag PCR
duplicates against FI's per-contig coordinate space (the upstream BAM's
0x400 flags are tied to genome coordinates and are not meaningful after
re-alignment to short fusion contigs).

Both mates of every kept pair are extracted via samtools view -N. Outputs
gzipped paired FASTQs plus the filtered fusion-list TSV.
"""
from __future__ import annotations

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Artifact regex set, applied to gene names from STAR-Fusion / Arriba TSVs
# and the consolidated fusion list. Drop a candidate if EITHER partner
# matches any pattern, or if the partners are identical (self-fusion).
ARTIFACT_PATTERNS = [
    re.compile(r"^MT-"),                  # mitochondrial
    re.compile(r"^IG[HKL][@JCV]"),        # IG@-ext, IGHJ*, IGLC*, IGKV*, ...
    re.compile(r"^IGH-?@"),
    re.compile(r"^IGL-?@"),
    re.compile(r"^IGK-?@"),
    re.compile(r"^TR[ABDG][JCV]"),        # T-cell receptor family
    re.compile(r"^RNU\d"),                # snRNA
    re.compile(r"^SNORD\d"),              # snoRNA C/D
    re.compile(r"^SNORA\d"),              # snoRNA H/ACA
    re.compile(r"^SCARNA\d"),             # small Cajal-body RNA
    re.compile(r"^RNA5S"),                # 5S rRNA family (incl. RNA5SP*)
    re.compile(r"^RNA5-8S"),              # 5.8S rRNA
    re.compile(r"^RNA28S"),               # 28S rRNA
    re.compile(r"^MIR\d"),                # microRNAs
    re.compile(r"\dP\d+$"),               # pseudogenes (digit-anchored to
                                          # avoid false-positives like ARHGAP26)
    re.compile(r"-AS\d+$"),               # antisense lncRNAs
    re.compile(r"^ENSG\d"),               # unannotated ENSG* loci
    re.compile(r"^LINC\d"),               # generic lincRNAs
    re.compile(r"^AC\d{6}\.\d"),          # uncharacterized clones
    re.compile(r"^AL\d{6}\.\d"),
    re.compile(r"^AP\d{6}\.\d"),
    re.compile(r"^FP\d{6}\.\d"),
]


def is_artifact(gene: str) -> bool:
    return any(p.search(gene) for p in ARTIFACT_PATTERNS)


def pair_keep(a: str, b: str) -> bool:
    return bool(a) and bool(b) and a != b and not is_artifact(a) and not is_artifact(b)


def run(cmd, **kwargs):
    """Run a command and stream its output."""
    if isinstance(cmd, str):
        logging.info("running: %s", cmd)
    else:
        logging.info("running: %s", " ".join(str(c) for c in cmd))
    return subprocess.run(cmd, check=True, **kwargs)


def ensure_sorted_indexed(in_bam: Path, work: Path, threads: int) -> Path:
    """Sort and index the BAM if it isn't already coord-sorted+indexed.
    Returns the path to the indexed BAM (either in_bam or work/sorted.bam)."""
    bai_a = Path(str(in_bam) + ".bai")
    bai_b = in_bam.with_suffix(".bai")
    if bai_a.exists() or bai_b.exists():
        logging.info("input already indexed at %s; skipping sort", in_bam)
        return in_bam

    header = subprocess.run(
        ["samtools", "view", "-H", str(in_bam)],
        check=True, capture_output=True, text=True,
    ).stdout
    so = "unknown"
    for line in header.splitlines():
        if line.startswith("@HD"):
            for f in line.split("\t"):
                if f.startswith("SO:"):
                    so = f[3:].strip()
            break

    sorted_path = work / "sorted.bam"
    if so == "coordinate":
        logging.info("input is coord-sorted but unindexed; symlinking and indexing")
        sorted_path.symlink_to(in_bam.resolve())
    else:
        logging.info("sorting input BAM (header SO=%s)", so)
        run([
            "samtools", "sort",
            "-@", str(threads), "-m", "2G",
            "-T", str(work / "sort"),
            "-o", str(sorted_path),
            str(in_bam),
        ])
    run(["samtools", "index", "-@", str(threads), str(sorted_path)])
    return sorted_path


def parse_starfusion_breakpoints(tsv: Path):
    """Yield (chrom, pos) for each breakpoint of every surviving fusion."""
    if not tsv.exists() or os.path.getsize(tsv) == 0:
        return
    with open(tsv) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        try:
            iName = header.index("FusionName")
            iL = header.index("LeftBreakpoint")
            iR = header.index("RightBreakpoint")
        except ValueError as e:
            logging.warning("STAR-Fusion TSV missing column: %s", e)
            return
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(iL, iR):
                continue
            name = f[iName]
            if "--" not in name:
                continue
            a, b = name.split("--", 1)
            if not pair_keep(a, b):
                continue
            for col in (iL, iR):
                parts = f[col].split(":")
                if len(parts) < 2:
                    continue
                try:
                    yield parts[0], int(parts[1])
                except ValueError:
                    continue


def parse_arriba_breakpoints(tsv: Path):
    if not tsv.exists() or os.path.getsize(tsv) == 0:
        return
    with open(tsv) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        try:
            i1 = next(i for i, c in enumerate(header) if c.lstrip("#") == "gene1")
            i2 = next(i for i, c in enumerate(header) if c == "gene2")
            iB1 = header.index("breakpoint1")
            iB2 = header.index("breakpoint2")
        except (ValueError, StopIteration) as e:
            logging.warning("Arriba TSV missing column: %s", e)
            return
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(iB1, iB2):
                continue
            g1 = f[i1] if i1 < len(f) else ""
            g2 = f[i2] if i2 < len(f) else ""
            if not pair_keep(g1, g2):
                continue
            for col in (iB1, iB2):
                parts = f[col].split(":")
                if len(parts) < 2:
                    continue
                try:
                    yield parts[0], int(parts[1])
                except ValueError:
                    continue


def filter_fusion_list(in_path: Path, out_path: Path) -> tuple[int, int]:
    """Filter the consolidated geneA--geneB list. Returns (kept, dropped)."""
    kept = dropped = 0
    with open(in_path) as src, open(out_path, "w") as dst:
        for line in src:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            name = stripped.split("\t", 1)[0]
            if "--" not in name:
                continue
            a, b = name.split("--", 1)
            if pair_keep(a, b):
                dst.write(stripped + "\n")
                kept += 1
            else:
                dropped += 1
    return kept, dropped


def build_breakpoint_bed(
    starfusion_tsv: Path, arriba_tsv: Path, win: int, out_path: Path
) -> tuple[int, int]:
    """Sort and merge ±win bp intervals around all surviving fusion breakpoints.

    win <= 0 short-circuits to an empty BED, signalling chimeric-only mode
    (collect_read_names then skips the loci-scan branch and the read set
    collapses to supplementary alignments only).
    """
    if win <= 0:
        out_path.write_text("")
        return 0, 0
    raw: list[tuple[str, int, int]] = []
    for chrom, pos in parse_starfusion_breakpoints(starfusion_tsv):
        raw.append((chrom, max(0, pos - win), pos + win))
    for chrom, pos in parse_arriba_breakpoints(arriba_tsv):
        raw.append((chrom, max(0, pos - win), pos + win))
    raw.sort()
    merged: list[list] = []
    for chrom, s, e in raw:
        if merged and merged[-1][0] == chrom and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([chrom, s, e])
    with open(out_path, "w") as fh:
        for chrom, s, e in merged:
            fh.write(f"{chrom}\t{s}\t{e}\n")
    return len(raw), len(merged)


def chrom_sanity_check(bed_path: Path, indexed_bam: Path) -> None:
    bed_chroms: set[str] = set()
    with open(bed_path) as fh:
        for line in fh:
            if line.strip():
                bed_chroms.add(line.split("\t", 1)[0])
    if not bed_chroms:
        return
    header = subprocess.run(
        ["samtools", "view", "-H", str(indexed_bam)],
        check=True, capture_output=True, text=True,
    ).stdout
    bam_chroms: set[str] = set()
    for line in header.splitlines():
        if line.startswith("@SQ"):
            for f in line.split("\t"):
                if f.startswith("SN:"):
                    bam_chroms.add(f[3:])
    if not bed_chroms & bam_chroms:
        logging.warning(
            "chrom mismatch: BED uses %s but BAM has %s; spanning-fragment filter will return 0 reads",
            sorted(bed_chroms)[:5], sorted(bam_chroms)[:5],
        )


def collect_read_names(
    indexed_bam: Path, bed_path: Path, work: Path, threads: int
) -> tuple[Path, int, int, int]:
    """Build the union of supplementary-alignment + breakpoint-locus read names.
    Returns (names_file, n_supp, n_loci, n_total)."""
    supp = work / "names_supp.txt"
    loci = work / "names_loci.txt"
    out = work / "read_names.txt"

    # No -F 1024: dup-flagged reads are intentionally retained so post-FI
    # samtools markdup can re-flag them against FI's per-contig coordinate
    # space. Excluding them here would leave us with no way to recover
    # dup-flagged status after re-alignment.
    run(
        f"samtools view -@ {threads} -f 2048 {indexed_bam} | cut -f1 | sort -u > {supp}",
        shell=True,
    )
    if bed_path.exists() and os.path.getsize(bed_path) > 0:
        run(
            f"samtools view -@ {threads} -L {bed_path} -M {indexed_bam} | cut -f1 | sort -u > {loci}",
            shell=True,
        )
    else:
        loci.write_text("")
    run(f"sort -m -u {supp} {loci} > {out}", shell=True)
    n_supp = sum(1 for _ in open(supp))
    n_loci = sum(1 for _ in open(loci))
    n_total = sum(1 for _ in open(out))
    return out, n_supp, n_loci, n_total


def extract_pairs_to_fastq(
    indexed_bam: Path, names_file: Path, out_r1: Path, out_r2: Path,
    work: Path, threads: int,
) -> None:
    namesort = work / "filtered_namesort.bam"
    # Pulls all records (including dup-flagged) for each kept read name —
    # FI re-aligns to per-contig fusion sequences and the post-FI markdup
    # pass re-flags dups in that new coordinate space.
    run(
        f"samtools view -@ {threads} -h -N {names_file} -b {indexed_bam} | "
        f"samtools sort -@ {threads} -n -m 2G -T {work}/sortn -o {namesort} -",
        shell=True,
    )
    run([
        "samtools", "fastq",
        "-@", str(threads),
        "-1", str(out_r1), "-2", str(out_r2),
        "-0", "/dev/null", "-s", "/dev/null",
        "-N", str(namesort),
    ])
    namesort.unlink(missing_ok=True)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in-bam", required=True, type=Path)
    p.add_argument("--starfusion-tsv", type=Path, default=None,
                   help="STAR-Fusion fusion-calls TSV (LeftBreakpoint/RightBreakpoint cols).")
    p.add_argument("--arriba-tsv", type=Path, default=None,
                   help="Arriba fusion-calls TSV (breakpoint1/breakpoint2 cols).")
    p.add_argument("--fusion-list-tsv", required=True, type=Path,
                   help="Consolidated geneA--geneB list (one per line).")
    p.add_argument("--out-r1", required=True, type=Path)
    p.add_argument("--out-r2", required=True, type=Path)
    p.add_argument("--out-filtered-list", required=True, type=Path,
                   help="Output path for the artifact-filtered fusion list.")
    p.add_argument("--read-count-out", default=None, type=Path,
                   help="If given, write the total read-name count to this file.")
    p.add_argument("--breakpoint-window", type=int, default=0,
                   help="±bp window around each fusion breakpoint to collect "
                        "spanning fragments. 0 disables the BED-based scan and "
                        "the script returns chimeric (supplementary) reads only.")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--work", type=Path, default=Path("/tmp/extract"))
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    args.work.mkdir(parents=True, exist_ok=True)
    args.out_r1.parent.mkdir(parents=True, exist_ok=True)
    args.out_filtered_list.parent.mkdir(parents=True, exist_ok=True)

    # Ignore "." sentinel paths from Argo (used for optional inputs).
    sf_tsv = args.starfusion_tsv if args.starfusion_tsv and str(args.starfusion_tsv) != "." else Path("/dev/null")
    ar_tsv = args.arriba_tsv if args.arriba_tsv and str(args.arriba_tsv) != "." else Path("/dev/null")

    # 1. Filter the consolidated fusion list.
    kept, dropped = filter_fusion_list(args.fusion_list_tsv, args.out_filtered_list)
    logging.info("Fusion list: %d kept, %d dropped (artifact / self-fusion filter)",
                 kept, dropped)
    if kept == 0:
        logging.warning("All upstream candidates filtered out — emitting empty FASTQs and exiting.")
        for p_ in (args.out_r1, args.out_r2):
            run(["bash", "-c", f"echo -n '' | gzip -c > {p_}"])
        if args.read_count_out:
            args.read_count_out.write_text("0\n")
        return 0

    # 2. Sort + index the input BAM if not already.
    indexed_bam = ensure_sorted_indexed(args.in_bam, args.work, args.threads)

    # 3. Build the merged breakpoint BED from filtered caller TSVs.
    raw_n, merged_n = build_breakpoint_bed(
        sf_tsv, ar_tsv, args.breakpoint_window, args.work / "breakpoints.bed",
    )
    logging.info("Breakpoint windows: %d raw -> %d after merge (±%d bp)",
                 raw_n, merged_n, args.breakpoint_window)

    # 4. Sanity-check chrom-name compat.
    chrom_sanity_check(args.work / "breakpoints.bed", indexed_bam)

    # 5. Collect read names (supplementary + breakpoint-window).
    names_file, n_supp, n_loci, n_total = collect_read_names(
        indexed_bam, args.work / "breakpoints.bed", args.work, args.threads,
    )
    logging.info("Read names: %d supp + %d breakpoint-window -> %d unique",
                 n_supp, n_loci, n_total)
    if args.read_count_out:
        args.read_count_out.write_text(f"{n_total}\n")
    if n_total == 0:
        logging.error("Zero read names after filter; nothing for FI to do.")
        return 1

    # 6. Pull both mates of each pair and convert to gzipped FASTQs.
    extract_pairs_to_fastq(
        indexed_bam, names_file, args.out_r1, args.out_r2,
        args.work, args.threads,
    )
    logging.info("Filtered FASTQs ready: %s, %s", args.out_r1, args.out_r2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
