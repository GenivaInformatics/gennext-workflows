#!/usr/bin/env python3
"""Propagate the BAM_FDUP (0x400) flag from a genome-coord MarkDuplicates BAM
onto a FusionInspector consolidated BAM by QNAME.

Coordinate-based dedup on small (~30-80 kb) fusion contigs over-flags
fusion-supporting reads because every supporter aligns at the breakpoint
position by construction; running Picard or samtools markdup on FI's BAM
collapses legitimate fusion evidence into "duplicates" (88% dup rate
observed empirically on AGI-S5).

The right view for IGV — and the one arriba and star-fusion actually
consumed — is the dedup verdict produced at *genome* reference
coordinates by the upstream Picard MarkDuplicates pass. This script
collects the QNAMEs flagged 0x400 in the upstream BAM and OR-s 0x400
onto matching records in the FI consolidated BAM. Records not in the
upstream dup set get 0x400 cleared (idempotent on re-runs).

Inputs/outputs:
  --upstream-bam  Genome-aligned, MarkDup'd BAM (e.g. Agi-S5_RNA.dedup.bam).
  --fi-bam        FusionInspector consolidated BAM (no dup flags expected).
  --out-bam       Path to write the dup-flag-propagated BAM. The .bai index
                  is created alongside.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys

import pysam


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--upstream-bam", required=True,
                   help="Genome-coord MarkDuplicates BAM whose 0x400 flags "
                        "are the source of truth.")
    p.add_argument("--fi-bam", required=True,
                   help="FI consolidated BAM (input).")
    p.add_argument("--out-bam", required=True,
                   help="Output BAM path; .bai index is created alongside.")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    logging.info("Building dup-QNAME set from %s", args.upstream_bam)
    dup: set[str] = set()
    with pysam.AlignmentFile(args.upstream_bam, "rb") as up:
        for r in up:
            if r.is_duplicate:
                dup.add(r.query_name)
    logging.info("Dup QNAMEs: %d", len(dup))

    # In-place rewrite is the common case here (Argo mounts fi-bam and
    # out-bam to the same hostPath via different mountPaths, e.g.
    # /mnt/fi/X vs /mnt/out/X). os.path.realpath can't see the
    # hostPath aliasing, so unconditionally stream to a `.tmp` next to
    # out-bam and rename into place. If the dirs are actually distinct,
    # this just adds one rename — cheap and always safe.
    write_path = args.out_bam + ".tmp"

    logging.info("Rewriting %s -> %s (via %s)",
                 args.fi_bam, args.out_bam, write_path)
    n_total = n_set = n_clear = 0
    with pysam.AlignmentFile(args.fi_bam, "rb") as fi:
        with pysam.AlignmentFile(write_path, "wb", template=fi) as out:
            for r in fi:
                n_total += 1
                if r.query_name in dup:
                    if not r.is_duplicate:
                        r.is_duplicate = True
                        n_set += 1
                else:
                    if r.is_duplicate:
                        r.is_duplicate = False
                        n_clear += 1
                out.write(r)
    logging.info(
        "Records: total=%d, set 0x400=%d, cleared 0x400=%d",
        n_total, n_set, n_clear,
    )
    os.replace(write_path, args.out_bam)

    logging.info("Indexing %s", args.out_bam)
    pysam.index(args.out_bam)
    logging.info("Done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
