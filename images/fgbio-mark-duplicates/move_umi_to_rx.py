#!/usr/bin/env python3
"""Move the UMI from the QNAME suffix into the RX:Z BAM tag.

fgbio's GroupReadsByUmi (and most UMI-aware tools that key off a tag) reads
the UMI from a configurable BAM tag, defaulting to RX. fastp's UMI extraction
puts the UMI into the QNAME instead — appended as `:UMI_<R1bases>_<R2bases>`
on each read name. This pre-step bridges those conventions:

    Input QNAME : R213040190036:20260204212056:V350396806:1:674046:2:58:GT:UMI_GGCTG_ATGTC
    Output QNAME: R213040190036:20260204212056:V350396806:1:674046:2:58:GT
    Output RX  : Z:GGCTG_ATGTC

The composite UMI token includes the underscore between R1 and R2 UMIs;
fgbio's Adjacency strategy treats the whole 11-character string as one UMI
(opaque clustering, edit-distance 1 by default), which is correct for
paired-end UMI inference.

The match is anchored to the literal `:UMI_` suffix produced by fastp's
`--umi_prefix=UMI` setting; if that ever changes, update SUFFIX_RE.

Records that don't match the expected suffix are passed through unchanged
and counted in --no-umi-out (default: warn at end).
"""
from __future__ import annotations

import argparse
import logging
import re
import sys

import pysam


SUFFIX_RE = re.compile(r"^(?P<head>.*?):UMI_(?P<umi>[ACGTN_]+)$")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-bam", required=True,
                   help="Input BAM (sorted+indexed); QNAMEs must end in "
                        "`:UMI_<R1>_<R2>` (the layout fastp produces with "
                        "--umi --umi_prefix=UMI).")
    p.add_argument("--out-bam", required=True,
                   help="Output BAM. Same sort order as input. RX tag set on "
                        "every record whose QNAME matched the suffix; QNAME "
                        "rewritten to drop the suffix.")
    p.add_argument("--rx-tag", default="RX",
                   help="BAM tag to write the captured UMI to (default: RX).")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    n_total = n_set = n_no_umi = 0
    sample_qname = None
    sample_rx = None

    with pysam.AlignmentFile(args.in_bam, "rb") as inb:
        with pysam.AlignmentFile(args.out_bam, "wb", template=inb) as outb:
            for r in inb:
                n_total += 1
                m = SUFFIX_RE.match(r.query_name)
                if m is None:
                    n_no_umi += 1
                    outb.write(r)
                    continue
                head = m.group("head")
                umi = m.group("umi")
                r.query_name = head
                r.set_tag(args.rx_tag, umi, value_type="Z")
                n_set += 1
                if sample_qname is None:
                    sample_qname = head
                    sample_rx = umi
                outb.write(r)

    logging.info("Records: %d total, %d UMI tagged, %d no UMI suffix",
                 n_total, n_set, n_no_umi)
    if sample_qname is not None:
        logging.info("Sample mapping: head=%r %s=Z:%s",
                     sample_qname, args.rx_tag, sample_rx)
    if n_set == 0:
        logging.error("No records matched the UMI suffix pattern; "
                      "check fastp's --umi_prefix is 'UMI' and that the "
                      "delimiter before UMI is ':'.")
        return 2
    if n_no_umi > 0:
        logging.warning("%d records had no UMI suffix and were passed through "
                        "unchanged. They will be ignored by fgbio's UMI "
                        "grouping (no RX tag set).", n_no_umi)
    return 0


if __name__ == "__main__":
    sys.exit(main())
