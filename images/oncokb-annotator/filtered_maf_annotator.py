#!/usr/bin/env python3
"""Annotate every MAF row with OncoKB (Genome kit: skip OncoKB), retaining input text."""

import argparse
import csv
import logging
import os
from pathlib import Path
import sys
import tempfile

import AnnotatorCore
from AnnotatorCore import MutationStatus, get_oncokb_annotation_column_headers


csv.field_size_limit(sys.maxsize)

GENOME_KIT_REGIONS_BASENAME = "hg38_canonical_chrom_contig_regions.bed"
GENOME_KIT_VALIDATED_BASENAME = "hg38_canonical_chrom_contig_regions.valid.bed"


def is_genome_kit_regions_file(regions_file):
    return Path(regions_file).name in {
        GENOME_KIT_REGIONS_BASENAME, GENOME_KIT_VALIDATED_BASENAME,
    }


def maf_annotation_header():
    """The OncoKB columns MafAnnotator.py appends to a MAF.

    Since upstream v4.0, MafAnnotator.py always annotates in mixed
    somatic/germline mode (processalterationevents_mixed, combined_mode=True).
    That adds PATHOGENIC (Germline), PENETRANCE (Germline) and
    GENOMIC_INDICATOR (Germline) after ONCOGENIC. Our MAFs have no
    Mutation_Status column, so every row is annotated as somatic and those three
    columns stay empty.
    """
    return get_oncokb_annotation_column_headers(
        False, True, MutationStatus.SOMATIC, combined_mode=True)


# Rows per upstream call. Upstream v3 posted genomic-change queries to OncoKB
# in batches of POST_QUERIES_THRESHOLD_GC_HGVSG (100). v4's mixed mode sends a
# whole file's queries in one POST, which OncoKB answers with HTTP 500 for a
# MAF of tens of thousands of rows. Annotating the MAF in chunks of that size
# restores v3's request size.
CHUNK_ROWS = AnnotatorCore.POST_QUERIES_THRESHOLD_GC_HGVSG


def annotate_maf(input_path, output_path, token, reference, tumor_type=".",
                 api_url=None, chunk_rows=CHUNK_ROWS):
    """MafAnnotator.py -q Genomic_Change, run chunk by chunk. Returns the rows sent.

    Does what MafAnnotator.main does for our options, once: set the API URL and
    token, validate the token. Then it calls upstream's processalterationevents
    on consecutive chunks of the input and joins their output. Each chunk is
    the input's header line plus a run of its data lines, copied byte for
    byte, so upstream reads the same text as from the whole file. The token
    is never put on a command line.
    """
    if api_url:
        AnnotatorCore.setoncokbbaseurl(api_url)
    AnnotatorCore.setoncokbapitoken(token)
    AnnotatorCore.validate_oncokb_token()
    query_type = AnnotatorCore.QueryType.GENOMIC_CHANGE
    reference_genome = AnnotatorCore.ReferenceGenome[reference.upper()]
    cancer_type = "" if tumor_type == "." else tumor_type

    with tempfile.TemporaryDirectory(prefix="chunks-", dir=output_path.parent) as temp_dir, \
            input_path.open(newline="") as source, output_path.open("w", newline="") as target:
        chunk_in = Path(temp_dir) / "in.maf"
        chunk_out = Path(temp_dir) / "out.maf"
        header = None
        state = {"sent": 0, "chunks": 0, "header_written": False}

        def annotate(chunk):
            with chunk_in.open("w", newline="") as part:
                part.write(header)
                part.writelines(chunk)
            AnnotatorCore.processalterationevents(
                str(chunk_in), str(chunk_out), "", cancer_type, {}, False,
                query_type, reference_genome, False)
            with chunk_out.open(newline="") as part:
                part_header = part.readline()
                if not state["header_written"]:
                    target.write(part_header)
                    state["header_written"] = True
                target.writelines(part)
            state["sent"] += len(chunk)
            state["chunks"] += 1
            if state["chunks"] % 20 == 0:
                print(f"OncoKB: {state['sent']} MAF rows annotated", flush=True)

        chunk = []
        for line in source:
            if header is None:
                if line.strip() and not line.startswith("#"):
                    header = line if line.endswith("\n") else line + "\n"
                continue
            chunk.append(line if line.endswith("\n") else line + "\n")
            if len(chunk) == chunk_rows:
                annotate(chunk)
                chunk = []
        if chunk:
            annotate(chunk)
        return state["sent"]


def oncokb_echo(row):
    # AnnotatorCore.append_annotation_to_file writes every row through
    # encode('ascii', 'ignore'), so its echo of the input drops non-ASCII
    # characters (ClinVar CLNDN "Muir-Torré" comes back as "Muir-Torr").
    return [value.encode("ascii", "ignore").decode("ascii") for value in row]


def read_header(reader):
    for row in reader:
        if row and not row[0].startswith("#"):
            return row
    raise ValueError("MAF has no header")


def count_maf(input_path):
    total = 0
    with input_path.open(newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = read_header(reader)
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"MAF row {total + 1} has {len(row)} columns, expected {len(header)}")
            total += 1
    return header, total


def merge_maf(input_path, annotated_path, output_path, input_header, skip_all=False):
    annotation_header = maf_annotation_header()
    temp_output = output_path.with_name(output_path.name + ".tmp")
    try:
        with input_path.open(newline="") as source, temp_output.open("w", newline="") as target:
            reader = csv.reader(source, delimiter="\t")
            if read_header(reader) != input_header:
                raise ValueError("Input MAF changed after filtering")
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            writer.writerow(input_header + annotation_header)

            annotated_file = annotated_reader = None
            if not skip_all:
                annotated_file = annotated_path.open(newline="")
                annotated_reader = csv.reader(annotated_file, delimiter="\t")
                header = read_header(annotated_reader)
                if header != input_header + annotation_header:
                    raise ValueError("OncoKB output header differs from expected MAF schema")

            try:
                seen = used = 0
                for row in reader:
                    if len(row) != len(input_header):
                        raise ValueError(f"Input MAF row {seen + 1} changed")
                    seen += 1
                    suffix = [""] * len(annotation_header)
                    if not skip_all:
                        annotated = next(annotated_reader, None)
                        if annotated is None or annotated[:len(input_header)] != oncokb_echo(row):
                            raise ValueError(f"OncoKB output row {used + 1} missing or out of order")
                        suffix = annotated[len(input_header):]
                        if len(suffix) != len(annotation_header):
                            raise ValueError(f"OncoKB output row {used + 1} has wrong column count")
                        used += 1
                    writer.writerow(row + suffix)
                if annotated_reader is not None and next(annotated_reader, None) is not None:
                    raise ValueError("OncoKB output has more rows than the input MAF")
            finally:
                if annotated_file is not None:
                    annotated_file.close()
        os.replace(temp_output, output_path)
        return seen
    finally:
        temp_output.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--tumor-type", default=".")
    parser.add_argument("--regions-file", default=".")
    # Only for tests against demo.oncokb.org; production uses upstream's default.
    parser.add_argument("--api-url", default=None)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if is_genome_kit_regions_file(args.regions_file):
        with args.input.open(newline="") as source:
            header = read_header(csv.reader(source, delimiter="\t"))
        written = merge_maf(args.input, None, args.output, header, skip_all=True)
        print(f"Genome kit regions file: {args.regions_file}; skipped OncoKB for all {written} MAF rows", flush=True)
        return

    with tempfile.TemporaryDirectory(prefix="oncokb-maf-", dir=args.output.parent) as temp_dir:
        annotated = Path(temp_dir) / "annotated.maf"
        header, total = count_maf(args.input)
        print(f"MAF rows: {total}; sending all to OncoKB", flush=True)
        if total:
            sent = annotate_maf(args.input, annotated, os.environ["ONCOKB_API_TOKEN"].strip(),
                                args.reference, args.tumor_type, args.api_url)
            if sent != total:
                # A record spanning lines would split across chunks.
                raise ValueError(f"MAF has {total} records but {sent} data lines")
        written = merge_maf(args.input, annotated, args.output, header, skip_all=not total)
        if written != total:
            raise ValueError("Merged MAF row count differs from input")
        print(f"Wrote {written} MAF rows with OncoKB annotations", flush=True)

if __name__ == "__main__":
    main()
