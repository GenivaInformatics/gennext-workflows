#!/usr/bin/env python3
"""Annotate the rows a filtered OncoKB run skipped, in place of a full re-run.

Between 2026-09-27 and 2026-09-30, filtered_maf_annotator.py sent only rows whose
Func.refGene was exonic/splicing/UTR/upstream (QUERY_CLASSES below) to OncoKB and
left the other rows' OncoKB columns empty. This takes such an *_oncokb.maf,
sends only the skipped rows to MafAnnotator, and fills their OncoKB columns.
Rows that were already annotated are left untouched. The input is replaced
atomically; the original is kept as <name>.filtered.bak.
"""

import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from AnnotatorCore import get_oncokb_annotation_column_headers
from filtered_maf_annotator import oncokb_echo, read_header


csv.field_size_limit(sys.maxsize)

# The filter's selection, as it was: rows matching none of these were skipped.
QUERY_CLASSES = {
    "exonic", "splicing", "utr5", "utr3", "upstream",
    "ncrna_exonic", "ncrna_splicing", "ncrna_utr5", "ncrna_utr3",
}


def was_skipped(row, func_index):
    return not any(value.strip().lower() in QUERY_CLASSES
                   for value in row[func_index].split(";"))


def split_skipped(maf_path, query_path):
    annotation_header = get_oncokb_annotation_column_headers(False, True)
    with maf_path.open(newline="") as source, query_path.open("w", newline="") as target:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        header = read_header(reader)
        n_input = len(header) - len(annotation_header)
        if header[n_input:] != annotation_header:
            raise ValueError("MAF does not end with the OncoKB annotation columns")
        func_index = header.index("Func.refGene")
        writer.writerow(header[:n_input])
        total = skipped = 0
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"MAF row {total + 1} has {len(row)} columns, expected {len(header)}")
            total += 1
            if was_skipped(row, func_index):
                if any(row[n_input:]):
                    raise ValueError(f"Row {total} matches the filter's skip rule but has OncoKB values")
                writer.writerow(row[:n_input])
                skipped += 1
    return header, n_input, func_index, total, skipped


def merge_backfill(maf_path, annotated_path, output_path, header, n_input, func_index):
    filled = 0
    with maf_path.open(newline="") as source, annotated_path.open(newline="") as annotated_file, \
            output_path.open("w", newline="") as target:
        reader = csv.reader(source, delimiter="\t")
        annotated_reader = csv.reader(annotated_file, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        if read_header(reader) != header:
            raise ValueError("MAF changed during backfill")
        if read_header(annotated_reader) != header:
            raise ValueError("OncoKB output header differs from expected MAF schema")
        writer.writerow(header)
        for row in reader:
            if was_skipped(row, func_index):
                annotated = next(annotated_reader, None)
                if annotated is None or annotated[:n_input] != oncokb_echo(row[:n_input]):
                    raise ValueError(f"OncoKB output row {filled + 1} missing or out of order")
                if len(annotated) != len(header):
                    raise ValueError(f"OncoKB output row {filled + 1} has wrong column count")
                row = row[:n_input] + annotated[n_input:]
                filled += 1
            writer.writerow(row)
        if next(annotated_reader, None) is not None:
            raise ValueError("OncoKB output has more rows than skipped input rows")
    return filled


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maf", required=True, type=Path, help="filtered *_oncokb.maf, updated in place")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--tumor-type", default=".")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="oncokb-backfill-", dir=args.maf.parent) as temp_dir:
        query = Path(temp_dir) / "skipped.maf"
        annotated = Path(temp_dir) / "annotated.maf"
        merged = Path(temp_dir) / "merged.maf"
        header, n_input, func_index, total, skipped = split_skipped(args.maf, query)
        print(f"MAF rows: {total}; previously skipped: {skipped}", flush=True)
        if not skipped:
            return
        command = [sys.executable, "/app/oncokb/MafAnnotator.py", "-i", str(query),
                   "-o", str(annotated), "-q", "Genomic_Change", "-r", args.reference,
                   "-b", os.environ["ONCOKB_API_TOKEN"]]
        if args.tumor_type != ".":
            command.extend(["-t", args.tumor_type])
        subprocess.run(command, check=True)
        filled = merge_backfill(args.maf, annotated, merged, header, n_input, func_index)
        if filled != skipped:
            raise ValueError(f"Filled {filled} rows, expected {skipped}")
        backup = args.maf.with_name(args.maf.name + ".filtered.bak")
        if not backup.exists():
            shutil.copy2(args.maf, backup)
        os.replace(merged, args.maf)
        print(f"Filled OncoKB columns for {filled} of {total} rows", flush=True)


if __name__ == "__main__":
    main()
