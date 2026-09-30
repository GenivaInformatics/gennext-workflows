#!/usr/bin/env python3
"""Annotate gene-relevant MAF rows with OncoKB, retaining every input row."""

import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from AnnotatorCore import get_oncokb_annotation_column_headers


csv.field_size_limit(sys.maxsize)

# ANNOVAR's Func.refGene is more precise here than maftools' classification:
# maftools can leave non-exonic classifications blank. Include promoter/UTR and
# splice events, but avoid querying millions of intergenic/intronic WGS calls.
QUERY_CLASSES = {
    "exonic", "splicing", "utr5", "utr3", "upstream",
    "ncrna_exonic", "ncrna_splicing", "ncrna_utr5", "ncrna_utr3",
}
GENOME_KIT_REGIONS_BASENAME = "hg38_canonical_chrom_contig_regions.bed"
GENOME_KIT_VALIDATED_BASENAME = "hg38_canonical_chrom_contig_regions.valid.bed"


def is_genome_kit_regions_file(regions_file):
    return Path(regions_file).name in {
        GENOME_KIT_REGIONS_BASENAME, GENOME_KIT_VALIDATED_BASENAME,
    }


def should_query(row, func_index):
    return any(value.strip().lower() in QUERY_CLASSES
               for value in row[func_index].split(";"))


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


def split_maf(input_path, query_path):
    total = selected = 0
    with input_path.open(newline="") as source, query_path.open("w", newline="") as target:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        header = read_header(reader)
        if "Func.refGene" not in header:
            raise ValueError("MAF lacks Func.refGene; refusing to send unfiltered variants")
        func_index = header.index("Func.refGene")
        writer.writerow(header)
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"MAF row {total + 1} has {len(row)} columns, expected {len(header)}")
            total += 1
            if should_query(row, func_index):
                writer.writerow(row)
                selected += 1
    return header, func_index, total, selected


def merge_maf(input_path, annotated_path, output_path, input_header, func_index, selected,
              skip_all=False):
    annotation_header = get_oncokb_annotation_column_headers(False, True)
    temp_output = output_path.with_name(output_path.name + ".tmp")
    try:
        with input_path.open(newline="") as source, temp_output.open("w", newline="") as target:
            reader = csv.reader(source, delimiter="\t")
            if read_header(reader) != input_header:
                raise ValueError("Input MAF changed after filtering")
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            writer.writerow(input_header + annotation_header)

            annotated_file = annotated_reader = None
            if selected:
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
                    if not skip_all and should_query(row, func_index):
                        annotated = next(annotated_reader, None)
                        if annotated is None or annotated[:len(input_header)] != oncokb_echo(row):
                            raise ValueError(f"OncoKB output row {used + 1} missing or out of order")
                        suffix = annotated[len(input_header):]
                        if len(suffix) != len(annotation_header):
                            raise ValueError(f"OncoKB output row {used + 1} has wrong column count")
                        used += 1
                    writer.writerow(row + suffix)
                if used != selected or (annotated_reader is not None and next(annotated_reader, None) is not None):
                    raise ValueError("OncoKB output row count differs from selected input")
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
    args = parser.parse_args()
    if is_genome_kit_regions_file(args.regions_file):
        with args.input.open(newline="") as source:
            header = read_header(csv.reader(source, delimiter="\t"))
        written = merge_maf(args.input, None, args.output, header, 0, 0, skip_all=True)
        print(f"Genome kit regions file: {args.regions_file}; skipped OncoKB for all {written} MAF rows", flush=True)
        return

    with tempfile.TemporaryDirectory(prefix="oncokb-maf-", dir=args.output.parent) as temp_dir:
        query = Path(temp_dir) / "query.maf"
        annotated = Path(temp_dir) / "annotated.maf"
        header, func_index, total, selected = split_maf(args.input, query)
        print(f"MAF rows: {total}; OncoKB candidates: {selected}; skipped: {total - selected}", flush=True)
        if selected:
            token = os.environ["ONCOKB_API_TOKEN"]
            command = [sys.executable, "/app/oncokb/MafAnnotator.py", "-i", str(query),
                       "-o", str(annotated), "-q", "Genomic_Change", "-r", args.reference,
                       "-b", token]
            if args.tumor_type != ".":
                command.extend(["-t", args.tumor_type])
            subprocess.run(command, check=True)
        written = merge_maf(args.input, annotated, args.output, header, func_index, selected)
        if written != total:
            raise ValueError("Merged MAF row count differs from input")
        print(f"Wrote {written} MAF rows with {selected} OncoKB annotations", flush=True)


if __name__ == "__main__":
    main()
