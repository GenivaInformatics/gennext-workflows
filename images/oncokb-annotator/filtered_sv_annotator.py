#!/usr/bin/env python3
"""Keep OncoKB SV/CNV schema without API calls for the Genome kit."""

import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys

from AnnotatorCore import get_oncokb_annotation_column_headers
from filtered_maf_annotator import is_genome_kit_regions_file, read_header


csv.field_size_limit(sys.maxsize)


def write_empty_annotations(input_path, output_path):
    suffix = get_oncokb_annotation_column_headers(False, False)
    temporary = output_path.with_name(output_path.name + ".tmp")
    try:
        with input_path.open(newline="") as source, temporary.open("w", newline="") as target:
            reader = csv.reader(source, delimiter="\t")
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            header = read_header(reader)
            writer.writerow(header + suffix)
            total = 0
            for row in reader:
                if len(row) != len(header):
                    raise ValueError(f"SV/CNV row {total + 1} has the wrong column count")
                writer.writerow(row + [""] * len(suffix))
                total += 1
        os.replace(temporary, output_path)
        return total
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--regions-file", default=".")
    parser.add_argument("--tumor-type", default=".")
    args = parser.parse_args()
    if is_genome_kit_regions_file(args.regions_file):
        total = write_empty_annotations(args.input, args.output)
        print(f"Genome kit regions file: {args.regions_file}; skipped OncoKB for all {total} SV/CNV rows", flush=True)
        return
    command = [sys.executable, "/app/oncokb/StructuralVariantAnnotator.py", "-i", str(args.input),
               "-o", str(args.output), "-b", os.environ["ONCOKB_API_TOKEN"]]
    if args.tumor_type != ".":
        command.extend(["-t", args.tumor_type])
    subprocess.run(command, check=True)
    if not args.output.is_file() or args.output.stat().st_size == 0:
        raise RuntimeError("OncoKB SV/CNV annotator produced no output")


if __name__ == "__main__":
    main()
