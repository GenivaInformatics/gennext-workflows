import csv
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from AnnotatorCore import get_oncokb_annotation_column_headers
from filtered_maf_annotator import (is_genome_kit_regions_file, main, merge_maf,
                                    should_query, split_maf)
from filtered_sv_annotator import main as sv_main


HEADER = ["Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position",
          "Reference_Allele", "Tumor_Seq_Allele2", "Hugo_Symbol",
          "Variant_Classification", "Func.refGene", "Otherinfo10"]
ROWS = [
    ["S1", "chr1", "1", "1", "A", "G", "Unknown", "IGR", "intergenic", "PASS"],
    ["S1", "chr1", "2", "2", "A", "G", "GENE1", "Intron", "intronic", "PASS"],
    ["S1", "chr1", "3", "3", "A", "G", "GENE1", "Missense_Mutation", "exonic", "PASS"],
    ["S1", "chr1", "4", "4", "A", "G", "GENE2", "Splice_Site", "splicing", "PASS"],
    ["S1", "chr1", "5", "5", "A", "G", "TERT", "5'Flank", "upstream", "PASS"],
]


def write_tsv(path, header, rows):
    with path.open("w", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


class FilteredMafTest(unittest.TestCase):
    def test_genome_kit_writes_empty_oncokb_columns_without_token_or_api(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = [Path(temp) / name for name in ("in.maf", "out.maf")]
            write_tsv(source, HEADER, ROWS)
            argv = ["filtered_maf_annotator.py", "--input", str(source), "--output",
                    str(output), "--reference", "GRCh38", "--regions-file",
                    "/sample/hg38_canonical_chrom_contig_regions.bed"]
            with patch.object(sys, "argv", argv), patch.dict(os.environ, {}, clear=True), \
                    patch("filtered_maf_annotator.subprocess.run") as run:
                main()
                run.assert_not_called()
            with output.open(newline="") as file:
                rows = list(csv.reader(file, delimiter="\t"))
            annotation_header = get_oncokb_annotation_column_headers(False, True)
            self.assertEqual(rows[0], HEADER + annotation_header)
            self.assertEqual(rows[1:], [row + [""] * len(annotation_header) for row in ROWS])
            self.assertEqual(is_genome_kit_regions_file("/sample/panel.bed"), False)
            self.assertTrue(is_genome_kit_regions_file(
                "/sample/hg38_canonical_chrom_contig_regions.valid.bed"))

    def test_genome_kit_sv_decoy_preserves_rows_without_token_or_api(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = [Path(temp) / name for name in ("sv.tsv", "annot.tsv")]
            header = ["Sample", "GeneA", "GeneB", "SVType"]
            data = [["S1", "GENE1", "GENE2", "DELETION"]]
            write_tsv(source, header, data)
            argv = ["filtered_sv_annotator.py", "--input", str(source), "--output",
                    str(output), "--regions-file",
                    "/sample/hg38_canonical_chrom_contig_regions.valid.bed"]
            with patch.object(sys, "argv", argv), patch.dict(os.environ, {}, clear=True), \
                    patch("filtered_sv_annotator.subprocess.run") as run:
                sv_main()
                run.assert_not_called()
            with output.open(newline="") as file:
                rows = list(csv.reader(file, delimiter="\t"))
            suffix = get_oncokb_annotation_column_headers(False, False)
            self.assertEqual(rows, [header + suffix, data[0] + [""] * len(suffix)])

    def test_selects_gene_relevant_rows_and_preserves_full_maf(self):
        with tempfile.TemporaryDirectory() as temp:
            source, query, annotated, output = [Path(temp) / name for name in
                                                ("input.maf", "query.maf", "annotated.maf", "out.maf")]
            write_tsv(source, HEADER, ROWS)
            header, func_index, total, selected = split_maf(source, query)
            self.assertEqual((total, selected), (5, 3))
            with query.open(newline="") as file:
                query_rows = list(csv.reader(file, delimiter="\t"))
            self.assertEqual(query_rows[1:], ROWS[2:])

            annotation_header = get_oncokb_annotation_column_headers(False, True)
            annotations = []
            for row in query_rows[1:]:
                suffix = [""] * len(annotation_header)
                suffix[annotation_header.index("ONCOGENIC")] = "Oncogenic"
                annotations.append(row + suffix)
            write_tsv(annotated, HEADER + annotation_header, annotations)
            self.assertEqual(merge_maf(source, annotated, output, header, func_index, selected), 5)
            with output.open(newline="") as file:
                rows = list(csv.reader(file, delimiter="\t"))
            self.assertEqual(rows[0], HEADER + annotation_header)
            self.assertEqual([row[:len(HEADER)] for row in rows[1:]], ROWS)
            self.assertEqual([row[rows[0].index("ONCOGENIC")] for row in rows[1:]],
                             ["", "", "Oncogenic", "Oncogenic", "Oncogenic"])

    def test_empty_candidate_maf_skips_api_and_keeps_rows(self):
        with tempfile.TemporaryDirectory() as temp:
            source, query, output = [Path(temp) / name for name in ("in.maf", "query.maf", "out.maf")]
            write_tsv(source, HEADER, ROWS[:2])
            header, func_index, total, selected = split_maf(source, query)
            self.assertEqual((total, selected), (2, 0))
            merge_maf(source, Path(temp) / "missing.maf", output, header, func_index, selected)
            with output.open(newline="") as file:
                rows = list(csv.reader(file, delimiter="\t"))
            self.assertEqual(len(rows), 3)
            self.assertEqual([row[:len(HEADER)] for row in rows[1:]], ROWS[:2])

    def test_refuses_missing_classification_column(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "in.maf"
            write_tsv(source, HEADER[:-2], [row[:-2] for row in ROWS])
            with self.assertRaisesRegex(ValueError, "Func.refGene"):
                split_maf(source, Path(temp) / "query.maf")

    def test_refuses_mismatched_oncokb_row(self):
        with tempfile.TemporaryDirectory() as temp:
            source, query, annotated, output = [Path(temp) / name for name in
                                                ("in.maf", "query.maf", "annotated.maf", "out.maf")]
            write_tsv(source, HEADER, ROWS)
            header, func_index, _, selected = split_maf(source, query)
            annotation_header = get_oncokb_annotation_column_headers(False, True)
            bad = ROWS[2].copy()
            bad[2] = "999"
            write_tsv(annotated, HEADER + annotation_header,
                      [bad + [""] * len(annotation_header)] * selected)
            with self.assertRaisesRegex(ValueError, "out of order"):
                merge_maf(source, annotated, output, header, func_index, selected)
            self.assertFalse(output.exists())

    def test_accepts_oncokb_ascii_stripped_echo_and_keeps_original_text(self):
        with tempfile.TemporaryDirectory() as temp:
            source, query, annotated, output = [Path(temp) / name for name in
                                                ("in.maf", "query.maf", "annotated.maf", "out.maf")]
            rows = [row.copy() for row in ROWS]
            rows[2][-1] = "Muir-Torré_syndrome"
            write_tsv(source, HEADER, rows)
            header, func_index, _, selected = split_maf(source, query)
            annotation_header = get_oncokb_annotation_column_headers(False, True)
            with query.open(newline="") as file:
                query_rows = list(csv.reader(file, delimiter="\t"))[1:]
            # Mirror AnnotatorCore.append_annotation_to_file.
            with annotated.open("w") as file:
                file.write("\t".join(HEADER + annotation_header) + "\n")
                for row in query_rows:
                    line = "\t".join(row + [""] * len(annotation_header))
                    file.write(line.encode("ascii", "ignore").decode("ascii") + "\n")
            self.assertEqual(merge_maf(source, annotated, output, header, func_index, selected), 5)
            with output.open(newline="") as file:
                out_rows = list(csv.reader(file, delimiter="\t"))
            self.assertEqual([row[:len(HEADER)] for row in out_rows[1:]], rows)


if __name__ == "__main__":
    unittest.main()
