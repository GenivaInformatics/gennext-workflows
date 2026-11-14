import csv
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import io

import AnnotatorCore
from AnnotatorCore import MutationStatus, append_headers, get_oncokb_annotation_column_headers
from filtered_maf_annotator import (
    CHUNK_ROWS, annotate_maf, count_maf, is_genome_kit_regions_file, main,
    maf_annotation_header, merge_maf)
from filtered_sv_annotator import main as sv_main
from backfill_filtered_maf import main as backfill_main


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
                    patch("AnnotatorCore.validate_oncokb_token") as validate, \
                    patch("AnnotatorCore.processalterationevents") as run:
                main()
                run.assert_not_called()
                validate.assert_not_called()
            with output.open(newline="") as file:
                rows = list(csv.reader(file, delimiter="\t"))
            annotation_header = maf_annotation_header()
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

    def fake_processalterationevents(self, sent):
        """Stands in for AnnotatorCore.processalterationevents: echoes the rows
        it is given (ASCII-folded, as upstream does) with ONCOGENIC=Unknown."""
        annotation_header = maf_annotation_header()

        def fake(eventfile, outfile, *args):
            with open(eventfile, newline="") as file:
                text = file.read()
            rows = list(csv.reader(text.splitlines(), delimiter="\t"))
            sent.append(text)
            suffix = [""] * len(annotation_header)
            suffix[annotation_header.index("ONCOGENIC")] = "Unknown"
            with open(outfile, "w") as file:
                for row in [rows[0] + annotation_header] + [r + suffix for r in rows[1:]]:
                    file.write("\t".join(row).encode("ascii", "ignore").decode("ascii") + "\n")
        return fake

    def test_panel_kit_sends_every_row_to_oncokb(self):
        annotation_header = maf_annotation_header()
        sent = []
        fake_annotator = self.fake_processalterationevents(sent)

        with tempfile.TemporaryDirectory() as temp:
            source, output = [Path(temp) / name for name in ("in.maf", "out.maf")]
            write_tsv(source, HEADER, ROWS)
            argv = ["filtered_maf_annotator.py", "--input", str(source), "--output",
                    str(output), "--reference", "GRCh38", "--regions-file", "/sample/panel.bed"]
            with patch.object(sys, "argv", argv), \
                    patch.dict(os.environ, {"ONCOKB_API_TOKEN": "t "}, clear=True), \
                    patch("AnnotatorCore.validate_oncokb_token") as validate, \
                    patch("AnnotatorCore.setoncokbapitoken") as set_token, \
                    patch("AnnotatorCore.processalterationevents", side_effect=fake_annotator) as run:
                main()
                validate.assert_called_once()
                set_token.assert_called_once_with("t")
                self.assertEqual(run.call_count, 1)
                # Genomic change, GRCh38, no tumour type, no hotspots/descriptions.
                args = run.call_args.args
                self.assertEqual(args[2:], ("", "", {}, False, AnnotatorCore.QueryType.GENOMIC_CHANGE,
                                            AnnotatorCore.ReferenceGenome.GRCH38, False))
            with source.open(newline="") as file:
                self.assertEqual(sent, [file.read()])
            with output.open(newline="") as file:
                rows = list(csv.reader(file, delimiter="\t"))
            self.assertEqual(rows[0], HEADER + annotation_header)
            self.assertEqual([row[:len(HEADER)] for row in rows[1:]], ROWS)
            self.assertEqual([row[rows[0].index("ONCOGENIC")] for row in rows[1:]], ["Unknown"] * 5)

    def test_count_maf_rejects_ragged_rows(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "in.maf"
            write_tsv(source, HEADER, ROWS + [ROWS[0][:-1]])
            with self.assertRaisesRegex(ValueError, "columns"):
                count_maf(source)

    def test_refuses_mismatched_oncokb_row(self):
        with tempfile.TemporaryDirectory() as temp:
            source, annotated, output = [Path(temp) / name for name in
                                         ("in.maf", "annotated.maf", "out.maf")]
            write_tsv(source, HEADER, ROWS)
            header, total = count_maf(source)
            annotation_header = maf_annotation_header()
            bad = ROWS[2].copy()
            bad[2] = "999"
            write_tsv(annotated, HEADER + annotation_header,
                      [bad + [""] * len(annotation_header)] * total)
            with self.assertRaisesRegex(ValueError, "out of order"):
                merge_maf(source, annotated, output, header)
            self.assertFalse(output.exists())

    def test_refuses_extra_oncokb_rows(self):
        with tempfile.TemporaryDirectory() as temp:
            source, annotated, output = [Path(temp) / name for name in
                                         ("in.maf", "annotated.maf", "out.maf")]
            write_tsv(source, HEADER, ROWS)
            header, _ = count_maf(source)
            annotation_header = maf_annotation_header()
            write_tsv(annotated, HEADER + annotation_header,
                      [row + [""] * len(annotation_header) for row in ROWS + ROWS[:1]])
            with self.assertRaisesRegex(ValueError, "more rows"):
                merge_maf(source, annotated, output, header)
            self.assertFalse(output.exists())

    def test_accepts_oncokb_ascii_stripped_echo_and_keeps_original_text(self):
        with tempfile.TemporaryDirectory() as temp:
            source, annotated, output = [Path(temp) / name for name in
                                         ("in.maf", "annotated.maf", "out.maf")]
            rows = [row.copy() for row in ROWS]
            rows[2][-1] = "Muir-Torré_syndrome"
            write_tsv(source, HEADER, rows)
            header, _ = count_maf(source)
            annotation_header = maf_annotation_header()
            # Mirror AnnotatorCore.append_annotation_to_file.
            with annotated.open("w") as file:
                file.write("\t".join(HEADER + annotation_header) + "\n")
                for row in rows:
                    line = "\t".join(row + [""] * len(annotation_header))
                    file.write(line.encode("ascii", "ignore").decode("ascii") + "\n")
            self.assertEqual(merge_maf(source, annotated, output, header), 5)
            with output.open(newline="") as file:
                out_rows = list(csv.reader(file, delimiter="\t"))
            self.assertEqual([row[:len(HEADER)] for row in out_rows[1:]], rows)

    def test_backfill_annotates_only_previously_skipped_rows(self):
        annotation_header = maf_annotation_header()
        oncogenic = annotation_header.index("ONCOGENIC")
        sent = []

        def fake_annotator(command, check):
            source = Path(command[command.index("-i") + 1])
            target = Path(command[command.index("-o") + 1])
            with source.open(newline="") as file:
                rows = list(csv.reader(file, delimiter="\t"))
            sent.extend(rows[1:])
            suffix = [""] * len(annotation_header)
            suffix[oncogenic] = "Unknown"
            # Mirror AnnotatorCore's ASCII stripping of the echoed input.
            with target.open("w") as file:
                for row in [rows[0] + annotation_header] + [r + suffix for r in rows[1:]]:
                    file.write("\t".join(row).encode("ascii", "ignore").decode("ascii") + "\n")

        rows = [row.copy() for row in ROWS]
        rows[1][-1] = "Café"
        filtered = []
        for row in rows:
            suffix = [""] * len(annotation_header)
            if row[HEADER.index("Func.refGene")] not in ("intergenic", "intronic"):
                suffix[oncogenic] = "Oncogenic"
            filtered.append(row + suffix)
        with tempfile.TemporaryDirectory() as temp:
            maf = Path(temp) / "S1_oncokb.maf"
            write_tsv(maf, HEADER + annotation_header, filtered)
            argv = ["backfill_filtered_maf.py", "--maf", str(maf), "--reference", "GRCh38"]
            with patch.object(sys, "argv", argv), \
                    patch.dict(os.environ, {"ONCOKB_API_TOKEN": "t"}, clear=True), \
                    patch("backfill_filtered_maf.subprocess.run", side_effect=fake_annotator):
                backfill_main()
            self.assertEqual(sent, rows[:2])
            with maf.open(newline="") as file:
                out = list(csv.reader(file, delimiter="\t"))
            self.assertEqual([row[:len(HEADER)] for row in out[1:]], rows)
            self.assertEqual([row[len(HEADER) + oncogenic] for row in out[1:]],
                             ["Unknown", "Unknown", "Oncogenic", "Oncogenic", "Oncogenic"])
            with maf.with_name(maf.name + ".filtered.bak").open(newline="") as file:
                self.assertEqual(list(csv.reader(file, delimiter="\t"))[1:], filtered)

    def test_maf_header_is_what_mafannotator_writes(self):
        # MafAnnotator.py -> processalterationevents_mixed writes its header
        # with append_headers(..., True, SOMATIC, combined_mode=True).
        out = io.StringIO()
        append_headers(out, 0, False, True, MutationStatus.SOMATIC, combined_mode=True)
        self.assertEqual(out.getvalue().rstrip("\n").split("\t"), maf_annotation_header())

    def test_maf_header_adds_germline_columns_after_oncogenic(self):
        v4 = maf_annotation_header()
        v3 = get_oncokb_annotation_column_headers(False, True)
        germline = ["PATHOGENIC (Germline)", "PENETRANCE (Germline)",
                    "GENOMIC_INDICATOR (Germline)"]
        at = v4.index("ONCOGENIC") + 1
        self.assertEqual(v4[at:at + 3], germline)
        self.assertEqual(v4[:at] + v4[at + 3:], v3)

    def test_backfill_refuses_a_v3_layout_maf(self):
        v3_header = get_oncokb_annotation_column_headers(False, True)
        with tempfile.TemporaryDirectory() as temp:
            maf = Path(temp) / "S1_oncokb.maf"
            write_tsv(maf, HEADER + v3_header,
                      [row + [""] * len(v3_header) for row in ROWS])
            argv = ["backfill_filtered_maf.py", "--maf", str(maf), "--reference", "GRCh38"]
            with patch.object(sys, "argv", argv), \
                    patch.dict(os.environ, {"ONCOKB_API_TOKEN": "t"}, clear=True), \
                    patch("backfill_filtered_maf.subprocess.run") as run:
                with self.assertRaisesRegex(ValueError, "OncoKB annotation columns"):
                    backfill_main()
                run.assert_not_called()


    def test_annotates_in_chunks_of_v3_post_size(self):
        self.assertEqual(CHUNK_ROWS, 100)
        sent = []
        rows = [row.copy() for row in ROWS]
        rows[1][-1] = "Café"
        with tempfile.TemporaryDirectory() as temp:
            source, annotated, output = [Path(temp) / n for n in ("in.maf", "a.maf", "out.maf")]
            write_tsv(source, HEADER, rows)
            with source.open(newline="") as file:
                lines = file.read().splitlines(keepends=True)
            with patch("AnnotatorCore.validate_oncokb_token"), \
                    patch("AnnotatorCore.processalterationevents",
                          side_effect=self.fake_processalterationevents(sent)) as run:
                self.assertEqual(annotate_maf(source, annotated, "t", "GRCh38", "MEL",
                                              chunk_rows=2), 5)
            self.assertEqual(run.call_count, 3)
            self.assertEqual(run.call_args.args[3], "MEL")
            # Each chunk is the header plus the next input lines, byte for byte.
            self.assertEqual(sent, [lines[0] + "".join(lines[1:3]), lines[0] + "".join(lines[3:5]),
                                    lines[0] + lines[5]])
            with annotated.open(newline="") as file:
                out = list(csv.reader(file, delimiter="\t"))
            self.assertEqual(out[0], HEADER + maf_annotation_header())
            self.assertEqual(len(out), 6)
            header, _ = count_maf(source)
            self.assertEqual(merge_maf(source, annotated, output, header), 5)
            with output.open(newline="") as file:
                merged = list(csv.reader(file, delimiter="\t"))
            self.assertEqual([row[:len(HEADER)] for row in merged[1:]], rows)

    def test_refuses_a_record_spanning_lines(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = [Path(temp) / name for name in ("in.maf", "out.maf")]
            rows = [row.copy() for row in ROWS]
            rows[0][-1] = "PASS\nsecond line"
            write_tsv(source, HEADER, rows)
            argv = ["filtered_maf_annotator.py", "--input", str(source), "--output",
                    str(output), "--reference", "GRCh38", "--regions-file", "/sample/panel.bed"]
            with patch.object(sys, "argv", argv), \
                    patch.dict(os.environ, {"ONCOKB_API_TOKEN": "t"}, clear=True), \
                    patch("AnnotatorCore.validate_oncokb_token"), \
                    patch("AnnotatorCore.processalterationevents",
                          side_effect=self.fake_processalterationevents([])):
                with self.assertRaisesRegex(ValueError, "records but"):
                    main()
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
