import csv
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import contextlib
import io
import subprocess

import AnnotatorCore
from AnnotatorCore import MutationStatus, append_headers, get_oncokb_annotation_column_headers
import filtered_maf_annotator
from filtered_maf_annotator import (
    CHUNK_ROWS, DEFAULT_REFHUB_URL, REFHUB_BAD, REFHUB_RETRY, OncoKBServerError, RefhubCache, _response,
    annotate_maf, count_maf, fail_fast_on_500, refhub_url_for, route_through_refhub,
    is_genome_kit_regions_file, main,
    maf_annotation_header, mask_output, merge_maf, run_masked)
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
                    patch("filtered_sv_annotator.run_masked") as run:
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

        def fake_annotator(command, check=False):
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
            return subprocess.CompletedProcess(command, 0)

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
                    patch("backfill_filtered_maf.run_masked", side_effect=lambda c: fake_annotator(c).returncode):
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
                    patch("backfill_filtered_maf.run_masked") as run:
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


    def test_sv_wrapper_reports_upstream_failure_without_the_token(self):
        with tempfile.TemporaryDirectory() as temp:
            source, output = [Path(temp) / name for name in ("sv.tsv", "annot.tsv")]
            write_tsv(source, ["Sample", "GeneA", "GeneB", "SVType"], [["S1", "BRAF", "BRAF", "DELETION"]])
            argv = ["filtered_sv_annotator.py", "--input", str(source), "--output", str(output)]
            err = io.StringIO()
            with patch.object(sys, "argv", argv), \
                    patch.dict(os.environ, {"ONCOKB_API_TOKEN": "secret-token"}, clear=True), \
                    patch("filtered_sv_annotator.run_masked", return_value=3), \
                    contextlib.redirect_stderr(err), contextlib.redirect_stdout(err):
                with self.assertRaises(SystemExit) as exit_:
                    sv_main()
            self.assertEqual(exit_.exception.code, 3)
            self.assertIn("StructuralVariantAnnotator.py exited with status 3", err.getvalue())
            self.assertNotIn("secret-token", err.getvalue())


    def test_upstream_output_is_masked(self):
        # Upstream prints the token-validation URL in connection errors and
        # urllib3 retry warnings.
        script = ("import sys; print('GET /api/v1/tokens/secret-token failed');"
                  "print('Traceback: secret-token', file=sys.stderr); sys.exit(4)")
        out = io.StringIO()
        with patch.dict(os.environ, {"ONCOKB_API_TOKEN": "secret-token "}), \
                patch.object(sys, "stdout", out):
            status = run_masked([sys.executable, "-c", script])
        self.assertEqual(status, 4)
        self.assertNotIn("secret-token", out.getvalue())
        self.assertIn("/api/v1/tokens/*** failed", out.getvalue())
        self.assertIn("Traceback: ***", out.getvalue())

    def test_own_output_is_masked(self):
        out, err = io.StringIO(), io.StringIO()
        with patch.dict(os.environ, {"ONCOKB_API_TOKEN": "secret-token"}), \
                patch.object(sys, "stdout", out), patch.object(sys, "stderr", err):
            mask_output()
            print("token secret-token in a message")
            print("and secret-token in a traceback", file=sys.stderr)
        self.assertEqual(out.getvalue(), "token *** in a message\n")
        self.assertEqual(err.getvalue(), "and *** in a traceback\n")

    def test_logging_is_masked(self):
        # AnnotatorCore configures logging on import, before mask_output runs.
        import logging
        err = io.StringIO()
        handler = logging.StreamHandler(err)
        logging.root.addHandler(handler)
        try:
            with patch.dict(os.environ, {"ONCOKB_API_TOKEN": "secret-token"}), \
                    patch.object(sys, "stdout", io.StringIO()), patch.object(sys, "stderr", err):
                mask_output()
                logging.getLogger("urllib3.connectionpool").warning(
                    "Retrying after connection broken: /api/v1/tokens/secret-token")
        finally:
            logging.root.removeHandler(handler)
        self.assertNotIn("secret-token", err.getvalue())
        self.assertIn("/api/v1/tokens/***", err.getvalue())

    def test_rejected_token_fails_the_maf_step(self):
        # validate_oncokb_token calls exit() on a rejected token: SystemExit(None).
        with tempfile.TemporaryDirectory() as temp:
            source, output = [Path(temp) / name for name in ("in.maf", "out.maf")]
            write_tsv(source, HEADER, ROWS)
            argv = ["filtered_maf_annotator.py", "--input", str(source), "--output",
                    str(output), "--reference", "GRCh38", "--regions-file", "/sample/panel.bed"]
            err = io.StringIO()
            with patch.object(sys, "argv", argv), \
                    patch.dict(os.environ, {"ONCOKB_API_TOKEN": "t"}, clear=True), \
                    patch("AnnotatorCore.validate_oncokb_token", side_effect=SystemExit(None)), \
                    contextlib.redirect_stderr(err):
                with self.assertRaises(SystemExit) as stop:
                    main()
            self.assertEqual(stop.exception.code, 1)
            self.assertIn("token missing, rejected or expired", err.getvalue())
            self.assertFalse(output.exists())


    @staticmethod
    def retry_error(status):
        """What upstream raises when OncoKB keeps answering `status`."""
        import requests
        from urllib3.exceptions import MaxRetryError, ResponseError
        reason = ResponseError(ResponseError.SPECIFIC_ERROR.format(status_code=status))
        return requests.exceptions.RetryError(MaxRetryError(None, "/annotate", reason))

    def poisoned_annotator(self, sent, poison, status=500):
        """fake_processalterationevents, except that a request holding a row
        whose Start_Position is in `poison` fails as upstream does. `poison`
        maps a position to how many requests holding it fail (None: all)."""
        fake = self.fake_processalterationevents(sent)
        failures = dict(poison) if isinstance(poison, dict) else {p: None for p in poison}

        def run(eventfile, outfile, *args):
            with open(eventfile, newline="") as file:
                rows = list(csv.reader(file, delimiter="\t"))[1:]
            hit = [row[2] for row in rows if row[2] in failures]
            if any(failures[p] is None or failures[p] > 0 for p in hit):
                for p in hit:
                    if failures[p]:
                        failures[p] -= 1
                raise self.retry_error(status)
            fake(eventfile, outfile, *args)
        return run

    def test_rows_oncokb_answers_with_500_are_isolated_and_marked(self):
        rows = [row.copy() for row in ROWS] + [
            ["S1", "chr1", str(n), str(n), "TTTT", "-", "Unknown", "Intron", "intronic", "PASS"]
            for n in range(6, 10)]
        sent, waits = [], []
        with tempfile.TemporaryDirectory() as temp:
            source, annotated, output = [Path(temp) / n for n in ("in.maf", "a.maf", "out.maf")]
            write_tsv(source, HEADER, rows)
            out = io.StringIO()
            with patch("AnnotatorCore.validate_oncokb_token"), \
                    patch("AnnotatorCore.processalterationevents",
                          side_effect=self.poisoned_annotator(sent, {"1", "7"})), \
                    contextlib.redirect_stdout(out):
                self.assertEqual(annotate_maf(source, annotated, "t", "GRCh38", chunk_rows=4,
                                              sleep=waits.append), 9)
            header, _ = count_maf(source)
            self.assertEqual(merge_maf(source, annotated, output, header), 9)
            with output.open(newline="") as file:
                merged = list(csv.reader(file, delimiter="\t"))
        names = maf_annotation_header()
        self.assertEqual([row[:len(HEADER)] for row in merged[1:]], rows)
        annotation = {row[2]: dict(zip(names, row[len(HEADER):])) for row in merged[1:]}
        for position in ("1", "7"):
            self.assertEqual(annotation[position]["ANNOTATED"], "False")
            self.assertEqual(annotation[position]["GENE_IN_ONCOKB"], "False")
            self.assertEqual(annotation[position]["ONCOGENIC"], "Unknown")
        for position in ("2", "3", "4", "5", "6", "8", "9"):
            self.assertEqual(annotation[position]["ONCOGENIC"], "Unknown")
            self.assertEqual(annotation[position]["ANNOTATED"], "")  # the fake leaves it empty
        # Every answered row was sent in a request that succeeded.
        answered = [line.split("\t")[2] for text in sent for line in text.splitlines()[1:]]
        self.assertEqual(sorted(answered, key=int), ["2", "3", "4", "5", "6", "8", "9"])
        self.assertIn("2 written unannotated", out.getvalue())
        # Each poisoned row is tried alone once in its chunk's pass and once after the pause.
        self.assertEqual(waits, [60])
        # Chunk 1-4 and 5-8 fail and go row by row; chunk 9 succeeds whole.
        self.assertEqual(sorted(len(text.splitlines()) - 1 for text in sent), [1] * 7)

    def test_a_coding_row_oncokb_cannot_answer_fails_the_step(self):
        with tempfile.TemporaryDirectory() as temp:
            source, annotated = [Path(temp) / n for n in ("in.maf", "a.maf")]
            write_tsv(source, HEADER, ROWS)
            with patch("AnnotatorCore.validate_oncokb_token"), \
                    patch("AnnotatorCore.processalterationevents",
                          side_effect=self.poisoned_annotator([], {"3"})), \
                    contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(OncoKBServerError, "exonic row"):
                    annotate_maf(source, annotated, "t", "GRCh38", sleep=lambda s: None)

    def test_a_row_oncokb_never_answers_is_handled_like_a_500(self):
        import requests
        from urllib3.exceptions import MaxRetryError, ReadTimeoutError
        timeout = requests.exceptions.ConnectionError(
            MaxRetryError(None, "/annotate", ReadTimeoutError(None, "/annotate", "Read timed out.")))
        fake = self.fake_processalterationevents([])

        def run(eventfile, outfile, *args):
            with open(eventfile, newline="") as file:
                rows = list(csv.reader(file, delimiter="\t"))[1:]
            if any(row[2] == "1" for row in rows):
                raise timeout
            fake(eventfile, outfile, *args)
        with tempfile.TemporaryDirectory() as temp:
            source, annotated = [Path(temp) / n for n in ("in.maf", "a.maf")]
            write_tsv(source, HEADER, ROWS)
            with patch("AnnotatorCore.validate_oncokb_token"), \
                    patch("AnnotatorCore.processalterationevents", side_effect=run), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(annotate_maf(source, annotated, "t", "GRCh38", sleep=lambda s: None), 5)
            with annotated.open(newline="") as file:
                out = list(csv.reader(file, delimiter="\t"))
        names = maf_annotation_header()
        self.assertEqual(dict(zip(names, out[1][len(HEADER):]))["ANNOTATED"], "False")

    def test_many_failing_rows_are_an_outage(self):
        rows = [["S1", "chr1", str(n), str(n), "TT", "-", "Unknown", "Intron", "intronic", "PASS"]
                for n in range(1, 41)]
        for failing, raises in ((10, False), (11, True)):
            with tempfile.TemporaryDirectory() as temp:
                source, annotated = [Path(temp) / n for n in ("in.maf", "a.maf")]
                write_tsv(source, HEADER, rows)
                with patch("AnnotatorCore.validate_oncokb_token"), \
                        patch("AnnotatorCore.processalterationevents",
                              side_effect=self.poisoned_annotator([], {str(n) for n in range(1, failing + 1)})), \
                        contextlib.redirect_stdout(io.StringIO()):
                    if raises:
                        with self.assertRaisesRegex(OncoKBServerError, "outage"):
                            annotate_maf(source, annotated, "t", "GRCh38", sleep=lambda s: None)
                    else:
                        self.assertEqual(annotate_maf(source, annotated, "t", "GRCh38",
                                                      sleep=lambda s: None), 40)

    def test_an_outage_stops_before_trying_every_row(self):
        rows = [["S1", "chr1", str(n), str(n), "TT", "-", "Unknown", "Intron", "intronic", "PASS"]
                for n in range(1, 101)]
        calls = []
        with tempfile.TemporaryDirectory() as temp:
            source, annotated = [Path(temp) / n for n in ("in.maf", "a.maf")]
            write_tsv(source, HEADER, rows)
            with patch("AnnotatorCore.validate_oncokb_token"), \
                    patch("AnnotatorCore.processalterationevents",
                          side_effect=lambda *a: calls.append(1) or (_ for _ in ()).throw(OncoKBServerError("500"))), \
                    contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(OncoKBServerError, "outage"):
                    annotate_maf(source, annotated, "t", "GRCh38", sleep=lambda s: None)
        # One chunk request, then single rows until more than 2 x 10 failed.
        self.assertEqual(len(calls), 1 + 21)

    def test_other_oncokb_errors_still_fail_the_step(self):
        import requests
        for error in (self.retry_error(503), self.retry_error(429),
                      requests.exceptions.ConnectionError("network")):
            with tempfile.TemporaryDirectory() as temp:
                source, annotated = [Path(temp) / n for n in ("in.maf", "a.maf")]
                write_tsv(source, HEADER, ROWS)
                with patch("AnnotatorCore.validate_oncokb_token"), \
                        patch("AnnotatorCore.processalterationevents", side_effect=error) as run:
                    with self.assertRaises(type(error)):
                        annotate_maf(source, annotated, "t", "GRCh38")
                self.assertEqual(run.call_count, 1)


    def test_a_row_that_fails_only_next_to_others_is_annotated_on_retry(self):
        # Position 3 is exonic and fails twice (in its chunk and alone right
        # after), then succeeds after the first pause: annotated, step passes.
        sent, waits = [], []
        with tempfile.TemporaryDirectory() as temp:
            source, annotated, output = [Path(temp) / n for n in ("in.maf", "a.maf", "out.maf")]
            write_tsv(source, HEADER, ROWS)
            with patch("AnnotatorCore.validate_oncokb_token"), \
                    patch("AnnotatorCore.processalterationevents",
                          side_effect=self.poisoned_annotator(sent, {"3": 2})), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(annotate_maf(source, annotated, "t", "GRCh38", sleep=waits.append), 5)
            header, _ = count_maf(source)
            self.assertEqual(merge_maf(source, annotated, output, header), 5)
            with output.open(newline="") as file:
                merged = list(csv.reader(file, delimiter="\t"))
        self.assertEqual(waits, [60])
        names = maf_annotation_header()
        self.assertEqual([dict(zip(names, row[len(HEADER):]))["ANNOTATED"] for row in merged[1:]],
                         [""] * 5)

    def test_upstream_requests_raise_on_500_without_retrying_it(self):
        import requests
        fail_fast_on_500()
        fail_fast_on_500()
        session = AnnotatorCore.requests_retry_session(allowed_methods=["POST"])
        retry = session.get_adapter("https://www.oncokb.org").max_retries
        self.assertNotIn(500, retry.status_forcelist)
        self.assertIn(503, retry.status_forcelist)
        self.assertEqual(retry.read, 0)
        self.assertEqual(retry.connect, 3)
        self.assertEqual(AnnotatorCore.REQUEST_TIMEOUT, 90)
        self.assertEqual(len(session.hooks["response"]), 1)
        response = requests.models.Response()
        response.status_code = 500
        response.url = "https://www.oncokb.org/api/v1/annotate?x=1"
        response.request = requests.Request("POST", response.url).prepare()
        with self.assertRaisesRegex(OncoKBServerError, "POST https://www.oncokb.org/api/v1/annotate$"):
            session.hooks["response"][0](response)
        response.status_code = 200
        session.hooks["response"][0](response)


ONCOKB_POST = "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange"
REFHUB = "https://refhub.test/api/v0/refhub"


def oncokb_answer(query):
    """An OncoKB genomic-change answer, enough for upstream to write a row."""
    return {"geneExist": True, "variantExist": False, "oncogenic": "Unknown",
            "mutationEffect": {"knownEffect": "Unknown", "description": "",
                               "citations": {"pmids": [], "abstracts": []}},
            "treatments": [], "query": {"hugoSymbol": "GENE1", "alteration": "X1Y",
                                        "consequence": "missense_variant",
                                        "genomicLocation": query["genomicLocation"]},
            "dataVersion": "v6.0", "lastUpdate": "11/01/2026"}


class FakeRefhub:
    """A requests session standing in for refhub. `reply(queries)` returns
    (status, payload) or an exception to raise."""

    def __init__(self, reply):
        self.reply, self.calls = reply, []

    def post(self, url, json=None, timeout=None, headers=None):
        self.calls.append({"url": url, "queries": json, "headers": headers, "timeout": timeout})
        result = self.reply(json)
        if isinstance(result, Exception):
            raise result
        status, payload = result
        return _response(url, status, payload)


def refhub_answers(bad=(), failed=(), unavailable=()):
    """refhub's reply: an annotation per query, except queries whose
    genomicLocation starts with a position in `bad`, `failed`, `unavailable`."""
    def reply(queries):
        answers = []
        for query in queries:
            position = query["genomicLocation"].split(",")[1]
            if position in bad:
                answers.append({"status": "bad", "reason": "http_500", "source": "hit"})
            elif position in failed:
                answers.append({"status": "failed", "reason": "http_500"})
            elif position in unavailable:
                answers.append({"status": "unavailable", "reason": "http_503"})
            else:
                answers.append({"status": "annotation", "source": "hit", "annotation": oncokb_answer(query)})
        return 200, answers
    return reply


class RefhubCacheTest(unittest.TestCase):
    def setUp(self):
        self.log = []

    def tearDown(self):
        filtered_maf_annotator._refhub = None

    def cache(self, reply):
        session = FakeRefhub(reply)
        return RefhubCache(REFHUB, "tok", session, log=self.log.append), session

    @staticmethod
    def queries(*positions):
        return [AnnotatorCore.GenomicChangeQuery("chr1", p, p, "A", "G", "", AnnotatorCore.ReferenceGenome.GRCH38)
                for p in positions]

    def test_answers_come_back_as_oncokb_would_send_them(self):
        cache, refhub = self.cache(refhub_answers())
        cache.path = "genomic"
        response = cache.post(ONCOKB_POST, self.queries("5", "6"))
        self.assertEqual(response.status_code, 200)
        self.assertEqual([a["query"]["genomicLocation"] for a in response.json()],
                         ["1,5,5,A,G", "1,6,6,A,G"])
        call, = refhub.calls
        self.assertEqual(call["url"], REFHUB + "/oncokb/cache/annotate/mutations/byGenomicChange")
        self.assertEqual(call["headers"], {"Authorization": "Bearer tok"})
        # The query objects upstream would post to OncoKB, unchanged.
        self.assertEqual(call["queries"][0],
                         {"genomicLocation": "1,5,5,A,G", "tumorType": "", "referenceGenome": "GRCh38"})
        self.assertEqual(self.log, [])

    def test_germline_goes_to_its_own_endpoint_and_others_skip_refhub(self):
        cache, refhub = self.cache(refhub_answers())
        cache.post("https://www.oncokb.org/api/v1/annotate/germline/mutations/byGenomicChange", self.queries("5"))
        self.assertTrue(refhub.calls[0]["url"].endswith("/oncokb/cache/annotate/germline/mutations/byGenomicChange"))
        self.assertIsNone(cache.post("https://www.oncokb.org/api/v1/annotate/copyNumberAlterations", []))
        self.assertEqual(len(refhub.calls), 1)

    def test_unreachable_refhub_falls_back_and_is_dropped_after_three_failures(self):
        import requests
        cache, refhub = self.cache(lambda q: requests.exceptions.ConnectionError("no route"))
        for _ in range(4):
            self.assertIsNone(cache.post(ONCOKB_POST, self.queries("5")))
        self.assertEqual(len(refhub.calls), 3)
        self.assertFalse(cache.enabled)
        self.assertIn("WARNING: refhub cache unreachable (ConnectionError); asking OncoKB directly", self.log)
        self.assertIn("for the rest of this run", self.log[-1])

    def test_refhub_errors_fall_back_and_a_success_resets_the_count(self):
        replies = iter([(503, {"error": "OncoKB cache is not configured"}), (200, "not a list"),
                        (200, [{"status": "annotation", "annotation": {}, "source": "hit"}]),
                        (502, "Bad Gateway"), (504, "Gateway Timeout")])
        cache, refhub = self.cache(lambda q: next(replies))
        self.assertIsNone(cache.post(ONCOKB_POST, self.queries("5")))
        self.assertIsNone(cache.post(ONCOKB_POST, self.queries("5")))
        self.assertEqual(cache.post(ONCOKB_POST, self.queries("5")).status_code, 200)
        self.assertIsNone(cache.post(ONCOKB_POST, self.queries("5")))
        self.assertIsNone(cache.post(ONCOKB_POST, self.queries("5")))
        self.assertTrue(cache.enabled)
        self.assertIn("HTTP 503", self.log[0])

    def test_a_rejected_token_stops_asking_refhub(self):
        cache, refhub = self.cache(lambda q: (401, {"error": "send refhub's OncoKB token"}))
        self.assertIsNone(cache.post(ONCOKB_POST, self.queries("5")))
        self.assertIsNone(cache.post(ONCOKB_POST, self.queries("5")))
        self.assertEqual(len(refhub.calls), 1)

    def test_oncokb_unavailable_to_refhub_marks_rows_for_later_not_direct(self):
        cache, _ = self.cache(refhub_answers(unavailable={"5", "6"}))
        answers = cache.post(ONCOKB_POST, self.queries("5", "6")).json()
        self.assertEqual([a["query"]["consequence"] for a in answers], [REFHUB_RETRY, REFHUB_RETRY])
        self.assertIn("could not reach OncoKB (http_503); trying these rows again later", self.log[0])
        self.assertEqual(cache.counts["direct"], 0)

    def test_somatic_queries_ask_refhub_genomic_first_by_default(self):
        cache, refhub = self.cache(refhub_answers())
        with patch.dict(os.environ, {}, clear=True):
            RefhubCache(REFHUB, "tok", refhub).post(ONCOKB_POST, self.queries("5"))
            RefhubCache(REFHUB, "tok", refhub).post(
                "https://www.oncokb.org/api/v1/annotate/germline/mutations/byGenomicChange", self.queries("5"))
        with patch.dict(os.environ, {"REFHUB_ONCOKB_PATH": "genomic"}):
            RefhubCache(REFHUB, "tok", refhub).post(ONCOKB_POST, self.queries("5"))
        with patch.dict(os.environ, {"REFHUB_ONCOKB_PATH": "protein"}):
            RefhubCache(REFHUB, "tok", refhub).post(ONCOKB_POST, self.queries("5"))
        urls = [c["url"] for c in refhub.calls]
        self.assertEqual(urls, [
            REFHUB + "/oncokb/cache/annotate/mutations/byGenomicChange?path=genomic_first",
            REFHUB + "/oncokb/cache/annotate/germline/mutations/byGenomicChange",
            REFHUB + "/oncokb/cache/annotate/mutations/byGenomicChange",
            REFHUB + "/oncokb/cache/annotate/mutations/byGenomicChange?path=protein",
        ])

    def test_protein_route_answers_are_counted_for_the_log_and_passed_on_as_oncokbs(self):
        def answers(queries):
            out = []
            for q in queries:
                position = q["genomicLocation"].split(",")[1]
                annotation = {"geneExist": True, "oncogenic": "Likely Oncogenic", "treatments": [],
                              "query": {"genomicLocation": q["genomicLocation"], "consequence": "missense_variant"}}
                answer = {"status": "annotation", "annotation": annotation, "source": "miss"}
                if position == "6":
                    answer["route"] = "protein"
                out.append(answer)
            return 200, out
        cache, _ = self.cache(answers)
        reply = cache.post(ONCOKB_POST, self.queries("5", "6")).json()
        # What upstream gets is OncoKB's annotation only: no route anywhere.
        self.assertNotIn("route", json.dumps(reply))
        self.assertEqual([a["oncogenic"] for a in reply], ["Likely Oncogenic"] * 2)
        self.assertEqual(cache.counts["protein_route"], 1)
        self.assertIn("1 answered by protein change", cache.summary())

    def test_refhub_without_genomic_first_is_asked_by_genomic_change(self):
        replies = iter([(400, {"error": "path is genomic or protein"}),
                        (200, [{"status": "annotation", "annotation": {}, "source": "hit"}])])
        cache, refhub = self.cache(lambda q: next(replies))
        cache.path = "genomic_first"
        self.assertEqual(cache.post(ONCOKB_POST, self.queries("5")).status_code, 200)
        urls = [c["url"] for c in refhub.calls]
        self.assertEqual(urls, [
            REFHUB + "/oncokb/cache/annotate/mutations/byGenomicChange?path=genomic_first",
            REFHUB + "/oncokb/cache/annotate/mutations/byGenomicChange",
        ])
        self.assertEqual(cache.counts["direct"], 0)

    def test_rows_refhub_cannot_annotate_come_back_marked_in_one_answer(self):
        cache, _ = self.cache(refhub_answers(bad={"6"}, failed={"7"}, unavailable={"8"}))
        answers = cache.post(ONCOKB_POST, self.queries("5", "6", "7", "8")).json()
        self.assertEqual([a["query"]["consequence"] for a in answers],
                         ["missense_variant", REFHUB_BAD, REFHUB_RETRY, REFHUB_RETRY])
        self.assertEqual((cache.counts["hit"], cache.counts["bad"], cache.counts["failed"]), (1, 1, 2))

    def test_direct_fallback_sends_oncokb_at_most_100_queries_a_request(self):
        import requests
        route_through_refhub()
        cache, _ = self.cache(lambda q: requests.exceptions.ConnectionError("no route"))
        filtered_maf_annotator._refhub = cache
        sizes = []

        def direct(url, body):
            sizes.append(len(body))
            return _response(url, 200, [{"n": i} for i in range(len(body))])
        with patch("filtered_maf_annotator._direct_post", side_effect=direct):
            response = AnnotatorCore.makeoncokbpostrequest(ONCOKB_POST, self.queries(*map(str, range(250))))
        self.assertEqual(sizes, [100, 100, 50])
        self.assertEqual(len(response.json()), 250)

    def test_upstream_posts_fall_back_to_oncokb_when_refhub_cannot_answer(self):
        import requests
        route_through_refhub()
        route_through_refhub()
        cache, _ = self.cache(lambda q: requests.exceptions.ReadTimeout("slow"))
        filtered_maf_annotator._refhub = cache
        with patch("filtered_maf_annotator._direct_post", return_value="direct") as direct:
            self.assertEqual(AnnotatorCore.makeoncokbpostrequest(ONCOKB_POST, self.queries("5")), "direct")
            filtered_maf_annotator._refhub = None
            self.assertEqual(AnnotatorCore.makeoncokbpostrequest(ONCOKB_POST, []), "direct")
        self.assertEqual(direct.call_count, 2)

    def test_refhub_url(self):
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual(refhub_url_for(None), DEFAULT_REFHUB_URL)
            self.assertIsNone(refhub_url_for("https://demo.oncokb.org"))
        with patch.dict(os.environ, {"REFHUB_URL": " "}):
            self.assertIsNone(refhub_url_for(None))
        with patch.dict(os.environ, {"REFHUB_URL": REFHUB}):
            self.assertEqual(refhub_url_for(None), REFHUB)

    def run_upstream(self, rows, reply, sleep=None):
        """annotate_maf with upstream's real processalterationevents, refhub
        faked and OncoKB itself unreachable. Returns (rows by position, refhub
        calls, waits, out); refhub chunks of 3 rows."""
        session = FakeRefhub(reply)
        waits = []
        with tempfile.TemporaryDirectory() as temp:
            source, annotated, output = [Path(temp) / n for n in ("in.maf", "a.maf", "out.maf")]
            write_tsv(source, HEADER, rows)
            out = io.StringIO()
            with patch("AnnotatorCore.validate_oncokb_token"), \
                    patch("filtered_maf_annotator._direct_post",
                          side_effect=AssertionError("OncoKB asked directly")), \
                    contextlib.redirect_stdout(out):
                annotate_maf(source, annotated, "tok", "GRCh38", refhub_url=REFHUB,
                             refhub_session=session, refhub_chunk_rows=3, sleep=waits.append)
            header, _ = count_maf(source)
            merge_maf(source, annotated, output, header)
            with output.open(newline="") as file:
                merged = list(csv.reader(file, delimiter="\t"))
        names = maf_annotation_header()
        by_position = {row[2]: dict(zip(names, row[len(HEADER):])) for row in merged[1:]}
        return by_position, session.calls, waits, out.getvalue()

    def test_upstream_annotates_through_refhub_and_known_bad_rows_are_not_retried(self):
        rows = [row.copy() for row in ROWS] + [
            ["S1", "chr1", "119510721", "119510740", "TTTTTTTTTTTTTTTTTTTT", "-", "NOTCH2", "Intron",
             "intronic", "PASS"]]
        by_position, calls, waits, out = self.run_upstream(rows, refhub_answers(bad={"119510721"}))
        for position in ("1", "2", "3", "4", "5"):
            self.assertEqual(by_position[position]["ANNOTATED"], "True")
            self.assertEqual(by_position[position]["GENE_IN_ONCOKB"], "True")
        bad = by_position["119510721"]
        self.assertEqual((bad["ANNOTATED"], bad["GENE_IN_ONCOKB"], bad["ONCOGENIC"]), ("False", "False", "Unknown"))
        # Chunks 1-3 and 4-6, each once: the bad row needs no resend, no pause.
        self.assertEqual([len(c["queries"]) for c in calls], [3, 3])
        self.assertEqual(waits, [])
        self.assertIn("known bad in refhub's cache, written unannotated", out)
        self.assertIn("refhub cache: 2 requests; answers: 5 hits, 0 asked of OncoKB, 0 stale, 1 known bad, "
                      "0 not answered; 0 answered by protein change (OncoKB failed them by genomic "
                      "change); 0 requests sent to OncoKB directly", out)

    def test_known_bad_rows_are_no_outage(self):
        # 40 known-bad intronic rows of 45: more than max_unanswered (10).
        rows = [row.copy() for row in ROWS] + [
            ["S1", "chr1", str(n), str(n), "TTTT", "-", "NOTCH2", "Intron", "intronic", "PASS"]
            for n in range(100, 140)]
        by_position, calls, waits, out = self.run_upstream(rows, refhub_answers(bad={str(n) for n in range(100, 140)}))
        self.assertEqual(sum(1 for r in by_position.values() if r["ANNOTATED"] == "False"), 40)
        self.assertEqual(len(calls), 15)
        self.assertEqual(waits, [])

    def test_a_chunk_refhub_cannot_answer_is_retried_together_not_row_by_row(self):
        # OncoKB down for refhub on the first pass (all unavailable), then back.
        replies = {"n": 0}

        def reply(queries):
            replies["n"] += 1
            return refhub_answers(unavailable={r[2] for r in ROWS} if replies["n"] <= 2 else set())(queries)
        by_position, calls, waits, out = self.run_upstream(ROWS, reply)
        self.assertEqual([len(c["queries"]) for c in calls], [3, 2, 3, 2])
        self.assertEqual(waits, [60])
        self.assertTrue(all(r["ANNOTATED"] == "True" for r in by_position.values()))

    def test_a_failed_direct_fallback_keeps_the_rows_for_refhub(self):
        # refhub answers HTTP 502 once; the direct fallback fails (500).
        replies = {"n": 0}

        def reply(queries):
            replies["n"] += 1
            return (502, "Bad Gateway") if replies["n"] == 1 else refhub_answers()(queries)
        session = FakeRefhub(reply)
        waits = []
        with tempfile.TemporaryDirectory() as temp:
            source, annotated = [Path(temp) / n for n in ("in.maf", "a.maf")]
            write_tsv(source, HEADER, ROWS)
            with patch("AnnotatorCore.validate_oncokb_token"), \
                    patch("filtered_maf_annotator._direct_post", side_effect=OncoKBServerError("HTTP 500")) as direct, \
                    contextlib.redirect_stdout(io.StringIO()):
                annotate_maf(source, annotated, "tok", "GRCh38", refhub_url=REFHUB,
                             refhub_session=session, refhub_chunk_rows=3, sleep=waits.append)
        # One direct try for chunk 1-3; its rows went again with chunk 4-5's
        # leftovers after the pause, never one by one.
        self.assertEqual(direct.call_count, 1)
        self.assertEqual([len(c["queries"]) for c in session.calls], [3, 2, 3])
        self.assertEqual(waits, [60])

    def test_many_rows_refhub_could_not_answer_are_still_an_outage(self):
        rows = [["S1", "chr1", str(n), str(n), "TTTT", "-", "NOTCH2", "Intron", "intronic", "PASS"]
                for n in range(100, 140)]
        with self.assertRaisesRegex(OncoKBServerError, "outage"):
            self.run_upstream(rows, refhub_answers(failed={str(n) for n in range(100, 140)}))

    def test_a_coding_row_refhub_knows_is_bad_fails_the_step(self):
        with self.assertRaisesRegex(OncoKBServerError, "exonic row"):
            self.run_upstream(ROWS, refhub_answers(bad={"3"}))

    def test_a_row_refhub_could_not_confirm_is_tried_again_after_the_pause(self):
        replies = {"n": 0}

        def reply(queries):
            replies["n"] += 1
            return refhub_answers(failed={"2"} if replies["n"] <= 2 else set())(queries)
        by_position, calls, waits, _ = self.run_upstream(ROWS, reply)
        self.assertEqual(waits, [60])
        self.assertEqual(by_position["2"]["ANNOTATED"], "True")
        # Chunks 1-3 and 4-5; then row 2 again with the other unanswered rows.
        self.assertEqual([len(c["queries"]) for c in calls], [3, 2, 1])


if __name__ == "__main__":
    unittest.main()
