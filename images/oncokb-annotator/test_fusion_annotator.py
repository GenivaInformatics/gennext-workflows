"""fusion_annotator.py: upstream FusionAnnotator.py's output, unchanged, plus
ONCOKB_MATCHED_ALTERATIONS and ONCOKB_VARIANT_SUMMARY. OncoKB is faked;
upstream's real AnnotatorCore and FusionAnnotator.py run."""
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import AnnotatorCore
import fusion_annotator
from filtered_maf_annotator import _response


def answer(oncogenic, variant_exist, summary, treatments=(), dx=()):
    return {
        "geneExist": True, "variantExist": variant_exist, "oncogenic": oncogenic,
        "mutationEffect": {"knownEffect": "Gain-of-function", "description": "",
                           "citations": {"pmids": [], "abstracts": []}},
        "treatments": [{"level": level, "drugs": [{"drugName": drug}], "alterations": list(alts),
                        "pmids": ["1"], "abstracts": []} for level, drug, alts in treatments],
        "diagnosticImplications": [{"levelOfEvidence": level, "alterations": list(alts),
                                    "tumorType": {"name": "Lung Adenocarcinoma", "code": "LUAD",
                                                  "mainType": {"name": "Non-Small Cell Lung Cancer"}},
                                    "pmids": [], "abstracts": []} for level, alts in dx],
        "prognosticImplications": [],
        "highestDiagnosticImplicationLevel": dx[0][0] if dx else None,
        "highestPrognosticImplicationLevel": None,
        "variantSummary": summary,
    }


ANSWERS = {
    ("EML4", "ALK"): answer("Oncogenic", True, "The EML4-ALK fusion is known to be oncogenic.",
                            treatments=[("LEVEL_1", "Alectinib", ["EML4-ALK Fusion"]),
                                        ("LEVEL_3A", "Neladalkib", ["Fusions", "EML4-ALK Fusion"])]),
    ("NTRK2", "MYO9A"): answer("Likely Oncogenic", False,
                               "The NTRK2-MYO9A fusion has not specifically been reviewed by the OncoKB team.\n"
                               "However, NTRK2 fusions are likely oncogenic.",
                               treatments=[("LEVEL_1", "Larotrectinib", ["Fusions"])]),
    ("MSMB", "NCOA4"): answer("Unknown", False, "", dx=[("LEVEL_Dx2", ["Fusions"])]),
    ("MSH2", "KCNK12"): answer("Likely Oncogenic", True,
                               "The MSH2-KCNK12 fusion is likely oncogenic."),
}


def fake_oncokb(url, body):
    queries = json.loads(json.dumps(body, default=lambda o: o.__dict__))
    out = []
    for q in queries:
        pair = (q["geneA"]["hugoSymbol"], q["geneB"]["hugoSymbol"])
        out.append(ANSWERS.get(pair))
    return _response(url, 200, out)


def fake_get(url):
    # The per-query GET after a failed batch: not used here.
    return _response(url, 500, {})


class FusionAnnotatorTest(unittest.TestCase):
    def annotate(self, with_columns, rows, post=fake_oncokb):
        with tempfile.TemporaryDirectory() as temp:
            source, target = Path(temp) / "in.tsv", Path(temp) / "out.tsv"
            with source.open("w") as f:
                f.write("Tumor_Sample_Barcode\tFusion\n")
                for fusion in rows:
                    f.write(f"S1\t{fusion}\n")
            argv = ["-i", str(source), "-o", str(target), "-b", "tok", "-t", "LUAD"]
            with patch("AnnotatorCore.validate_oncokb_token"), \
                    patch("AnnotatorCore.makeoncokbpostrequest", side_effect=post), \
                    patch("AnnotatorCore.makeoncokbgetrequest", side_effect=fake_get):
                if with_columns:
                    fusion_annotator.run(argv)
                else:
                    import runpy
                    import sys
                    script = str(Path(AnnotatorCore.__file__).with_name("FusionAnnotator.py"))
                    saved, sys.argv = sys.argv, [script] + argv
                    try:
                        runpy.run_path(script, run_name="__main__")
                    finally:
                        sys.argv = saved
            with target.open(newline="") as f:
                return list(csv.reader(f, delimiter="\t"))

    def test_upstream_columns_are_unchanged_and_two_are_added(self):
        rows = ["EML4-ALK", "NTRK2-MYO9A", "MSMB-NCOA4", "MSH2-KCNK12"]
        upstream = self.annotate(False, rows)
        ours = self.annotate(True, rows)
        self.assertEqual(ours[0], upstream[0] + fusion_annotator.EXTRA_HEADERS)
        width = len(upstream[0])
        for a, b in zip(upstream[1:], ours[1:]):
            self.assertEqual(b[:width], a)
        by_fusion = {r[1]: dict(zip(ours[0], r)) for r in ours[1:]}
        self.assertEqual(by_fusion["EML4-ALK"]["ONCOKB_MATCHED_ALTERATIONS"], "EML4-ALK Fusion;Fusions")
        self.assertEqual(by_fusion["NTRK2-MYO9A"]["ONCOKB_MATCHED_ALTERATIONS"], "Fusions")
        self.assertEqual(by_fusion["MSMB-NCOA4"]["ONCOKB_MATCHED_ALTERATIONS"], "Fusions")
        # No treatment or implication: OncoKB names no matched alteration.
        self.assertEqual(by_fusion["MSH2-KCNK12"]["ONCOKB_MATCHED_ALTERATIONS"], "")
        self.assertEqual(by_fusion["NTRK2-MYO9A"]["ONCOKB_VARIANT_SUMMARY"],
                         "The NTRK2-MYO9A fusion has not specifically been reviewed by the OncoKB team. "
                         "However, NTRK2 fusions are likely oncogenic.")
        self.assertEqual(by_fusion["EML4-ALK"]["HIGHEST_LEVEL"], "LEVEL_1")

    def test_an_unanswered_query_keeps_its_columns_in_place(self):
        # The batch fails, then each GET fails too: upstream writes ANNOTATED=False
        # and pads; the extra columns are empty and in their own place.
        def failing(url, body):
            return _response(url, 500, {})
        ours = self.annotate(True, ["EML4-ALK"], post=failing)
        row = dict(zip(ours[0], ours[1]))
        self.assertEqual(len(ours[1]), len(ours[0]))
        self.assertEqual(row["ANNOTATED"], "False")
        self.assertEqual((row["ONCOKB_MATCHED_ALTERATIONS"], row["ONCOKB_VARIANT_SUMMARY"]), ("", ""))

    def test_the_patch_is_undone(self):
        headers = AnnotatorCore.get_oncokb_annotation_column_headers
        pull = AnnotatorCore.pull_structural_variant_info
        self.annotate(True, ["EML4-ALK"])
        self.assertIs(AnnotatorCore.get_oncokb_annotation_column_headers, headers)
        self.assertIs(AnnotatorCore.pull_structural_variant_info, pull)

    def test_matched_alterations_reads_every_implication_kind(self):
        a = {"treatments": [{"alterations": ["B"]}], "diagnosticImplications": [{"alterations": ["A", "B"]}],
             "prognosticImplications": [{"alterations": [{"alteration": "C"}]}]}
        self.assertEqual(fusion_annotator.matched_alterations(a), "A;B;C")
        self.assertEqual(fusion_annotator.matched_alterations(None), "")
        self.assertEqual(fusion_annotator.matched_alterations({"treatments": None}), "")


if __name__ == "__main__":
    unittest.main()
