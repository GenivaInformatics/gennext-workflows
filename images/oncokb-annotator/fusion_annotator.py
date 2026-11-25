#!/usr/bin/env python3
"""Upstream FusionAnnotator.py, unchanged, plus two columns from the same
OncoKB answers (no extra calls):

- MATCHED_ALTERATIONS: the alterations OncoKB matched, as it names
  them in treatments[].alterations, diagnosticImplications[].alterations and
  prognosticImplications[].alterations, e.g. "EML4-ALK Fusion" for a curated
  pair or "Fusions" for a gene-level match; ";"-joined, sorted, empty when
  the answer has no treatment or implication (OncoKB names no matched
  alteration then).
- VARIANT_SUMMARY: OncoKB's variantSummary sentence, on one line. With
  -d (include descriptions) upstream writes its own VARIANT_SUMMARY, so only
  MATCHED_ALTERATIONS is added then.

Both come after upstream's own columns; those, and every row's values in
them, are what FusionAnnotator.py writes. Same arguments as
FusionAnnotator.py.
"""
import os
import runpy
import sys

import AnnotatorCore

EXTRA_HEADERS = ["MATCHED_ALTERATIONS", "VARIANT_SUMMARY"]


def extra_headers(include_descriptions):
    """The columns added: upstream's -d already writes VARIANT_SUMMARY."""
    return EXTRA_HEADERS[:1] if include_descriptions else EXTRA_HEADERS


_IMPLICATIONS = ("treatments", "diagnosticImplications", "prognosticImplications")


def matched_alterations(answer):
    """OncoKB's names for the alterations it matched, ";"-joined."""
    if not isinstance(answer, dict):
        return ""
    names = set()
    for key in _IMPLICATIONS:
        for item in answer.get(key) or []:
            for alteration in (item or {}).get("alterations") or []:
                if isinstance(alteration, dict):
                    alteration = alteration.get("alteration") or alteration.get("name") or ""
                alteration = str(alteration).strip()
                if alteration:
                    names.add(alteration)
    return ";".join(sorted(names))


def variant_summary(answer):
    if not isinstance(answer, dict):
        return ""
    text = answer.get("variantSummary") or ""
    return " ".join(str(text).split())


def install():
    """Adds the two columns to upstream's fusion output. Returns a function
    that undoes it."""
    headers = AnnotatorCore.get_oncokb_annotation_column_headers
    pull = AnnotatorCore.pull_structural_variant_info

    def headers_with_extra(include_descriptions, genomic_change_annotation, *args, **kwargs):
        return (headers(include_descriptions, genomic_change_annotation, *args, **kwargs)
                + extra_headers(include_descriptions))

    def pull_with_extra(queries, include_descriptions):
        raw = []
        process = AnnotatorCore.process_oncokb_annotation

        def capture(annotation, *args, **kwargs):
            raw.append(annotation)
            return process(annotation, *args, **kwargs)

        AnnotatorCore.process_oncokb_annotation = capture
        try:
            processed = pull(queries, include_descriptions)
        finally:
            AnnotatorCore.process_oncokb_annotation = process
        # An unanswered query gives a short list (upstream pads it when
        # writing): pad it first, so the extra columns land in their place.
        width = len(headers(include_descriptions, False))
        extra = len(extra_headers(include_descriptions))
        return [list(columns) + [""] * (width - len(columns))
                + [matched_alterations(answer), variant_summary(answer)][:extra]
                for columns, answer in zip(processed, raw)]

    AnnotatorCore.get_oncokb_annotation_column_headers = headers_with_extra
    AnnotatorCore.pull_structural_variant_info = pull_with_extra

    def restore():
        AnnotatorCore.get_oncokb_annotation_column_headers = headers
        AnnotatorCore.pull_structural_variant_info = pull
    return restore


def run(argv):
    """FusionAnnotator.py with these arguments, with the two columns."""
    script = os.path.join(os.path.dirname(os.path.abspath(AnnotatorCore.__file__)), "FusionAnnotator.py")
    restore = install()
    saved = sys.argv
    sys.argv = [script] + list(argv)
    try:
        runpy.run_path(script, run_name="__main__")
    finally:
        sys.argv = saved
        restore()


if __name__ == "__main__":
    run(sys.argv[1:])
