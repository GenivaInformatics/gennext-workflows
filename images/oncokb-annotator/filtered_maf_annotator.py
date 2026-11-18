#!/usr/bin/env python3
"""Annotate every MAF row with OncoKB (Genome kit: skip OncoKB), retaining input text."""

import argparse
import csv
import logging
import os
from pathlib import Path
import subprocess
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

# What a row OncoKB cannot answer gets. ANNOTATED=False marks it; the other
# values are the ones OncoKB returns for a variant it has no knowledge of, so
# the columns gennext-api-edge requires (GENE_IN_ONCOKB, ONCOGENIC) are set.
UNANSWERED = {
    "ANNOTATED": "False",
    "GENE_IN_ONCOKB": "False",
    "VARIANT_IN_ONCOKB": "False",
    "MUTATION_EFFECT": "Unknown",
    "ONCOGENIC": "Unknown",
}

# Func.refGene values OncoKB must answer. Its knowledge is about these rows, so
# the step fails rather than leave one unannotated.
MUST_ANNOTATE_FUNC = {"exonic", "splicing"}


class OncoKBServerError(Exception):
    """OncoKB answered a request with HTTP 500 after upstream's retries."""


def is_oncokb_500(error):
    """True when upstream gave up on HTTP 500 responses.

    Only 500: OncoKB answers 500 to some queries it cannot process, such as
    indels in long homopolymers. 429, 502, 503 and 504 are the service, not
    the rows, and connection errors are the network; those still fail the step.
    """
    import requests
    from urllib3.exceptions import MaxRetryError, ResponseError

    if not isinstance(error, requests.exceptions.RetryError):
        return False
    cause = error.args[0] if error.args else None
    return (isinstance(cause, MaxRetryError) and isinstance(cause.reason, ResponseError)
            and ResponseError.SPECIFIC_ERROR.format(status_code=500) in str(cause.reason))


def unanswered_suffix():
    return [UNANSWERED.get(name, "") for name in maf_annotation_header()]


def annotate_maf(input_path, output_path, token, reference, tumor_type=".",
                 api_url=None, chunk_rows=CHUNK_ROWS):
    """MafAnnotator.py -q Genomic_Change, run chunk by chunk. Returns the rows sent.

    Does what MafAnnotator.main does for our options, once: set the API URL and
    token, validate the token. Then it calls upstream's processalterationevents
    on consecutive chunks of the input and joins their output. Each chunk is
    the input's header line plus a run of its data lines, copied byte for
    byte, so upstream reads the same text as from the whole file. The token
    is never put on a command line.

    One query OncoKB answers with HTTP 500 fails its whole request. A chunk
    that gets a 500 is split in halves until the rows that cause it stand
    alone; those get UNANSWERED values, every other row its OncoKB answer.
    A coding or splice-site row (MUST_ANNOTATE_FUNC) is never left
    unanswered: the step fails instead.
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
        header_fields = []
        state = {"sent": 0, "chunks": 0, "header_written": False, "unanswered": 0}

        def write_header(text):
            if not state["header_written"]:
                target.write(text)
                state["header_written"] = True

        def unanswered(line):
            fields = next(csv.reader([line], delimiter="\t"))
            func_column = header_fields.index("Func.refGene") if "Func.refGene" in header_fields else None
            func = fields[func_column].split(";")[0] if func_column is not None else ""
            location = "\t".join(fields[1:6])
            if func in MUST_ANNOTATE_FUNC:
                raise OncoKBServerError(
                    f"OncoKB answers HTTP 500 for a {func} row ({location}); not leaving it unannotated")
            print(f"OncoKB: HTTP 500 for one row, written unannotated: {location} ({func or 'no Func.refGene'})",
                  flush=True)
            write_header(header.rstrip("\r\n") + "\t" + "\t".join(maf_annotation_header()) + "\n")
            echo = line.rstrip("\r\n").encode("ascii", "ignore").decode("ascii")
            target.write(echo + "\t" + "\t".join(unanswered_suffix()) + "\n")
            state["unanswered"] += 1
            state["sent"] += 1

        def annotate(chunk):
            with chunk_in.open("w", newline="") as part:
                part.write(header)
                part.writelines(chunk)
            try:
                AnnotatorCore.processalterationevents(
                    str(chunk_in), str(chunk_out), "", cancer_type, {}, False,
                    query_type, reference_genome, False)
            except Exception as error:
                if not is_oncokb_500(error):
                    raise
                if len(chunk) == 1:
                    unanswered(chunk[0])
                else:
                    half = len(chunk) // 2
                    annotate(chunk[:half])
                    annotate(chunk[half:])
                return
            with chunk_out.open(newline="") as part:
                write_header(part.readline())
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
                    header_fields = next(csv.reader([header.rstrip("\r\n")], delimiter="\t"))
                continue
            chunk.append(line if line.endswith("\n") else line + "\n")
            if len(chunk) == chunk_rows:
                annotate(chunk)
                chunk = []
        if chunk:
            annotate(chunk)
        if state["unanswered"]:
            print(f"OncoKB: {state['unanswered']} of {state['sent']} rows got HTTP 500 and are written "
                  "unannotated (ANNOTATED=False)", flush=True)
        return state["sent"]


TOKEN_MASK = "***"


def _token():
    return os.environ.get("ONCOKB_API_TOKEN", "").strip()


def mask_token(text, token=None):
    """The text with the OncoKB token replaced by ***."""
    token = _token() if token is None else token
    return text.replace(token, TOKEN_MASK) if token else text


class _MaskingStream:
    def __init__(self, stream, token):
        self._stream, self._token = stream, token

    def write(self, text):
        return self._stream.write(mask_token(text, self._token))

    def __getattr__(self, name):
        return getattr(self._stream, name)


def mask_output():
    """Masks the token in everything this process prints from here on.

    The token reaches output in more ways than our own code: upstream validates
    it with GET /tokens/<token>, so a connection error, urllib3's retry
    warnings and a traceback all print that URL. Logging handlers already
    writing to stdout or stderr are moved to the masked streams.
    """
    token = _token()
    if not token or isinstance(sys.stdout, _MaskingStream):
        return
    original = {id(sys.stdout): None, id(sys.stderr): None}
    sys.stdout = _MaskingStream(sys.stdout, token)
    sys.stderr = _MaskingStream(sys.stderr, token)
    original[id(sys.stdout._stream)] = sys.stdout
    original[id(sys.stderr._stream)] = sys.stderr
    # AnnotatorCore calls logging.basicConfig when imported, so the root
    # handler already holds the unmasked stream: point it at the masked one.
    for handler in logging.root.handlers:
        stream = getattr(handler, "stream", None)
        if isinstance(handler, logging.StreamHandler) and original.get(id(stream)):
            handler.setStream(original[id(stream)])


def run_masked(command):
    """Runs an upstream script with its output, masked, on ours. Returns its exit status.

    Upstream's scripts take the token as -b, so the command itself is never
    printed either.
    """
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, bufsize=1) as process:
        for line in process.stdout:
            sys.stdout.write(mask_token(line))
            sys.stdout.flush()
    return process.returncode


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
    mask_output()
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
            try:
                sent = annotate_maf(args.input, annotated, os.environ["ONCOKB_API_TOKEN"].strip(),
                                    args.reference, args.tumor_type, args.api_url)
            except SystemExit as stop:
                # Upstream's validate_oncokb_token calls exit() when it rejects
                # the token, which is status 0 when called in-process.
                if stop.code in (None, 0):
                    print("OncoKB: upstream stopped (token missing, rejected or expired); failing",
                          file=sys.stderr)
                    sys.exit(1)
                raise
            if sent != total:
                # A record spanning lines would split across chunks.
                raise ValueError(f"MAF has {total} records but {sent} data lines")
        written = merge_maf(args.input, annotated, args.output, header, skip_all=not total)
        if written != total:
            raise ValueError("Merged MAF row count differs from input")
        print(f"Wrote {written} MAF rows with OncoKB annotations", flush=True)

if __name__ == "__main__":
    main()
