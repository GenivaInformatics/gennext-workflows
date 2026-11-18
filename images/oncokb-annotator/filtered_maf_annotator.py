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
import time

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


# Pauses before each extra pass over the rows that failed. A query OncoKB
# cannot process takes it about 31 s to answer with HTTP 500, or it never
# answers; requests sent right after one can fail too, and succeed alone a
# little later.
RETRY_WAITS = (60,)

# Seconds to wait for OncoKB's answer to one annotation request (upstream: 240,
# retried 3 times on a read timeout). A healthy 100-row request takes seconds.
REQUEST_TIMEOUT = 90

# More rows than this left unanswered is an outage, not a few queries OncoKB
# cannot process: the step fails.
def max_unanswered(rows):
    return max(10, rows // 20)


class OncoKBServerError(Exception):
    """OncoKB answered a request with HTTP 500."""


def is_row_failure(error):
    """True when OncoKB failed on what was sent rather than being unavailable.

    That is an HTTP 500 (raised by our hook, or by upstream's retries) or no
    answer within REQUEST_TIMEOUT: OncoKB does both for queries it cannot
    process, such as some indels in long homopolymers. 429, 502, 503, 504 and
    connection errors are the service or the network; those fail the step.
    """
    import requests
    from urllib3.exceptions import MaxRetryError, ReadTimeoutError, ResponseError

    if isinstance(error, (OncoKBServerError, requests.exceptions.ReadTimeout)):
        return True
    if isinstance(error, requests.exceptions.ConnectionError):
        cause = error.args[0] if error.args else None
        return isinstance(cause, ReadTimeoutError) or (
            isinstance(cause, MaxRetryError) and isinstance(cause.reason, ReadTimeoutError))
    if not isinstance(error, requests.exceptions.RetryError):
        return False
    cause = error.args[0] if error.args else None
    return (isinstance(cause, MaxRetryError) and isinstance(cause.reason, ResponseError)
            and ResponseError.SPECIFIC_ERROR.format(status_code=500) in str(cause.reason))


def raise_on_500(response, *args, **kwargs):
    if response.status_code == 500:
        raise OncoKBServerError(f"HTTP 500 from {response.request.method} {response.url.split('?')[0]}")


def fail_fast_on_500():
    """Make upstream's requests fail fast on queries OncoKB cannot process.

    - HTTP 500 raises OncoKBServerError, unretried. Upstream retries 500 three
      times, about two minutes for each such query. Without the raise, a POST
      answered with 500 would reach upstream's fallback, which re-sends each
      row as a GET without the reference genome (OncoKB then assumes GRCh37).
    - Requests time out after REQUEST_TIMEOUT and a read timeout is not
      retried; upstream waits 240 s and retries 3 times.
    429, 502, 503, 504 and connection failures keep upstream's retries.
    Idempotent.
    """
    original = AnnotatorCore.requests_retry_session
    if getattr(original, "_fails_fast_on_500", False):
        return

    def session(*args, **kwargs):
        if len(args) >= 3:
            args = list(args)
            args[2] = [code for code in args[2] if code != 500]
        else:
            kwargs["status_forcelist"] = [
                code for code in kwargs.get("status_forcelist",
                                            AnnotatorCore.API_REQUEST_RETRY_STATUS_FORCELIST)
                if code != 500]
        result = original(*args, **kwargs)
        for prefix, adapter in list(result.adapters.items()):
            adapter.max_retries = adapter.max_retries.new(read=0)
        result.hooks["response"].append(raise_on_500)
        return result

    session._fails_fast_on_500 = True
    AnnotatorCore.requests_retry_session = session
    AnnotatorCore.REQUEST_TIMEOUT = REQUEST_TIMEOUT


def unanswered_suffix():
    return [UNANSWERED.get(name, "") for name in maf_annotation_header()]


def annotate_maf(input_path, output_path, token, reference, tumor_type=".",
                 api_url=None, chunk_rows=CHUNK_ROWS, retry_waits=RETRY_WAITS, sleep=time.sleep):
    """MafAnnotator.py -q Genomic_Change, run chunk by chunk. Returns the rows written.

    Does what MafAnnotator.main does for our options, once: set the API URL and
    token, validate the token. Then it calls upstream's processalterationevents
    on consecutive chunks of the input and joins their output in input order.
    Each chunk is the input's header line plus data lines, copied byte for
    byte, so upstream reads the same text as from the whole file. The token is
    never put on a command line.

    One query OncoKB cannot process fails its whole request (HTTP 500 or no
    answer, see is_row_failure). Such a chunk is sent again one row per
    request. Rows that still fail are tried again, alone, after each pause in
    retry_waits. Rows that fail every time are written with UNANSWERED
    values; every other row keeps its OncoKB answer. The step fails instead
    when one of them is a coding or splice-site row (MUST_ANNOTATE_FUNC), or
    when there are more than max_unanswered of them.
    """
    fail_fast_on_500()
    if api_url:
        AnnotatorCore.setoncokbbaseurl(api_url)
    AnnotatorCore.setoncokbapitoken(token)
    AnnotatorCore.validate_oncokb_token()
    query_type = AnnotatorCore.QueryType.GENOMIC_CHANGE
    reference_genome = AnnotatorCore.ReferenceGenome[reference.upper()]
    cancer_type = "" if tumor_type == "." else tumor_type

    header = None
    lines = []
    with input_path.open(newline="") as source:
        for line in source:
            if header is None:
                if line.strip() and not line.startswith("#"):
                    header = line if line.endswith("\n") else line + "\n"
                continue
            lines.append(line if line.endswith("\n") else line + "\n")
    if header is None:
        raise ValueError("MAF has no header")
    records = sum(1 for _ in csv.reader(lines, delimiter="\t"))
    if records != len(lines):
        # Upstream is sent whole lines, so a record spanning lines would split.
        raise ValueError(f"MAF has {records} records but {len(lines)} data lines")
    header_fields = next(csv.reader([header.rstrip("\r\n")], delimiter="\t"))
    func_column = header_fields.index("Func.refGene") if "Func.refGene" in header_fields else None

    results = [None] * len(lines)
    state = {"header": None, "requests": 0}

    with tempfile.TemporaryDirectory(prefix="chunks-", dir=output_path.parent) as temp_dir:
        chunk_in = Path(temp_dir) / "in.maf"
        chunk_out = Path(temp_dir) / "out.maf"

        def run(indices):
            """Annotate these lines in one upstream call. False on a row failure."""
            with chunk_in.open("w", newline="") as part:
                part.write(header)
                part.writelines(lines[i] for i in indices)
            state["requests"] += 1
            try:
                AnnotatorCore.processalterationevents(
                    str(chunk_in), str(chunk_out), "", cancer_type, {}, False,
                    query_type, reference_genome, False)
            except Exception as error:
                if is_row_failure(error):
                    return False
                raise
            with chunk_out.open(newline="") as part:
                part_header = part.readline()
                out = part.readlines()
            if len(out) != len(indices):
                raise ValueError(f"OncoKB output has {len(out)} rows for {len(indices)} sent")
            if state["header"] is None:
                state["header"] = part_header
            for i, line in zip(indices, out):
                results[i] = line
            return True

        limit = max_unanswered(len(lines))

        def one_by_one(indices, earlier=0):
            still = []
            for i in indices:
                if not run([i]):
                    still.append(i)
                    if earlier + len(still) > 2 * limit:
                        # Stop early: an outage would take REQUEST_TIMEOUT per row.
                        raise OncoKBServerError(
                            f"OncoKB failed on more than {2 * limit} of {len(lines)} rows; "
                            "treated as an outage")
            return still

        failed = []
        starts = range(0, len(lines), chunk_rows)
        for number, start in enumerate(starts, 1):
            indices = list(range(start, min(start + chunk_rows, len(lines))))
            if not run(indices):
                failed += one_by_one(indices, len(failed))
            if number % 5 == 0 or number == len(starts):
                print(f"OncoKB: {indices[-1] + 1} of {len(lines)} MAF rows done, "
                      f"{len(failed)} failed so far", flush=True)
        for wait in retry_waits:
            if not failed:
                break
            print(f"OncoKB: {len(failed)} rows failed; trying them again alone in {wait} s", flush=True)
            sleep(wait)
            failed = one_by_one(failed)

        if len(failed) > limit:
            raise OncoKBServerError(
                f"OncoKB failed on {len(failed)} of {len(lines)} rows every time; "
                f"more than {limit} is treated as an outage")
        for i in failed:
            fields = next(csv.reader([lines[i].rstrip("\r\n")], delimiter="\t"))
            func = fields[func_column].split(";")[0] if func_column is not None else ""
            location = "\t".join(fields[1:6])
            if func in MUST_ANNOTATE_FUNC:
                raise OncoKBServerError(
                    f"OncoKB fails on a {func} row ({location}) every time; "
                    "not leaving it unannotated")
            print(f"OncoKB: failed every time, written unannotated: {location} ({func or 'no Func.refGene'})",
                  flush=True)
            echo = lines[i].rstrip("\r\n").encode("ascii", "ignore").decode("ascii")
            results[i] = echo + "\t" + "\t".join(unanswered_suffix()) + "\n"

    with output_path.open("w", newline="") as target:
        target.write(state["header"] or
                     header.rstrip("\r\n") + "\t" + "\t".join(maf_annotation_header()) + "\n")
        target.writelines(results)
    print(f"OncoKB: {len(lines)} rows in {state['requests']} requests; {len(failed)} written unannotated "
          "(ANNOTATED=False)", flush=True)
    return len(lines)


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
