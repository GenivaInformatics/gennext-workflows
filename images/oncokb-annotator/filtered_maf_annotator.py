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


# ---------------------------------------------------------------------------
# refhub's OncoKB cache
# ---------------------------------------------------------------------------

# Where upstream's genomic-change POSTs go first (REFHUB_URL overrides; empty
# turns the cache off). refhub answers from its cache, asks OncoKB for what
# it lacks, and knows the locations OncoKB fails on (see refhub's README,
# "OncoKB cache").
DEFAULT_REFHUB_URL = "https://api.gennext.bio/api/v0/refhub"

# Seconds to connect to refhub, and to wait for its answer. On a cold cache
# refhub asks OncoKB for every row, in parallel, and bisects batches OncoKB
# fails on, so a chunk can take minutes.
REFHUB_TIMEOUT = (10, 1800)

# MAF rows per upstream call when refhub answers: refhub batches and
# parallelises its OncoKB requests itself, and takes up to 1000 queries.
REFHUB_CHUNK_ROWS = 500

# What a query refhub could not annotate becomes in upstream's answer list,
# marked in ONCOKB_CONSEQUENCE so the wrapper finds the row after upstream
# wrote it: known bad (filled at once) or not answered (tried again).
REFHUB_BAD = "refhub:known_bad"
REFHUB_RETRY = "refhub:not_answered"

# After this many refhub failures in a row, this run stops asking it.
REFHUB_MAX_FAILURES = 3

# How refhub asks OncoKB for somatic queries:
# - "genomic_first" (the default, refhub 0.9+): by genomic change, OncoKB's
#   own answer as it is; only the rows OncoKB fails on that way are answered
#   by protein change (counted in this step's log, nowhere else);
# - "protein": refhub maps the genomic change to a protein change itself and
#   asks byProteinChange;
# - "genomic": byGenomicChange only.
# REFHUB_ONCOKB_PATH overrides. refhub before 0.8 ignores the setting; refhub
# 0.8 refuses genomic_first, and the cache then asks by genomic change.
DEFAULT_REFHUB_PATH = "genomic_first"

# refhub 0.9.0's answer to a query carrying proteinFallback (HTTP 400, serde's
# "unknown field `proteinFallback`"). Only that clears the fallbacks.
REFHUB_090_FALLBACK_REFUSAL = "unknown field"

# refhub 0.8's answer to ?path=genomic_first (HTTP 400). Any other 400 is
# handled as a refhub failure, and genomic first stays on.
REFHUB_08_PATH_REFUSAL = "path is genomic or protein"

# The ANNOVAR step's oncokbGene facts (maftools-rs annovar-to-maf 0.2+):
# the protein change on OncoKB's own isoform, for refhub's strict
# fallback when OncoKB fails a row by genomic change. Absent or NA: none.
ONCOKB_GENE_COLUMNS = ("oncokb_isoform", "oncokb_isoform_aa", "oncokb_isoform_ref_ok",
                       "oncokb_isoform_exon_edge", "oncokb_gene_db")
_NA = {"", "NA", ".", "None", "nan"}


def protein_fallbacks(header_fields, lines):
    """genomicLocation (as upstream builds it) -> refhub's proteinFallback,
    for rows whose oncokbGene facts make them eligible: isoform and protein
    change set, residue check TRUE, not at an exon edge, the DB named.
    refhub checks again and applies the gene policy (curated isoform, not
    SRC/FGFR1). A location whose rows disagree gets none."""
    if not all(c in header_fields for c in ONCOKB_GENE_COLUMNS):
        return {}
    col = {name: i for i, name in enumerate(header_fields)}
    needed = ("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")
    if not all(c in col for c in needed):
        return {}
    out = {}
    for row in csv.reader(lines, delimiter="\t"):
        def cell(name):
            i = col.get(name)
            return row[i].strip() if i is not None and i < len(row) else ""
        isoform, alteration = cell("oncokb_isoform"), cell("oncokb_isoform_aa")
        eligible = (isoform not in _NA and alteration not in _NA
                    and cell("oncokb_isoform_ref_ok").upper() == "TRUE"
                    and cell("oncokb_isoform_exon_edge").upper() == "FALSE"
                    and cell("oncokb_gene_db") not in _NA)
        ref = cell("Reference_Allele")
        var = AnnotatorCore.get_var_allele(ref, cell("Tumor_Seq_Allele1"), cell("Tumor_Seq_Allele2"))
        location = AnnotatorCore.GenomicChangeQuery(
            cell("Chromosome"), cell("Start_Position"), cell("End_Position"), ref, var, "").genomicLocation
        fallback = None
        if eligible:
            fallback = {"transcript": isoform, "alteration": alteration}
            # OncoKB's hugoSymbol of the isoform (maftools-rs 0.2.0; it can
            # differ from Hugo_Symbol when refGene names another gene first).
            gene = cell("oncokb_isoform_gene") if "oncokb_isoform_gene" in col else cell("Hugo_Symbol")
            if gene not in _NA:
                fallback["hugoSymbol"] = gene
        if location in out and out[location] != fallback:
            out[location] = None
        else:
            out.setdefault(location, fallback)
    return {loc: fb for loc, fb in out.items() if fb is not None}


# OncoKB endpoints refhub caches, by the path upstream posts to.
REFHUB_ENDPOINTS = {
    "/annotate/mutations/byGenomicChange": "/oncokb/cache/annotate/mutations/byGenomicChange",
    "/annotate/germline/mutations/byGenomicChange":
        "/oncokb/cache/annotate/germline/mutations/byGenomicChange",
}


def _marker_annotation(marker):
    """An OncoKB-shaped answer upstream writes as a row whose
    ONCOKB_CONSEQUENCE is `marker`."""
    return {"geneExist": False, "variantExist": False, "oncogenic": "Unknown", "treatments": [],
            "query": {"hugoSymbol": "", "alteration": "", "consequence": marker}}


def _response(url, status, payload):
    import json
    import requests
    response = requests.models.Response()
    response.status_code = status
    response.url = url
    response.headers["Content-Type"] = "application/json"
    response._content = json.dumps(payload).encode()
    return response


class RefhubCache:
    """Posts upstream's genomic-change queries to refhub instead of OncoKB.

    post() returns what upstream expects from OncoKB: a 200 response holding
    one answer per query. A query refhub knows OncoKB fails on gets a marker
    answer (REFHUB_BAD), one it could not get answered another
    (REFHUB_RETRY); the wrapper finds those rows in upstream's output. Or it
    returns None: then the request goes to OncoKB directly, as without the
    cache. refhub unreachable, refusing the request, or reporting OncoKB
    unavailable for every query all give None, with a warning. The cache
    never fails a step.
    """

    def __init__(self, base_url, token, session=None, log=None, path=None):
        import requests
        self.base_url = base_url.rstrip("/")
        self.path = path or os.environ.get("REFHUB_ONCOKB_PATH", DEFAULT_REFHUB_PATH).strip() or "genomic"
        self.token = token
        self.session = session or requests.Session()
        self.log = log or (lambda message: print(f"OncoKB: {message}", flush=True))
        self.failures = 0
        self.enabled = True
        self.counts = {"requests": 0, "direct": 0, "hit": 0, "miss": 0, "refreshed": 0,
                       "stale": 0, "bad": 0, "failed": 0, "protein_route": 0, "fallback_offered": 0}
        # genomicLocation -> proteinFallback, for the MAF being annotated.
        self.fallbacks = {}

    def _give_up(self, why):
        self.failures += 1
        self.counts["direct"] += 1
        if self.failures >= REFHUB_MAX_FAILURES and self.enabled:
            self.enabled = False
            self.log(f"WARNING: refhub cache {why}; {self.failures} times in a row, "
                     "asking OncoKB directly for the rest of this run")
        else:
            self.log(f"WARNING: refhub cache {why}; asking OncoKB directly")
        return None

    def post(self, url, body):
        import json
        import requests
        endpoint = next((cache for path, cache in REFHUB_ENDPOINTS.items() if url.endswith(path)), None)
        if endpoint is None or not self.enabled:
            return None
        # Serialized as upstream serializes its POST body.
        queries = json.loads(json.dumps(body, default=lambda o: o.__dict__))
        self.counts["requests"] += 1
        target = self.base_url + endpoint
        if self.path != "genomic" and "/germline/" not in endpoint:
            target += f"?path={self.path}"
        offered = 0
        if self.path == "genomic_first" and "/germline/" not in endpoint and self.fallbacks:
            for query in queries:
                fallback = self.fallbacks.get(query.get("genomicLocation"))
                if fallback:
                    query["proteinFallback"] = fallback
                    offered += 1
        try:
            reply = self.session.post(
                target, json=queries, timeout=REFHUB_TIMEOUT,
                headers={"Authorization": f"Bearer {self.token}"})
        except requests.exceptions.RequestException as error:
            return self._give_up(f"unreachable ({type(error).__name__})")
        if (reply.status_code == 400 and offered
                and REFHUB_090_FALLBACK_REFUSAL in reply.text and "proteinFallback" in reply.text):
            # refhub 0.9.0 doesn't take the oncokbGene fallback: ask without it.
            self.log("refhub cache does not take proteinFallback; asking genomic first without it")
            self.fallbacks = {}
            self.counts["requests"] -= 1
            return self.post(url, body)
        if (reply.status_code == 400 and self.path == "genomic_first"
                and REFHUB_08_PATH_REFUSAL in reply.text):
            # refhub 0.8 doesn't know genomic first: ask it by genomic change.
            self.log("refhub cache does not offer genomic first; asking it by genomic change")
            self.path = "genomic"
            self.counts["requests"] -= 1
            return self.post(url, body)
        if reply.status_code != 200:
            if reply.status_code == 401:
                self.enabled = False
            return self._give_up(f"answered HTTP {reply.status_code}: {reply.text[:200]}")
        try:
            answers = reply.json()
            statuses = [answer["status"] for answer in answers]
        except (ValueError, KeyError, TypeError):
            return self._give_up("sent an unreadable answer")
        if len(answers) != len(queries):
            return self._give_up(f"sent {len(answers)} answers for {len(queries)} queries")
        self.failures = 0
        # Counted once refhub has taken the request (not for refused tries).
        self.counts["fallback_offered"] += offered

        if all(status == "unavailable" for status in statuses):
            # refhub could not reach OncoKB for any of these rows: ask OncoKB
            # directly, as without the cache (the cache never fails a step).
            # refhub itself answered, so this is not a refhub failure.
            reasons = sorted({a.get("reason", "") for a in answers})
            self.log(f"refhub cache could not reach OncoKB ({', '.join(reasons)}); "
                     "asking OncoKB directly")
            self.counts["direct"] += 1
            return None
        reply_answers = []
        for answer in answers:
            status = answer["status"]
            if status == "annotation":
                self.counts[answer.get("source")] = self.counts.get(answer.get("source"), 0) + 1
                if answer.get("route") == "protein":
                    # OncoKB failed the row by genomic change; refhub answered
                    # it by protein change. Logged here only, never written.
                    self.counts["protein_route"] += 1
                reply_answers.append(answer["annotation"])
            elif status == "bad":
                self.counts["bad"] += 1
                reply_answers.append(_marker_annotation(REFHUB_BAD))
            else:
                self.counts["failed"] += 1
                reply_answers.append(_marker_annotation(REFHUB_RETRY))
        return _response(url, 200, reply_answers)

    def summary(self):
        c = self.counts
        return (f"refhub cache: {c['requests']} requests; answers: {c['hit']} hits, "
                f"{c['miss'] + c['refreshed']} asked of OncoKB, {c['stale']} stale, "
                f"{c['bad']} known bad, {c['failed']} not answered; "
                f"{c['protein_route']} answered by protein change (OncoKB failed them "
                f"by genomic change; {c['fallback_offered']} rows sent with an oncokbGene "
                f"protein change); {c['direct']} requests sent to OncoKB directly")


# The RefhubCache upstream's POSTs go through, or None; set by annotate_maf.
_refhub = None
# Upstream's own makeoncokbpostrequest.
_direct_post = None


def _post(url, body):
    if _refhub is not None:
        response = _refhub.post(url, body)
        if response is not None:
            return response
        return _direct_in_parts(url, body)
    return _direct_post(url, body)


def _direct_in_parts(url, body):
    """Upstream's own POST, at most POST_QUERIES_THRESHOLD_GC_HGVSG queries at
    a time: refhub's chunks are larger than OncoKB takes in one request."""
    size = AnnotatorCore.POST_QUERIES_THRESHOLD_GC_HGVSG
    if len(body) <= size:
        return _direct_post(url, body)
    answers = []
    for start in range(0, len(body), size):
        response = _direct_post(url, body[start:start + size])
        if response.status_code != 200:
            return response
        answers += response.json()
    return _response(url, 200, answers)


def route_through_refhub():
    """Make upstream's POSTs ask _refhub first. Idempotent."""
    global _direct_post
    if AnnotatorCore.makeoncokbpostrequest is not _post:
        _direct_post = AnnotatorCore.makeoncokbpostrequest
        AnnotatorCore.makeoncokbpostrequest = _post


def refhub_url_for(api_url):
    """The refhub URL to use, or None. Not for --api-url: refhub caches what
    production OncoKB answers, not demo.oncokb.org."""
    if api_url:
        return None
    return os.environ.get("REFHUB_URL", DEFAULT_REFHUB_URL).strip() or None


def annotate_maf(input_path, output_path, token, reference, tumor_type=".",
                 api_url=None, chunk_rows=CHUNK_ROWS, retry_waits=RETRY_WAITS, sleep=time.sleep,
                 refhub_url=None, refhub_session=None, refhub_chunk_rows=REFHUB_CHUNK_ROWS):
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

    With refhub_url, upstream's genomic-change POSTs ask refhub's OncoKB
    cache first (RefhubCache), refhub_chunk_rows rows per upstream call.
    refhub answers row by row, so a chunk is not sent again row by row: a row
    refhub knows OncoKB fails on is filled at once (not retried, and not
    counted towards max_unanswered: it is no outage); a row refhub could not
    get answered is tried again after the pauses, with the other such rows.
    """
    global _refhub
    fail_fast_on_500()
    route_through_refhub()
    if api_url:
        AnnotatorCore.setoncokbbaseurl(api_url)
    AnnotatorCore.setoncokbapitoken(token)
    AnnotatorCore.validate_oncokb_token()
    _refhub = RefhubCache(refhub_url, token, refhub_session) if refhub_url else None
    if _refhub:
        print(f"OncoKB: asking refhub's cache first ({refhub_url})", flush=True)
        chunk_rows = refhub_chunk_rows
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
    if _refhub:
        _refhub.fallbacks = protein_fallbacks(header_fields, lines)
        if _refhub.fallbacks:
            print(f"OncoKB: {len(_refhub.fallbacks)} locations carry an oncokbGene protein change "
                  "(refhub's fallback if OncoKB fails them by genomic change)", flush=True)

    results = [None] * len(lines)
    state = {"header": None, "requests": 0}
    known_bad = set()
    annotation_columns = len(maf_annotation_header())
    consequence = maf_annotation_header().index("ONCOKB_CONSEQUENCE") - annotation_columns

    with tempfile.TemporaryDirectory(prefix="chunks-", dir=output_path.parent) as temp_dir:
        chunk_in = Path(temp_dir) / "in.maf"
        chunk_out = Path(temp_dir) / "out.maf"

        def run(indices):
            """Annotate these lines in one upstream call. Returns the rows left
            unanswered (refhub's markers), or None on a row failure."""
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
                    return None
                raise
            with chunk_out.open(newline="") as part:
                part_header = part.readline()
                out = part.readlines()
            if len(out) != len(indices):
                raise ValueError(f"OncoKB output has {len(out)} rows for {len(indices)} sent")
            if state["header"] is None:
                state["header"] = part_header
            left = []
            for i, line in zip(indices, out):
                marker = line.rstrip("\r\n").rsplit("\t", annotation_columns)[consequence]
                if marker == REFHUB_BAD:
                    known_bad.add(i)
                    left.append(i)
                elif marker == REFHUB_RETRY:
                    left.append(i)
                else:
                    results[i] = line
            return left

        limit = max_unanswered(len(lines))

        def one_by_one(indices, earlier=0):
            still = []
            for i in indices:
                if run([i]) != []:
                    still.append(i)
                    if earlier + len([j for j in still if j not in known_bad]) > 2 * limit:
                        # Stop early: an outage would take REQUEST_TIMEOUT per row.
                        raise OncoKBServerError(
                            f"OncoKB failed on more than {2 * limit} of {len(lines)} rows; "
                            "treated as an outage")
            return still

        def unknown(rows):
            """Rows OncoKB left unanswered that refhub doesn't know are bad."""
            return len([i for i in rows if i not in known_bad])

        def through_refhub():
            return bool(_refhub and _refhub.enabled)

        def again(rows, last):
            """Rows still unanswered after another try."""
            if not through_refhub():
                return one_by_one(rows)
            still = []
            for start in range(0, len(rows), chunk_rows):
                part = rows[start:start + chunk_rows]
                left = run(part)
                if left is None:
                    # The direct fallback failed. With refhub still on and a
                    # pass to come, keep the rows for it; on the last pass,
                    # send them one by one, so one failing row doesn't take
                    # its chunk with it.
                    if through_refhub() and not last:
                        left = part
                    else:
                        left = one_by_one(part, unknown(still))
                        if last:
                            # Sent alone only now, straight after the
                            # chunk's 500: they get one more pause.
                            solo_late.update(left)
                still += left
            return still

        failed = []
        # Rows first sent alone on the last pass (their chunk failed): as
        # 0.3.7 did for every row, they are tried alone once more after a pause.
        solo_late = set()
        starts = range(0, len(lines), chunk_rows)
        for number, start in enumerate(starts, 1):
            indices = list(range(start, min(start + chunk_rows, len(lines))))
            left = run(indices)
            if left is None:
                # A row failure from asking OncoKB directly. With refhub still
                # on and a retry pass to come, the rows are tried again
                # together after the pause (a 500-row chunk one by one takes
                # about 15 min); with no retry pass, one by one now.
                if through_refhub() and retry_waits:
                    left = indices
                else:
                    left = one_by_one(indices, unknown(failed))
            failed += left
            if number % 5 == 0 or number == len(starts):
                print(f"OncoKB: {indices[-1] + 1} of {len(lines)} MAF rows done, "
                      f"{unknown(failed)} failed so far, {len(failed) - unknown(failed)} known bad",
                      flush=True)
        for number, wait in enumerate(retry_waits, 1):
            retry = [i for i in failed if i not in known_bad]
            if not retry:
                break
            print(f"OncoKB: {len(retry)} rows failed; trying them again alone in {wait} s", flush=True)
            sleep(wait)
            failed = [i for i in failed if i in known_bad] + again(retry, number == len(retry_waits))
        late = [i for i in failed if i in solo_late and i not in known_bad]
        if late and retry_waits:
            wait = retry_waits[-1]
            print(f"OncoKB: {len(late)} rows failed alone right after their chunk; "
                  f"trying them alone again in {wait} s", flush=True)
            sleep(wait)
            failed = [i for i in failed if i not in late] + one_by_one(late)

        if unknown(failed) > limit:
            raise OncoKBServerError(
                f"OncoKB failed on {unknown(failed)} of {len(lines)} rows every time; "
                f"more than {limit} is treated as an outage")
        for i in failed:
            fields = next(csv.reader([lines[i].rstrip("\r\n")], delimiter="\t"))
            func = fields[func_column].split(";")[0] if func_column is not None else ""
            location = "\t".join(fields[1:6])
            if func in MUST_ANNOTATE_FUNC:
                raise OncoKBServerError(
                    f"OncoKB fails on a {func} row ({location}) every time; "
                    "not leaving it unannotated")
            how = "known bad in refhub's cache" if i in known_bad else "failed every time"
            print(f"OncoKB: {how}, written unannotated: {location} ({func or 'no Func.refGene'})",
                  flush=True)
            echo = lines[i].rstrip("\r\n").encode("ascii", "ignore").decode("ascii")
            results[i] = echo + "\t" + "\t".join(unanswered_suffix()) + "\n"

    with output_path.open("w", newline="") as target:
        target.write(state["header"] or
                     header.rstrip("\r\n") + "\t" + "\t".join(maf_annotation_header()) + "\n")
        target.writelines(results)
    print(f"OncoKB: {len(lines)} rows in {state['requests']} requests; {len(failed)} written unannotated "
          "(ANNOTATED=False)", flush=True)
    if _refhub:
        print(f"OncoKB: {_refhub.summary()}", flush=True)
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
                                    args.reference, args.tumor_type, args.api_url,
                                    refhub_url=refhub_url_for(args.api_url))
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
