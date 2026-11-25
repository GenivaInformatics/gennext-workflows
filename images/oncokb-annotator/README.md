# OncoKB annotator image

`bergun/oncokb-annotator:0.3.15-upstream-v4.1.0` is upstream
[oncokb-annotator](https://github.com/oncokb/oncokb-annotator) at tag `v4.1.0`
(commit `77521cd`), unmodified, plus the three wrapper scripts in this
directory. The previous image, `0.2.0-upstream-v3.4.2`, is upstream `v3.4.2`
(commit `1cbdfb2`) with the same wrappers. OncoKB publishes no container image, so the Dockerfile builds
upstream's source at that pinned commit on a digest-pinned
`python:3.12.12-alpine3.22`, installs upstream's own requirement files, and
copies the wrappers in last.

Differences from upstream, all of them:

| What | Kind |
| --- | --- |
| `filtered_maf_annotator.py`, `filtered_sv_annotator.py`, `backfill_filtered_maf.py` | our wrappers, added to `/app/oncokb` |
| `constraints.txt` holds `certifi`, `charset-normalizer`, `idna` at the versions `0.1.7-unfiltered` shipped; upstream's own files pin `requests` 2.32.4, `urllib3` 2.2.3 and `kiwisolver` 1.4.7 | packaging only |
| Python from the official `python` image instead of Alpine's `apk` package (3.12.12 in both) | packaging only |

No upstream file is patched: `/app/oncokb` is byte-identical to the upstream
tag, plus the three wrappers. (`0.2.0` was also byte-identical to
`0.1.7-unfiltered` for every shared file, with the same Python packages.) The unused `custom/` helpers and the image-level
`Dockerfile`/`.dockerignore` of the old base image are gone; upstream's
`.github/` and `.gitignore` are now present.

The same image serves `oncokb-templates` (SNV, SV/CNV, fusions) and
`oncokb-fusion-annotation-template`.

## SNV annotation wrapper

`filtered_maf_annotator.py` sends every MAF row to OncoKB through upstream's
`AnnotatorCore`, then merges the result back onto the original input rows.

It does what `MafAnnotator.py` does for our options (`-q Genomic_Change -r
<ref> [-t <tumour type>]`): set the token, validate it, call
`processalterationevents`. But it calls `processalterationevents` once per
chunk of 100 rows, because v4 sends all of a file's queries in one POST.
Upstream v3 posted 100 genomic-change queries at a time
(`POST_QUERIES_THRESHOLD_GC_HGVSG`). v4's mixed mode posts a whole file at
once, and OncoKB answers a MAF of tens of thousands of rows with HTTP 500
(upstream issue oncokb/oncokb-annotator#238 asks for request limits; there is
no fix as of v4.1.0).

Each chunk is the header line plus the next 100 data lines, copied byte for
byte. The wrapper refuses a MAF whose records span lines. The token is no
longer passed on a command line, so it doesn't show in a failed step's log.
`0.3.0` lacks the chunking and isn't used.

**Rows OncoKB cannot process.** One query OncoKB cannot process fails its
whole request. OncoKB works on it for about 31 s and answers HTTP 500, or
does not answer at all; requests sent right after such a query can fail too,
and succeed when sent alone a little later. The NanOnCT panel produces such
queries: insertions and deletions in long homopolymers at its MSI and intronic
sites (for example the poly-T run at chr1:119,510,721-119,510,740).

From `0.3.4`, with the handling below from `0.3.5`:

- upstream no longer retries a 500 or a read timeout, and waits 90 s for an
  answer instead of 240 s (`fail_fast_on_500()`). A 500 must raise: upstream's
  fallback for a failed POST re-sends each row as a GET without the reference
  genome, which OncoKB reads as GRCh37;
- a 100-row chunk that fails is sent again one row per request;
- rows that still fail are tried again, alone, after 60 s;
- a row that fails every time is written with `ANNOTATED=False`,
  `GENE_IN_ONCOKB=False`, `VARIANT_IN_ONCOKB=False`, `MUTATION_EFFECT=Unknown`
  and `ONCOGENIC=Unknown`, the rest empty, so the columns gennext-api-edge
  requires are set. The step log lists it. In DuckDB it looks like a variant
  OncoKB has no knowledge of; `ANNOTATED` is not carried into DuckDB.

The step fails instead when such a row is `exonic` or `splicing`
(`Func.refGene`), or when more than max(10, 5% of rows) are left (an outage).
429, 502, 503, 504 and connection errors still fail the step.

**refhub's OncoKB cache.** From `0.3.6`, upstream's genomic-change POSTs go
to gennext-refhub's OncoKB cache first (`RefhubCache`, refhub's README
"OncoKB cache"), with the same query objects and the same token. refhub
answers from its cache (annotations kept 1 month), asks OncoKB for the rest,
and knows the locations OncoKB fails on (kept 3 months, after a confirming
second try alone after a quiet pause), so each run stops spending 31-90 s
twice on each of them and stops sending OncoKB the same queries again.

- refhub URL: `REFHUB_URL`, default `https://api.gennext.bio/api/v0/refhub`;
  empty turns the cache off. `--api-url` (demo.oncokb.org tests) turns it off
  too: refhub caches production OncoKB.
- From `0.3.7`, upstream is called on 500 rows at a time and refhub answers
  every row in one go (it batches, parallelises and bisects its own OncoKB
  requests). A row refhub reports `bad` comes back marked
  (`ONCOKB_CONSEQUENCE=refhub:known_bad`) and is written with the
  `ANNOTATED=False` values above; its chunk is not sent again and there is no
  60 s retry. Known-bad rows don't count towards the outage cap (a NanOnCT
  MAF holds hundreds); an `exonic` or `splicing` one still fails the step. A
  row refhub reports `failed` or `unavailable` (`refhub:not_answered`) is
  tried again after 60 s, with the other such rows, and is then handled as
  above (filled, or the step fails).
- From `0.3.14`, a MAF row with the ANNOVAR step's oncokbGene facts
  (`oncokb_isoform`, `oncokb_isoform_aa`, `oncokb_isoform_ref_ok`,
  `oncokb_isoform_exon_edge`, `oncokb_gene_db`; maftools-rs annovar-to-maf)
  sends, with its genomic query, `proteinFallback: {transcript, alteration,
  hugoSymbol}` (hugoSymbol: `oncokb_isoform_gene`, maftools-rs 0.2.0, else
  `Hugo_Symbol`) when isoform and protein change are set, the residue check is
  TRUE, the row is not at an exon edge and the DB is named. refhub 0.9.1 uses
  it only when OncoKB fails the row by genomic change, after its own checks
  (curated isoform, not SRC/FGFR1, residues, exon edge). MAFs without the
  columns send none. refhub 0.9.0 refusing the field is asked without it.
  The step log counts rows sent with one; the output is unchanged.
- From `0.3.9`, somatic queries ask refhub genomic first
  (`?path=genomic_first`, refhub 0.9; answers must be OncoKB's own).
  refhub asks OncoKB by genomic change and passes on its answer
  as it is; only rows OncoKB fails on that way (known bad, or failed) are
  answered by protein change. The step log counts those
  (`N answered by protein change`); nothing in the output says so, and the
  MAF format is unchanged. `REFHUB_ONCOKB_PATH=protein` or `genomic` picks
  the older paths. refhub 0.8 refuses genomic first, and the wrapper then
  asks it by genomic change.
- From `0.3.8`, somatic queries ask refhub by protein change (`?path=protein`;
  `REFHUB_ONCOKB_PATH=genomic` for byGenomicChange): refhub maps the genomic
  change to a protein change on OncoKB's canonical transcript and asks
  OncoKB's byProteinChange, so Genome Nexus, where the queries OncoKB fails on
  hang, is out of the path. refhub before 0.8 ignores it. If refhub fails and
  the direct fallback fails too, the rows are tried again together through
  refhub after 60 s instead of one by one (a 500-row chunk takes about 15 min
  one by one). On the last pass they go one by one, so one failing row
  doesn't take its chunk with it, and rows that still fail there are tried
  alone once more after another 60 s (0.3.13). Rows refhub reports
  `unavailable` for a whole request go to OncoKB directly (0.3.11).
- The cache never fails a step. refhub unreachable, slow (30 min), answering
  anything but 200, or reporting OncoKB `unavailable` for every query sends
  that request to OncoKB directly, 100 queries a request, with a warning;
  after 3 refhub failures in a row the rest of the run asks OncoKB directly.
- The step log ends with refhub's counts (hits, asked of OncoKB, known bad,
  sent directly).
- A row refhub answers by protein change (genomic first, OncoKB failed it by
  genomic change) gets OncoKB's answer to that byProteinChange query, so
  `ONCOKB_PROTEIN_CHANGE`, `ONCOKB_CONSEQUENCE` and the other query-derived
  columns are that query's (gene, alteration, consequence), not Genome
  Nexus's. Nothing in the row marks the route; the step log counts it.

Rows are not filtered by gene: OncoKB maps intronic rows our MAF calls
`Unknown` to their gene, MET exon 14 skipping variants can be intronic, and
the TERT promoter rows are `upstream`.

**Fusions** (`fusion_annotator.py`, from `0.3.15`; used by
`oncokb-fusion-annotation-template`). It runs upstream's `FusionAnnotator.py`
unchanged, with the same arguments, and appends two columns built from the
same `/annotate/structuralVariants` answers, with no extra calls:

- `ONCOKB_MATCHED_ALTERATIONS`: the alterations OncoKB matched, as it names
  them in `treatments[].alterations`, `diagnosticImplications[].alterations`
  and `prognosticImplications[].alterations`. Examples: `EML4-ALK Fusion` for
  a curated pair, `Fusions` for a gene-level match. They are `;`-joined and
  sorted. The column is empty when the answer has no treatment or
  implication, since OncoKB names no matched alteration then.
- `ONCOKB_VARIANT_SUMMARY`: OncoKB's `variantSummary` sentence, on one line.

Upstream's columns, and their values, are unchanged. The template copies the
two columns into `fusions_with_oncokb.json`'s `OncoKB` object along with the
others.

**The token in logs.** The token reaches output in more ways than a command
line: upstream validates it with `GET /tokens/<token>`, so a connection error,
urllib3's retry warnings and a traceback all print that URL. From `0.3.3`:

- every wrapper masks the token (as `***`) in everything it prints, including
  logging and tracebacks (`mask_output()`);
- the SV wrapper and `backfill_filtered_maf.py`, which run upstream's
  command-line scripts with `-b <token>`, pass upstream's output through the
  same mask and report a failure by exit status only, never the command
  (`run_masked()`);
- the task scripts that call upstream directly (`oncokb-cnv`, `oncokb-fusions`,
  `oncokb-fusion-annotation`) use the same helpers or a `sed` mask.

Upstream's `validate_oncokb_token` calls `exit()` when it rejects a token,
which is status 0 when the MAF wrapper calls it in-process; the wrapper turns
that into exit 1. The live `0.2.0` tasks print the token on any such failure.
Upstream writes each output row through `encode('ascii', 'ignore')`, which
drops non-ASCII characters from echoed input columns (ClinVar `CLNDN`
"Muir-Torré"); the merge checks the echo against the ASCII-folded input and
keeps the original text. The script fails if the annotator changes the row
order, row count, or output schema.

Since upstream v4.0, `MafAnnotator.py` always annotates in mixed
somatic/germline mode. Its header has three more columns after `ONCOGENIC`:
`PATHOGENIC (Germline)`, `PENETRANCE (Germline)` and
`GENOMIC_INDICATOR (Germline)`. `maf_annotation_header()` asks `AnnotatorCore`
for exactly that header. Rows go to the germline endpoints only when a
`Mutation_Status` column says `germline`; our MAFs have no such column, so
every row is annotated as somatic and the three columns stay empty. Every
reader of the MAF (`duckdb-snv`, `determine-otherinfo-column`) selects columns
by name, so MAFs of both layouts load. One more v4 change: a row without
genomic coordinates is no longer queried (v3 sent it and got "Unknown"); our
ANNOVAR-built MAFs always have coordinates.

The SV/CNV, CNA and fusion annotators write the same columns as in v3.

`backfill_filtered_maf.py` only handles MAFs from the v3 annotator's era, so
run it from `0.2.0-upstream-v3.4.2`; with this image it refuses them.

An earlier version sent only exonic, splice, UTR and upstream rows (by
`Func.refGene`) to reduce WGS runtime. It applied to every kit, left panel
samples' intronic and intergenic rows with empty OncoKB fields, and has been
removed.

For the Genome kit, identified by the exact regions-file basename
`hg38_canonical_chrom_contig_regions.bed` (or its validated `.valid.bed`
copy), the SNV and SV/CNV wrappers do not call OncoKB. They retain every input
row and write the expected OncoKB columns empty. The somatic DNA v2 FASTQ,
BAM, and VCF workflows forward the original `dna-input-files.regionsFile`
value. Per-sample and batch CNV workflows forward their regions file to the
same SV/CNV wrapper. Other kits retain annotation.

Both `oncokb-snv` and `oncokb-snv-on-node` in
`tasks/base_tasks/oncokb.yaml` use this image. A running Argo workflow may
have cached the old template in `status.storedTemplates`, so deploying the
WorkflowTemplate alone does not update already-running workflows.

Run tests in the built image before pushing:

```sh
docker buildx build --load -t oncokb-annotator:test images/oncokb-annotator
docker run --rm --entrypoint sh \
  -v "$PWD/images/oncokb-annotator:/work:ro" \
  oncokb-annotator:test \
  -c 'cd /tmp && PYTHONPATH=/work:/app/oncokb python3 -m unittest -v test_filtered_maf_annotator'
```
