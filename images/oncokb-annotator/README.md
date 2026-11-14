# OncoKB annotator image

`bergun/oncokb-annotator:0.3.2-upstream-v4.1.0` is upstream
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

The SV wrapper and `backfill_filtered_maf.py` still run upstream's command-line
scripts, which take the token as `-b`. When upstream exits non-zero they print
only the script name and exit status, never the command (from `0.3.2`; before,
an uncaught `CalledProcessError` printed it). The live `0.2.0` tasks still
print it on such a failure.
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
