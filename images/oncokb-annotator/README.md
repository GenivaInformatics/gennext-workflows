# OncoKB annotator image

`bergun/oncokb-annotator:0.2.0-upstream-v3.4.2` is upstream
[oncokb-annotator](https://github.com/oncokb/oncokb-annotator) at tag `v3.4.2`
(commit `1cbdfb2`), unmodified, plus the three wrapper scripts in this
directory. OncoKB publishes no container image, so the Dockerfile builds
upstream's source at that pinned commit on a digest-pinned
`python:3.12.12-alpine3.22`, installs upstream's own requirement files, and
copies the wrappers in last.

Differences from upstream, all of them:

| What | Kind |
| --- | --- |
| `filtered_maf_annotator.py`, `filtered_sv_annotator.py`, `backfill_filtered_maf.py` | our wrappers, added to `/app/oncokb` |
| `constraints.txt` holds `certifi`, `charset-normalizer`, `idna` at the versions `0.1.7-unfiltered` shipped | packaging only |
| Python from the official `python` image instead of Alpine's `apk` package (3.12.12 in both) | packaging only |

No upstream file is patched. Compared with `0.1.7-unfiltered`, every shared
file under `/app/oncokb` is byte-identical and the installed Python packages
are the same versions. The unused `custom/` helpers and the image-level
`Dockerfile`/`.dockerignore` of the old base image are gone; upstream's
`.github/` and `.gitignore` are now present.

The same image serves `oncokb-templates` (SNV, SV/CNV, fusions) and
`oncokb-fusion-annotation-template`.

## SNV annotation wrapper

`filtered_maf_annotator.py` sends every MAF row to the upstream
`MafAnnotator.py`, then merges the result back onto the original input rows.
Upstream writes each output row through `encode('ascii', 'ignore')`, which
drops non-ASCII characters from echoed input columns (ClinVar `CLNDN`
"Muir-Torré"); the merge checks the echo against the ASCII-folded input and
keeps the original text. The script fails if the annotator changes the row
order, row count, or output schema.

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
