# OncoKB SNV annotation wrapper

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

Run tests in the upstream image before building:

```sh
docker run --rm --entrypoint sh \
  -v "$PWD/images/oncokb-annotator:/work:ro" \
  bergun/oncokb-annotator:v0.1.2 \
  -c 'PYTHONPATH=/work:/app/oncokb python -m unittest -v test_filtered_maf_annotator'
```
