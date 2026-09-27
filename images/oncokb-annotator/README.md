# Filtered OncoKB SNV annotation

Whole-genome ANNOVAR MAFs contain millions of intergenic and intronic calls.
Sending every call to the OncoKB API makes the somatic post-processing step
impractically long.

`filtered_maf_annotator.py` selects rows by `Func.refGene` before calling the
upstream `MafAnnotator.py`. It sends exonic, splice, UTR, upstream/promoter,
and corresponding noncoding RNA classes. Intergenic, intronic, and downstream
rows are not sent. The output still contains every input row in the original
order and all 30 OncoKB columns; skipped rows have empty OncoKB values. The
script fails if the MAF lacks `Func.refGene` or if the annotator changes the
selected row order or output schema.

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
