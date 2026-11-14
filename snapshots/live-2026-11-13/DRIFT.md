# Live WorkflowTemplates vs origin/main

Snapshot of namespace `argo` taken 2026-11-13 from the live cluster (196 templates; 9 templates of a separate product are excluded from this snapshot and its counts), compared by `spec` with WorkflowTemplate documents on `origin/main` (c8c7807).

- matches-main: 125
- drift (repo has the name, spec differs): 20
- not-in-repo: 42

| Template | Status | Created | Repo file(s) |
|---|---|---|---|
| `alter-bam-contig-names-template` | matches-main | 2024-03-27 | tasks/base_tasks/alter_bam_contig_names.yaml |
| `annotsv-template` | matches-main | 2024-01-22 | tasks/base_tasks/annotsv.yaml |
| `annotsv-to-oncokb-template` | drift | 2024-01-22 | tasks/base_tasks/annotsv-to-oncokb.yaml |
| `annovar-template` | matches-main | 2023-12-24 | tasks/base_tasks/annovar.yaml |
| `annovar-to-maf-template` | matches-main | 2023-12-24 | tasks/base_tasks/annovar-to-maf.yaml |
| `arriba-fusion-detection` | not-in-repo | 2026-02-10 |  |
| `arriba-fusion-detection-template` | matches-main | 2026-04-11 | tasks/base_tasks/arriba.yaml |
| `bcftools-template` | matches-main | 2024-07-21 | tasks/base_tasks/bcftools.yaml |
| `bed-split-template` | matches-main | 2023-09-24 | tasks/bed_split.yaml |
| `bed-split-workflow` | drift | 2023-05-31 | tasks/bed_split_v1.yaml |
| `bed-validator-template` | matches-main | 2025-11-19 | tasks/base_tasks/bed_validator.yaml |
| `bgzip-tabix-template` | matches-main | 2024-03-26 | tasks/base_tasks/bgzip-tabix.yaml |
| `bwa-mem2-template` | matches-main | 2023-05-24 | tasks/base_tasks/bwa-mem2.yaml |
| `caller-workflow` | not-in-repo | 2023-09-20 |  |
| `cancervar-template` | matches-main | 2023-12-24 | tasks/base_tasks/cancervar.yaml |
| `cat-json` | matches-main | 2024-07-21 | tasks/test_tasks/cat_json.yaml |
| `ces-workflow` | matches-main | 2023-07-01 | workflows/ces-workflow.yaml |
| `change-run-state-template` | matches-main | 2023-10-28 | tasks/change_run_state.yaml |
| `change-sample-state-template` | matches-main | 2026-03-29 | tasks/change_sample_state.yaml |
| `check-file-path` | matches-main | 2024-09-10 | tasks/check-file-path.yaml |
| `cnv-calling` | not-in-repo | 2026-02-10 |  |
| `cnv-fusion-combined` | not-in-repo | 2026-02-11 |  |
| `cnv-workflow-template` | matches-main | 2025-07-23 | workflows/cnv_workflow.yaml |
| `cnvkit-analysis` | not-in-repo | 2026-02-10 |  |
| `cnvkit-annotate` | not-in-repo | 2026-02-10 |  |
| `cnvkit-export` | not-in-repo | 2026-02-10 |  |
| `cnvkit-filter` | not-in-repo | 2026-02-10 |  |
| `cnvpytor-template` | matches-main | 2025-03-03 | tasks/base_tasks/cnvpytor.yaml |
| `collect-qc-data-templates` | not-in-repo | 2023-11-23 |  |
| `collect-qc-data-v2-templates` | matches-main | 2025-09-10 | tasks/collect_qc_data.yaml |
| `collect-somatic-data-templates` | matches-main | 2024-02-18 | tasks/collect_somatic_data.yaml |
| `collect-tmb-data-template` | matches-main | 2025-07-18 | tasks/collect_tmb_data.yaml |
| `complete-workflow-test` | not-in-repo | 2023-06-26 |  |
| `compose-read-group-str` | matches-main | 2023-07-02 | tasks/compose-read-group-str.yaml |
| `copy-inputs-to-scratch` | drift | 2023-11-20 | tasks/copy_inputs_to_scratch.yaml |
| `copy-inputs-to-temp` | matches-main | 2023-07-02 | tasks/copy_inputs_to_temp.yaml |
| `copy-temp-to-output` | drift | 2023-11-20 | tasks/copy_temp_to_output.yaml |
| `cravat-templates` | matches-main | 2023-10-27 | tasks/cravat/cravat.yaml |
| `cravat-templates-nodeselect` | matches-main | 2024-03-05 | tasks/cravat/cravat_nodeSelector.yaml |
| `cravat-templates-trio` | matches-main | 2025-11-06 | tasks/cravat/cravat-trio.yaml |
| `cuda-vector-add-test` | matches-main | 2024-07-13 | tasks/test_tasks/gpu_task.yaml |
| `deepvariant-template` | matches-main | 2023-09-24 | tasks/callers/deepvariant.yaml |
| `dir-clean-up` | matches-main | 2023-11-20 | tasks/base_tasks/dir_cleanup.yaml |
| `duckdb-template` | drift | 2023-12-24 | tasks/base_tasks/duckdb.yaml |
| `elprep-template` | matches-main | 2023-05-25 | tasks/base_tasks/elprep.yaml |
| `example-template-3ffd` | not-in-repo | 2025-10-08 |  |
| `example-template-3fv` | not-in-repo | 2025-10-06 |  |
| `exomedepth-template` | drift | 2025-04-13 | tasks/base_tasks/exomedepth.yaml |
| `extract-reads-for-fusion-inspector-template` | drift | 2026-05-06 | tasks/base_tasks/extract-reads-for-fusion-inspector.yaml |
| `fan-out-in-params-workflow` | not-in-repo | 2023-05-26 |  |
| `fastp-template` | matches-main | 2026-04-29 | tasks/base_tasks/fastp.yaml |
| `fastqc-template` | drift | 2023-05-17 | tasks/base_tasks/fastqc.yaml |
| `fastqc-v2-template` | not-in-repo | 2023-06-05 |  |
| `fgbio-mark-duplicates-template` | matches-main | 2026-05-07 | tasks/base_tasks/fgbio-mark-duplicates.yaml |
| `fgbio-umi-dna-template` | matches-main | 2026-10-25 | tasks/base_tasks/fgbio-umi-dna.yaml |
| `filter-and-annot-workflow` | not-in-repo | 2024-07-16 |  |
| `filter-annotsv-template` | matches-main | 2024-01-23 | tasks/base_tasks/filter-annotsv.yaml |
| `filter-artifact-fusions-template` | matches-main | 2026-05-07 | tasks/base_tasks/filter-artifact-fusions.yaml |
| `filter-chimeric-junction-by-bam-template` | matches-main | 2026-05-07 | tasks/base_tasks/filter-chimeric-junction-by-bam.yaml |
| `filter-indelible-template` | matches-main | 2023-12-24 | tasks/base_tasks/filter-indelible.yaml |
| `filter-msisensor-pro-template` | drift | 2025-08-31 | tasks/base_tasks/filter-msisensor-pro.yaml |
| `fusion-detection-ffpe` | not-in-repo | 2026-02-10 |  |
| `fusion-inspector-template` | matches-main | 2026-04-29 | tasks/base_tasks/fusion-inspector.yaml |
| `fusion-report-generation` | not-in-repo | 2026-02-10 |  |
| `fusion-report-template` | matches-main | 2026-04-11 | tasks/base_tasks/fusion-report.yaml |
| `generate-run-info` | matches-main | 2024-09-05 | tasks/generate_run_info.yaml |
| `gennext-prs-matcher-template` | matches-main | 2025-06-04 | tasks/prs/gennext_prs_matcher.yaml |
| `gennext-prs-template` | matches-main | 2024-12-15 | tasks/prs/gennext_prs.yaml |
| `genotype-conditional-workflow-template` | drift | 2025-04-24 | tasks/prs/genotype_conditional_workflow.yaml |
| `germline-rna-gpu-workflow` | matches-main | 2025-05-23 | workflows/germline/germline-rna-gpu-workflow.yaml |
| `germline-tps-gpu-workflow` | matches-main | 2025-05-18 | workflows/germline/germline-tps-gpu-workflow.yaml |
| `germline-wes-gpu-workflow` | matches-main | 2024-07-25 | workflows/germline/germline-wes-gpu-workflow.yaml |
| `germline-wgs-gpu-workflow` | matches-main | 2024-07-18 | workflows/germline/germline-wgs-gpu-workflow.yaml |
| `get-abs-ref-dirs` | matches-main | 2023-07-03 | tasks/get_abs_ref_dirs.yaml |
| `get-base-node-config` | matches-main | 2023-11-20 | tasks/get-base-node-config.yaml |
| `get-node-base-dirs` | matches-main | 2023-06-16 | tasks/get_node_base_dirs.yaml |
| `get-node-config` | matches-main | 2023-07-03 | tasks/get-node-config.yaml |
| `get-node-rel-dirs` | not-in-repo | 2023-06-16 |  |
| `get-rel-dirs` | not-in-repo | 2023-06-19 |  |
| `get-rel-ref-dirs` | matches-main | 2023-06-19 | tasks/get_rel_ref_dirs.yaml |
| `get-run-data` | matches-main | 2023-07-01 | tasks/get_run_data.yaml |
| `gridss-template` | not-in-repo | 2025-03-18 |  |
| `haplotypecaller-template` | matches-main | 2023-09-19 | tasks/callers/haplotypecaller.yaml |
| `haplotypecaller-wrapper-template` | matches-main | 2023-09-24 | workflows/haplotypecaller-wrapper-workflow.yaml |
| `hostpath-example` | not-in-repo | 2023-06-26 |  |
| `indelible-template` | matches-main | 2023-12-24 | tasks/base_tasks/indelible.yaml |
| `load-and-prepare` | matches-main | 2026-04-02 | tasks/load_and_prepare.yaml |
| `maf-to-db-template` | matches-main | 2024-07-22 | tasks/base_tasks/maf-to-db.yaml |
| `manta-calls` | not-in-repo | 2024-01-28 |  |
| `manta-template` | drift | 2024-01-22 | tasks/base_tasks/manta.yaml |
| `mosdepth-template` | matches-main | 2023-11-14 | tasks/base_tasks/mosdepth.yaml |
| `msi-tmb-coverage-workflow-template` | drift | 2025-07-29 | workflows/msi_tmb_coverage_workflow.yaml |
| `msisensor-pro-template` | matches-main | 2024-01-02 | tasks/base_tasks/msisensor-pro.yaml |
| `multiqc-template` | matches-main | 2023-11-14 | tasks/base_tasks/multiqc.yaml |
| `multiqc-template-test` | not-in-repo | 2023-11-17 |  |
| `mutect2-sharded-template` | matches-main | 2026-04-28 | tasks/base_tasks/mutect2-sharded.yaml |
| `mutect2-template` | matches-main | 2023-12-24 | tasks/base_tasks/mutect2.yaml |
| `no-filter-mutect2-calls` | not-in-repo | 2024-01-27 |  |
| `oncokb-fusion-annotation` | not-in-repo | 2026-02-10 |  |
| `oncokb-fusion-annotation-template` | matches-main | 2026-04-11 | tasks/base_tasks/oncokb-fusion-annotation.yaml |
| `oncokb-templates` | matches-main | 2023-12-24 | tasks/base_tasks/oncokb.yaml |
| `parabricks-bam2fastq-template` | not-in-repo | 2024-07-20 |  |
| `parabricks-bam2fq-template` | matches-main | 2024-07-20 | tasks/base_tasks/gpu_tasks/parabricks-bam2fq.yaml |
| `parabricks-bamsort-template` | matches-main | 2026-04-29 | tasks/base_tasks/gpu_tasks/parabricks-bamsort.yaml |
| `parabricks-bqsr-template` | drift | 2024-07-18 | tasks/base_tasks/gpu_tasks/parabricks-bqsr.yaml |
| `parabricks-fastq2bam-template` | not-in-repo | 2024-07-19 |  |
| `parabricks-fq2bam-template` | matches-main | 2024-07-19 | tasks/base_tasks/gpu_tasks/parabricks-fq2bam.yaml |
| `parabricks-mutect2-template` | matches-main | 2024-07-22 | tasks/base_tasks/gpu_tasks/parabricks-mutect2.yaml |
| `parabricks-rnafq2bam-template` | matches-main | 2024-07-29 | tasks/base_tasks/gpu_tasks/parabricks-rnafq2bam.yaml |
| `parabricks-somatic-pipeline-template` | matches-main | 2024-07-19 | tasks/base_tasks/gpu_tasks/parabricks-somatic.yaml |
| `parabricks-starfusion-template` | matches-main | 2024-07-29 | tasks/base_tasks/gpu_tasks/parabricks-starfusion.yaml |
| `parabricks-template` | not-in-repo | 2024-07-14 |  |
| `per-sample-cnv-analysis` | matches-main | 2026-04-08 | workflows/per-sample-cnv-analysis.yaml |
| `per-sample-cnv-cohort-analysis` | matches-main | 2026-05-03 | workflows/per-sample-cnv-cohort-analysis.yaml |
| `per-sample-germline` | matches-main | 2026-07-19 | workflows/per-sample-germline.yaml |
| `per-sample-somatic-dna` | matches-main | 2026-03-29 | workflows/per-sample-somatic-dna.yaml |
| `per-sample-somatic-dna-rna` | matches-main | 2026-03-29 | workflows/per-sample-somatic-dna-rna.yaml |
| `per-sample-somatic-fusion-gpu` | matches-main | 2026-04-29 | workflows/per-sample-somatic-fusion-gpu.yaml |
| `per-sample-somatic-rna` | matches-main | 2026-03-29 | workflows/per-sample-somatic-rna.yaml |
| `pgs-workflow` | not-in-repo | 2024-12-15 |  |
| `picard-collect-rnaseq-metrics-template` | matches-main | 2026-04-29 | tasks/base_tasks/picard-collect-rnaseq-metrics.yaml |
| `picard-markduplicates-template` | matches-main | 2026-04-29 | tasks/base_tasks/picard-markduplicates.yaml |
| `post-process-cnv` | drift | 2026-04-06 | workflows/post-process-cnv.yaml |
| `post-process-variant-counter` | matches-main | 2026-04-08 | workflows/post-process-variant-counter.yaml |
| `preprocess-fastq` | matches-main | 2023-07-02 | workflows/preprocess-fastq.yaml |
| `propagate-fi-dup-flags-template` | matches-main | 2026-05-07 | tasks/base_tasks/propagate-fi-dup-flags.yaml |
| `prs-workflow` | matches-main | 2025-04-24 | workflows/prs/prs_workflow.yaml |
| `qorts-template` | matches-main | 2024-07-30 | tasks/base_tasks/qorts.yaml |
| `recompute-cds-coverage-workflow-template` | not-in-repo | 2026-02-12 |  |
| `remove-temp-fqs-template` | matches-main | 2024-09-18 | tasks/base_tasks/remove-temp-fqs.yaml |
| `rustqc-template` | matches-main | 2026-04-11 | tasks/base_tasks/rustqc.yaml |
| `samtools-template` | matches-main | 2023-11-20 | tasks/base_tasks/samtools.yaml |
| `scarhrd` | not-in-repo | 2026-02-11 |  |
| `seqkit-bedtools-template` | matches-main | 2024-07-08 | tasks/base_tasks/seqkit_bedtools.yaml |
| `snap-template` | matches-main | 2023-07-02 | tasks/base_tasks/snap.yaml |
| `somatic-dna-bam-v2` | matches-main | 2026-04-03 | workflows/somatic/v2/somatic-dna-bam.yaml |
| `somatic-dna-fastq-v2` | matches-main | 2026-04-03 | workflows/somatic/v2/somatic-dna-fastq.yaml |
| `somatic-dna-padding-repair-v1` | matches-main | 2026-09-01 | workflows/somatic/v2/somatic-dna-padding-repair.yaml |
| `somatic-dna-vcf-v2` | matches-main | 2026-04-03 | workflows/somatic/v2/somatic-dna-vcf.yaml |
| `somatic-fusion-fastq-gpu-v2` | drift | 2026-04-29 | workflows/somatic/v2/somatic-fusion-fastq-gpu.yaml |
| `somatic-fusion-fastq-v2` | drift | 2026-04-29 | workflows/somatic/v2/somatic-fusion-fastq.yaml |
| `somatic-rna-bam-v2` | drift | 2026-04-03 | workflows/somatic/v2/somatic-rna-bam.yaml |
| `somatic-rna-fastq-v2` | drift | 2026-04-03 | workflows/somatic/v2/somatic-rna-fastq.yaml |
| `somatic-testing-dna-gpu-workflow` | matches-main | 2024-07-19 | workflows/somatic/somatic-testing-dna-gpu-workflow.yaml |
| `somatic-testing-rna-gpu-workflow` | matches-main | 2024-07-31 | workflows/somatic/somatic-testing-rna-gpu-workflow.yaml |
| `somatic-testing-workflow` | matches-main | 2023-12-24 | workflows/somatic-testing-workflow.yaml |
| `sortmerna-template` | matches-main | 2026-04-29 | tasks/base_tasks/sortmerna.yaml |
| `sprig-os-base-test` | not-in-repo | 2023-07-03 |  |
| `star-align-template` | matches-main | 2026-04-29 | tasks/base_tasks/star-align.yaml |
| `star-fusion-from-junction-template` | matches-main | 2026-04-29 | tasks/base_tasks/star-fusion-from-junction.yaml |
| `star-fusion-template` | matches-main | 2024-01-02 | tasks/base_tasks/star-fusion.yaml |
| `star-fusion-to-oncokb-template` | matches-main | 2024-01-02 | tasks/base_tasks/star-fusion-to-oncokb.yaml |
| `starfusion-detection` | not-in-repo | 2026-02-10 |  |
| `structural-variant-analysis` | not-in-repo | 2025-03-19 |  |
| `subtract-test` | not-in-repo | 2023-07-02 |  |
| `survivor` | not-in-repo | 2025-03-18 |  |
| `test-workflow` | not-in-repo | 2023-06-04 |  |
| `test-workflow-v2` | not-in-repo | 2023-06-05 |  |
| `tiddit` | not-in-repo | 2025-03-18 |  |
| `tps-workflow` | matches-main | 2024-03-26 | workflows/tps-workflow.yaml |
| `trim-galore-template` | drift | 2023-07-02 | tasks/base_tasks/trim-galore.yaml |
| `trim-galore-test` | not-in-repo | 2022-11-17 |  |
| `umi-detect-template` | matches-main | 2026-10-25 | tasks/base_tasks/umi-detect.yaml |
| `umi-tools-dedup-template` | matches-main | 2026-05-07 | tasks/base_tasks/umi-tools-dedup.yaml |
| `update-genotype-snps-template` | matches-main | 2024-12-15 | tasks/prs/update_genotype_snps.yaml |
| `update-run-status-template` | not-in-repo | 2023-10-28 |  |
| `vcf-harmonizer-template` | matches-main | 2026-01-19 | tasks/base_tasks/vcf-harmonizer.yaml |
| `verifybamid-template` | matches-main | 2023-11-14 | tasks/base_tasks/verifybamid.yaml |
| `verifybamid-template-test` | not-in-repo | 2023-11-16 |  |
| `vlod-template` | matches-main | 2023-12-24 | tasks/base_tasks/vlod.yaml |
| `wes-workflow` | matches-main | 2023-11-14 | workflows/wes-workflow.yaml |
| `wgs-workflow` | matches-main | 2024-03-25 | workflows/wgs-workflow.yaml |
| `wrapper-annotation` | matches-main | 2023-11-21 | workflows/wrapper-annotation.yaml |
| `wrapper-batch-somatic-dna` | matches-main | 2026-03-29 | workflows/wrappers/wrapper-batch-somatic-dna.yaml |
| `wrapper-batch-somatic-dna-rna-rna` | matches-main | 2026-03-29 | workflows/wrappers/wrapper-batch-somatic-dna-rna.yaml |
| `wrapper-batch-somatic-fusion-gpu` | matches-main | 2026-04-29 | workflows/wrappers/wrapper-batch-somatic-fusion-gpu.yaml |
| `wrapper-batch-somatic-rna` | matches-main | 2026-03-29 | workflows/wrappers/wrapper-batch-somatic-rna.yaml |
| `wrapper-rna` | matches-main | 2026-03-05 | workflows/wrappers/wrapper-rna.yaml |
| `wrapper-somatic-dna` | matches-main | 2026-03-05 | workflows/wrappers/wrapper-somatic-dna.yaml |
| `wrapper-somatic-dna-rna` | matches-main | 2026-03-05 | workflows/wrappers/wrapper-somatic-dna-rna.yaml |
| `wrapper-somatic-rna` | matches-main | 2026-03-05 | workflows/wrappers/wrapper-somatic-rna.yaml |
| `wrapper-tps` | matches-main | 2026-03-05 | workflows/wrappers/wrapper-tps.yaml |
| `wrapper-wes` | matches-main | 2026-03-05 | workflows/wrappers/wrapper-wes.yaml |
| `wrapper-wgs` | matches-main | 2026-03-05 | workflows/wrappers/wrapper-wgs.yaml |
| `wrapper-workflow` | matches-main | 2023-10-27 | workflows/wrapper-workflow.yaml |
| `wrapper-workflow-prs` | drift | 2024-12-15 | workflows/wrapper-workflow-prs.yaml |
| `wrapper-workflow-somatic` | matches-main | 2024-07-09 | workflows/wrapper-workflow-somatic.yaml |
