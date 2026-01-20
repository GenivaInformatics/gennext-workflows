# Germline Variant Processing Pipeline

## Overview

This directory contains Argo Workflows for germline variant calling. The pipelines support multiple input types (FASTQ, BAM, VCF) and sequencing assays (WGS, WES, TPS, RNA).

## Workflow Files

| File | Description |
|------|-------------|
| `germline-wgs-gpu-workflow.yaml` | Whole Genome Sequencing with GPU acceleration |
| `germline-wes-gpu-workflow.yaml` | Whole Exome Sequencing with GPU acceleration |
| `germline-tps-gpu-workflow.yaml` | Targeted Panel Sequencing with GPU acceleration |
| `germline-rna-gpu-workflow.yaml` | RNA-seq workflow with GPU acceleration |
| `germline-wgs-workflow.yaml` | CPU-only WGS workflow variant |

## Main Pipeline Flow (WGS)

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              INPUT                                          │
│                          (FASTQ files)                                      │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           trim-galore                                       │
│                      (Quality trimming)                                     │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                     pb-fq2bam (Parabricks GPU)                              │
│              ┌────────────────────────────────────────┐                     │
│              │  bwa-mem2 alignment + MarkDuplicates   │                     │
│              └────────────────────────────────────────┘                     │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      pb-bqsr → pb-applybqsr                                 │
│              ┌────────────────────────────────────────┐                     │
│              │   Base Quality Score Recalibration     │                     │
│              └────────────────────────────────────────┘                     │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                    ┌───────────────┼───────────────┐
                    ▼               ▼               ▼
         ┌──────────────┐  ┌──────────────┐  ┌──────────────┐
         │ DeepVariant  │  │    Manta     │  │  CNVpytor    │
         │  (GPU)       │  │              │  │              │
         │ SNVs/Indels  │  │ Structural   │  │    CNVs      │
         │              │  │  Variants    │  │              │
         └──────────────┘  └──────────────┘  └──────────────┘
                    │               │               │
                    │               ▼               │
                    │      ┌──────────────┐         │
                    │      │   AnnotSV    │         │
                    │      │ SV Annotation│         │
                    │      └──────────────┘         │
                    │               │               │
                    └───────────────┼───────────────┘
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                              CRAVAT                                         │
│     (Variant Annotation: ClinVar, CADD, gnomAD4, AlphaMissense, etc.)       │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                              OUTPUT                                         │
│               VCF + gVCF + Annotated variants + QC metrics                  │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Simplified Linear View

```
FASTQ → trim-galore → bwa-mem2/MarkDups → BQSR → DeepVariant → CRAVAT → VCF
                         (pb-fq2bam)      (pb-bqsr)
```

## Pipeline Steps

| Step | Tool | Purpose |
|------|------|---------|
| Trimming | trim-galore | Remove low-quality bases and adapters |
| Alignment + Dedup | pb-fq2bam | GPU-accelerated bwa-mem2 + MarkDuplicates |
| BQSR | pb-bqsr/applybqsr | Recalibrate base quality scores using known variants |
| SNV/Indel Calling | DeepVariant | Deep learning-based variant caller (GPU) |
| SV Calling | Manta | Structural variant detection |
| CNV Calling | CNVpytor | Copy number variant detection |
| Annotation | CRAVAT | Functional annotation with multiple databases |

## Input Types

The workflows automatically route based on input data type:

- **FASTQ**: Full pipeline (trim → align → BQSR → variant calling → annotation)
- **BAM**: Converts to FASTQ first, then follows full pipeline
- **VCF**: Skips to annotation only (CRAVAT)

## QC Checkpoints

Quality control is performed at multiple stages:

- **FastQC**: Read quality assessment
- **VerifyBamID**: Contamination detection
- **Mosdepth**: Coverage depth analysis
- **MultiQC**: Aggregated QC report
