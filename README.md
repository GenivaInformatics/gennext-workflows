# Gennext Workflows

This repository contains the Argo Workflow templates used for the Gennext genomic data processing pipeline.

## Overview

Gennext Workflows is a collection of bioinformatics pipelines implemented using Argo Workflows for Kubernetes. The pipelines are designed to process various types of genomic data, including:

- Whole Exome Sequencing (WES)
- Whole Genome Sequencing (WGS)
- Targeted Panel Sequencing (TPS/CES)
- Somatic Variant Analysis (Tumor/RNA/Tumor-RNA)
- Polygenic Risk Scores (PRS)

## Workflow Components

### Wrapper Workflows

These serve as entry points for the pipeline execution:

- `wrapper-workflow.yaml` - Main wrapper for germline analysis pipelines
- `wrapper-workflow-somatic.yaml` - Wrapper for somatic analysis pipelines
- `wrapper-workflow-prs.yaml` - Wrapper for PRS analysis

### Germline Workflows

- `germline-wes-gpu-workflow.yaml` - Whole Exome Sequencing pipeline with GPU acceleration
- `germline-wgs-gpu-workflow.yaml` - Whole Genome Sequencing pipeline with GPU acceleration
- `wes-workflow.yaml` - Original Whole Exome Sequencing pipeline
- `wgs-workflow.yaml` - Original Whole Genome Sequencing pipeline
- `tps-workflow.yaml` - Targeted Panel Sequencing pipeline

### Somatic Workflows

- `somatic-testing-workflow.yaml` - Main somatic analysis workflow
- `somatic-testing-dna-gpu-workflow.yaml` - GPU-accelerated somatic DNA analysis
- `somatic-testing-rna-gpu-workflow.yaml` - GPU-accelerated somatic RNA analysis

### Annotation Workflows

- `wrapper-annotation.yaml` - Annotation workflow for solo, duo, and trio samples
- `prs-workflow.yaml` - Polygenic Risk Score calculation

### Helper Components

- `preprocess-fastq.yaml` - FASTQ preprocessing steps
- `caller-dv-workflow.yaml`, `caller-hc-workflow.yaml` - Variant calling workflows
- `haplotypecaller-wrapper-workflow.yaml` - Wrapper for GATK HaplotypeCaller

## Dockerfiles

The repository includes Dockerfiles for building container images used in the workflows:

- `dockerfiles/seGMM/Dockerfile` - Container for segmm tools
- `dockerfiles/seqkit_bedtools/Dockerfile` - Container with seqkit and bedtools

## Architecture

The workflows follow a DAG (Directed Acyclic Graph) structure where tasks are organized in a dependency tree. Common patterns include:

1. Configuration and resource allocation
2. Reference data retrieval
3. Input data preparation
4. Alignment and preprocessing
5. Variant calling
6. Post-processing and annotation
7. Quality control
8. Results handling

The pipelines support different input types (FASTQ, BAM, VCF) and automatically select the appropriate workflow based on the input data type.

## GPU Acceleration

Many workflows have GPU-accelerated versions that leverage NVIDIA Parabricks for faster processing, particularly for compute-intensive tasks like alignment and variant calling.

## Requirements

- Kubernetes cluster
- Argo Workflows
- Storage for input data, reference data, and results
- Container images for bioinformatics tools
- GPU resources (for accelerated workflows)

## Development

When modifying these workflows, ensure all dependency requirements are correctly specified and consider the resource requirements of each step.
