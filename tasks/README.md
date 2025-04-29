# Gennext Workflows Tasks Documentation

## Overview

The `tasks` directory contains the core building blocks of the Gennext genomic data processing pipeline. These are Argo Workflow templates that define individual steps in the genomic analysis pipelines, from basic operations like file handling to complex bioinformatics analyses.

## Task Categories

### Base Tasks

The `base_tasks` directory contains fundamental operations used across various workflows:

- **File Operations**
  - `alter_bam_contig_names.yaml`: Updates contig names in BAM files (e.g., adding "chr" prefix)
  - `bgzip-tabix.yaml`: Compresses files using bgzip and indexes with tabix
  - `dir_cleanup.yaml`: Cleans up temporary directories
  - `samtools.yaml`: Provides templates for common SAMtools operations

- **Quality Control**
  - `fastqc.yaml`: Performs quality control on FASTQ files
  - `mosdepth.yaml`: Calculates sequencing depth and coverage metrics
  - `multiqc.yaml`: Aggregates QC reports into a comprehensive dashboard
  - `qorts.yaml`: Quality control for RNA-seq data
  - `verifybamid.yaml`: Estimates sample contamination

- **Alignment**
  - `bwa-mem2.yaml`: Aligns reads to a reference genome using BWA-MEM2
  - `elprep.yaml`: Performs preprocessing of aligned reads (marking duplicates, etc.)
  - `snap.yaml`: Provides an alternative aligner
  - `seqkit_bedtools.yaml`: Utilities for sequence manipulation

- **Variant Calling & Analysis**
  - `bcftools.yaml`: Operations for VCF file manipulation
  - `freebayes.yaml`: Bayesian variant calling
  - `haplotypecaller.yaml` (in callers/): GATK HaplotypeCaller implementation
  - `mutect2.yaml`: Somatic variant calling
  - `calc-tmb.yaml`: Calculates tumor mutational burden

- **Annotation**
  - `annovar.yaml`, `annovar-to-maf.yaml`: Performs genetic variant annotation
  - `annotsv.yaml`, `annotsv-to-oncokb.yaml`: Annotates structural variants
  - `cancervar.yaml`: Clinical annotation of cancer variants
  - `duckdb.yaml`: Database operations for variant storage and retrieval
  - `vlod.yaml`: Variant likelihood calculation

- **RNA Analysis**
  - `arriba.yaml`, `arriba-to-oncokb.yaml`: Fusion gene detection tools
  - `star-fusion.yaml`, `star-fusion-to-oncokb.yaml`: Alternative fusion detection

- **Structural Variant Analysis**
  - `cnvpytor.yaml`: Copy number variant detection
  - `exomedepth.yaml`: CNV detection for exome data
  - `filter-annotsv.yaml`, `filter-indelible.yaml`, `filter-manta.yaml`: SV filtering
  - `indelible.yaml`: Structural variant detection
  - `manta.yaml`: Structural variant caller
  - `msisensor-pro.yaml`: Microsatellite instability analysis

### GPU-Accelerated Tasks

The `gpu_tasks` directory contains GPU-optimized versions of various tasks:

- `parabricks-bam2fq.yaml`: Convert BAM to FASTQ using GPU acceleration
- `parabricks-bqsr.yaml`: Base quality score recalibration with GPU
- `parabricks-fq2bam.yaml`: Alignment with GPU acceleration
- `parabricks-germline.yaml`: Complete germline variant calling pipeline
- `parabricks-mutect2.yaml`: Somatic variant calling with GPU
- `parabricks-rnafq2bam.yaml`: RNA-seq alignment with GPU
- `parabricks-somatic.yaml`: Somatic variant analysis pipeline
- `parabricks-starfusion.yaml`: Fusion detection with GPU acceleration

### CRAVAT Annotation

The `cravat` directory provides templates for running the Clinical Relevance of Variants Annotation Tool:

- `cravat-solo.yaml`: Annotation for a single sample
- `cravat-duo.yaml`: Annotation for a duo (e.g., tumor-normal pair)
- `cravat-trio.yaml`: Annotation for a trio (e.g., child-mother-father)
- `cravat.yaml`: Combined template with multiple capabilities
- `cravat_nodeSelector.yaml`: Node-specific implementations

### Polygenic Risk Score (PRS) Analysis

The `prs` directory contains templates for calculating polygenic risk scores:

- `gennext_prs.yaml`: Main PRS calculation template
- `genotype_conditional_workflow.yaml`: Workflow for genotype processing
- `update_genotype_snps.yaml`: SNP update for genotype data

### Utility Tasks

- `bed_split.yaml`: Splits BED files for parallel processing
- `compose-read-group-str.yaml`: Creates read group strings for sequencing data
- `copy_inputs_to_scratch.yaml`, `copy_inputs_to_temp.yaml`, `copy_temp_to_output.yaml`: File transfer utilities
- `generate_run_info.yaml`: Generates run metadata
- `get-node-config.yaml`, `get-base-node-config.yaml`: Retrieve node configurations
- `get_abs_ref_dirs.yaml`, `get_rel_ref_dirs.yaml`: Reference directory handling
- `get_run_data.yaml`: Retrieves run data from storage
- `change_run_state.yaml`: Updates the state of a workflow run
- `check-file-path.yaml`: Validates file paths
- `collect_qc_data.yaml`, `collect_somatic_data.yaml`: Data collection utilities

## Architecture

The task templates follow a consistent architecture:

1. **Inputs/Outputs**: Clearly defined parameters for inputs and outputs
2. **Volumes**: Mount points for data files and reference materials
3. **Scripts**: The actual execution logic, typically in bash, Python, or R
4. **Resource Specifications**: CPU, memory, and GPU requests

## Common Patterns

- **File Path Handling**: Most tasks use `sprig` template functions for path manipulation
- **Conditional Execution**: Tasks often include logic to run different steps based on input characteristics
- **GPU Acceleration**: Many compute-intensive tasks have GPU-accelerated alternatives
- **Node Selection**: Some tasks include node selection to target specific hardware

## Integration with Workflows

These tasks are designed to be composable building blocks that can be combined in different ways to create complete analysis pipelines. The main workflows (in the `workflows` directory) reference these tasks to create end-to-end pipelines for different types of genomic analyses.
