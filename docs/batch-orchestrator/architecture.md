# Architecture

## Platform Overview

Gennext is a genomics analysis platform with 5 repositories:

| Repo | Role | Language |
|------|------|----------|
| `gennext-backend` | Vertex logic, DB, run/sample CRUD | Rust (Axum, SeaORM, PostgreSQL) |
| `gennext-api-edge` | Edge node data API (variants, coverage, files) | Rust (Axum, DuckDB) |
| `gennext-ft` | P2P/relay file transfer | Go (libp2p) |
| `gennext-k8s-deployments` | K8s deployments, Argo Workflow templates | YAML |
| `gennext-frontend` | UI and file transfer client | React/Vite |

All repos have branch `feature/batch-sample-processor` with the batch changes.

## Naming Convention (Critical)

Two different sample identifiers — never confuse them:

| Field | Type | Source | Used For |
|-------|------|--------|----------|
| `sample.id` / `sampleUuid` | UUID | DB-assigned | Directory paths, internal routing, API addressing |
| `sample.sampleId` | String | User-defined | Filenames, tool CLI args (--sample-name), display |

## Directory Structure

```
{base}/{user_id}/{run_id}/                          ← run level
  {run_id}.json                                      ← run metadata (all samples)
  {sample.id UUID}/                                  ← per-sample level
    {sample.sampleId}_R1.fastq.gz                    ← input files use sampleId in name
    results/
      db/{sample.sampleId}_results.duckdb
      snv/{sample.sampleId}_mutect2_vlod.vcf.gz
      sv/{sample.sampleId}_tumorSV.vcf.gz
      cnv/{sample.sampleId}_recal.cns
      metrics/msi/{sample.sampleId}_msi_filtered.tsv
      fusions/fusion-report/fusions.json
    bam/{sample.sampleId}_recal.bam
    metrics/mosdepth/{sample.sampleId}.regions.bed.gz
```

## Run Data JSON

Located at `{input_dir}/{user_id}/{run_id}/{run_id}.json`. Created by the frontend during file upload. Contains all samples:

```json
{
  "userId": "3c41bacd-...",
  "runId": "76b2eb95-...",
  "processType": "som-test-tumor",
  "dataType": "fastq",
  "refVersion": "hg38",
  "enrichmentKit": {"company": "Agilent", "title": "...", "filename": "A3416642"},
  "samples": [
    {
      "id": "c9961cf4-...",
      "sampleId": "BATCH-S1",
      "files": {
        "fqR1": {"name": "SRR100910_R1.fastq.gz"},
        "fqR2": {"name": "SRR100910_R2.fastq.gz"}
      }
    },
    {
      "id": "0b830101-...",
      "sampleId": "BATCH-S2",
      "files": {
        "fqR1": {"name": "USB81218_R1.fastq.gz"},
        "fqR2": {"name": "USB81218_R2.fastq.gz"}
      }
    }
  ]
}
```

## Workflow Flow

### Current (working, but limited)

```
Sensor event → wrapper-batch-somatic-dna (orchestrator)
  ├── get-node-config        → reads K8s node annotations for paths + concurrency
  ├── get-run-data           → reads {run_id}.json, patches all samples
  ├── update-semaphore       → writes max-concurrent to ConfigMap
  ├── extract-samples        → outputs lightweight [{id, sampleId}, ...] array
  ├── fan-out-samples        → withParam over samples
  │   └── submit-sample-workflow (resource:create + semaphore)
  │       └── per-sample-somatic-dna (independent child Workflow)
  │           ├── set-sample-running    → POST /runs/{rid}/samples/{sid}/state
  │           ├── load-and-prepare      → reads JSON, finds sample, builds paths
  │           ├── run-somatic-dna       → existing pipeline (unchanged)
  │           ├── set-sample-completed
  │           └── onExit: set-sample-error
  └── determine-run-status
```

### Desired (custom orchestrator)

```
Sensor event → custom-orchestrator
  ├── Read run data JSON
  ├── For each sample:
  │   ├── Submit child Workflow (Argo API or kubectl)
  │   ├── Track child state continuously (poll or watch)
  │   ├── Handle manual resume/retry (re-sync state)
  │   └── Report per-sample status to backend
  ├── Respect concurrency limit (from node annotation)
  └── Set aggregate run status when all children finish
```

## Sample Status Model

Samples have their own `status` column in the database (same enum as runs):
- `submitted` → `running` → `completed` or `error`

Run status is the aggregate:
- All submitted → `submitted`
- Any running → `running`
- All completed → `completed`
- All error → `error`
- Mix of completed + error → `partial`

### Backend Endpoint

```
POST /api/v0/runs/{run_id}/samples/{sample_uuid}/state
Body: { "token": "<cluster_token>", "state": "running" | "completed" | "error" }
```

Authenticated via `GENNEXT_CLUSTER_TOKEN` (same as run state endpoint). Updates the sample status and recomputes the aggregate run status automatically.

## Concurrency Control

- Node annotation: `gennext.bio/max-concurrent-samples` (default: `1`)
- Read by `get-node-config` task
- Currently written to ConfigMap `batch-sample-semaphore` for Argo's semaphore mechanism
- The custom orchestrator should read this directly and manage its own concurrency

## Existing Per-Sample Workflow Templates

These work and don't need changes:

| Template | Process Type | Calls |
|----------|-------------|-------|
| `per-sample-somatic-dna` | `som-test-tumor` | `somatic-testing-workflow/main-dna` |
| `per-sample-somatic-rna` | `som-test-rna` | `somatic-testing-workflow/main-rna` |
| `per-sample-somatic-dna-rna` | `som-test-tumor-rna` | `somatic-testing-workflow/main-dna-rna` |

Each accepts simple string params: `sample-uuid`, `sample-id`, `run-id`, `user-id`, `input-rel-dir`, `output-rel-dir`, `temp-rel-dir`, `scratch-rel-dir`, `ref-base-dir`, `thread-count`.

The `load-and-prepare` step inside each child:
1. Mounts `input-rel-dir` as hostPath
2. Reads `{run-id}.json`
3. Finds the sample matching `sample-uuid`
4. Builds single-sample `run-data` with `samples = [matched_sample]`
5. Appends `/{sample-uuid}/` to all directory paths

This means existing pipeline templates (`somatic-testing-workflow`, etc.) remain completely unchanged — they still reference `$.samples[0]` and it works because the child workflow's `load-and-prepare` reconstructs the data with only the target sample.
