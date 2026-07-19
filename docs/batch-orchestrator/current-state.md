# Current State

## What's Been Built and Deployed

### Backend (`gennext-backend` — branch `feature/batch-sample-processor`)

**DB Changes (applied to production):**
- Added `partial` to `run_status` PostgreSQL enum
- Added `status` column to `samples` table (type `run_status`, default `submitted`)
- Backfilled all existing samples to match their parent run's status

**Code Changes:**
- `CreateSample` accepts optional `data_type` per sample
- API responses include both `id` and `sampleUuid` fields (backwards compat)
- `Sample` and `SummarySample` include `status` field
- New endpoint: `POST /api/v0/runs/{run_id}/samples/{sample_uuid}/state`
- `compute_aggregate_run_status()` derives run status from sample statuses
- `RunStatus::from_str()` handles `partial` and `waiting_for_data`
- All `RunStatus` match statements handle `Partial` variant

**Deployed:** `bergun/gennext-backend:0.18.0` to both `gennext-dev` and `gennext-prod`

### Edge API (`gennext-api-edge` — branch `feature/batch-sample-processor`)

**Code Changes:**
- `src/utils/paths.rs` — `resolve_output_path()` and `build_output_dir()` try `{run_id}/{sample_uuid}/{path}` first, fall back to legacy `{run_id}/{path}`
- All request/query models accept optional `sample_uuid`
- ~20 path constructions updated across `queries/files.rs`, `routes/variants.rs`, `routes/somatic.rs`, `routes/coverage.rs`, `routes/qualities.rs`
- Fusions endpoint accepts `sample_uuid` query param
- `TinySample` in cluster model uses `#[serde(alias = "id")] sample_uuid`

**Deployed:** `bergun/gennext-api-edge:0.50.0-beta` to `gennext-01-gennext-api-edge` (gennext-prod)

### File Transfer (`gennext-ft` — branch `feature/batch-sample-processor`)

**Code Changes:**
- `FileInfo`, `UploadV2Metadata`, `PartialUpload` — added `SampleID` field
- `pathKey()`, `InputFilePath()`, `GetFileByPath()`, `SaveFileToInput()` — accept `sampleID` param
- Upload handler passes `sampleID` through to storage
- Discovery recursion handles sample UUID subdirectories
- HTTP upload reads optional `sampleId` form value

**Deployed:** `bergun/gennext-ft-edge:2.4.0-beta` to `gennext-ft-v2-edge` (gennext-prod, all nodes)

### Workflows (`gennext-k8s-deployments` — branch `feature/batch-sample-processor`)

**Modified Tasks:**
- `tasks/get_run_data.yaml` — patches ALL samples (not just `samples[0]`)
- `tasks/get-node-config.yaml` — outputs `max-concurrent-samples` from node annotation

**New Tasks:**
- `tasks/change_sample_state.yaml` — calls backend per-sample state endpoint

**New WorkflowTemplates:**
- `workflows/per-sample-somatic-dna.yaml` — child workflow for DNA samples
- `workflows/per-sample-somatic-rna.yaml` — child workflow for RNA samples
- `workflows/per-sample-somatic-dna-rna.yaml` — child workflow for DNA+RNA samples
- `workflows/wrappers/wrapper-batch-somatic-dna.yaml` — orchestrator (DNA)
- `workflows/wrappers/wrapper-batch-somatic-rna.yaml` — orchestrator (RNA)
- `workflows/wrappers/wrapper-batch-somatic-dna-rna.yaml` — orchestrator (DNA+RNA)

**New Resources:**
- `deployments/batch-sample-semaphore.yaml` — ConfigMap for concurrency control
- RBAC Role + RoleBinding for `gennext-workflows` SA (create/get/watch/list workflows + create/get/update/patch configmaps)

**All applied to `argo` namespace.**

### Frontend (`gennext-frontend` — branch `feature/batch-sample-processor`)

- Service layer sends `sampleUuid` in all edge API calls
- Navigation carries `sampleUuid` through URLs
- Detail pages resolve selected sample via `findSelectedSample()` instead of `samples[0]`
- Pushed to remote. Frontend team has Linear issue GNV-299 for review.

### Infrastructure

- Node annotation `gennext.bio/max-concurrent-samples=2` set on `gennext-01`
- Sensor triggers NOT yet switched — old wrappers still active for production

## Test Data

### Test batch run (existing data, symlinked)

- Run ID: `9fefe802-4fac-4699-a4ff-ece89de21281`
- Samples: BRG-B1-S2 (`59372b2f-...`) and BRG-B1-S5 (`1a910ede-...`)
- Data: symlinked from original single-sample runs
- Status: completed (both samples)
- Used to verify edge API path resolution with `sampleUuid`

### Test batch run (real processing)

- Run ID: `76b2eb95-544a-4d0a-9530-bff3dd75f30c`
- Samples: BATCH-S1 (`c9961cf4-...`) and BATCH-S2 (`0b830101-...`)
- Data: real FASTQ files in per-sample UUID subdirectories
- Used to test the orchestrator workflow submission
- Child workflows were created and ran the real pipeline

## Known Issues

1. **Orchestrator state sync**: Argo's `resource` template is one-shot — once it sees a child fail, it records failure even if the child is manually resumed. This is why we need a custom orchestrator.

2. **Sensor triggers not switched**: Production still uses old single-sample wrappers. Switch requires changing 3 entries in `events/gennext-workflow-from-sensor.yaml`.

3. **Frontend file transfer**: Upload metadata doesn't yet carry `sampleId` (UUID) for per-sample directory routing. This is tracked in GNV-299.
