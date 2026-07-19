# Implementation Guide: Custom Batch Orchestrator

## Goal

Replace the Argo `resource` template-based orchestrator with a custom orchestrator that continuously tracks child workflow state, supports manual interventions, and maintains sync between the orchestrator and child workflows.

## Approach Options

### Option A: Argo Workflow with Polling Loop

An Argo Workflow template that:
1. Submits child workflows via `kubectl`
2. Runs a polling loop that checks child status every N seconds
3. Updates sample states in the backend as children progress
4. Handles completion, failure, and manual resume

**Pros:** Stays within the Argo ecosystem. No new infrastructure.
**Cons:** Polling is wasteful. Long-running pods for the orchestrator.

### Option B: Standalone Controller (Pod/Deployment)

A lightweight service (Go or Python) deployed as a K8s Deployment that:
1. Watches for orchestrator trigger events (webhook or CRD)
2. Submits child workflows via K8s API
3. Watches child workflow status via K8s watch API (not polling)
4. Updates backend sample states in real-time
5. Handles resume/retry by re-watching resumed workflows

**Pros:** Event-driven (no polling). Clean separation. Can run permanently.
**Cons:** New service to maintain. More infrastructure.

### Option C: Argo Workflow + Kubernetes Watch Script

Similar to Option A but instead of polling, uses a script step with `kubectl get --watch`:
1. Submits all child workflows
2. Watches their status changes via streaming watch
3. Reacts to state transitions immediately

**Pros:** Near real-time. Still Argo-native.
**Cons:** Complex shell scripting. Watch can disconnect.

## Recommended: Option B (Standalone Controller)

The cleanest approach for our needs. It can be a simple Go program using the K8s client-go library.

## How Child Workflows Are Submitted

The child workflow templates (`per-sample-somatic-dna`, etc.) accept simple string params:

```yaml
apiVersion: argoproj.io/v1alpha1
kind: Workflow
metadata:
  generateName: "sample-dna-{sample-uuid}-"
  namespace: argo
  labels:
    gennext.bio/orchestrator: "{orchestrator-run-id}"
    gennext.bio/sample-id: "BATCH-S1"
    gennext.bio/run-id: "{run-id}"
  annotations:
    gennext.bio/sample-uuid: "{sample-uuid}"
    gennext.bio/sample-id: "BATCH-S1"
spec:
  serviceAccountName: gennext-workflows
  nodeSelector:
    kubernetes.io/hostname: "gennext-01"
  podGC:
    strategy: OnPodSuccess
  entrypoint: main
  arguments:
    parameters:
      - name: sample-uuid
        value: "c9961cf4-4d90-44eb-88d9-cb03a6943bd4"
      - name: sample-id
        value: "BATCH-S1"
      - name: run-id
        value: "76b2eb95-544a-4d0a-9530-bff3dd75f30c"
      - name: user-id
        value: "3c41bacd-22a6-4639-ab15-2fc9141a8d4f"
      - name: input-rel-dir
        value: "/runspace/gennext/data/input/3c41bacd-.../76b2eb95-.../"
      - name: output-rel-dir
        value: "/runspace/gennext/data/output/3c41bacd-.../76b2eb95-.../"
      - name: temp-rel-dir
        value: "/runspace/gennext/data/temp/3c41bacd-.../76b2eb95-..._VInHgs/"
      - name: scratch-rel-dir
        value: "/runspace/gennext/data/scratch/3c41bacd-.../76b2eb95-..._VInHgs/"
      - name: ref-base-dir
        value: "/runspace/gennext/data/ref/"
      - name: thread-count
        value: "20"
  workflowTemplateRef:
    name: per-sample-somatic-dna
```

The child workflow internally:
1. Calls `change-sample-state` to set status to `running`
2. Runs `load-and-prepare` which reads `{run-id}.json` from disk, finds its sample by UUID, builds per-sample paths
3. Runs the actual pipeline (`somatic-testing-workflow/main-dna`)
4. Calls `change-sample-state` to set status to `completed`
5. On failure (via `onExit`), sets status to `error`

## What the Orchestrator Must Do

### 1. Receive Trigger

Currently triggered by Argo Events sensor. The orchestrator receives:
- `node-name`: Which K8s node to run on
- `user-id`: User UUID
- `job-id`: Run UUID
- `process-type`: `som-test-tumor` | `som-test-rna` | `som-test-tumor-rna`

### 2. Read Configuration

- Read node config (same as `get-node-config` task):
  - Directory paths from node annotations
  - `max-concurrent-samples` from node annotation (default: 1)
- Read run data JSON from `{input-dir}/{user-id}/{run-id}/{run-id}.json`
- Patch samples (same as `get-run-data` task): parse files strings, default regionsFile, somaticInfo

### 3. Submit Child Workflows

For each sample in the JSON:
- Create a Workflow resource via K8s API
- Use the `per-sample-somatic-{dna|rna|dna-rna}` WorkflowTemplate (based on process-type)
- Pass the simple string parameters listed above
- Label with `gennext.bio/orchestrator`, `gennext.bio/run-id`, `gennext.bio/sample-id`
- Respect concurrency limit: only N children running at a time

### 4. Monitor Children

- Watch child Workflow status changes via K8s watch API
- On state change:
  - `Running` → sample is processing (state already set by child itself)
  - `Succeeded` → sample completed (state already set by child itself)
  - `Failed` → sample errored (state already set by child's onExit)
- If a child is manually resumed (status changes from `Failed` back to `Running`):
  - Detect this and update tracking
  - Wait for the resumed child to reach a terminal state
- When a semaphore slot frees up, submit the next pending sample

### 5. Report Aggregate Status

After all children reach terminal state (Succeeded/Failed):
- The backend already computes aggregate run status from sample statuses
- The orchestrator can optionally call `POST /runs/{run-id}/state` as a safety net

## Key APIs

### Backend API

```
# Per-sample state (called by child workflows, but orchestrator can also call)
POST /api/v0/runs/{run_id}/samples/{sample_uuid}/state
Body: { "token": "<GENNEXT_CLUSTER_TOKEN>", "state": "running|completed|error" }

# Run state (backup, normally computed by aggregate)
POST /api/v0/runs/{run_id}/state
Body: { "token": "<GENNEXT_CLUSTER_TOKEN>", "state": "completed|error|partial" }
```

### Kubernetes API

```
# Submit child workflow
kubectl create -f workflow.yaml

# Watch child workflows
kubectl get workflows -n argo -l gennext.bio/run-id={run-id} --watch

# Get child workflow status
kubectl get workflow {name} -n argo -o jsonpath='{.status.phase}'
```

## Node Configuration

Directory paths and concurrency come from K8s node annotations:

```
gennext.bio/input-base-dir=/runspace/gennext/data/input/
gennext.bio/output-base-dir=/runspace/gennext/data/output/
gennext.bio/ref-base-dir=/runspace/gennext/data/ref/
gennext.bio/temp-base-dir=/runspace/gennext/data/temp/
gennext.bio/scratch-base-dir=/runspace/gennext/data/scratch/
gennext.bio/thread-count=20
gennext.bio/max-concurrent-samples=2
```

The input/output dirs are constructed as: `{base-dir}/{user-id}/{run-id}/`
The temp/scratch dirs add a random suffix: `{base-dir}/{user-id}/{run-id}_{random}/`

## RBAC

The orchestrator's service account needs:

```yaml
rules:
  - apiGroups: ["argoproj.io"]
    resources: ["workflows"]
    verbs: ["create", "get", "watch", "list"]
  - apiGroups: [""]
    resources: ["configmaps"]
    verbs: ["create", "get", "update", "patch"]
  - apiGroups: [""]
    resources: ["nodes"]
    verbs: ["get"]
```

## Environment

- Cluster: RKE2 v1.28
- Argo Workflows: v3.4.11
- Namespace: `argo` (workflows), `gennext-dev`/`gennext-prod` (services)
- Backend URL (from within cluster): `gennext-backend.gennext-dev.svc.cluster.local:3030`
- Cluster token: stored in secret `gennext-cluster-secrets` key `gennext_cluster_token`

## Testing

### Manual Test Submission

```bash
argo submit --from workflowtemplate/wrapper-batch-somatic-dna \
  -n argo \
  -p node-name=gennext-01 \
  -p user-id=3c41bacd-22a6-4639-ab15-2fc9141a8d4f \
  -p job-id=76b2eb95-544a-4d0a-9530-bff3dd75f30c \
  -p process-type=som-test-tumor
```

### Verify Child Workflows

```bash
# List children of an orchestrator
argo list -n argo -l "gennext.bio/run-id=76b2eb95-544a-4d0a-9530-bff3dd75f30c"

# Check sample annotations
kubectl get workflow {name} -n argo -o jsonpath='{.metadata.annotations}'
```

### Test Run Data

Run ID `76b2eb95-544a-4d0a-9530-bff3dd75f30c` has 2 samples with real FASTQ data on gennext-01. The JSON and input files are in place. Use this for end-to-end testing.
