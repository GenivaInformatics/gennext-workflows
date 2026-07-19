# Post-Processing Workflow Templates

Post-processing workflows run batch-level analysis (e.g. CNV) across groups of samples that share the same enrichment kit, after all sample workflows in a run have completed.

## Current templates

| Template | WorkflowTemplate name | Status |
|---|---|---|
| `workflows/post-process-cnv.yaml` | `post-process-cnv` | Dummy (sleep steps, replace with real CNV) |

## Template contract

Every PP workflow template must follow this contract for the orchestrator and backend to manage it correctly.

### Required parameters

| Parameter | Description | Example |
|---|---|---|
| `run-id` | Run UUID | `76b2eb95-544a-4d0a-9530-bff3dd75f30c` |
| `user-id` | User UUID | `3c41bacd-22a6-4639-ab15-2fc9141a8d4f` |
| `enrichment-kit` | Kit display name | `Agilent - SureSelect Cancer CGP DNA+RNA Assay` |
| `node-name` | Compute node hostname | `gennext-01` |
| `ref-base-dir` | Reference data base path | `/runspace/gennext/data/ref/` |
| `thread-count` | Thread count for tools | `20` |
| `samples` | JSON string array of sample refs | See below |

### Samples parameter format

```json
[
  {"sample-uuid": "0b830101-...", "sample-id": "BATCH-S1", "output-dir": "/runspace/.../output/user-id/run-id/"},
  {"sample-uuid": "c9961cf4-...", "sample-id": "BATCH-S2", "output-dir": "/runspace/.../output/user-id/run-id/"}
]
```

Note: `samples` is passed as a **JSON string**, not a nested object. Parse with `json.loads()` / equivalent inside workflow steps.

### Required callbacks

The workflow must call the backend on completion and failure:

**On success** — a DAG task after all processing:
```
POST /api/v0/runs/{run-id}/postprocess-complete
Body: {"token": "$GENNEXT_CLUSTER_TOKEN", "kit": "{enrichment-kit}", "step": "{pp-type}"}
```

This should also set affected samples back to `completed`:
```
POST /api/v0/runs/{run-id}/samples/{sample-uuid}/state
Body: {"token": "$GENNEXT_CLUSTER_TOKEN", "state": "completed"}
```

**On failure** — via `onExit` handler with `when: "{{workflow.status}} == Failed"`:
```
POST /api/v0/runs/{run-id}/postprocess-error
Body: {"token": "$GENNEXT_CLUSTER_TOKEN", "kit": "{enrichment-kit}", "step": "{pp-type}"}
```

This should also set affected samples back to `completed` (their per-sample data is valid):
```
POST /api/v0/runs/{run-id}/samples/{sample-uuid}/state
Body: {"token": "$GENNEXT_CLUSTER_TOKEN", "state": "completed"}
```

### Required annotations

```yaml
metadata:
  annotations:
    gennext.bio/workflow-type: post-process
    gennext.bio/post-process-step: cnv  # or whatever the step name is
```

These are informational — the orchestrator sets its own labels when submitting. But they help identify the template's purpose.

### Cluster token

Access via environment variable from the shared secret:
```yaml
env:
  - name: GENNEXT_CLUSTER_TOKEN
    valueFrom:
      secretKeyRef:
        name: gennext-cluster-secrets
        key: gennext_cluster_token
```

The backend address for callbacks is `gennext-backend.gennext-prod.svc.cluster.local:3030` (adjust namespace as needed).

## DAG structure pattern

```
set-samples-postprocessing (optional — backend already sets this)
  → gather/prepare step
    → actual analysis (CNV, etc.)
      → write results
        → set-postprocess-completed (calls backend)

onExit:
  → set-postprocess-error (calls backend, when: Failed)
```

The `set-samples-postprocessing` step in the dummy template is redundant — the backend's `batch_complete` handler already sets samples to `PostProcessing` before submitting the workflow. It can be removed in the real implementation or kept as a safety net.

## Adding a new PP template

1. Create the WorkflowTemplate YAML following the contract above
2. Apply to cluster: `kubectl apply -f workflows/post-process-{type}.yaml -n argo`
3. Add to orchestrator's `postProcessTemplateMap` in `internal/argo/submit.go`:
   ```go
   var postProcessTemplateMap = map[string]string{
       "cnv": "post-process-cnv",
       "new-type": "post-process-new-type",  // add here
   }
   ```
4. Add the step name to the enrichment kit's `post_processing_steps` column in the DB:
   ```sql
   UPDATE enrichment_kits SET post_processing_steps = '{cnv,new-type}' WHERE name = '...';
   ```
5. Steps are executed sequentially per kit group — order in the array matters

## Auto-retry

The orchestrator auto-retries failed PP workflows up to 3 times (same as sample workflows). The `onExit` error handler only fires after the final retry is exhausted. Retry count is tracked via the `gennext.bio/retry-count` workflow annotation.
