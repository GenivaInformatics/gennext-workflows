# Fix Orphaned GPU Pods on Workflow Stop

## Problem

When an Argo workflow is stopped (via `Stop` strategy), pods stuck in the init container phase are not terminated. This leaves GPU resources reserved indefinitely, blocking subsequent workflows from scheduling.

### Observed incident (2026-03-06)

- Workflow `somatic-dna-from-sensor-699cc` was stopped.
- Pod `parabricks-applybqsr-3753418387` was stuck in `Init:0/1` because a hostPath volume mount failed (`tre432.recal is not a file`).
- The pod was never cleaned up, holding 1x NVIDIA RTX A4000 GPU on `gennext-03`.
- The next workflow (`somatic-dna-from-sensor-vl4nb`) could not schedule its GPU step for the entire duration.

### Root cause

Argo's `Stop` strategy sends shutdown signals to running containers. Pods in init/pending phase have no running container to signal, so the kubelet keeps retrying the mount indefinitely. The workflow controller marks the node as stopped in its status but does not force-delete the pod.

## Proposed mitigations

Three independent options, ordered by recommendation. They can be combined.

### Option 1: `activeDeadlineSeconds` on GPU pod specs

Add a hard timeout to GPU-consuming pod templates. If a pod (including init) doesn't complete within the deadline, the kubelet kills it and frees resources.

**Where:** Each Parabricks template that requests `nvidia.com/gpu`.

**Example change (in the WorkflowTemplate pod spec):**

```yaml
podSpecPatch: |
  activeDeadlineSeconds: 7200   # 2 hours — adjust per step
```

**Pros:**
- Simple, built-in Kubernetes mechanism
- Catches any stuck state (init, pending mount, hung process)
- No cluster-wide config change needed

**Cons:**
- Requires choosing a timeout per step (too short = false kills on slow runs)
- Needs to be added to each GPU template individually

### Option 2: Workflow-level pod GC policy

Configure Argo's pod garbage collection to clean up pods from completed/failed/stopped workflows.

**Where:** Workflow defaults in the Argo ConfigMap or on each workflow spec.

**Example (workflow spec level):**

```yaml
spec:
  podGC:
    strategy: OnWorkflowCompletion
    # Optional: keep pods for debugging for a short window
    # deleteDelayDuration: 300s
```

**Example (cluster default in `workflow-controller-configmap`):**

```yaml
podGC: |
  strategy: OnWorkflowCompletion
```

**Pros:**
- Covers all pods in a workflow, not just GPU ones
- Single config change if set at cluster level
- `deleteDelayDuration` allows a window for log collection

**Cons:**
- Pod logs are lost after deletion (mitigate with Argo artifact log archiving)
- Applies to all pods — may affect debugging if no delay is set
- Does not help if the workflow itself is still "Running" and only one step is stuck

### Option 3: Label-based CronJob cleanup

A CronJob that periodically finds and deletes pods from terminated workflows that are still in Init/Pending state.

**Where:** New CronJob in the `argo` namespace.

**Script logic:**

```sh
# Find pods belonging to workflows that are no longer Running
# and are stuck in Init or Pending for > N minutes
kubectl get pods -n argo \
  -l workflows.argoproj.io/workflow \
  --field-selector=status.phase=Pending \
  -o json | \
jq -r '.items[] |
  select(
    (.status.initContainerStatuses[]?.state.waiting != null) and
    (now - (.status.startTime | fromdateiso8601) > 1800)
  ) | .metadata.name' | \
xargs -r kubectl delete pod -n argo
```

**Pros:**
- Catches edge cases the other options miss
- Does not require changes to workflow templates

**Cons:**
- Additional moving part to maintain
- Blunt instrument — needs careful filtering to avoid killing legitimately waiting pods

## Recommendation

**Implement Options 1 + 2 together:**

1. Add `activeDeadlineSeconds` to all GPU-requesting templates as a safety net (Option 1).
2. Set `podGC.strategy: OnWorkflowCompletion` at the cluster level with a `deleteDelayDuration` of 5 minutes (Option 2).

This combination ensures:
- Stuck init pods on active workflows get killed after the deadline (Option 1)
- All pods from stopped/completed workflows get cleaned up promptly (Option 2)
- A 5-minute delay window remains for checking logs before pod deletion

## Files to modify

| Change | File(s) |
|---|---|
| `activeDeadlineSeconds` on GPU templates | `tasks/parabricks_*.yaml` (each GPU-requesting template) |
| Cluster-wide pod GC | Argo `workflow-controller-configmap` in `argo` namespace, or `workflows/*/` workflow specs |

## Verification

1. Submit a GPU workflow and stop it while a pod is in init phase.
2. Confirm the pod is deleted within `activeDeadlineSeconds` or `deleteDelayDuration`.
3. Confirm the GPU becomes available for the next workflow without manual intervention.
