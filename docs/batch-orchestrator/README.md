# Batch Sample Orchestrator

## What This Is

A batch processing system for the Gennext genomics platform that allows multiple samples to be processed in a single run. Previously, each run processed exactly one sample. Now, an orchestrator workflow spawns independent child workflows per sample with concurrency control.

## Why We Need a Custom Orchestrator

### The Problem with Argo's Built-in Patterns

We initially used Argo's `resource` template with `action: create` to spawn child Workflow resources. This worked for the basic case but has fundamental limitations:

1. **One-shot state tracking**: Once the orchestrator's `resource` template sees a child workflow reach `Failed`, it marks the step as failed and moves on. If someone manually resumes the child workflow (e.g., after fixing a transient issue), the orchestrator has already recorded the failure and won't update.

2. **No state sync**: The orchestrator and child workflows are decoupled after creation. There's no mechanism to re-check child status, handle manual interventions, or coordinate retries.

3. **YAML manifest limitations**: Passing complex JSON data (like sample metadata) through the `resource` template's inline YAML manifest causes parsing errors due to nested quotes. We solved this by having children read the JSON file directly, but it added complexity.

### What We Want Instead

A custom orchestrator that:
- Submits child workflows as independent Workflow resources (each visible in Argo UI)
- Continuously monitors their state (not one-shot)
- Supports manual resume/retry of individual samples without losing sync
- Reports aggregate status back to the backend API
- Respects concurrency limits per node

## Architecture Overview

See [architecture.md](./architecture.md) for the full system design.

## Current State

See [current-state.md](./current-state.md) for what's been built and deployed.

## Implementation Guide

See [implementation-guide.md](./implementation-guide.md) for how to build the custom orchestrator.
