#!/usr/bin/env bash
# Copy the oncokbGene tarball to one node's ANNOVAR database and install it there,
# through a short-lived pod pinned to that node (no ssh needed).
#
# usage: copy-to-node.sh --context <kube-context> <node> <annovardb-host-path> <tarball> <tarball-sha256> [namespace]
#        copy-to-node.sh --context <kube-context> --remove <node> <annovardb-host-path> [namespace]
#
#   --context            the kubeconfig context to use (required, so the cluster is explicit)
#   annovardb-host-path  the node's directory that annovar-template mounts
#                        (ref-paths $.other.annovar): absolute, ending in /annovardb
#   namespace            where the pod runs (default gennext-dev)
#
# The pod mounts only that directory, read-write, has no service-account token, runs
# with all capabilities dropped but DAC_OVERRIDE (to write a directory it may not
# own), a read-only root filesystem and a 256Mi memory limit, and is deleted at the
# end. Runs install.sh from this folder; see it for what changes. An interrupted run
# is safe to repeat: install.sh clears its own leftovers and completes the switch.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
image="alpine:3.20@sha256:d9e853e87e55526f6b2917df91a2115c36dd7c696a35be12163d44e6e2a4b6bc"
usage() { sed -n '5,6p' "$0" | sed 's/^# *//' >&2; exit 2; }

[ "${1:-}" = "--context" ] && [ -n "${2:-}" ] || usage
ctx=$2
shift 2
remove=
if [ "${1:-}" = "--remove" ]; then remove=1; shift; fi
node=${1:-}
dbpath=${2:-}
[ -n "$node" ] && [ -n "$dbpath" ] || usage
if [ -n "$remove" ]; then ns=${3:-gennext-dev}; else
    tarball=${3:-}; sha=${4:-}; ns=${5:-gennext-dev}
    [ -n "$tarball" ] && [ -n "$sha" ] || usage
fi

# These go into the pod spec: keep them to plain names and paths.
[[ $node =~ ^[a-z0-9]([a-z0-9.-]*[a-z0-9])?$ ]] || { echo "bad node name: $node" >&2; exit 2; }
[[ $ns =~ ^[a-z0-9]([a-z0-9-]*[a-z0-9])?$ ]] || { echo "bad namespace: $ns" >&2; exit 2; }
dbpath=${dbpath%/}
[[ $dbpath =~ ^/[A-Za-z0-9._/-]+/annovardb$ ]] && [[ $dbpath != *..* ]] \
    || { echo "annovardb-host-path must be absolute and end in /annovardb: $dbpath" >&2; exit 2; }

if [ -z "$remove" ]; then
    got=$(sha256sum "$tarball" | cut -d' ' -f1)
    [ "$got" = "$sha" ] || { echo "local tarball sha256 $got, expected $sha" >&2; exit 1; }
fi

k() { kubectl --context "$ctx" -n "$ns" "$@"; }
kubectl config get-contexts "$ctx" >/dev/null || { echo "no kube context $ctx" >&2; exit 2; }
echo "context $ctx, namespace $ns, node $node, $dbpath"

pod="oncokbgene-install-$(date +%Y%m%d%H%M%S)-$(od -An -N4 -tx4 /dev/urandom | tr -d ' ')"
trap 'k delete pod "$pod" --ignore-not-found --wait=false >/dev/null' EXIT
k create -f - <<EOF
apiVersion: v1
kind: Pod
metadata:
  name: "$pod"
  labels: {app: oncokbgene-install}
spec:
  nodeName: "$node"
  restartPolicy: Never
  automountServiceAccountToken: false
  enableServiceLinks: false
  securityContext:
    seccompProfile: {type: RuntimeDefault}
  containers:
    - name: install
      image: "$image"
      command: [sleep, "1800"]
      securityContext:
        allowPrivilegeEscalation: false
        readOnlyRootFilesystem: true
        capabilities: {drop: [ALL], add: [DAC_OVERRIDE]}
      resources:
        requests: {cpu: 100m, memory: 64Mi}
        limits: {memory: 256Mi}
      volumeMounts:
        - {name: annovardb, mountPath: /mnt/annovardb}
        - {name: tmp, mountPath: /tmp}
  volumes:
    - name: annovardb
      hostPath: {path: "$dbpath", type: Directory}
    - name: tmp
      emptyDir: {sizeLimit: 128Mi}
EOF
k wait --for=condition=Ready "pod/$pod" --timeout=180s

k cp "$here/install.sh" "$pod:/tmp/install.sh"
if [ -n "$remove" ]; then
    k exec "$pod" -- sh /tmp/install.sh --remove /mnt/annovardb
else
    k cp "$tarball" "$pod:/tmp/oncokbGene.tar.gz"
    k exec "$pod" -- sh /tmp/install.sh /tmp/oncokbGene.tar.gz /mnt/annovardb "$sha"
fi
