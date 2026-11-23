#!/usr/bin/env bash
# Copy the oncokbGene tarball to one node's ANNOVAR database and install it there,
# through a short-lived pod pinned to that node (no ssh needed).
#
# usage: copy-to-node.sh <node> <annovardb-host-path> <tarball> <tarball-sha256> [namespace]
#        copy-to-node.sh --remove <node> <annovardb-host-path> [namespace]
#
#   annovardb-host-path  the node's directory that annovar-template mounts
#                        (ref-paths $.other.annovar), e.g. <ref-root>/hg38/other/annovardb
#   namespace            where the pod runs (default gennext-dev)
#
# The pod mounts only that directory, read-write, with a 256Mi memory limit, and is
# deleted at the end. Runs install.sh from this folder; see it for what changes.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
image="alpine:3.20@sha256:d9e853e87e55526f6b2917df91a2115c36dd7c696a35be12163d44e6e2a4b6bc"

remove=
if [ "${1:-}" = "--remove" ]; then remove=1; shift; fi
node=${1:?node}
dbpath=${2:?annovardb host path}
if [ -n "$remove" ]; then ns=${3:-gennext-dev}; else
    tarball=${3:?tarball}; sha=${4:?tarball sha256}; ns=${5:-gennext-dev}
    got=$(sha256sum "$tarball" | cut -d' ' -f1)
    [ "$got" = "$sha" ] || { echo "local tarball sha256 $got, expected $sha" >&2; exit 1; }
fi
pod="oncokbgene-install-$(date +%s)"

kubectl -n "$ns" apply -f - <<EOF
apiVersion: v1
kind: Pod
metadata:
  name: $pod
  labels: {app: oncokbgene-install}
spec:
  nodeName: $node
  restartPolicy: Never
  containers:
    - name: install
      image: $image
      command: [sleep, "1800"]
      resources:
        requests: {cpu: 100m, memory: 64Mi}
        limits: {memory: 256Mi}
      volumeMounts:
        - {name: annovardb, mountPath: /mnt/annovardb}
  volumes:
    - name: annovardb
      hostPath: {path: "$dbpath", type: Directory}
EOF
trap 'kubectl -n "$ns" delete pod "$pod" --wait=false >/dev/null' EXIT
kubectl -n "$ns" wait --for=condition=Ready "pod/$pod" --timeout=180s

kubectl -n "$ns" cp "$here/install.sh" "$pod:/tmp/install.sh"
if [ -n "$remove" ]; then
    kubectl -n "$ns" exec "$pod" -- sh /tmp/install.sh --remove /mnt/annovardb
else
    kubectl -n "$ns" cp "$tarball" "$pod:/tmp/oncokbGene.tar.gz"
    kubectl -n "$ns" exec "$pod" -- sh /tmp/install.sh /tmp/oncokbGene.tar.gz /mnt/annovardb "$sha"
fi
