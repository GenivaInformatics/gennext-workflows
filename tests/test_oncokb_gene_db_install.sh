#!/usr/bin/env bash
# ref-data/oncokb-gene-db/install.sh in the alpine image copy-to-node.sh uses (busybox),
# on tarballs made from the fixture: checksum checks, the link layout, a re-run, a
# version switch, refusing to touch real files, and --remove.
set -euo pipefail
cd "$(dirname "$0")/.."
image="alpine:3.20@sha256:d9e853e87e55526f6b2917df91a2115c36dd7c696a35be12163d44e6e2a4b6bc"
fix="$PWD/tests/fixtures/oncokb-gene-db"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

tarball() { # version -> $work/<version>.tar.gz
    local v=$1 d="$work/src/$1"
    mkdir -p "$d"
    cp "$fix/hg38_oncokbGene.txt" "$d/hg38_oncokbGene.$v.txt"
    cp "$fix/hg38_oncokbGeneMrna.fa" "$d/hg38_oncokbGeneMrna.$v.fa"
    cp "$fix/hg38_oncokbGene.isoforms.tsv" "$d/hg38_oncokbGene.$v.isoforms.tsv"
    printf 'name=oncokbGene\nversion=%s\n' "$v" > "$d/hg38_oncokbGene.$v.version"
    (cd "$d" && sha256sum hg38_oncokbGene* > SHA256SUMS)
    tar -C "$work/src" -czf "$work/$v.tar.gz" "$v"
}
tarball v1
tarball v2
mkdir -p "$work/db" && touch "$work/db/hg38_refGene.txt"
sha() { sha256sum "$work/$1.tar.gz" | cut -d' ' -f1; }

inst() { # args...
    docker run --rm --network none --memory=256m -u "$(id -u):$(id -g)" \
        -v "$PWD/ref-data/oncokb-gene-db/install.sh:/install.sh:ro" -v "$work:/w" "$image" \
        sh /install.sh "$@" >> "$work/log" 2>&1
}
fail=0
check() { if eval "$2"; then echo "ok: $1"; else echo "FAIL: $1"; tail -5 "$work/log"; fail=1; fi; }
link() { readlink "$work/db/$1"; }

check "wrong sha256 is refused" "! inst /w/v1.tar.gz /w/db 0000"
check "nothing installed after a refusal" "[ \$(ls $work/db | wc -l) -eq 1 ]"
check "install v1" "inst /w/v1.tar.gz /w/db $(sha v1)"
check "links point at v1" "[ \"\$(link hg38_oncokbGene.txt)\" = hg38_oncokbGene.v1.txt ] && [ \"\$(link hg38_oncokbGeneMrna.fa)\" = hg38_oncokbGeneMrna.v1.fa ] && [ \"\$(link hg38_oncokbGene.isoforms.tsv)\" = hg38_oncokbGene.v1.isoforms.tsv ] && [ \"\$(link hg38_oncokbGene.version)\" = hg38_oncokbGene.v1.version ]"
check "content through the links" "cmp -s $work/db/hg38_oncokbGene.txt $fix/hg38_oncokbGene.txt"
check "re-run is a no-op" "inst /w/v1.tar.gz /w/db $(sha v1) && [ \"\$(link hg38_oncokbGene.txt)\" = hg38_oncokbGene.v1.txt ]"
check "no staging left" "! ls -a $work/db | grep -q oncokbGene-install"
check "install v2 switches the links, keeps v1" "inst /w/v2.tar.gz /w/db $(sha v2) && [ \"\$(link hg38_oncokbGene.txt)\" = hg38_oncokbGene.v2.txt ] && [ -f $work/db/hg38_oncokbGene.v1.txt ]"
# an interrupted v1 -> v2 switch: three links on v2, hg38_oncokbGene.txt still on v1,
# a stale temp link and a staging dir left behind. A re-run completes it.
inst /w/v1.tar.gz /w/db "$(sha v1)" >/dev/null 2>&1 || true
for l in hg38_oncokbGeneMrna.fa:hg38_oncokbGeneMrna.v2.fa hg38_oncokbGene.isoforms.tsv:hg38_oncokbGene.v2.isoforms.tsv hg38_oncokbGene.version:hg38_oncokbGene.v2.version; do
    ln -sfn "${l#*:}" "$work/db/${l%%:*}"
done
ln -sfn hg38_oncokbGene.v2.txt "$work/db/.hg38_oncokbGene.txt.new"
mkdir -p "$work/db/.oncokbGene-install.abc123/v2"
check "re-run after an interrupted switch completes it" \
    "inst /w/v2.tar.gz /w/db $(sha v2) && for l in hg38_oncokbGene.txt hg38_oncokbGeneMrna.fa hg38_oncokbGene.isoforms.tsv hg38_oncokbGene.version; do readlink $work/db/\$l | grep -q '\.v2\.' || exit 1; done"
check "leftovers of the interrupted run are gone" "! ls -A $work/db | grep -q -e '^\.oncokbGene-install' -e '\.new\$'"
check "remove takes the links only" "inst --remove /w/db && [ ! -e $work/db/hg38_oncokbGene.txt ] && [ -f $work/db/hg38_oncokbGene.v2.txt ]"
cp "$fix/hg38_oncokbGene.txt" "$work/db/hg38_oncokbGene.txt"
check "a real hg38_oncokbGene.txt is not replaced" "! inst /w/v1.tar.gz /w/db $(sha v1) && [ ! -L $work/db/hg38_oncokbGene.txt ]"
exit $fail
