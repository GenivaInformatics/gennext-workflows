#!/bin/sh
# Install a versioned oncokbGene tarball into an ANNOVAR database directory.
#
# usage: install.sh <oncokbGene-<version>.tar.gz> <annovardb-dir> <tarball-sha256>
#        install.sh --remove <annovardb-dir>
#
# The tarball holds versioned files (hg38_oncokbGene.<version>.txt, ...). They are
# checked against the tarball's SHA256SUMS and moved into <annovardb-dir>; existing
# files are never overwritten. Then the names ANNOVAR and annovar-template read
# become relative symlinks to them:
#   hg38_oncokbGeneMrna.fa, hg38_oncokbGene.isoforms.tsv, hg38_oncokbGene.version,
#   and hg38_oncokbGene.txt last
# annovar-template runs the oncokbGene protocol only when all four links point at
# one version's files, so a sample annotating during the install, or on a node left
# half-switched, sees either no oncokbGene or a complete one. Each link is replaced with a rename, never removed first. Leftovers of an
# interrupted run (a staging dir, a temp link) are removed first, so a re-run
# completes a half-finished switch.
#
# --remove takes the links away (hg38_oncokbGene.txt first) and keeps the versioned
# files, so later samples skip the protocol again.
# POSIX sh + busybox: runs on a node or in an alpine pod.
set -eu

links="hg38_oncokbGene.txt hg38_oncokbGeneMrna.fa hg38_oncokbGene.isoforms.tsv hg38_oncokbGene.version"

# One install (or --remove) at a time per directory: a second one exits instead of
# racing the first. The kernel drops the lock when this process ends, however it ends.
lock() {
    exec 9>"$1/.oncokbGene.lock"
    flock -n 9 || { echo "another oncokbGene install is running in $1; try again when it ends" >&2; exit 1; }
}

if [ "${1:-}" = "--remove" ]; then
    db=${2:?annovardb dir}
    [ -d "$db" ] || { echo "not a directory: $db" >&2; exit 1; }
    lock "$db"
    for l in $links; do
        if [ -L "$db/$l" ]; then rm -f "$db/$l"; echo "removed $db/$l"; fi
    done
    exit 0
fi

tarball=${1:?tarball}
db=${2:?annovardb dir}
want=${3:?tarball sha256}
[ -d "$db" ] || { echo "not a directory: $db" >&2; exit 1; }
lock "$db"
for l in $links; do
    if [ -e "$db/$l" ] && [ ! -L "$db/$l" ]; then
        echo "$db/$l exists and is not a symlink; not touching it" >&2
        exit 1
    fi
done

got=$(sha256sum "$tarball" | cut -d' ' -f1)
[ "$got" = "$want" ] || { echo "tarball sha256 $got, expected $want" >&2; exit 1; }

# leftovers of an interrupted run (staging dirs, temp links) go first
rm -rf "$db"/.oncokbGene-install.*
for l in $links; do rm -f "$db/.$l.new"; done
stage=$(mktemp -d "$db/.oncokbGene-install.XXXXXX")
trap 'rm -rf "$stage"' EXIT
tar -xzf "$tarball" -C "$stage"
set -- "$stage"/*/
[ "$#" -eq 1 ] && [ -d "$1" ] || { echo "tarball must hold one version directory" >&2; exit 1; }
src=${1%/}
version=$(basename "$src")
(cd "$src" && sha256sum -c SHA256SUMS >/dev/null) || { echo "SHA256SUMS check failed" >&2; exit 1; }

txt="hg38_oncokbGene.$version.txt"
mrna="hg38_oncokbGeneMrna.$version.fa"
iso="hg38_oncokbGene.$version.isoforms.tsv"
ver="hg38_oncokbGene.$version.version"
for f in "$txt" "$mrna" "$iso" "$ver"; do
    [ -f "$src/$f" ] || { echo "missing $f in the tarball" >&2; exit 1; }
    if [ -e "$db/$f" ]; then
        cmp -s "$src/$f" "$db/$f" || { echo "$db/$f exists with other content" >&2; exit 1; }
    else
        chmod 0444 "$src/$f"
        mv "$src/$f" "$db/$f"
    fi
done

link() { # target name: atomic replace via rename
    rm -f "$db/.$2.new"
    ln -s "$1" "$db/.$2.new"
    mv -f "$db/.$2.new" "$db/$2"
}
link "$mrna" hg38_oncokbGeneMrna.fa
link "$iso" hg38_oncokbGene.isoforms.tsv
link "$ver" hg38_oncokbGene.version
link "$txt" hg38_oncokbGene.txt

# all four must now point at this version (annovar-template checks the same)
for pair in "hg38_oncokbGene.txt:$txt" "hg38_oncokbGeneMrna.fa:$mrna" "hg38_oncokbGene.isoforms.tsv:$iso" "hg38_oncokbGene.version:$ver"; do
    [ "$(readlink "$db/${pair%%:*}")" = "${pair#*:}" ] || { echo "${pair%%:*} does not point at ${pair#*:}" >&2; exit 1; }
done
echo "installed $version in $db:"
ls -l "$db"/hg38_oncokbGene* | sed 's/^/  /'
