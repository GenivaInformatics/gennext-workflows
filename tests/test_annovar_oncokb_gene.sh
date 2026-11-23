#!/usr/bin/env bash
# annovar-template's optional oncokbGene protocol and annovar-to-maf, run in their
# production images on the synthetic fixture (tests/fixtures/oncokb-gene-db: MANE
# refGene/ensGene and OncoKB-isoform oncokbGene for 13 genes, 17 made-up variants).
#
# The baseline is tests/fixtures/oncokb-gene-db/expected: the multianno and MAF the
# templates wrote before oncokbGene existed (commit $PRE below). REGEN_GOLDEN=1
# rewrites them from that commit's templates (needs the commit in the clone).
#
#  1. database installed: the five *.oncokbGene columns come last before Otherinfo,
#     every other column is the baseline's, and the isoform table and version are
#     copied next to the multianno;
#  2. no database, oncokb-gene-db=off, or links that don't name one version: the
#     multianno is the baseline, byte for byte;
#  3. an unknown oncokb-gene-db value, and a database switched while ANNOVAR runs,
#     fail the step;
#  4. annovar-to-maf: without the database the MAF is the baseline; with it, the
#     OncoKB columns are appended;
#  5. CancerVar (--skip_annovar reads the same multianno) gives the same output,
#     when CANCERVARDB names CancerVar's data folder (mise's task fetches it);
#  6. the production annotator list gets oncokbGene, "g" and "-hgvs" appended
#     (table_annovar.pl stubbed).
#
# KEEP_MULTIANNO=<path> keeps the with-database multianno (maftools-rs's fixture).
set -euo pipefail
cd "$(dirname "$0")/.."
PRE=345413bbdf36c16d0327b5a74e5f495b2ff85625
image="bergun/cancervar@sha256:ffad71ec62c2abcea5656102ed528b3a9db940af63b768bf32f7e5d9ba125754"
fix="$PWD/tests/fixtures/oncokb-gene-db"
golden="$fix/expected"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

render() { # yaml-file out-script annotators annotator-opts oncokb-gene-db
    python3 - "$@" <<'EOF'
import sys, yaml
src, out, annotators, opts, okb = sys.argv[1:6]
t = yaml.safe_load(open(src))["spec"]["templates"][0]
p = {x["name"]: x.get("value") for x in t["inputs"]["parameters"]}
s = t["script"]["source"]
vals = {"sample-id": "S1", "annotators": annotators or p["annotators"],
        "annotator-opts": opts or p["annotator-opts"], "oncokb-gene-db": okb}
for k, v in vals.items():
    s = s.replace("{{inputs.parameters.%s}}" % k, v)
s = s.replace("{{=sprig.osBase(inputs.parameters['input-vcf-f'])}}", "fixture.vcf")
assert "{{" not in s, s
open(out, "w").write(s)
EOF
}

mkdir -p "$work/in" && cp "$fix/fixture.vcf" "$work/in/"
install_like() { # dir version...: the layout install.sh leaves, links to the first version
    local d=$1; shift
    mkdir -p "$d"
    cp "$fix"/hg38_refGene* "$fix"/hg38_ensGene* "$d/"
    for v in "$@"; do
        cp "$fix/hg38_oncokbGene.txt" "$d/hg38_oncokbGene.$v.txt"
        cp "$fix/hg38_oncokbGeneMrna.fa" "$d/hg38_oncokbGeneMrna.$v.fa"
        cp "$fix/hg38_oncokbGene.isoforms.tsv" "$d/hg38_oncokbGene.$v.isoforms.tsv"
        cp "$fix/hg38_oncokbGene.version" "$d/hg38_oncokbGene.$v.version"
    done
    ln -s "hg38_oncokbGene.$1.txt" "$d/hg38_oncokbGene.txt"
    ln -s "hg38_oncokbGeneMrna.$1.fa" "$d/hg38_oncokbGeneMrna.fa"
    ln -s "hg38_oncokbGene.$1.isoforms.tsv" "$d/hg38_oncokbGene.isoforms.tsv"
    ln -s "hg38_oncokbGene.$1.version" "$d/hg38_oncokbGene.version"
}
mkdir -p "$work/db-without" && cp "$fix"/hg38_refGene* "$fix"/hg38_ensGene* "$work/db-without/"
install_like "$work/db-with" fixture
install_like "$work/db-mixed" v1 v2
ln -sfn hg38_oncokbGeneMrna.v2.fa "$work/db-mixed/hg38_oncokbGeneMrna.fa" # half-switched

run() { # name yaml db okb [stub-dir]: ANNOVAR's output in $work/out-<name>, exit code in $work/<name>.rc
    local name=$1 yaml=$2 db=$3 okb=$4 stub=${5:-}
    mkdir -p "$work/out-$name"
    if [ -n "$stub" ]; then
        render "$yaml" "$work/$name.sh" "" "" "$okb"   # the production list
    else
        render "$yaml" "$work/$name.sh" "refGene,ensGene" "'-hgvs','-hgvs'" "$okb"
    fi
    local mode=ro wd=/tools/annovar/bin extra=()
    if [ -n "$stub" ]; then mode=rw; wd=/stub; extra=(-v "$stub:/stub:ro"); fi
    set +e
    docker run --rm --network none --memory=2g -u "$(id -u):$(id -g)" \
        -v "$work/in:/mnt/input:ro" -v "$work/out-$name:/mnt/output" \
        -v "$work/$db:/mnt/annovardb:$mode" -v "$work/$name.sh:/run.sh:ro" "${extra[@]}" \
        -w "$wd" "$image" sh /run.sh > "$work/$name.log" 2>&1
    echo $? > "$work/$name.rc"
    set -e
}
rc() { cat "$work/$1.rc"; }

fail=0
check() { if eval "$2"; then echo "ok: $1"; else echo "FAIL: $1"; fail=1; fi; }

maf() { # out-name yaml
    local name=$1 yaml=$2 img
    img=$(python3 -c "import sys,yaml;print(yaml.safe_load(open(sys.argv[1]))['spec']['templates'][0]['script']['image'])" "$yaml")
    python3 -c "
import sys, yaml
s = yaml.safe_load(open(sys.argv[1]))['spec']['templates'][0]['script']['source']
s = s.replace('{{inputs.parameters.sample-id}}', 'S1').replace(\"{{=sprig.osBase(inputs.parameters['annovar-out-f'])}}\", 'S1.hg38_multianno.txt')
assert '{{' not in s, s
open(sys.argv[2], 'w').write(s)" "$yaml" "$work/maf-$name.sh"
    docker run --rm --network none --memory=2g -u "$(id -u):$(id -g)" -v "$work/out-$name:/mnt/input" \
        -v "$work/maf-$name.sh:/run.sh:ro" -w /mnt/input "$img" sh /run.sh > "$work/maf-$name.log" 2>&1 \
        || { cat "$work/maf-$name.log"; echo "FAIL: annovar-to-maf $name"; exit 1; }
}

if [ "${REGEN_GOLDEN:-}" = 1 ]; then
    git show "$PRE:tasks/base_tasks/annovar.yaml" > "$work/pre.yaml"
    git show "$PRE:tasks/base_tasks/annovar-to-maf.yaml" > "$work/pre-maf.yaml"
    run pre "$work/pre.yaml" db-with auto
    [ "$(rc pre)" = 0 ] || { cat "$work/pre.log"; exit 1; }
    maf pre "$work/pre-maf.yaml"
    mkdir -p "$golden"
    cp "$work/out-pre/S1.hg38_multianno.txt" "$work/out-pre/S1.maf" "$golden/"
    echo "rewrote $golden from $PRE"
fi

t=tasks/base_tasks/annovar.yaml
run with "$t" db-with auto
run without "$t" db-without auto
run off "$t" db-with off
run mixed "$t" db-mixed auto
run badvalue "$t" db-with on
for n in with without off mixed; do
    [ "$(rc $n)" = 0 ] || { cat "$work/$n.log"; echo "FAIL: $n exited $(rc $n)"; exit 1; }
done

m="$golden/S1.hg38_multianno.txt"
check "no database: multianno is the baseline" "cmp -s $m $work/out-without/S1.hg38_multianno.txt"
check "oncokb-gene-db=off: multianno is the baseline" "cmp -s $m $work/out-off/S1.hg38_multianno.txt"
check "links naming two versions: skipped, multianno is the baseline" \
    "cmp -s $m $work/out-mixed/S1.hg38_multianno.txt && grep -q 'not one complete version' $work/mixed.log && [ ! -e $work/out-mixed/S1.oncokbGene.isoforms.tsv ]"
check "no database: no isoform table" "[ ! -e $work/out-without/S1.oncokbGene.isoforms.tsv ]"
check "unknown oncokb-gene-db value fails the step" "[ \"\$(rc badvalue)\" = 2 ] && grep -q 'must be auto or off' $work/badvalue.log"
check "with database: isoform table and version copied" \
    "cmp -s $fix/hg38_oncokbGene.isoforms.tsv $work/out-with/S1.oncokbGene.isoforms.tsv && cmp -s $fix/hg38_oncokbGene.version $work/out-with/S1.oncokbGene.version"
[ -z "${KEEP_MULTIANNO:-}" ] || cp "$work/out-with/S1.hg38_multianno.txt" "$KEEP_MULTIANNO"
python3 - "$m" "$work/out-with/S1.hg38_multianno.txt" <<'EOF' || fail=1
import csv, sys
a = list(csv.reader(open(sys.argv[1]), delimiter="\t"))
b = list(csv.reader(open(sys.argv[2]), delimiter="\t"))
new = [f"{c}.oncokbGene" for c in ("Func", "Gene", "GeneDetail", "ExonicFunc", "AAChange")]
i = a[0].index("Otherinfo1")
assert b[0] == a[0][:i] + new + a[0][i:], b[0]
assert len(a) == len(b)
for ra, rb in zip(a, b):
    assert ra == rb[:i] + rb[i + 5:], (ra, rb)
okb = " ".join(r[i + 4] for r in b[1:])
mane = " ".join(r[a[0].index("AAChange.refGene")] for r in a[1:])
assert "ROS1:ENST00000368508.7:" in okb and "ENST00000368508.7" not in mane, okb
assert "BRAF:ENST00000646891.2:exon15:c.1799T>A:p.V600E" in okb, okb
print("ok: with database: five oncokbGene columns before Otherinfo1, all other columns the baseline's")
EOF

# a database switched while ANNOVAR runs: table_annovar.pl stubbed to switch one link
install_like "$work/db-switch" v1 v2
mkdir -p "$work/stub-switch"
printf '%s\n' 'symlink("hg38_oncokbGeneMrna.v2.fa", "/mnt/annovardb/.m") or die; rename("/mnt/annovardb/.m", "/mnt/annovardb/hg38_oncokbGeneMrna.fa") or die;' \
    > "$work/stub-switch/table_annovar.pl"
run switch "$t" db-switch auto "$work/stub-switch"
check "database switched during the run fails the step" "[ \"\$(rc switch)\" = 1 ] && grep -q 'changed during the run (was v1)' $work/switch.log"

maf without tasks/base_tasks/annovar-to-maf.yaml
maf with tasks/base_tasks/annovar-to-maf.yaml
check "no database: MAF is the baseline" "cmp -s $golden/S1.maf $work/out-without/S1.maf"
python3 - "$work/out-with/S1.maf" <<'EOF' || fail=1
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
assert list(rows[0])[-10:] == ["tx_mane", "exon_mane", "txChange_mane", "aaChange_mane", "oncokb_isoform",
    "oncokb_isoform_gene", "oncokb_isoform_aa", "oncokb_isoform_ref_ok", "oncokb_isoform_exon_edge", "oncokb_gene_db"]
ros1 = next(r for r in rows if r["Hugo_Symbol"] == "ROS1")
assert (ros1["aaChange"], ros1["aaChange_mane"], ros1["oncokb_isoform_aa"]) == ("p.L1000V", "p.L995V", "L1000V"), ros1
assert all(r["oncokb_gene_db"] == "fixture" for r in rows)
print("ok: with database: the MAF shows ROS1 on OncoKB's isoform, keeps MANE, appends the OncoKB columns")
EOF

# CancerVar. PYTHONHASHSEED fixes the order of its Therap/Diag/Prog lists (Python sets).
if [ -n "${CANCERVARDB:-}" ]; then
    mkdir -p "$work/cv-pre" "$work/cv-with"
    cp "$m" "$work/cv-pre/"
    cp "$work/out-with/S1.hg38_multianno.txt" "$work/cv-with/"
    for name in pre with; do
        docker run --rm --network none --memory=2g -e PYTHONHASHSEED=0 -u "$(id -u):$(id -g)" \
            -v "$work/cv-$name:/mnt/output" -v "$CANCERVARDB:/mnt/annovardb:ro" --entrypoint sh \
            -w /tools/annovar/bin "$image" -c 'touch /mnt/output/S1.vcf && python ./CancerVar.py -b hg38 -i /mnt/output/S1.vcf --input_type=VCF -o /mnt/output/S1 -c ./config_edit.ini --skip_annovar' \
            > "$work/cv-$name.log" 2>&1
    done
    check "CancerVar: same output with the oncokbGene columns" \
        "[ -s $work/cv-pre/S1.hg38_multianno.txt.cancervar ] && cmp -s $work/cv-pre/S1.hg38_multianno.txt.cancervar $work/cv-with/S1.hg38_multianno.txt.cancervar"
else
    echo "skip: CancerVar (set CANCERVARDB to CancerVar's cancervardb folder)"
fi

# the production annotator list, table_annovar.pl stubbed to print its arguments
mkdir -p "$work/stub-args"
printf '%s\n' 'print join("\n", @ARGV), "\n";' > "$work/stub-args/table_annovar.pl"
cp -a "$work/db-with" "$work/db-args"
run args "$t" db-args auto "$work/stub-args"
args=$(tr '\n' ' ' < "$work/args.log")
check "production list: protocol gets oncokbGene last" \
    "grep -q -- '--protocol refGene,ensGene,esp6500siv2_all,1000g2015aug_all,exac03,avsnp150,dbnsfp42c,dbscsnv11,dbnsfp31a_interpro,clinvar_20221231,cosmic91,icgc28,gnomad_genome,oncokbGene ' <<< \"\$args\""
check "production list: operation g,g,f x11,g" "grep -q -- '--operation g,g,f,f,f,f,f,f,f,f,f,f,f,g ' <<< \"\$args\""
check "production list: --arg gets -hgvs for oncokbGene" "grep -q -- '--arg -hgvs,-hgvs,,,,,,,,,,,,-hgvs ' <<< \"\$args\""
exit $fail
