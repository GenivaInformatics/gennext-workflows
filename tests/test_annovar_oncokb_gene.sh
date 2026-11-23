#!/usr/bin/env bash
# annovar-template's optional oncokbGene protocol, run in the production ANNOVAR image
# on the synthetic fixture (tests/fixtures/oncokb-gene-db: MANE refGene/ensGene and
# OncoKB-isoform oncokbGene for 13 genes, 17 made-up variants).
#
#  1. database present: *.oncokbGene columns come last before Otherinfo, the isoform
#     table and version are copied next to the multianno, and every other column is
#     byte-identical to the template on main;
#  2. database absent, and oncokb-gene-db=off: the multianno is byte-identical to main's;
#  3. the production annotator list: oncokbGene is appended to --protocol, "g" to
#     --operation and "-hgvs" to --arg (table_annovar.pl stubbed);
#  4. annovar-to-maf on those outputs: without the database the MAF is
#     byte-identical to main's; with it, the OncoKB columns are appended.
#
# KEEP_MULTIANNO=<path> keeps the with-database multianno (maftools-rs's fixture).
set -euo pipefail
cd "$(dirname "$0")/.."
image="bergun/cancervar@sha256:ffad71ec62c2abcea5656102ed528b3a9db940af63b768bf32f7e5d9ba125754"
fix="$PWD/tests/fixtures/oncokb-gene-db"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

render() { # yaml-file out-script annotators annotator-opts oncokb-gene-db
    python3 - "$@" <<'EOF'
import sys, yaml
src, out, annotators, opts, okb = sys.argv[1:6]
t = yaml.safe_load(open(src))["spec"]["templates"][0]
s = t["script"]["source"]
for k, v in {"sample-id": "S1", "annotators": annotators, "annotator-opts": opts,
             "oncokb-gene-db": okb}.items():
    s = s.replace("{{inputs.parameters.%s}}" % k, v)
s = s.replace("{{=sprig.osBase(inputs.parameters['input-vcf-f'])}}", "fixture.vcf")
assert "{{" not in s, s
open(out, "w").write(s)
EOF
}

git show origin/main:tasks/base_tasks/annovar.yaml > "$work/main.yaml"
mkdir -p "$work/in" && cp "$fix/fixture.vcf" "$work/in/"
mkdir -p "$work/db-with" "$work/db-without"
cp "$fix"/hg38_refGene* "$fix"/hg38_ensGene* "$work/db-without/"
cp "$fix"/hg38_refGene* "$fix"/hg38_ensGene* "$work/db-with/"
# the layout install.sh leaves: versioned files and relative symlinks
for f in hg38_oncokbGene.txt:hg38_oncokbGene.fixture.txt hg38_oncokbGeneMrna.fa:hg38_oncokbGeneMrna.fixture.fa \
         hg38_oncokbGene.isoforms.tsv:hg38_oncokbGene.fixture.isoforms.tsv hg38_oncokbGene.version:hg38_oncokbGene.fixture.version; do
    cp "$fix/${f%%:*}" "$work/db-with/${f#*:}"
    ln -s "${f#*:}" "$work/db-with/${f%%:*}"
done

run() { # name yaml db okb
    local name=$1 yaml=$2 db=$3 okb=$4
    mkdir -p "$work/out-$name"
    render "$yaml" "$work/$name.sh" "refGene,ensGene" "'-hgvs','-hgvs'" "$okb"
    docker run --rm --network none --memory=2g -u "$(id -u):$(id -g)" \
        -v "$work/in:/mnt/input:ro" -v "$work/out-$name:/mnt/output" -v "$work/$db:/mnt/annovardb:ro" \
        -v "$work/$name.sh:/run.sh:ro" -w /tools/annovar/bin "$image" sh /run.sh > "$work/$name.log" 2>&1 \
        || { cat "$work/$name.log"; echo "FAIL: $name exited non-zero"; exit 1; }
}

run main "$work/main.yaml" db-with auto
run with tasks/base_tasks/annovar.yaml db-with auto
run without tasks/base_tasks/annovar.yaml db-without auto
run off tasks/base_tasks/annovar.yaml db-with off

fail=0
check() { if eval "$2"; then echo "ok: $1"; else echo "FAIL: $1"; fail=1; fi; }
m="$work/out-main/S1.hg38_multianno.txt"
check "no database: multianno identical to main" "cmp -s $m $work/out-without/S1.hg38_multianno.txt"
check "oncokb-gene-db=off: multianno identical to main" "cmp -s $m $work/out-off/S1.hg38_multianno.txt"
check "no database: no isoform table" "[ ! -e $work/out-without/S1.oncokbGene.isoforms.tsv ]"
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
print("ok: with database: five oncokbGene columns before Otherinfo1, all other columns unchanged")
EOF

# 4. annovar-to-maf on those outputs: main's template (maftools-rs 0.1.0) vs this one
maf() { # out-name yaml
    local name=$1 yaml=$2 image
    image=$(python3 -c "import sys,yaml;print(yaml.safe_load(open(sys.argv[1]))['spec']['templates'][0]['script']['image'])" "$yaml")
    python3 -c "
import sys, yaml
s = yaml.safe_load(open(sys.argv[1]))['spec']['templates'][0]['script']['source']
s = s.replace('{{inputs.parameters.sample-id}}', 'S1').replace(\"{{=sprig.osBase(inputs.parameters['annovar-out-f'])}}\", 'S1.hg38_multianno.txt')
assert '{{' not in s, s
open(sys.argv[2], 'w').write(s)" "$yaml" "$work/maf-$name.sh"
    docker run --rm --network none --memory=2g -u "$(id -u):$(id -g)" -v "$work/out-$name:/mnt/input" \
        -v "$work/maf-$name.sh:/run.sh:ro" -w /mnt/input "$image" sh /run.sh > "$work/maf-$name.log" 2>&1 \
        || { cat "$work/maf-$name.log"; echo "FAIL: annovar-to-maf $name"; exit 1; }
}
git show origin/main:tasks/base_tasks/annovar-to-maf.yaml > "$work/main-maf.yaml"
maf main "$work/main-maf.yaml"
maf without tasks/base_tasks/annovar-to-maf.yaml
maf with tasks/base_tasks/annovar-to-maf.yaml
check "no database: MAF identical to main's" "cmp -s $work/out-main/S1.maf $work/out-without/S1.maf"
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

# 5. CancerVar (--skip_annovar reads the same multianno): identical output with the
#    oncokbGene columns. Needs CancerVar's data files (the cancervardb folder of
#    github.com/WGLab/CancerVar) in CANCERVARDB; skipped without. PYTHONHASHSEED
#    fixes the order of its Therap/Diag/Prog lists, which come from Python sets.
if [ -n "${CANCERVARDB:-}" ]; then
    for name in main with; do
        mkdir -p "$work/cv-$name" && cp "$work/out-$name/S1.hg38_multianno.txt" "$work/cv-$name/"
        docker run --rm --network none --memory=2g -e PYTHONHASHSEED=0 -u "$(id -u):$(id -g)" \
            -v "$work/cv-$name:/mnt/output" -v "$CANCERVARDB:/mnt/annovardb:ro" --entrypoint sh \
            -w /tools/annovar/bin "$image" -c 'touch /mnt/output/S1.vcf && python ./CancerVar.py -b hg38 -i /mnt/output/S1.vcf --input_type=VCF -o /mnt/output/S1 -c ./config_edit.ini --skip_annovar' \
            > "$work/cv-$name.log" 2>&1
    done
    check "CancerVar: identical output with the oncokbGene columns" \
        "[ -s $work/cv-main/S1.hg38_multianno.txt.cancervar ] && cmp -s $work/cv-main/S1.hg38_multianno.txt.cancervar $work/cv-with/S1.hg38_multianno.txt.cancervar"
else
    echo "skip: CancerVar (set CANCERVARDB to CancerVar's cancervardb folder)"
fi

# 3. production annotators, table_annovar.pl stubbed to print its arguments
mkdir -p "$work/stub" "$work/out-stub"
printf '%s\n' 'print join("\n", @ARGV), "\n";' > "$work/stub/table_annovar.pl"
python3 - tasks/base_tasks/annovar.yaml "$work/stub.sh" <<'EOF'
import sys, yaml
t = yaml.safe_load(open(sys.argv[1]))["spec"]["templates"][0]
p = {x["name"]: x.get("value") for x in t["inputs"]["parameters"]}
s = t["script"]["source"]
for k in ("annotators", "annotator-opts", "oncokb-gene-db"):
    s = s.replace("{{inputs.parameters.%s}}" % k, p[k])
s = s.replace("{{inputs.parameters.sample-id}}", "S1").replace("{{=sprig.osBase(inputs.parameters['input-vcf-f'])}}", "fixture.vcf")
open(sys.argv[2], "w").write(s)
EOF
docker run --rm --network none --memory=512m -u "$(id -u):$(id -g)" -v "$work/out-stub:/mnt/output" \
    -v "$work/db-with:/mnt/annovardb:ro" -v "$work/stub:/stub:ro" -v "$work/stub.sh:/run.sh:ro" -w /stub \
    "$image" sh /run.sh > "$work/stub.log"
args=$(tr '\n' ' ' < "$work/stub.log")
check "production list: protocol gets oncokbGene last" \
    "grep -q -- '--protocol refGene,ensGene,esp6500siv2_all,1000g2015aug_all,exac03,avsnp150,dbnsfp42c,dbscsnv11,dbnsfp31a_interpro,clinvar_20221231,cosmic91,icgc28,gnomad_genome,oncokbGene ' <<< \"\$args\""
check "production list: operation g,g,f x11,g" "grep -q -- '--operation g,g,f,f,f,f,f,f,f,f,f,f,f,g ' <<< \"\$args\""
check "production list: --arg gets -hgvs for oncokbGene" "grep -q -- '--arg -hgvs,-hgvs,,,,,,,,,,,,-hgvs ' <<< \"\$args\""
exit $fail
