"""Exercise real CNVkit commands against a synthetic, known-gain WGS fixture."""
import importlib.util
import math
import os
from pathlib import Path
import subprocess
import tempfile

spec=importlib.util.spec_from_file_location('builder',os.environ.get('BUILDER_PATH','/usr/local/bin/cnv-cohort-builder'))
# Installed executable has no .py suffix.
from importlib.machinery import SourceFileLoader
b=SourceFileLoader('builder',os.environ.get('BUILDER_PATH','/usr/local/bin/cnv-cohort-builder')).load_module()
root=Path(tempfile.mkdtemp(prefix='wgs-smoke-'))
fa=root/'genome.fa'
# Include a long inaccessible gap and an unselected alternate contig.
fa.write_text('>chr1\n'+('ACGT'*100000)+'N'*20000+('ACGT'*145000)+'\n>chrUn\n'+'ACGT'*2500+'\n')
import pysam
pysam.faidx(str(fa))
bed=root/'genome.bed'; bed.write_text('chr1\t0\t1000000\n')
ref=root/'refFlat.txt'; ref.write_text('TEST\tNM_000001\tchr1\t+\t0\t1000000\t0\t1000000\t1\t0,\t1000000,\n')
kit=root/'kits'/'whole'; kit.mkdir(parents=True)
logs=kit/'logs'
assert b.resolve_mode(bed,fa,'auto')=='wgs'
access=b.ensure_access_bed(fa,kit,logs)
targets,anti=b.build_target_beds(bed,ref,access,kit,'wgs',logs,10000)
bins=b.load_bins(targets,None)
assert 90 <= len(bins) <= 110, len(bins)
assert all(c=='chr1' and (e<=400000 or s>=420000) for c,s,e,g in bins)
assert anti.stat().st_size==0
# Exercise D4 creation, multi-track extraction, reference build and cache reuse.
from argparse import Namespace
genome=root/'chrom.sizes'; genome.write_text('chr1\t1000000\nchrUn\t10000\n')
tracks=[]
for n in range(3):
    depth=root/f'normal{n}.bedgraph'
    depth.write_text(f'chr1\t0\t1000000\t{29+n}\nchrUn\t0\t10000\t{29+n}\n')
    d4=root/f'normal{n}.d4'
    subprocess.run(['d4tools','create','-g',str(genome),str(depth),str(d4)],check=True)
    tracks.append(str(d4))
cohort=root/'cohort.d4'
subprocess.run(['d4tools','merge',*tracks,str(cohort)],check=True,capture_output=True)
args=Namespace(kit_bed=str(bed),cohort_d4=str(cohort),ref_fa=str(fa),ref_flat=str(ref),
               kits_dir=str(kit.parent),kit_name='whole',mode='auto',workers=2,
               skip_existing=False,wgs_bin_size=10000)
assert b.cmd_ensure_kit(args)==0
assert b.cmd_ensure_kit(args)==0
reference=kit/'pooled_reference.cnn'
args.wgs_bin_size=20000
try:
    b.cmd_ensure_kit(args)
    raise AssertionError('incompatible cache was accepted')
except RuntimeError as error:
    assert 'Cached WGS reference' in str(error)
# First quarter is copy-gained; normals remain diploid throughout.
tumor=root/'tumor.targetcoverage.cnn'
b.write_cnn(tumor,bins,[60 if start<250000 else 30 for chrom,start,end,gene in bins])
empty=root/'empty.cnn'; b.write_cnn(empty,[],[])
cnr=root/'tumor.cnr'
b.cnvkit(['fix',str(tumor),str(empty),str(reference),'--no-edge','-o',str(cnr)])
import csv
rows=list(csv.DictReader(cnr.open(),delimiter='\t'))
gain=[float(r['log2']) for r in rows if int(r['start'])<250000]
base=[float(r['log2']) for r in rows if int(r['start'])>=500000]
assert gain and base
contrast=sum(gain)/len(gain)-sum(base)/len(base)
assert 0.85<contrast<1.15, contrast
b.cnvkit(['segment',str(cnr),'--method','cbs','-o',str(root/'tumor.cns')])
b.cnvkit(['call',str(root/'tumor.cns'),'-o',str(root/'tumor.call.cns')])
print(f'WGS_SMOKE_OK bins={len(bins)} gain_log2_contrast={contrast:.4f}',flush=True)
