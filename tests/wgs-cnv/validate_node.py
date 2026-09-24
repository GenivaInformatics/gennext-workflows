"""Read-only validation of real hg38 references and healthy-cohort D4 on a node."""
from importlib.machinery import SourceFileLoader
from pathlib import Path
import os
b=SourceFileLoader('builder',os.environ.get('BUILDER_PATH','/code/cnv_cohort_builder.py')).load_module()
base=Path('/reference/hg38')
fa=base/'genome/GCA_000001405.15_GRCh38_full_minus_alt_plus_hs38d1_plus_HLA_analysis_set/GCA_000001405.15_GRCh38_full_minus_alt_plus_hs38d1_plus_HLA_analysis_set.fa'
bed=base/'regions/hg38_canonical_chrom_contig_regions.bed'
ref=base/'regions/refFlat_hg38.txt'
assert b.resolve_mode(bed,fa,'auto')=='wgs'
from argparse import Namespace
out=Path('/tmp/real-wgs-reference')
b.cmd_prepare_bins(Namespace(kit_bed=str(bed),ref_fa=str(fa),ref_flat=str(ref),output_dir=str(out),mode='auto',wgs_bin_size=10000))
bins=b.load_bins(out/'targets.bed',None)
assert 200000<len(bins)<400000,len(bins)
allowed={x[0] for x in b.load_bins(bed,None)}
assert all(c in allowed for c,s,e,g in bins)
assert not (out/'antitargets.bed').stat().st_size
# Confirm production D4 contig naming and interval queries match generated bins.
cohort=base/'other/cohorts/cnv_cohort.d4'
tracks=b.d4_list_tracks(cohort)
subset=out/'probe.bed'
probe=bins[:10]+bins[len(bins)//2:len(bins)//2+10]
subset.write_text(''.join(f'{c}\t{s}\t{e}\n' for c,s,e,g in probe))
means=b.d4_region_means(cohort,tracks[0],subset)
assert len(means)==len(probe),(len(means),len(probe))
assert any(v>0 for v in means.values())
print(f'REAL_REFERENCE_OK bins={len(bins)} selected_contigs={len(set(x[0] for x in bins))} cohort_tracks={len(tracks)} d4_probe_intervals={len(means)}',flush=True)
