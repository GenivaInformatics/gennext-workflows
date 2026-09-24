import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('builder', Path(__file__).with_name('cnv_cohort_builder.py'))
b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b)

class WgsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fa = self.root / 'genome.fa'
        self.fa.write_text('>chr1\nACGT\n')
        Path(str(self.fa)+'.fai').write_text('chr1\t100000\t0\t0\t0\nchrUn\t4000\t0\t0\t0\n')
        self.bed = self.root / 'kit.bed'
        self.bed.write_text('chr1\t0\t100000\n')
    def test_auto_requires_full_indexed_contigs(self):
        self.assertEqual(b.resolve_mode(self.bed,self.fa,'auto'),'wgs')
        self.bed.write_text('chr1\t1\t100000\n')
        self.assertEqual(b.resolve_mode(self.bed,self.fa,'auto'),'amplicon')
        self.bed.write_text('chr2\t0\t100000\n')
        self.assertNotEqual(b.resolve_mode(self.bed,self.fa,'auto'),'wgs')
    def test_explicit_panel_mode_is_preserved(self):
        self.assertEqual(b.resolve_mode(self.bed,self.fa,'hybrid'),'hybrid')
    def test_wgs_excludes_gaps_and_unselected_contigs(self):
        access=self.root/'access.bed'
        access.write_text('chr1\t0\t40000\nchr1\t45000\t100000\nchrUn\t0\t4000\n')
        with patch.object(b,'cnvkit') as run:
            b.build_target_beds(self.bed,self.root/'genes.txt',access,self.root,'wgs',self.root,10000)
        self.assertEqual((self.root/'wgs-access.bed').read_text(),'chr1\t0\t40000\nchr1\t45000\t100000\n')
        self.assertEqual((self.root/'antitargets.bed').read_text(),'')
        self.assertEqual(run.call_count,1)
        self.assertIn('--avg-size',run.call_args.args[0])
    def test_nonoverlapping_reference_fails(self):
        access=self.root/'access.bed'; access.write_text('chrUn\t0\t4000\n')
        with self.assertRaisesRegex(RuntimeError,'no overlap'):
            b.intersect_access(self.bed,access,self.root/'out.bed')
    def test_wgs_reference_uses_no_edge_and_only_target_coverage(self):
        cov=self.root/'coverage'; cov.mkdir()
        (cov/'normal.targetcoverage.cnn').touch()
        (cov/'stale.antitargetcoverage.cnn').touch()
        with patch.object(b,'cnvkit') as run:
            b.build_pooled_reference(self.root,self.fa,self.root,'wgs')
        argv=run.call_args.args[0]
        self.assertIn('--no-edge',argv)
        self.assertNotIn(str(cov/'stale.antitargetcoverage.cnn'),argv)
    def test_panel_reference_retains_existing_corrections(self):
        cov=self.root/'coverage'; cov.mkdir(); (cov/'normal.targetcoverage.cnn').touch()
        with patch.object(b,'cnvkit') as run:
            b.build_pooled_reference(self.root,self.fa,self.root,'hybrid')
        self.assertNotIn('--no-edge',run.call_args.args[0])

if __name__=='__main__': unittest.main()
