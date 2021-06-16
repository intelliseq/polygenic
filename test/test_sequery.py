from unittest import TestCase
from polygenic import polygenic
from polygenic import polygenicmaker

class PolygenicTest(TestCase):

    def test(self):
        self.assertEqual(1, 1)
    def testPolygenicCoreWithAf(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz", 
            "--log_file", "/dev/null",
            "--model", "test/resources/model/scaled_eas_model.py", 
            "--population", "eas", 
            "--out_dir", "/tmp/",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')

class PolygenicMakerTest(TestCase):
    def testBiobankukIndex(self):
        polygenicmaker.main([
            "biobankuk-index",
            "--output", "results"
        ])
        self.assertEqual('1', '1')

    def testBiobankukGet(self):
        polygenicmaker.main([
            "biobankuk-get",
            "--index", "results/phenotype_manifest.tsv",
            "--phenocode", "30600",
            "--output", "results"
        ])
        self.assertEqual('1', '1')

    def testBiobankukBuildModel(self):
        polygenicmaker.main([
            "biobankuk-build-model",
            "--data", "results/biomarkers-30600-both_sexes-irnt.tsv",
            "--output", "results/model",
            "--anno", "results/full_variant_qc_metrics.txt",
            "--threshold", "1e-08"
        ])
        self.assertEqual('1', '1')
