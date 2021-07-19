from unittest import TestCase
from polygenic import polygenic
from polygenic import polygenicmaker
from polygenic import vcfstat

import tabix
import os
import configparser

class PolygenicTest(TestCase):

    def test(self):
        self.assertEqual(1, 1)
    def testPolygenicCoreWithAf(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz",
            "--sample", "yfyfy", 
            "--log-file", "/dev/null",
            "--model", "test/resources/model/scaled_eas_model.py", 
            "--population", "eas", 
            "--output-directory", "/tmp/polygenic",
            "--output-file-name", "bambala.json",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')
    
    def testPolygenicForBiobankModel(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz", 
            "--log_file", "/dev/null",
            "--model", "test/resources/model/biomarkers-30600-both_sexes-irnt.tsv.py", 
            "--population", "eas", 
            "--out_dir", "/tmp/",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')

    def testPolygenicForGbeModel(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz", 
            "--log_file", "/dev/null",
            "--model", "results/model/BIN1210.py", 
            "--population", "nfe", 
            "--out_dir", "/tmp/",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')

class PolygenicMakerTest(TestCase):

    def testTabix(self):
        config = configparser.ConfigParser()
        config.read(os.path.dirname(__file__) + "/../polygenic/polygenic.cfg")
        url = config['urls']['hg19-rsids']
        tb = tabix.open(url)
        print(type(tb).__name__)
        records = tb.query("16", 1650945, 1650946)
        for record in records:
            print(record[:3])

    def testGbeIndex(self):
        polygenicmaker.main([
            "gbe-index",
            "--output", "results"
        ])

    def testGbeGet(self):
        polygenicmaker.main([
            "gbe-get",
            "--code", "BIN1210",
            "--output", "results"
        ])

    def testGbePrepareModel(self):
        polygenicmaker.main([
            "gbe-prepare-model",
            "--data", "results/BIN1210",
            "--af", "/home/marpiech/data/af.vcf.gz",
            "--output", "results/model"
        ])

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

    def testPgsIndex(self):
        polygenicmaker.main([
            "pgs-index",
            "--output", "/tmp/polygenic"
        ])

    def testPgsGet(self):
        polygenicmaker.main([
            "pgs-get",
            "--code", "BIN1210",
            "--output", "results"
        ])

class VcfstatTest(TestCase):

    def testBaf(self):
        vcfstat.main([
            "baf",
            "--vcf", "/home/marpiech/data/clustered_204800980122_R01C02.vcf.gz",
            "--output-directory", "/tmp/baf"
        ])
        self.assertEqual('1', '1')

    def testZygosity(self):
        vcfstat.main([
            "zygosity",
            "--vcf", "/home/marpiech/data/clustered_204800980122_R01C02.vcf.gz",
            "--output-file", "/tmp/baf/stats.json"
        ])
        self.assertEqual('1', '1')

class Debug(TestCase):

    def testDebug(self):
        polygenic.main([
            "--vcf", "/tmp/clustered_204800980122_R01C02.vcf.gz",
            "--sample-name", "clustered_204800980122_R01C02",
            "--log-file", "/dev/null",
            "--output-name-appendix", "cancer",
            "--model", "/tmp/py/cancer1002.py",
            "--model", "/tmp/py/cancer1044.py", 
            "--population", "eas", 
            "--output-directory", "/tmp/polygenic/",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')