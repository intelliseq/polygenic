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

    def testPolygenicYml(self):
        polygenic.main([
            #'--vcf', '/home/marpiech/data/gtcgenome/clusterd.id.vcf.gz',
            '--vcf', '/home/marpiech/data/mygenome-bgi/bgi-genotype.nochr.vcf.gz',
            '--log-file', '/dev/null',
            #'--model', 'test/resources/model/gc_prodia.yml',  
            '--model', 'test/resources/model/hirisplex.yml',
            #'--model', 'test/resources/model/breast_prodia.yml',
            '--output-directory', '/tmp/polygenic',
            '--output-name-appendix', 'yml'])
            #"--vcf", "test/resources/vcf/my.vcf.gz",
            #"--sample", "yfyfy", 
            #"--log-file", "/dev/null",
            #"--model", "test/resources/model/variantset.yml",  
            #"--model", "test/resources/model/brest_prodia.yml",  
            #"--model", "test/resources/model/polygenicscore.yml",  
            #"--output-directory", "/tmp/polygenic",
            #"--output-name-appendix", "yml",
            #"--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')

    def testPolygenicGc(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz",
            "--sample", "yfyfy", 
            "--log-file", "/dev/null",
            #"--model", "test/resources/model/variantset.yml",  
            #"--model", "test/resources/model/gc_prodia.yml",
            "--model", "test/resources/model/breast_prodia.yml",
            #"--model", "test/resources/model/polygenicscore.yml",  
            "--output-directory", "/tmp/polygenic",
            "--output-name-appendix", "yml",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')

    def testPolygenicCoreWithAf(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz",
            "--sample", "yfyfy", 
            "--log-file", "/dev/null",
            "--model", "test/resources/model/scaled_eas_model.py", 
            "--population", "eas", 
            "--output-directory", "/tmp/polygenic",
            "--output-name-appendix", "bambala",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')
    
    def testPolygenicForBiobankModel(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz", 
            "--log-file", "/dev/null",
            "--model", "test/resources/model/biomarkers-30600-both_sexes-irnt.tsv.py", 
            "--population", "eas", 
            "--output-directory", "/tmp/",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')

    def testPolygenicForGbeModel(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz", 
            "--log-file", "/dev/null",
            "--model", "results/model/BIN1210.py", 
            "--population", "nfe", 
            "--output-directory", "/tmp/",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')

    def testPolygenicParameters(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz", 
            "--model", "test/resources/model/test-params.yml",
            "--parameters", "test/resources/json/test-params.json",
            "--output-directory", "/tmp/polygenic/",
            "--log-file", "/dev/null"])

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
            "--output", "/tmp/polygenic/results"
        ])

    def testGbeGet(self):
        polygenicmaker.main([
            "gbe-get",
            "--code", "HC710",
            "--output", "/tmp/polygenic/results"
        ])

    def testGbePrepare(self):
        polygenicmaker.main([
            "gbe-prepare",
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
            "--code", "PGS000004",
            "--output-path", "/tmp/polygenic/PGS000004.txt"
        ])

    def testPgsPrepare(self):
        polygenicmaker.main([
            "pgs-prepare",
            "--input", "/tmp/polygenic/PGS000004.txt",
            "--output-path", "/tmp/polygenic/PGS000004.py",
            "--af", "/home/marpiech/data/af.vcf.gz",
            "--origin-reference-vcf", "/tmp/dbsnp/grch37/00-common_all.vcf.gz",
            "--model-reference-vcf", "/tmp/dbsnp/grch38/00-common_all.vcf.gz"
        ])

    def testVcfIndex(self):
        polygenicmaker.main([
            "vcf-index",
            "--vcf", "/tmp/dbsnp/grch38/00-common_all.vcf.gz"
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