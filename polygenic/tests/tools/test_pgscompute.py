from unittest import TestCase
from polygenic.tools import pgscompute
from polygenic import pgstk

from pathlib import Path as path

import tabix
import os
import json

class PgsComputeTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(PgsComputeTest, self).__init__(*args, **kwargs)
        self.output_directory = "/tmp/polygenic/test"
        path(self.output_directory).mkdir(parents=True, exist_ok=True)

    # tests if can be executed directly
    def testPgsComputeDirectly(self):
        appendix = "directly"
        pgscompute.main([
            '--vcf', 'polygenic/tests/resources/vcf/test.sample.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test.model.yml',
            '--output-name-appendix', appendix,
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/testsample-test.model.yml-" + appendix + "-result.json", 'r') as output:
            results = json.load(output)
            self.assertTrue("score_model" in results)
            self.assertTrue("description" in results)
            self.assertTrue("randomentry" not in results)

    # test if can be executed through pgstk
    def testPgsComputePgstk(self):
        appendix = "pgstk"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test.sample.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test.model.yml',
            '--output-name-appendix', appendix,
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/testsample-test.model.yml-" + appendix + "-result.json", 'r') as output:
            results = json.load(output)
            self.assertTrue("score_model" in results)
            self.assertTrue("description" in results)
            self.assertTrue("randomentry" not in results)

    # test if af is used
    def testPgsComputeAf(self):
        appendix = "af"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test.sample.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test.model.yml',
            '--output-name-appendix', appendix,
            '--af', 'polygenic/tests/resources/vcf/test.af.vcf.gz',
            '--af-field', 'AF_nfe',
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/testsample-test.model.yml-" + appendix + "-result.json", 'r') as output:
            results = json.load(output)
            self.assertTrue("score_model" in results)
            self.assertEquals("af", results["genotypes"]["rs1570391830"]["genotype"]["source"])
            self.assertTrue("description" in results)
            self.assertTrue("randomentry" not in results)

    # test error
    def testPgsComputeError(self):
        appendix = "error"
        with self.assertRaises(SystemExit) as cm:
            pgstk.main([
                'pgs-compute',
                '--vcf', 'polygenic/tests/resources/vcf/test.sample.vcf.gz',
                '--model', 'polygenic/tests/resources/model/test.model.yml',
                '--output-name-appendix', appendix,
                '--af', 'polygenic/tests/resources/vcf/test.af.vcf.gz',
                '--af-field', 'bobodo',
                '--output-directory', self.output_directory])

        self.assertEqual(cm.exception.code, 1)

    # def testPolygenicMat(self):
    #     polygenic.main([
    #         "--vcf", "/tmp/marpiech/kenobi/polygenic/illu_merged-imputed.vcf.gz",
    #         "--model", "/tmp/marpiech/kenobi/polygenic/described-model.yml",
    #         "--output-directory", "/tmp/polygenic"])
    #     self.assertEqual('1', '1')

    # def testPolygenicGc(self):
    #     polygenic.main([
    #         "--vcf", "test/resources/vcf/my.vcf.gz",
    #         "--sample", "yfyfy", 
    #         "--log-file", "/dev/null",
    #         #"--model", "test/resources/model/variantset.yml",  
    #         #"--model", "test/resources/model/gc_prodia.yml",
    #         "--model", "test/resources/model/breast_prodia.yml",
    #         #"--model", "test/resources/model/polygenicscore.yml",  
    #         "--output-directory", "/tmp/polygenic",
    #         "--output-name-appendix", "yml",
    #         "--af", "test/resources/vcf/af.vcf.gz"])
    #     self.assertEqual('1', '1')

    # def testPolygenicCoreWithAf(self):
    #     polygenic.main([
    #         "--vcf", "test/resources/vcf/my.vcf.gz",
    #         "--sample", "yfyfy", 
    #         "--log-file", "/dev/null",
    #         "--model", "test/resources/model/scaled_eas_model.py", 
    #         "--population", "eas", 
    #         "--output-directory", "/tmp/polygenic",
    #         "--output-name-appendix", "bambala",
    #         "--af", "test/resources/vcf/af.vcf.gz"])
    #     self.assertEqual('1', '1')
    
    # def testPolygenicForBiobankModel(self):    def testPolygenicMat(self):
    #     polygenic.main([
    #         "--vcf", "/tmp/marpiech/kenobi/polygenic/illu_merged-imputed.vcf.gz",
    #         "--model", "/tmp/marpiech/kenobi/polygenic/described-model.yml",
    #         "--output-directory", "/tmp/polygenic"])
    #     self.assertEqual('1', '1')

    # def testPolygenicGc(self):
    #     polygenic.main([
    #         "--vcf", "test/resources/vcf/my.vcf.gz",
    #         "--sample", "yfyfy", 
    #         "--log-file", "/dev/null",
    #         #"--model", "test/resources/model/variantset.yml",  
    #         #"--model", "test/resources/model/gc_prodia.yml",
    #         "--model", "test/resources/model/breast_prodia.yml",
    #         #"--model", "test/resources/model/polygenicscore.yml",  
    #         "--output-directory", "/tmp/polygenic",
    #         "--output-name-appendix", "yml",
    #         "--af", "test/resources/vcf/af.vcf.gz"])
    #     self.assertEqual('1', '1')

    # def testPolygenicCoreWithAf(self):
    #     polygenic.main([
    #         "--vcf", "test/resources/vcf/my.vcf.gz",
    #         "--sample", "yfyfy", 
    #         "--log-file", "/dev/null",
    #         "--model", "test/resources/model/scaled_eas_model.py", 
    #         "--population", "eas", 
    #     polygenic.main([
    #         "--vcf", "test/resources/vcf/my.vcf.gz", 
    #         "--log-file", "/dev/null",
    #         "--model", "test/resources/model/biomarkers-30600-both_sexes-irnt.tsv.py", 
    #         "--population", "eas", 
    #         "--output-directory", "/tmp/",
    #         "--af", "test/resources/vcf/af.vcf.gz"])
    #     self.assertEqual('1', '1')

    # def testPolygenicForGbeModel(self):
    #     polygenic.main([
    #         "--vcf", "test/resources/vcf/my.vcf.gz", 
    #         "--log-file", "/dev/null",
    #         "--model", "results/model/BIN1210.py", 
    #         "--population", "nfe", 
    #         "--output-directory", "/tmp/",
    #         "--af", "test/resources/vcf/af.vcf.gz"])
    #     self.assertEqual('1', '1')

    # def testPolygenicParameters(self):
    #     polygenic.main([
    #         "--vcf", "test/resources/vcf/my.vcf.gz", 
    #         "--model", "test/resources/model/test-params.yml",
    #         "--parameters", "test/resources/json/test-params.json",
    #         "--output-directory", "/tmp/polygenic/",
    #         "--log-file", "/dev/null"])

    # def testPolygenicGeneSymbol(self):
    #     polygenic.main([
    #         "--vcf", "test/resources/vcf/my.vcf.gz", 
    #         "--model", "test/resources/model/test-gene-symbol.yml",
    #         "--output-directory", "/tmp/polygenic/",
    #         "--log-file", "/dev/null"])