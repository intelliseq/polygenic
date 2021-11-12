from unittest import TestCase
from polygenic.tools import pgscompute

from pathlib import Path as path

import tabix
import os
import configparser

class PgsComputeTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(PgsComputeTest, self).__init__(*args, **kwargs)
        self.output_directory = "/tmp/polygenic/test"
        path(self.output_directory).mkdir(parents=True, exist_ok=True)

    def testPgsCompute(self):
        pgscompute.main([
            '--vcf', 'polygenic/tests/resources/vcf/test.sample.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test.model.yml',
            '--output-directory', self.output_directory])

        # with open(result_path, 'r') as output:
        #     data = output.read()
        #     header = list(filter(lambda line: "score_model:" in line, data.split('\n')))
        #     self.assertEqual(1, len(header))
        #     print(str(data))

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