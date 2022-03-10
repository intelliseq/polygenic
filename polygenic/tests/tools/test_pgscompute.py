from unittest import TestCase
from polygenic.tools import pgscompute
from polygenic import pgstk

from pathlib import Path as path

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
            '--model', 'polygenic/tests/resources/model/test.model.yml', "polygenic/tests/resources/model/test.model.yml",
            '--output-name-appendix', appendix,
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/testsample-test.model.yml-" + appendix + "-result.json", 'r') as output:
            results = json.load(output)
            self.assertTrue("score_model" in results)
            self.assertTrue("description" in results)
            self.assertTrue("sample_name" in results["description"])
            self.assertTrue("model_name" in results["description"])
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

    # test if af is used
    def testPgsComputeMerge(self):
        appendix = "merge"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test.sample.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test.model.yml',
            '--output-name-appendix', appendix,
            '--merge-outputs',
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/testsample-test.model.yml-" + appendix + "-result.json", 'r') as output:
            results = json.load(output)
            self.assertTrue("test.model" in results)

    # test diplotype model categories
    def testPgsComputeDiplotype(self):
        appendix = "diplotype"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test.sample.vcf.gz',
            '--model', 'polygenic/tests/resources/model/diplotype_model.yml',
            '--output-name-appendix', appendix,
            '--af', 'polygenic/tests/resources/vcf/test.af.vcf.gz',
            '--merge-outputs',
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/testsample-diplotype_model.yml-" + appendix + "-result.json", 'r') as output:
            results = json.load(output)
            self.assertTrue("test.model" in results)

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

    # test haplotype
    def testPgsComputeHaplotype(self):
        appendix = "haplotype"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test.sample.vcf.gz',
            '--model', 'polygenic/tests/resources/model/bitter-taste.yml',
            '--output-name-appendix', appendix,
            '--af', 'polygenic/tests/resources/vcf/test.af.vcf.gz',
            '--af-field', 'AF_nfe',
            '--output-directory', self.output_directory])

    def testPgsComputeHaplotypePgx(self):
        appendix = "pgx"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test-pgs-compute-cyp-bug.vcf.gz',
            '--model', 'polygenic/models/pgx/cyp2d6-pharmvar.yml',
            '--output-name-appendix', appendix,
            '--af', '/tmp/marpiech/polygenic/gnomad.3.1.vcf.gz',
            '--af-field', 'AF_nfe',
            '--output-directory', self.output_directory,
            '--print'])

    # test if can be executed through pgstk
    def testPgsComputeWork(self):
        appendix = "work"
        pgstk.main([
            'pgs-compute',
            '--vcf', '/data/downloads/resources/prodia/samples/phased_203253100045_R01C01.vcf.gz',
            '--model', "/tmp/marpiech/prodia/nutrigx/alcohol.yml", "/tmp/marpiech/prodia/nutrigx/lactose.yml", "/tmp/marpiech/prodia/nutrigx/antioxidants.yml", "/tmp/marpiech/prodia/nutrigx/liver_detoxification.yml", "/tmp/marpiech/prodia/nutrigx/benefit_of_exercise.yml", "/tmp/marpiech/prodia/nutrigx/low_iron_status.yml", "/tmp/marpiech/prodia/nutrigx/benefit_of_exercise_in_blood_pressure.yml", "/tmp/marpiech/prodia/nutrigx/monounsaturated_fat.yml", "/tmp/marpiech/prodia/nutrigx/benefit_of_exercise_in_bmi.yml", "/tmp/marpiech/prodia/nutrigx/omega_profile.yml", "/tmp/marpiech/prodia/nutrigx/benefit_of_exercise_in_hdl_cholesterol.yml", "/tmp/marpiech/prodia/nutrigx/power.yml", "/tmp/marpiech/prodia/nutrigx/benefit_of_exercise_in_insulin_sensitivity.yml", "/tmp/marpiech/prodia/nutrigx/protein.yml", "/tmp/marpiech/prodia/nutrigx/bitter_taste.yml", "/tmp/marpiech/prodia/nutrigx/saturated_fat.yml", "/tmp/marpiech/prodia/nutrigx/caffeine.yml", "/tmp/marpiech/prodia/nutrigx/sodium.yml", "/tmp/marpiech/prodia/nutrigx/calcium.yml", "/tmp/marpiech/prodia/nutrigx/vitamin_a.yml", "/tmp/marpiech/prodia/nutrigx/choline.yml", "/tmp/marpiech/prodia/nutrigx/vitamin_b12.yml", "/tmp/marpiech/prodia/nutrigx/endurance.yml", "/tmp/marpiech/prodia/nutrigx/vitamin_b6.yml", "/tmp/marpiech/prodia/nutrigx/energy_balance.yml", "/tmp/marpiech/prodia/nutrigx/vitamin_c.yml", "/tmp/marpiech/prodia/nutrigx/folate.yml", "/tmp/marpiech/prodia/nutrigx/vitamin_d.yml", "/tmp/marpiech/prodia/nutrigx/gluten.yml", "/tmp/marpiech/prodia/nutrigx/vitamin_e.yml", "/tmp/marpiech/prodia/nutrigx/individual_eating_motivation.yml", "/tmp/marpiech/prodia/nutrigx/weight_balance.yml", "/tmp/marpiech/prodia/nutrigx/iron_overload.yml", "/tmp/marpiech/prodia/nutrigx/whole_grain.yml",
            '--output-name-appendix', appendix,
            '--merge-outputs',
            '--output-directory', self.output_directory])

    # test bug in haplotype model with missing variant
    def testPgcComputeCypBug(self):
        appendix = "cypbug"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test-pgs-compute-cyp-bug.vcf.gz',
            '--model', 'models/pgx/cyp2d6-pharmvar.yml',
            '--output-name-appendix', appendix,
            '--merge-outputs',
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/2824-" + appendix + "-result.json", 'r') as output:
            results = json.load(output)
            self.assertTrue("cyp2d6-pharmvar" in results)
            self.assertTrue("haplotype_model" in results["cyp2d6-pharmvar"])

    def testCyp2c19(self):
        appendix = "cyp2c19"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/2824.vcf.gz',
            '--model', 'models/pgx/cyp2c19-pharmvar-5.1.8.yml',
            '--output-name-appendix', appendix,
            '--af', 'polygenic/tests/resources/vcf/test.af.vcf.gz',
            '--af-field', 'AF_nfe',
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/2824-" + appendix + "-result.json", 'r') as output:
            results = json.load(output)
            self.assertTrue("cyp2c19-pharmvar" in results)
            self.assertTrue("haplotype_model" in results["cyp2c19-pharmvar"])

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