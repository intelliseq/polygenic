"""
Test PgsCompute
"""
import json

from unittest import TestCase
from pathlib import Path as path

from polygenic import pgstk

class PgsComputeTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(PgsComputeTest, self).__init__(*args, **kwargs)
        self.output_directory = "/tmp/polygenic/test"
        path(self.output_directory).mkdir(parents=True, exist_ok=True)

    # tests if can be executed directly
    def testPgsComputeDirectly(self):
        appendix = "directly"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test-vcf-general.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test-model-score.yml',
            '--output-name-appendix', appendix,
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/testsample-test-model-score-" + appendix + "-result.json", 'r') as output:
            results = json.load(output)
            self.assertTrue("score_model" in results)
            self.assertTrue("description" in results)
            self.assertTrue("randomentry" not in results)

    # test if can be executed through pgstk
    def testPgsComputePgstk(self):
        appendix = "pgstk"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test-vcf-general.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test-model-score.yml', "polygenic/tests/resources/model/test-model-score.yml",
            '--output-name-appendix', appendix,
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/testsample-test-model-score-" + appendix + "-result.json", 'r') as output:
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
            '--vcf', 'polygenic/tests/resources/vcf/test-vcf-general.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test-model-score.yml',
            '--output-name-appendix', appendix,
            '--af', 'polygenic/tests/resources/vcf/test-af.vcf.gz',
            '--af-field', 'AF_nfe',
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/testsample-test-model-score-" + appendix + "-result.json", 'r') as output:
            results = json.load(output)
            self.assertTrue("score_model" in results)
            self.assertEquals("af", results["genotypes"]["rs1570391830"]["genotype"]["source"])
            self.assertTrue("description" in results)
            self.assertTrue("randomentry" not in results)

    # test haplotype model with unphased
    def testPgsComputeUnphasedHaplotype(self):
        appendix = "unphased.haplotype"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test-vcf-haplotype-unphased.vcf.gz',
            '--model', 'models/pgx/cyp26a1-pharmvar-5.1.8.yml',
            '--output-name-appendix', appendix,
            '--output-directory', self.output_directory,
            '--print'])

    def testPgsComputePhasedHaplotype(self):
        appendix = "phased.haplotype"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test-vcf-cyp2d6-14-28.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test-model-cyp2d6.yml',
            '--output-name-appendix', appendix,
            '--output-directory', self.output_directory,
            '--print'])

    def testPgsComputeDiplotype(self):
        appendix = "diplotype"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test-vcf-general.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test-model-diplotype.yml',
            '--output-name-appendix', appendix,
            '--af', 'polygenic/tests/resources/vcf/test-af.vcf.gz',
            '--output-directory', self.output_directory])

        with open(self.output_directory + "/testsample-test-model-diplotype-" + appendix + "-result.json", 'r') as output:
            results = json.load(output)
            self.assertTrue("diplotype_model" in results)
            self.assertTrue("genotypes" in results)
            self.assertTrue("qc" in results["diplotype_model"])

    # test if af is used
    def test_pgs_compute_merge(self):
        """
        Test if can merge results
        """
        appendix = "merge"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test-vcf-general.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test-model-score.yml', 'polygenic/tests/resources/model/test-model-diplotype.yml',
            '--output-name-appendix', appendix,
            '--merge-outputs',
            '--output-directory', self.output_directory])

        output_file_name = "testsample-" + appendix + "-result.json"


        #with open(self.output_directory + "/" + output_file_name, mode='r', encoding="utf-8") as output:
        #    results = json.load(output)
        #    self.assertTrue("test-model-score" in results)

        appendix = "merge-array"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test-vcf-general.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test-model-score.yml', 'polygenic/tests/resources/model/test-model-diplotype.yml',
            '--output-name-appendix', appendix,
            '--merge-outputs',
            '--merge-as-array',
            '--output-directory', self.output_directory])

        output_file_name = "testsample-" + appendix + "-result.json"


    def test_pgs_compute_genotype_effect_allele(self):
        """
        Test if can output empty genotype effect allele
        """
        appendix = "effect_allele"
        pgstk.main([
            'pgs-compute',
            '--vcf', 'polygenic/tests/resources/vcf/test-vcf-general.vcf.gz',
            '--model', 'polygenic/tests/resources/model/test-model-noeff.yml',
            '--output-name-appendix', appendix,
            '--merge-outputs',
            '--output-directory', self.output_directory])

        output_file_name = "testsample-" + appendix + "-result.json"
        with open(self.output_directory + "/" + output_file_name, mode='r', encoding="utf-8") as output:
            results = json.load(output)
            self.assertTrue("genotypes" in results["test-model-noeff"])

    # # test diplotype model categories
    # def testPgsComputeDiplotype(self):
    #     appendix = "diplotype"
    #     pgstk.main([
    #         'pgs-compute',
    #         '--vcf', 'polygenic/tests/resources/vcf/test.sample.vcf.gz',
    #         '--model', 'polygenic/tests/resources/model/diplotype_model.yml',
    #         '--output-name-appendix', appendix,
    #         '--af', 'polygenic/tests/resources/vcf/test.af.vcf.gz',
    #         '--merge-outputs',
    #         '--output-directory', self.output_directory])

    #     with open(self.output_directory + "/testsample-diplotype_model.yml-" + appendix + "-result.json", 'r') as output:
    #         results = json.load(output)
    #         self.assertTrue("test.model" in results)

    # # test error
    # def testPgsComputeError(self):
    #     appendix = "error"
    #     with self.assertRaises(SystemExit) as cm:
    #         pgstk.main([
    #             'pgs-compute',
    #             '--vcf', 'polygenic/tests/resources/vcf/test.sample.vcf.gz',
    #             '--model', 'polygenic/tests/resources/model/test.model.yml',
    #             '--output-name-appendix', appendix,
    #             '--af', 'polygenic/tests/resources/vcf/test.af.vcf.gz',
    #             '--af-field', 'bobodo',
    #             '--output-directory', self.output_directory])

    #     self.assertEqual(cm.exception.code, 1)



    # # test bug in haplotype model with missing variant
    # def testPgcComputeCypBug(self):
    #     appendix = "cypbug"
    #     pgstk.main([
    #         'pgs-compute',
    #         '--vcf', 'polygenic/tests/resources/vcf/test-pgs-compute-cyp-bug.vcf.gz',
    #         '--model', 'models/pgx/cyp2d6-pharmvar.yml',
    #         '--output-name-appendix', appendix,
    #         '--merge-outputs',
    #         '--output-directory', self.output_directory])

    #     with open(self.output_directory + "/2824-" + appendix + "-result.json", 'r') as output:
    #         results = json.load(output)
    #         self.assertTrue("cyp2d6-pharmvar" in results)
    #         self.assertTrue("haplotype_model" in results["cyp2d6-pharmvar"])

    # def testCyp2c19(self):
    #     appendix = "cyp2c19"
    #     pgstk.main([
    #         'pgs-compute',
    #         '--vcf', 'polygenic/tests/resources/vcf/2824.vcf.gz',
    #         '--model', 'models/pgx/cyp2c19-pharmvar-5.1.8.yml',
    #         '--output-name-appendix', appendix,
    #         '--af', 'polygenic/tests/resources/vcf/test.af.vcf.gz',
    #         '--af-field', 'AF_nfe',
    #         '--output-directory', self.output_directory])

    #     with open(self.output_directory + "/2824-cyp2c19-pharmvar-5.1.8.yml-cyp2c19-result", 'r') as output:
    #         results = json.load(output)
    #         self.assertTrue("cyp2c19-pharmvar" in results)
    #         self.assertTrue("haplotype_model" in results["cyp2c19-pharmvar"])