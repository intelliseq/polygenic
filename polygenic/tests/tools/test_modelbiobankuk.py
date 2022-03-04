from unittest import TestCase
from polygenic.tools import pgscompute
from polygenic import pgstk

from pathlib import Path as path

import tabix
import os
import configparser

class ModelBiobankukTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(ModelBiobankukTest, self).__init__(*args, **kwargs)
        self.output_directory = "/tmp/polygenic/test"
        path(self.output_directory).mkdir(parents=True, exist_ok=True)

    def testModelBiobankukPgstk(self):
        pgstk.main([
            "model-biobankuk",
            "--code", "2395",
            "--sex", "both_sexes",
            "--coding", "4", 
            "--output-directory", self.output_directory,
            "--pvalue-threshold", "1e-08",
            "--clumping-vcf", "polygenic/tests/resources/largefiles/eur.phase3.biobank.set.vcf.gz",
            "--source-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch37.norm.vcf.gz",
            "--target-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch38.norm.vcf.gz",
            "--gene-positions", "polygenic/tests/resources/largefiles/ensembl-genes.104.tsv",
            "--test"
        ])

        result_path = self.output_directory + "/biobankuk-2395-both_sexes-4-hair_balding_pattern-EUR-1e-08.yml"

        with open(result_path, 'r') as output:
            data = output.read()
            header = list(filter(lambda line: "score_model:" in line, data.split('\n')))
            self.assertEqual(1, len(header))

    def testModelBiobankukCode(self):
        pgstk.main([
            "model-biobankuk",
            "--code", "thiazolidinedione|diabetes",
            "--sex", "both_sexes",
            "--output-directory", self.output_directory,
            "--pvalue-threshold", "1e-05",
            "--clumping-vcf", "polygenic/tests/resources/largefiles/eur.phase3.biobank.set.vcf.gz",
            "--source-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch37.norm.vcf.gz",
            "--target-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch38.norm.vcf.gz",
            "--gene-positions", "polygenic/tests/resources/largefiles/ensembl-genes.104.tsv",
            "--test", "true"
        ])

        result_path = self.output_directory + "/biobankuk-2395-both_sexes-4-hair_balding_pattern-EUR-1e-08.yml"

        with open(result_path, 'r') as output:
            data = output.read()
            header = list(filter(lambda line: "score_model:" in line, data.split('\n')))
            self.assertEqual(1, len(header))

    def testModel22046Bug(self):
        pgstk.main([
            "model-biobankuk",
            "--code", "22406",
            "--sex", "both_sexes",
            "--output-directory", self.output_directory,
            "--pvalue-threshold", "1e-12",
            "--clumping-vcf", "polygenic/tests/resources/largefiles/eur.phase3.biobank.set.vcf.gz",
            "--source-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch37.norm.vcf.gz",
            "--target-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch38.norm.vcf.gz",
            "--gene-positions", "polygenic/tests/resources/largefiles/ensembl-genes.104.tsv",
            "--test", "true"
        ])

    def testModelWhitespaceBug(self):
        pgstk.main([
            "model-biobankuk",
            "--code", "vitamin D derivative",
            "--sex", "both_sexes",
            "--output-directory", self.output_directory,
            "--pvalue-threshold", "1e-12",
            "--clumping-vcf", "polygenic/tests/resources/largefiles/eur.phase3.biobank.set.vcf.gz",
            "--source-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch37.norm.vcf.gz",
            "--target-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch38.norm.vcf.gz",
            "--gene-positions", "polygenic/tests/resources/largefiles/ensembl-genes.104.tsv",
            "--test", "true"
        ])

        result_path = self.output_directory + "/biobankuk-vitamin_d_derivative-both_sexes--na-EUR-1e-12.yml"

        # check if result path exists
        self.assertTrue(os.path.exists(result_path))

    def testModelApostropheBug(self):
        pgstk.main([
            "model-biobankuk",
            "--code", "acetylcholinesterase inhibitor|Alzheimer's",
            "--sex", "both_sexes",
            "--output-directory", self.output_directory,
            "--pvalue-threshold", "1e-12",
            "--clumping-vcf", "polygenic/tests/resources/largefiles/eur.phase3.biobank.set.vcf.gz",
            "--source-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch37.norm.vcf.gz",
            "--target-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch38.norm.vcf.gz",
            "--gene-positions", "polygenic/tests/resources/largefiles/ensembl-genes.104.tsv",
            "--test", "true"
        ])

        result_path = self.output_directory + "/acetylcholinesterase_inhibitor_alzheimer_s-both_sexes--na-EUR-1e-12.yml"

        # check if result path exists
        self.assertTrue(os.path.exists(result_path))

    def testFilteringByPBug(self):
        pgstk.main([
            "model-biobankuk",
            "--code", "20115",
            "--sex", "both_sexes",
            "--coding", "605",
            "--output-directory", self.output_directory,
            "--pvalue-threshold", "1e-8",
            "--clumping-vcf", "polygenic/tests/resources/largefiles/eur.phase3.biobank.set.vcf.gz",
            "--source-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch37.norm.vcf.gz",
            "--target-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch38.norm.vcf.gz",
            "--gene-positions", "polygenic/tests/resources/largefiles/ensembl-genes.104.tsv"
        ])

        result_path = self.output_directory + "/acetylcholinesterase_inhibitor_alzheimer_s-both_sexes--na-EUR-1e-12.yml"

        # check if result path exists
        self.assertTrue(os.path.exists(result_path))
