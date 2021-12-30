from unittest import TestCase
from polygenic.tools import pgscompute
from polygenic import pgstk

from pathlib import Path as path

import tabix
import os
import configparser

class ModelPgscatTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(ModelPgscatTest, self).__init__(*args, **kwargs)
        self.output_directory = "/tmp/polygenic/test"
        path(self.output_directory).mkdir(parents=True, exist_ok=True)

    def testModelPgscatGeneric(self):
        pgstk.main([
            "model-pgscat",
            "--code", "PGS000011",
            "--output-directory", self.output_directory,
            "--af", "polygenic/tests/resources/largefiles/gnomad.3.1.vcf.gz",
            "--origin-genome-build", "Grch38",
            "--source-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch37.norm.vcf.gz",
            "--target-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch38.norm.vcf.gz",
            "--gene-positions", "polygenic/tests/resources/largefiles/ensembl-genes.104.tsv"
        ])

        result_path = self.output_directory + "/pgscat-coronary_artery_disease-PGS000011-EUR.yml"

        with open(result_path, 'r') as output:
            data = output.read()
            header = list(filter(lambda line: "score_model:" in line, data.split('\n')))
            self.assertEqual(1, len(header))

    def testModelPgscatMissingPosition(self):
        pgstk.main([
            "model-pgscat",
            "--code", "PGS000004",
            "--output-directory", self.output_directory,
            "--af", "polygenic/tests/resources/largefiles/gnomad.3.1.vcf.gz",
            "--source-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch37.norm.vcf.gz",
            "--target-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch38.norm.vcf.gz",
            "--gene-positions", "polygenic/tests/resources/largefiles/ensembl-genes.104.tsv"
        ])

        result_path = self.output_directory + "/pgscat-breast_cancer-PGS000004-EUR.yml"

        with open(result_path, 'r') as output:
            data = output.read()
            header = list(filter(lambda line: "score_model:" in line, data.split('\n')))
            self.assertEqual(1, len(header))

    def testModelPgscatCreateModel(self):
        pgstk.main([
            "model-pgscat",
            "--code", "PGS000005",
            "--output-directory", self.output_directory,
            "--af", "polygenic/tests/resources/largefiles/gnomad.3.1.vcf.gz",
            "--af-field", "AF_nfe",
            "--origin-genome-build", "Grch37",
            "--source-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch37.norm.vcf.gz",
            "--target-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch38.norm.vcf.gz",
            "--gene-positions", "polygenic/tests/resources/largefiles/ensembl-genes.104.tsv"
        ])

        result_path = self.output_directory + "/pgscat-breast_cancer-PGS000004-EUR.yml"

        with open(result_path, 'r') as output:
            data = output.read()
            header = list(filter(lambda line: "score_model:" in line, data.split('\n')))
            self.assertEqual(1, len(header))