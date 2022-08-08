from unittest import TestCase
from polygenic import pgstk
from pathlib import Path as path

class ModelGwasFileTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(ModelGwasFileTest, self).__init__(*args, **kwargs)

    def test_vcf_index(self):

        ### run vcf-index
        pgstk.main([
            "--log-stdout",
            "--log-level", "DEBUG",
            "model-gwas-file",
            "--gwas-file", "/home/marpiech/downloads/EA4_additive_excl_23andMe_2.txt.gz",
            "--output", "/home/marpiech/downloads/EA4_additive_excl_23andMe.yml"
        ])

        ### define output path
        #result_path = "polygenic/tests/resources/vcf/test.vcf.gz.idx.db"

        ### check if index file exists
        #self.assertTrue(path(result_path).exists())