from unittest import TestCase
from polygenic import pgstk
from pathlib import Path as path

class ModelGwasFileTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(ModelGwasFileTest, self).__init__(*args, **kwargs)
        self.output_directory = "/tmp/polygenic/test"
        path(self.output_directory).mkdir(parents=True, exist_ok=True)


    def test_create_gwas_file_from_biobankuk(self):

        pgstk.main([
            "--log-stdout",
            "--log-level", "DEBUG",
            "model-gwas-file",
            "--gwas-file", "polygenic/tests/resources/largefiles/dat/continuous-3063-both_sexes-irnt.tsv.gz",
            "--output", self.output_directory + "/gwas_model.yml"
        ])

        ### define output path
        #result_path = "polygenic/tests/resources/vcf/test.vcf.gz.idx.db"

        ### check if index file exists
        #self.assertTrue(path(result_path).exists())