import subprocess
from unittest import TestCase
from polygenic import pgstk
from pathlib import Path as path

class GwasFileCreateTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(GwasFileCreateTest, self).__init__(*args, **kwargs)
        self.output_directory = "/tmp/polygenic/test"
        path(self.output_directory).mkdir(parents=True, exist_ok=True)

    def test_gwas_file_create_from_biobankuk(self):



        pgstk.main([
            '--log-stdout',
            '--log-level', 'DEBUG',
            'gwas-file-create',
            '--input', 'polygenic/tests/resources/tsv/biobankuk-volume-10klines.tsv',
            '--keyword', 'meta_hq',
            '--output', self.output_directory + '/gwas_file.dat',
        ])