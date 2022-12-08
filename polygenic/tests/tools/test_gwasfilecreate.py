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
            '--output', self.output_directory + '/gwas_file.dat',
            '--chromosome-column-name', 'chr',
            '--position-column-name', 'pos',
            '--ref-allele-column-name', 'ref',
            '--alt-allele-column-name', 'alt',
            '--effect-allele-column-name', 'alt',
            '--pvalue-column-name', 'pval_meta',
            '--beta-column-name', 'beta_meta',
            '--rsid-column-name', 'rsid',
            '--af-column-name', 'af_meta'
        ])

        ### check if index file exists
        #self.assertTrue(path(result_path).exists())