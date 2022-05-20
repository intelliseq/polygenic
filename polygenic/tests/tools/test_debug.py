from unittest import TestCase
from polygenic import pgstk
from pathlib import Path as path

class DebugTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(DebugTest, self).__init__(*args, **kwargs)

    def test_debug(self):
        appendix = "af"
        pgstk.main([
            'pgs-compute',
            '--vcf', '/tmp/tmp/inputs/1912153359/202504550232_R11C02.vcf.gz',
            '--model', '/tmp/tmp/inputs/371721200/alcohol.yml',
            '--output-name-appendix', appendix,
            '--af', '/home/marpiech/data/vcf/gnomad.3.1.vcf.gz',
            '--af-field', 'nfe',
            '--output-directory', '/tmp/polygenic/test'])
            