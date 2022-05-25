from unittest import TestCase
from polygenic import pgstk
from pathlib import Path as path

class DebugTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(DebugTest, self).__init__(*args, **kwargs)

    def test_debug(self):
        appendix = "debug"
        pgstk.main([
            'pgs-compute',
            '--vcf', '/home/marpiech/data/vcf/sportsmen-control.vcf.gz',
            '--model', 'polygenic/tests/resources/model/biobankuk-23121-both_sexes--arm_fat_free_mass_right_-EUR-1e-08.yml',
            '--sample-name', 'B539',
            '--output-name-appendix', appendix,
            '--af', '/home/marpiech/data/vcf/gnomad.3.1.vcf.gz',
            '--af-field', 'AF_nfe',
            '--output-directory', '/tmp/polygenic/test'])
            