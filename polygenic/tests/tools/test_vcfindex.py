from unittest import TestCase
from polygenic import pgstk
from pathlib import Path as path

class VcfIndexTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(VcfIndexTest, self).__init__(*args, **kwargs)

    def testVcfIndex(self):

        ### run vcf-index
        pgstk.main([
            "vcf-index",
            "--vcf", "polygenic/tests/resources/vcf/test.vcf.gz"
        ])

        ### define output path
        result_path = "polygenic/tests/resources/vcf/test.vcf.gz.idx.db"

        ### check if index file exists
        self.assertTrue(path(result_path).exists())