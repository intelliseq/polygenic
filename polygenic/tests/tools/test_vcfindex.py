"""
Test VcfIndex
"""
import os

from unittest import TestCase
from pathlib import Path as path

from polygenic import pgstk


class VcfIndexTest(TestCase):
    """
    Test VcfIndex
    """

    def __init__(self, *args, **kwargs):
        super(VcfIndexTest, self).__init__(*args, **kwargs)

    def test_vcf_index(self):
        """
        Test VcfIndex
        """

        ### delete index if exists
        path_to_vcf = "polygenic/tests/resources/vcf/test-vcf-general.vcf.gz"
        if os.path.exists(path_to_vcf + ".idx.db"):
            os.remove(path_to_vcf + ".idx.db")


        ### run vcf-index
        pgstk.main([
            "--log-stdout",
            "--log-level", "DEBUG",
            "vcf-index",
            "--vcf", path_to_vcf
        ])

        ### define output path
        result_path = path_to_vcf + ".idx.db"

        ### check if index file exists
        self.assertTrue(path(result_path).exists())