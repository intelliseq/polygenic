"""
Test vcfstatzygosity
"""
import os

from unittest import TestCase
from pathlib import Path as path

from polygenic import pgstk


class VcfStatZygosityTest(TestCase):
    """
    Test vcfstatzygosity
    """

    def __init__(self, *args, **kwargs):
        super(VcfStatZygosityTest, self).__init__(*args, **kwargs)

    def test_vcf_stat_baf(self):
        """
        Test VcfZygosity
        """

        path_to_vcf = "polygenic/tests/resources/vcf/test-vcf-general.vcf.gz"

        pgstk.main([
            "--log-stdout",
            "--log-level", "DEBUG",
            "vcf-stat-zygosity",
            "--vcf", path_to_vcf,
            "--output", "/tmp/polygenic/tests/vsfstatzygosity/zygosity.json"
        ])

        ### define output path
        result_path = "/tmp/polygenic/tests/vsfstatzygosity/"

        ### check if index file exists
        self.assertTrue(path(result_path).exists())