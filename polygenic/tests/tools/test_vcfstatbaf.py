"""
Test vcfstatbaf
"""
import os

from unittest import TestCase
from pathlib import Path as path

from polygenic import pgstk


class VcfStatBafTest(TestCase):
    """
    Test vcfstatbaf
    """

    def __init__(self, *args, **kwargs):
        super(VcfStatBafTest, self).__init__(*args, **kwargs)

    def test_vcf_stat_baf(self):
        """
        Test VcfIndex
        """

        ### delete index if exists
        path_to_vcf = "polygenic/tests/resources/vcf/test-vcf-general.vcf.gz"

        ### run vcf-index
        pgstk.main([
            "--log-stdout",
            "--log-level", "DEBUG",
            "vcf-stat-baf",
            "--vcf", path_to_vcf,
            "--output-directory", "/tmp/polygenic/tests/vcfstatbaf/"
        ])

        ### define output path
        result_path = path_to_vcf + ".idx.db"

        ### check if index file exists
        self.assertTrue(path(result_path).exists())