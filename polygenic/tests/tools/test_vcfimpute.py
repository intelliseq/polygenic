"""
Test VcfImpute
"""
import os
import argparse
import pysam

from unittest import TestCase
import polygenic.tools.vcfimpute as vcfimpute

class VcfImputeTest(TestCase):
    """
    Test VcfImpute
    """

    def __init__(self, *args, **kwargs):
        super(VcfImputeTest, self).__init__(*args, **kwargs)

    def test_vcf_impute(self):
        """
        Test VcfImpute
        """

        args = argparse.Namespace()
        #args.reference = "polygenic/tests/resources/vcf/test-impute-ref-eas.vcf.gz"
        args.reference = "/tmp/chr18.38.vcf.gz"
        #args.vcf = "polygenic/tests/resources/vcf/test-impute-target.vcf.gz"
        args.vcf = "/tmp/18_phased.reheader.vcf.gz"
        #args.region = "chr16:53553311-53841786"
        args.region = "18:15000001-16000001"
        args.out = "polygenic/tests/resources/vcf/test-impute-out.vcf.gz"
        args.window = 60
        args.missing_only = True
        args.ld_threshold = 0.95
        args.output = "/tmp/pgstk-tests/test-impute-out.vcf.gz"
        
        vcfimpute.run(args)

