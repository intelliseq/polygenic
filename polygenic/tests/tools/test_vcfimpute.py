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
        args.reference = "/tmp/polygenic/10.38.vcf.gz"
        #args.vcf = "polygenic/tests/resources/vcf/test-impute-target.vcf.gz"
        args.vcf = "/tmp/polygenic/chr10.vcf.gz"
        #args.region = "chr16:53553311-53841786"
        args.region = "chr10:14583-1777224"
        args.out = "polygenic/tests/resources/vcf/test-impute-out.vcf.gz"
        args.window = 60
        args.missing_only = True
        args.ld_threshold = 0.98
        args.output = "/tmp/pgstk-tests/test-impute-out.vcf.gz"
        
        vcfimpute.run(args)

