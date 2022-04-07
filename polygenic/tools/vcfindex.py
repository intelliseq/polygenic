import logging
from polygenic.data.vcf_accessor import VcfAccessor

def run(args):
    VcfAccessor(args.vcf)
    return
