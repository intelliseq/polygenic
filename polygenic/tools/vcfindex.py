from polygenic.data.vcf_accessor import VcfAccessor as vcf_accessor

def run(args):
    vcf_accessor.create_rsidx_index(args.vcf)
    return 0
