"""
vcf-impute
"""

import pysam

def read_vcf(vcf_file, region):
    with pysam.VariantFile(vcf_file, "r") as vcf:
        if region:
            vcf_iterator = vcf.fetch(region=region)
        return vcf_iterator

def get_header(vcf_file):
    with pysam.VariantFile(vcf_file, "r") as vcf:
        return vcf.header

def run(args):
    print("Running VcfImpute")
    ref = read_vcf(args.reference, args.region)
    target = read_vcf(args.vcf, args.region)

    samples = len(get_header(args.reference).samples)
    size = sum(1 for _ in ref)
    

    array_pos = [0] * size
    row_array = [[0] * samples] * size
    column_array = [[0] * size] * samples

#    for record in ref:
#        #print(record)
#        print(str(record.chrom) + " " + str(record.pos))
#        break