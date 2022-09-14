"""
vcf-impute
"""

import pysam
import numpy
from pathlib import Path
from collections import deque as Deque

def read_vcf(vcf_file, region):
    vcf = pysam.VariantFile(vcf_file, "r")
    if not region:
        return None
    vcf_iterator = vcf.fetch(region=region)
    return vcf_iterator

def get_header(vcf_file):
    with pysam.VariantFile(vcf_file, "r") as vcf:
        return vcf.header

def get_target_genotypes_as_dictionary(target_vcf, region):
    target = read_vcf(target_vcf, region)
    target_dict = {}
    for record in target:
        target_dict[record.pos] = record.samples.values()[0]['GT']
    return target_dict

def get_samples_count_in_vcf(vcf_file):
    with pysam.VariantFile(vcf_file, "r") as vcf:
        return len(vcf.header.samples)

def impute(index: int, ld_threshold: float, ref_row_array: list, target_row_array: list):
    nrow = len(ref_row_array)
    ncol = len(ref_row_array[0])
    index_row = numpy.array(ref_row_array[index]).astype(int)
    ldproxy = [0] * nrow

    # calculate correlations between all rows and the row in the middle of the window
    for i in range(nrow):
        if i == index:
            continue
        row = numpy.array(ref_row_array[i]).astype(int)
        row_sum = numpy.sum(row)
        ### print numer of row elements that are Nan
        if row_sum == 0 or row_sum == ncol:
            continue
        ldproxy[i] = numpy.corrcoef(row, index_row)[0,1]
    # if correlation is above threshold, impute
    ldproxy_order = (-(numpy.abs(ldproxy))).argsort()
    for ldproxy_index in ldproxy_order:
        if (not (abs(ldproxy[ldproxy_index]) > ld_threshold)) or (target_row_array[ldproxy_index][0] is None):
            continue
        proxy_invert = 1 if ldproxy[ldproxy_index] > 0 else -1
        return [target_row_array[ldproxy_index][0] * proxy_invert, target_row_array[ldproxy_index][0] * proxy_invert, abs(ldproxy[ldproxy_index])]
    return [None, None, None]

def run(args):

    # set comments
    numpy.seterr(divide='ignore', invalid='ignore')

    # prepare files
    ref_vcf = read_vcf(args.reference, args.region)
    ref_alleles_count = get_samples_count_in_vcf(args.reference) * 2
    
    # init arrays
    ref_pos_array = Deque([0] * args.window, maxlen=args.window)
    ref_row_array = Deque([[0] * (ref_alleles_count)] * args.window, maxlen=args.window)
    ref_column_array = [Deque([0] * args.window, maxlen = args.window)] * (ref_alleles_count)
    target_row_array = Deque([[0] * 2] * args.window, maxlen=args.window)
    target_dict = get_target_genotypes_as_dictionary(args.vcf, args.region)
    target_result = {}

    counter = 0
    middle_index = int(args.window / 2)

    with pysam.VariantFile(args.vcf, "r") as target_vcf:
        for record in target_vcf.fetch(region=args.region):
            target_result[record.pos] = record

    for record in ref_vcf:

        # include only records ./. for imputing
        if args.missing_only: 
            if record.pos not in target_dict:
                continue

        # add next position to array
        ref_pos_array.append(record.pos)

        # add next row to row array
        row = [genotype for sample in record.samples.values() for genotype in sample['GT']]
        ref_row_array.append(row)

        # add next row to column array
        for i in range(ref_alleles_count):
            (ref_column_array[i]).append(row[i])

        # add next haplotypes to target haplotypes
        target_row_array.append(list(target_dict.get(record.pos, [None,None])))

        # increment counter
        counter = counter + 1

        if counter < args.window:
            continue
        if counter == args.window:
            for i in range(middle_index):
                pass
        if counter >= args.window:
            # impute if enough data
            if target_row_array[middle_index] == [None, None]:
                imputed = impute(middle_index, args.ld_threshold, ref_row_array, target_row_array)
            
    target_vcf = pysam.VariantFile(args.vcf, "r")
    Path(args.output).parent.absolute().mkdir(parents=True, exist_ok=True)
    output_vcf = pysam.VariantFile(args.output, "w", header=target_vcf.header)
    output_vcf.header.info.add("IMP_PROB","1","Float","Imputation correctness probability")
    for record in target_result.values():
        output_vcf.write(record)
    

        
