"""
vcf-impute
"""


from pathlib import Path
from collections import deque as Deque

import numpy
from pysam import VariantFile

def read_vcf(vcf_file, region):
    """
    vcf file path to record iterator
    """
    vcf = VariantFile(vcf_file, "r")
    if not region:
        return None
    vcf_iterator = vcf.fetch(region=region)
    return vcf_iterator

def get_header(vcf_file):
    with VariantFile(vcf_file, "r") as vcf:
        return vcf.header

def get_target_genotypes_as_dictionary(target_vcf, region):
    target = read_vcf(target_vcf, region)
    target_dict = {}
    for record in target:
        target_dict[record.pos] = [genotype for sample in record.samples.values() for genotype in sample['GT']]
    return target_dict

def get_samples_count_in_vcf(vcf_file):
    with VariantFile(vcf_file, "r") as vcf:
        return len(vcf.header.samples)

def invert_genotype(genotype):
    if genotype is None:
        return None
    if genotype == 0:
        return 1
    if genotype == 1:
        return 0
    return 1

def __impute(index: int, ld_threshold: float, ref_row_array: list, target_row_array: list):
    ref_nrow = len(ref_row_array)
    ref_ncol = len(ref_row_array[0])
    target_ncol = len(target_row_array[0])

    index_row = numpy.array(ref_row_array[index]).astype(int)
    ldproxy = [0] * ref_nrow

    output_genotype = [None] * target_ncol
    output_probability = None

    # calculate correlations between all rows and the row in the middle of the window
    for i in range(ref_nrow):
        if i == index:
            continue
        row = numpy.array(ref_row_array[i]).astype(int)
        row_sum = numpy.sum(row)
        ### print numer of row elements that are Nan
        if row_sum == 0 or row_sum == ref_ncol:
            continue
        ldproxy[i] = numpy.corrcoef(row, index_row)[0,1]
    # if correlation is above threshold, impute
    ldproxy_order = (-(numpy.abs(ldproxy))).argsort()
    for sample_index in range(target_ncol):
        if target_row_array[index][sample_index] is not None:
            continue
        for ldproxy_index in ldproxy_order:
            proxy_ld = ldproxy[ldproxy_index]
            
            proxy_genotype = target_row_array[ldproxy_index][sample_index]
            if (not (abs(proxy_ld) > ld_threshold)) or (proxy_genotype is None):
                continue
            if proxy_ld < ld_threshold:
                break
            output_probability = abs(proxy_ld)
            output_genotype[sample_index] = proxy_genotype if proxy_ld > 0 else invert_genotype(proxy_genotype)
            break 
    return (output_genotype, output_probability)


def __update_record(target_result, position, genotype, probability):
    for sample_index in range(int(len(genotype) / 2)):
        target_result[position].samples.values()[sample_index]['GT'] = (genotype[sample_index * 2], genotype[sample_index * 2 + 1])
    target_result[position].info['IMP_PROB'] = probability

def impute(index: int, ld_threshold: float, ref_pos_array: list, ref_row_array: list, target_row_array: list, target_result: dict):
    if None not in target_row_array[index]:
        return
    result = __impute(index, ld_threshold, ref_row_array, target_row_array)
    if None not in result:
        __update_record(target_result, ref_pos_array[index], result[:-1], result[-1])

def __phase_by_proxy(proxy_genotype, proxy_ld, target_genotype):
    if proxy_genotype is None:
        return None
    if proxy_ld > 0:
        if target_genotype[0] == proxy_genotype[0]:
            return target_genotype
        if target_genotype[0] == proxy_genotype[1]:
            return (target_genotype[1], target_genotype[0])
        return None
    else:
        if target_genotype[0] != proxy_genotype[0]:
            return target_genotype
        if target_genotype[0] != proxy_genotype[1]:
            return (target_genotype[1], target_genotype[0])
        return None

def phase(index: int, ld_threshold: float, ref_pos_array: list, ref_row_array: list, ref_col_array:list, target_row_array: list, target_result: dict):

    ld_threshold = 0.5

    # check if there are any phased proxies in window
    if sum([1 for record in target_result.values() if record.samples[0].phased]) == 0:
        for sample_index in range(int(len(target_row_array[index]) / 2)):
            target_result[ref_pos_array[index]].samples[sample_index].phased = True
            target_result[ref_pos_array[index]].info.update({'PHASE_PROB': 0.5})

    ref_nrow = len(ref_row_array)
    ref_ncol = len(ref_row_array[0])
    target_ncol = len(target_row_array[0])
    target_nsamples = int(target_ncol / 2)
    ldproxy = [0] * ref_nrow
    index_row = numpy.array(ref_row_array[index]).astype(int)

    # calculate ld proxy with other snps in the window
    for i in range(ref_nrow):
        if i == index:
            continue
        row = numpy.array(ref_row_array[i]).astype(int)
        row_sum = numpy.sum(row)
        if row_sum == 0 or row_sum == ref_ncol:
            continue
        ldproxy[i] = numpy.corrcoef(row, index_row)[0,1]

    # if correlation is above threshold, phase
    ldproxy_order = (-(numpy.abs(ldproxy))).argsort()
    for sample_index in range(target_nsamples):
        sample_slice = slice(sample_index * 2, sample_index * 2 + 2)
        for ldproxy_index in ldproxy_order:
            proxy_ld = ldproxy[ldproxy_index]
            if abs(proxy_ld) < ld_threshold:
                break
            is_index_phased = target_result[ref_pos_array[index]].samples[sample_index].phased
            is_proxy_phased = target_result[ref_pos_array[ldproxy_index]].samples[sample_index].phased
            if  not is_index_phased and not is_proxy_phased:
                continue
            if is_index_phased and is_proxy_phased:
                continue
            if not is_index_phased:
                unphased_index = index
                phased_index = ldproxy_index
            if is_index_phased:
                unphased_index = ldproxy_index
                phased_index = index
            proxy_genotype = target_row_array[phased_index][sample_slice]
            phased_genotype = __phase_by_proxy(proxy_genotype, proxy_ld, target_row_array[unphased_index][sample_slice])
            if phased_genotype is None:
                continue
            target_result[ref_pos_array[unphased_index]].samples[sample_index]['GT'] = phased_genotype
            target_result[ref_pos_array[unphased_index]].samples[sample_index].phased = True
            target_result[ref_pos_array[unphased_index]].info.update({'PHASE_PROB': abs(proxy_ld)})

    for sample_index in range(target_nsamples):
        sample_cols = [[row[sample_index * 2] for row in target_row_array], [row[sample_index * 2 + 1] for row in target_row_array]]
        haplotype_cor = [0] * (ref_ncol * 2)
        for i in range(ref_ncol):
            if numpy.sum(ref_col_array) == 0 or numpy.sum(ref_col_array) == ref_nrow:
                continue
            print(sample_cols[0])
            print(ref_col_array[i])
            haplotype_cor[i] = numpy.corrcoef(sample_cols[0], ref_col_array[i])[0,1]
            haplotype_cor[i + ref_ncol] = numpy.corrcoef(sample_cols[1], ref_col_array[i])[0,1]
        print(haplotype_cor)



def run(args):

    # silence numpy warnings
    numpy.seterr(divide='ignore', invalid='ignore')

    # prepare files
    ref_vcf = read_vcf(args.reference, args.region)
    ref_alleles_count = get_samples_count_in_vcf(args.reference) * 2
    target_alleles_count = get_samples_count_in_vcf(args.vcf) * 2

    # init arrays
    ref_pos_array = Deque([0] * args.window, maxlen=args.window)
    ref_row_array = Deque([[0] * (ref_alleles_count)] * args.window, maxlen=args.window)
    ref_col_array = [Deque([0] * args.window, maxlen = args.window)] * (ref_alleles_count)
    target_row_array = Deque([[0] * target_alleles_count] * args.window, maxlen=args.window)
    target_dict = get_target_genotypes_as_dictionary(args.vcf, args.region)
    target_result = {}

    counter = 0
    middle_index = int(args.window / 2)

    with VariantFile(args.vcf, "r") as target_vcf:
        if "IMP_PROB" not in target_vcf.header.info:
            target_vcf.header.info.add("IMP_PROB","1","Float","Imputation correctness probability")
        target_vcf.header.info.add("PHASE_PROB","1","Float","Phasing correctness probability")
        for record in target_vcf.fetch(region = args.region):
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
            (ref_col_array[i]).append(row[i])

        # add next haplotypes to target haplotypes
        target_row_array.append(list(target_dict.get(record.pos, [None] * target_alleles_count)))

        # increment counter
        counter = counter + 1

        if counter < args.window:
            continue

        if counter == args.window:
            for i in range(middle_index):
                #phase(i, args.ld_threshold, ref_pos_array, ref_row_array, ref_col_array, target_row_array, target_result)
                impute(i, args.ld_threshold, ref_pos_array, ref_row_array, target_row_array, target_result)

        if counter >= args.window:
            #phase(middle_index, args.ld_threshold, ref_pos_array, ref_row_array, ref_col_array, target_row_array, target_result)
            impute(middle_index, args.ld_threshold, ref_pos_array, ref_row_array, target_row_array, target_result)

    if counter < args.window:
        for i in range(counter):
            impute(i, args.ld_threshold, ref_pos_array, ref_row_array, target_row_array, target_result)
    else:
        for i in range(middle_index + 1, args.window):
            impute(i, args.ld_threshold, ref_pos_array, ref_row_array, target_row_array, target_result)

    target_vcf = VariantFile(args.vcf, "r")
    Path(args.output).parent.absolute().mkdir(parents=True, exist_ok=True)
    output_vcf = VariantFile(args.output, "w", header=target_vcf.header)
    if "IMP_PROB" not in output_vcf.header.info:
        output_vcf.header.info.add("IMP_PROB","1","Float","Imputation correctness probability")
    #output_vcf.header.info.add("PHASE_PROB","1","Float","Phasing correctness probability")

    for position, record in target_result.items():
        # print(str(position) + " " + str(record.pos))
        output_vcf.write(record)
