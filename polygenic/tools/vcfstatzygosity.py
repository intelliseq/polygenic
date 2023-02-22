import os

# zygosity
import json

# baf
from pysam import VariantFile
import pandas as pandas
import plotly as plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from importlib_resources import files

def run(parsed_args):

    chrom_sizes = pandas.read_csv(
        files('polygenic.resources.chromsizes').joinpath('hg38.chrom.sizes'),
        delimiter='\t',
        index_col=0)
    chromosomes = list(chrom_sizes.index)
    vcf = VariantFile(parsed_args.vcf)

    # get output directory from file path and create if not exists
    output_directory = os.path.dirname(parsed_args.output)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    out = []

    for sample_index in range(len(vcf.header.samples)):
        sample_name = vcf.header.samples[sample_index]
        out.append({})
        out[sample_index]["id"] = sample_name
        out[sample_index]["qc"] = {}
        for chrom in ['chrX', 'chrY', 'chrM']:
            size = chrom_sizes.loc[chrom]['size']
            df = pandas.DataFrame()
            count = 0
            het = 0
            alt = 0
            ref = 0

            ### check if chrom exists in vcf
            ### if not, try prefixing with 'chr'
            if chrom not in vcf.header.contigs:
                chrom = 'chr' + chrom
            ### check if chrom exists in vcf
            ### if not, skip
            if chrom not in vcf.header.contigs:
                continue

            for record in vcf.fetch(chrom, 1, size):
                if 'BAF' in record.format:
                    gt = record.samples.values()[sample_index]['GT']
                    if gt[0] == 0 and gt[1] == 0:
                        count = count + 1; ref = ref + 1
                    if gt[0] == 1 and gt[1] == 1:
                        count = count + 1; alt = alt + 1
                    if (gt[0] == 0 and gt[1] == 1) or (gt[0] == 1 and gt[1] == 0):
                        count = count + 1; het = het + 1
            out[sample_index]["qc"][chrom + "heterozygosity"] = {}
            out[sample_index]["qc"][chrom + "heterozygosity"]["value"] = het / count
            out[sample_index]["qc"][chrom + "heterozygosity"]["pass"] = "passed"
            out[sample_index]["qc"][chrom + "heterozygosity"]["plot"] = sample_name + "-" + chrom + ".jpeg"

        out[sample_index]["qc"]["baflrr"] = {}
        out[sample_index]["qc"]["baflrr"]["value"] = "Click to view plots"
        out[sample_index]["qc"]["baflrr"]["pass"] = "passed"
        out[sample_index]["qc"]["baflrr"]["plot"] = sample_name + "-" + chrom + ".pdf"

    with open(parsed_args.output, 'w') as qcfile:
        json.dump(out, qcfile)
    
    return

