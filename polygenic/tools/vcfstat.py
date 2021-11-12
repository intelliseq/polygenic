import argparse
import logging

import sys
import os
import urllib.request

# zygosity
import json

# baf
from pysam import VariantFile
import pandas as pandas
import plotly as plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def zygosity(args):
    parser = argparse.ArgumentParser(description='vcfstat zygosity draws baf images for chromosomes')  # todo dodać opis
    parser.add_argument('-v', '--vcf', type=str, required=True, help='input vcf.gz')
    parser.add_argument('-o', '--output-file', type=str, default='stats.json', help='output directory')
    parsed_args = parser.parse_args(args)
    chrom_sizes = pandas.read_csv(
        os.path.dirname(__file__) + '/resources/chromsizes/hg38.chrom.sizes', 
        delimiter='\t',
        index_col=0)
    chromosomes = list(chrom_sizes.index)
    vcf = VariantFile(parsed_args.vcf)

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
            out[sample_index]["qc"][chrom + "heterozygosity"]["plot"] = sample_name + "-" + chrom + ".jpg"

        out[sample_index]["qc"]["baflrr"] = {}
        out[sample_index]["qc"]["baflrr"]["value"] = "Click to view plots"
        out[sample_index]["qc"]["baflrr"]["pass"] = "passed"
        out[sample_index]["qc"]["baflrr"]["plot"] = sample_name + "-" + chrom + ".pdf"

    with open(parsed_args.output_file, 'w') as qcfile:
        json.dump(out, qcfile)
    
    return

def baf(args):
    parser = argparse.ArgumentParser(description='vcfstat baf draws baf images for chromosomes')  # todo dodać opis
    parser.add_argument('-v', '--vcf', type=str, required=True, help='input vcf.gz')
    parser.add_argument('-o', '--output-directory', type=str, default='', help='output directory')
    parsed_args = parser.parse_args(args)
    chrom_sizes = pandas.read_csv(
        os.path.dirname(__file__) + '/resources/chromsizes/hg38.chrom.sizes', 
        delimiter='\t',
        index_col=0)
    chromosomes = list(chrom_sizes.index)
    vcf = VariantFile(parsed_args.vcf)

    if not os.path.exists(parsed_args.output_directory):
        os.mkdir(parsed_args.output_directory)

    for sample_index in range(len(vcf.header.samples)):
        sample_name = vcf.header.samples[sample_index]
        for chrom in chromosomes:
            size = chrom_sizes.loc[chrom]['size']
            df = pandas.DataFrame()
            position = []
            baf = []
            lrr = []

            for record in vcf.fetch(chrom, 1, size):
                if 'BAF' in record.format:
                    position.append(record.pos)
                    baf.append(record.samples.values()[sample_index]['BAF'])
                    lrr.append(record.samples.values()[sample_index]['LRR'])
            df = pandas.DataFrame({'position': position, 'baf': baf, 'lrr': lrr})
        
            fig = make_subplots(
                rows=2, 
                cols=1,
                shared_xaxes=True,
                vertical_spacing=0.02)

            fig.add_trace(
                go.Scatter(
                    x = list(df['position']), 
                    y = list(df['baf']),
                    mode = "markers",
                    
                ),
                row=1, col=1
            )

            fig.add_trace(
                go.Scatter(
                    x = list(df['position']), 
                    y = list(df['lrr']),
                    mode = "markers"
                ),
                row=2, col=1
            )

            fig.update_traces(
                marker=dict(size=0.8)
            )

            fig.update_layout(
                height=600, 
                width=800, 
                title_text=sample_name + " " + chrom,
                template='simple_white',
                showlegend = False
                )

            fig.update_yaxes(
                row=1, col=1,
                title_text = "BAF",
                range = [-0.1, 1.1]
            )

            fig.update_yaxes(
                row=2, col=1,
                title_text = "LRR",
                range = [-1, 1]
            )

            fig.write_image(parsed_args.output_directory + "/" + sample_name + "-" + chrom + ".jpg")

    return

def main(args = sys.argv[1:]):
    #try:
        if args[0] == 'baf':
            baf(args[1:])
        elif args[0] == 'zygosity':
            zygosity(args[1:])
        else:
            raise Exception()
    # except Exception as e:
    #     print(e)
    #     print("ERROR " + str(e))
    #     print("""
    #     Program: vcfstat (calculates stats for vcfs)
    #     Contact: Marcin Piechota <piechota.marcin@gmail.com>

    #     Usage:   vcfstat <command> [options]

    #     Command: 
    #     baf                downloads pan biobankuk index of gwas results
    #     zygosity           downloads gwas results for given phenocode
    #     """)

if __name__ == '__main__':
    main(sys.argv[1:])
