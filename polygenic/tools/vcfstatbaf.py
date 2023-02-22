import os
import pandas as pandas
import plotly as plotly
import plotly.graph_objects as go
from pysam import VariantFile
from plotly.subplots import make_subplots
from importlib_resources import files

def run(parsed_args):

    chrom_sizes = pandas.read_csv(
        files('polygenic.resources.chromsizes').joinpath('hg38.chrom.sizes'),
        delimiter = '\t',
        index_col = 0
    )
    chromosomes = list(chrom_sizes.index)
    vcf = VariantFile(parsed_args.vcf)

    if not os.path.exists(parsed_args.output_directory):
        os.makedirs(parsed_args.output_directory)

    for sample_index in range(len(vcf.header.samples)):
        sample_name = vcf.header.samples[sample_index]
        for chrom in chromosomes:
            size = chrom_sizes.loc[chrom]['size']
            df = pandas.DataFrame()
            position = []
            baf = []
            lrr = []

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

            fig.write_image(parsed_args.output_directory + "/" + sample_name + "-" + chrom + ".jpeg")

    return
