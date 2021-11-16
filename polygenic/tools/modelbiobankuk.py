import argparse
import sys
import os

def parse_args(args):
    parser = argparse.ArgumentParser(description='pgstk model-biobankuk prepares polygenic score model based on p value data')
    parser.add_argument('--code', '--phenocode', type=str, required=True, help='phenocode of phenotype form Uk Biobank')
    parser.add_argument('--sex', '--pheno_sex', type=str, default="both_sexes", help='pheno_sex of phenotype form Uk Biobank')
    parser.add_argument('--coding', type=str, default="", help='additional coding of phenotype form Uk Biobank')
    parser.add_argument('--index', type=str, required=True, help='path to Index file from PAN UKBiobank. It can be downloaded using gbe-get')
    parser.add_argument('--output-directory', type=str, default='', help='output directory')
    parser.add_argument('--variant-metrics', type=str, required=True, help='path to annotation file. It can be downloaded from https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz')
    parser.add_argument('--threshold', type=float, default=1e-08, help='significance cut-off threshold')
    parser.add_argument('--population', type=str, default='EUR', help='population: meta, AFR, AMR, CSA, EAS, EUR, MID')
    parser.add_argument('--clumping-vcf', type=str, default='eur.phase3.biobank.set.vcf.gz', help='')
    parser.add_argument('--source-ref-vcf', type=str, default='', help='')
    parser.add_argument('--target-ref-vcf', type=str, default='', help='')
    parsed_args = parser.parse_args(args)
    return parsed_args

def run(args):
    pass

def main(args = sys.argv[1:]):

    args = parse_args(args) 
    setup_logger(args.log_file) if args.log_file else setup_logger(args.output_directory + "/pgstk.log")

    try:
        run(args)
    except PolygenicException as e:
        error_exit(e)

if __name__ == '__main__':
    main(sys.argv[1:])