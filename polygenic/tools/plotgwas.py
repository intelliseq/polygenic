import argparse
import sys
import os

from polygenic.version import __version__ as version
from polygenic.tools.utils import error_exit
from polygenic.tools.utils import setup_logger

def parse_args(args):
    parser = argparse.ArgumentParser(description='plot-gwas plots manhattan plot')
    
    parser.add_argument('-i', '--input', required=True, help='tsv.gz file with gwas data')
    parser.add_argument('-g', '--genome-version', default="GRCh38", help="genome version GRCh37 or GRCh38 (default: GRCh38)")
    parser.add_argument('-c', '--chromosome-column', default="chr", help="column name for chromosome (default: chr)")
    parser.add_argument('-c', '--position-column', default="pos", help="column name for position (default: pos)")
    parser.add_argument('-c', '--pvalue-column', default="pval_meta", help="column name for position (default: pos)")

    parser.add_argument('-o', '--output-path', type=str, default='plot-gwas.svg', help='output path')
    parser.add_argument('-l', '--log-file', type=str, help='path to log file')
    parser.add_argument('--af', type=str, help='vcf file containing allele freq data')
    parser.add_argument('--af-field', type=str, default='AF',help='name of the INFO field to be used as allele frequency')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)

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
