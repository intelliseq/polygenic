"""
toolkit for polygenic trait analysis
"""

import argparse
import os
import logging
import sys
from datetime import datetime 

import polygenic.tools as tools
from polygenic.version import __version__ as version
from polygenic.error.polygenic_exception import PolygenicException

def main(args=None):
    """
    pgstk argument parser
    """

    parser = argparse.ArgumentParser(description='pgstk - the polygenic score toolkit')
    parser.add_argument('--log-level', type=str, default='INFO', help='logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default: INFO)')
    parser.add_argument('--log-stdout', action='store_true', default=False, help='log to stdout (default: False)')
    parser.add_argument('--log-file', type=str, default='~/.pgstk/log/pgstk.log', help='path to log file (default: $HOME/.pgstk/log/pgstk.log)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    subparsers = parser.add_subparsers(dest = 'tool')

    ### compute ###
    # pgs-compute
    pgs_compute_parser = subparsers.add_parser('pgs-compute', description='pgs-compute computes polygenic scores for genotyped sample in vcf format')
    pgs_compute_parser.add_argument('-i', '--vcf', required=True, help='vcf.gz file with genotypes')
    pgs_compute_parser.add_argument('-m', '--model', nargs='+', help='path to .yml model (can be specified multiple times with space as separator)')
    pgs_compute_parser.add_argument('--merge-outputs', default=False, action='store_true', help='combine outputs for multiple models into one file (default: False)')
    pgs_compute_parser.add_argument('-p', '--parameters', type=str, help="parameters json (to be used in formula models)")
    pgs_compute_parser.add_argument('-s', '--sample-name', type=str, help='sample name in vcf.gz to calculate')
    pgs_compute_parser.add_argument('-o', '--output-directory', type=str, default='.', help='output directory (default: .)')
    pgs_compute_parser.add_argument('-n', '--output-name-appendix', type=str, help='appendix for output file names')
    pgs_compute_parser.add_argument('--af', type=str, help='vcf file containing allele freq data')
    pgs_compute_parser.add_argument('--af-field', type=str, default='AF',help='name of the INFO field to be used as allele frequency')
    pgs_compute_parser.add_argument('--print', default=False, action='store_true', help='print output to stdout')

    ### build model ###
    # model-gwas-file
    model_gwas_file_parser = subparsers.add_parser('model-gwas-file', description='model-gwas-file builds a model from a gwas results')
    model_gwas_file_parser.add_argument('-i', '--gwas-file', required=True, help='gwas results file')
    model_gwas_file_parser.add_argument('-p', '--pvalue-column-name', type=str, default='PVALUE', help='name of the column containing p-values')
    model_gwas_file_parser.add_argument('-c', '--chromosome-column-name', type=str, default='CHROM', help='name of the column containing chromosome')
    model_gwas_file_parser.add_argument('-s', '--position-column-name', type=str, default='POS', help='name of the column containing position')
    model_gwas_file_parser.add_argument('-r', '--ref-allele-column-name', type=str, default='REF', help='name of the column containing reference allele')
    model_gwas_file_parser.add_argument('-a', '--alt-allele-column-name', type=str, default='ALT', help='name of the column containing alternate allele')
    model_gwas_file_parser.add_argument('-e', '--effect-allele-column-name', type=str, default='EFFECT', help='name of the column containing effect allele')
    model_gwas_file_parser.add_argument('-b', '--beta-column-name', type=str, default='BETA', help='name of the column containing beta')
    model_gwas_file_parser.add_argument('-d', '--rsid-column-name', type=str, default='RSID', help='name of the column containing rsid')
    model_gwas_file_parser.add_argument('-o', '--output', required=True, help='output file')
    model_gwas_file_parser.add_argument('--print', default=False, action='store_true', help='print output to stdout')

    ### utils ###
    # vcf-index
    vcf_index_parser = subparsers.add_parser('vcf-index', description='vcf-index creates index for vcf file')
    vcf_index_parser.add_argument('-i', '--vcf', required=True, help='path to vcf.gz')

    ### plots ###
    # plot-manhattan
    plot_manhattan_parser = subparsers.add_parser('plot-manhattan',
        description='plot-manhattan draws manhattan plot')
    plot_manhattan_parser.add_argument('-i', '--tsv', 
        required=True, help='tsv or tsv.gz file with gwas data')
    plot_manhattan_parser.add_argument('-d', '--delimiter', default='\t', help="tsv delimiter (default: '\\t')")
    plot_manhattan_parser.add_argument('-g', '--genome-version', default="GRCh38", choices=['GRCh37', 'GRCh38'], help="genome version GRCh37 or GRCh38 (default: GRCh38)")
    plot_manhattan_parser.add_argument('-c', '--chromosome-column', default="chr", help="column name for chromosome (default: chr)")
    plot_manhattan_parser.add_argument('-s', '--position-column', default="pos", help="column name for position (default: pos)")
    plot_manhattan_parser.add_argument('-p', '--pvalue-column', default="pval_meta", help="column name for pvalue (default: pos)")
    plot_manhattan_parser.add_argument('-f', '--output-format', default="pdf", help="output format {png, pdf} (default: png)")
    plot_manhattan_parser.add_argument('-o', '--output', help="output (default: {tsv}.{format}})")

    parsed_args = parser.parse_args(args)

    # configure logging
    logger = logging.getLogger()
    # set logger level based on argaprse
    logger.setLevel(parsed_args.log_level)
    #set logger format
    formatter = logging.Formatter("%(asctime)s [%(levelname)8s] %(message)s (%(filename)s:%(lineno)s)")    
    # get handlers for logger
    handlers = logger.handlers
    # remove all handlers
    for handler in handlers:
        logger.removeHandler(handler)
    # add file handler
    if parsed_args.log_file:
        path = os.path.abspath(os.path.expanduser(parsed_args.log_file))
        logging.info(path)
        log_directory = os.path.dirname(path)
        if log_directory and not os.path.exists(log_directory): os.makedirs(log_directory)
        file_handler = logging.FileHandler(path)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    # add stdout handler
    if parsed_args.log_stdout:
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    logging.debug("running %s", parsed_args.tool)

    try:
        if parsed_args.tool == 'pgs-compute':
            tools.pgscompute.run(parsed_args)
        elif parsed_args.tool == 'model-gwas-file':
            tools.modelgwasfile.run(parsed_args)
        elif parsed_args.tool == 'vcf-index':
            tools.vcfindex.run(parsed_args)
        elif parsed_args.tool == 'plot-manhattan':
            tools.plotmanhattan.run(parsed_args)
    except PolygenicException as exception:
        error_exit(exception)
    except RuntimeError as exception:
        error_exit(exception)
    return 0

def error_print(*args, **kwargs):
    """
    Prints error message to stderr
    """
    print(*args, file=sys.stderr, **kwargs)

def error_exit(exception):
    """
    Prints error message to stderr and exits with error code 1
    """
    time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    error_print("")
    error_print("  polygenic ERROR ")
    error_print("  version: " + version)
    error_print("  time: " + time)
    error_print("  command: pgstk " + (" ").join(sys.argv))
    error_print("  message: ")
    error_print("")
    error_print("  " + str(exception))
    error_print("")
    exit(1)

if __name__ == '__main__':
    main(sys.argv[1:])
