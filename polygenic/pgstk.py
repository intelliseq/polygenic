import argparse
import os
import logging
import sys
from datetime import datetime 

import polygenic.tools as tools
from polygenic.version import __version__ as version
from polygenic.error.polygenic_exception import PolygenicException

def main(args=sys.argv[1:]):

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

    ### utils ###
    # vcf-index
    vcf_index_parser = subparsers.add_parser('vcf-index', description='vcf-index creates index for vcf file')
    vcf_index_parser.add_argument('-i', '--vcf', required=True, help='path to vcf.gz')

    ### plots ###
    # plot-manhattan
    plot_manhattan_parser = subparsers.add_parser('plot-manhattan', description='plot-manhattan draws manhattan plot')
    plot_manhattan_parser.add_argument('-i', '--tsv', required=True, help='tsv or tsv.gz file with gwas data')
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

    logging.debug("running " + parsed_args.tool)

    try:    
        if parsed_args.tool == 'pgs-compute':
            tools.pgscompute.run(parsed_args)
        elif parsed_args.tool == 'vcf-index':
            tools.vcfindex.run(parsed_args)
        elif parsed_args.tool == 'plot-manhattan':
            tools.plotmanhattan.run(parsed_args)
    except PolygenicException as e:
        error_exit(e)
    except RuntimeError as e:
        error_exit(e)

    
    # args[0] == 'pgs-compute':
    #         tools.pgscompute.main(args[1:])
    #     elif args[0] == 'model-biobankuk':
    #         tools.modelbiobankuk.main(args[1:])
    #     elif args[0] == 'model-pgscat':
    #         tools.modelpgscat.main(args[1:])
    #     elif args[0] == 'vcf-index':
    #         tools.vcfindex.main(args[1:])
    #     else:
    # #         print('ERROR: Please select proper tool name"')
    #         print("""
    #         Program: polygenic toolkit (downloads gwas data, builds and computes polygenic scores)
    #         Contact: Marcin Piechota <piechota@intelliseq.com>
    #         Usage:   pgstk <command> [options]

    #         Commands:
    #         pgs-compute             computes pgs score for vcf file
    #         model-biobankuk         prepare polygenic score model based on gwas results from biobankuk
    #         model-pgscat            prepare polygenic score model based on gwas results from PGS Catalogue
    #         vcf-index               prepare rsidx for vcf
    #         plot-gwas               manhattan plot for gwas results
    #         """)
            


    return 0

def error_print(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def error_exit(e):
    time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    error_print("")
    error_print("  polygenic ERROR ")
    error_print("  version: " + version)
    error_print("  time: " + time)
    error_print("  command: pgstk " + (" ").join(sys.argv))
    error_print("  message: ")
    error_print("")
    error_print("  " + str(e))
    error_print("")
    exit(1)

if __name__ == '__main__':
    main(sys.argv[1:])
