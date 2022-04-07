import argparse
import os
import logging
import sys

import polygenic.tools as tools
from polygenic.version import __version__ as version
from polygenic.error.polygenic_exception import PolygenicException

def main(args=sys.argv[1:]):

    tools.utils.init_logger()

    parser = argparse.ArgumentParser(description='pgstk - the polygenic score toolkit')
    parser.add_argument('--log-level', type=str, default='INFO', help='logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default: INFO)')
    parser.add_argument('--log-stdout', action='store_true', default=False, help='log to stdout (default: False)')
    parser.add_argument('--log-file', type=str, default='~/.pgstk/log/pgstk.log', help='path to log file (default: $HOME/.pgstk/log/pgstk.log)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    subparsers = parser.add_subparsers(dest = 'tool')

    ### utils ###
    # vcf-index
    vcf_index_parser = subparsers.add_parser('vcf-index', description='vcf-index creates index for vcf file')
    vcf_index_parser.add_argument('-i', '--vcf', required=True, help='path to vcf.gz')

    ### plots ###
    # plot-manhattan
    plot_manhattan_parser = subparsers.add_parser('plot-manhattan', description='plot-manhattan draws manhattan plot')
    plot_manhattan_parser.add_argument('-i', '--tsv', required=True, help='tsv or tsv.gz file with gwas data')
    plot_manhattan_parser.add_argument('-g', '--genome-version', default="GRCh38", help="genome version GRCh37 or GRCh38 (default: GRCh38)")
    plot_manhattan_parser.add_argument('-c', '--chromosome-column', default="chr", help="column name for chromosome (default: chr)")
    plot_manhattan_parser.add_argument('-s', '--position-column', default="pos", help="column name for position (default: pos)")
    plot_manhattan_parser.add_argument('-p', '--pvalue-column', default="pval_meta", help="column name for pvalue (default: pos)")
    plot_manhattan_parser.add_argument('-f', '--format', default="png", help="output format {png, pdf} (default: png)")
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
        if parsed_args.tool == 'vcf-index':
            tools.vcfindex.run(parsed_args)
        elif parsed_args.tool == 'plot-manhattan':
            tools.plotmanhattan.run(parsed_args)
    except PolygenicException as e:
        tools.utils.error_exit(e)
    except RuntimeError as e:
        tools.utils.error_exit(e)

    
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

if __name__ == '__main__':
    main(sys.argv[1:])
