"""
pgstk command line tool for polygenic trait analysis
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
    """pgstk command line tool wrapper

    Args:
        args (list[str], optional): Arguments list. Defaults to None.

    Returns:
        int: exit code
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
    pgs_compute_parser.add_argument('--merge-as-array', default=False, action='store_true', help='combine outputs for multiple models into one array (default: False). Works only with --merge-outputs')
    pgs_compute_parser.add_argument('-p', '--parameters', type=str, help="parameters json (to be used in formula models)")
    pgs_compute_parser.add_argument('-s', '--sample-name', type=str, help='sample name in vcf.gz to calculate')
    pgs_compute_parser.add_argument('-o', '--output-directory', type=str, default='.', help='output directory (default: .)')
    pgs_compute_parser.add_argument('-n', '--output-name-appendix', type=str, help='appendix for output file names')
    pgs_compute_parser.add_argument('--af', type=str, help='vcf file containing allele freq data')
    pgs_compute_parser.add_argument('--af-field', type=str, default='AF',help='name of the INFO field to be used as allele frequency')
    pgs_compute_parser.add_argument('--print', default=False, action='store_true', help='print output to stdout')

    ### gwas file ###
    # gwas-file-create
    model_gwas_file_parser = subparsers.add_parser('gwas-file-create', description='model-gwas-file builds a model from a gwas results')
    model_gwas_file_parser.add_argument('-i', '--input', required=True, help='gwas results file')
    model_gwas_file_parser.add_argument('-k', '--keyword', required=False, help='keyword to be used for column names selection priority')
    model_gwas_file_parser.add_argument('-o', '--output', required=True, help='output file')


    ### build model ###
    # model-gwas-file
    model_gwas_file_parser = subparsers.add_parser('model-gwas-file', description='model-gwas-file builds a model from a gwas results')
    model_gwas_file_parser.add_argument('-i', '--gwas-file', required=True, help='gwas results file')
    model_gwas_file_parser.add_argument('-c', '--chromosome-column-name', type=str, default='CHROM', help='name of the column containing chromosome')
    model_gwas_file_parser.add_argument('-s', '--position-column-name', type=str, default='POS', help='name of the column containing position')
    model_gwas_file_parser.add_argument('-r', '--ref-allele-column-name', type=str, default='REF', help='name of the column containing reference allele')
    model_gwas_file_parser.add_argument('-a', '--alt-allele-column-name', type=str, default='ALT', help='name of the column containing alternate allele')
    model_gwas_file_parser.add_argument('-e', '--effect-allele-column-name', type=str, default='EFFECT', help='name of the column containing effect allele')
    model_gwas_file_parser.add_argument('-p', '--pvalue-column-name', type=str, default='PVALUE', help='name of the column containing p-values')
    model_gwas_file_parser.add_argument('-b', '--beta-column-name', type=str, default='BETA', help='name of the column containing beta')
    model_gwas_file_parser.add_argument('-d', '--rsid-column-name', type=str, default='RSID', help='name of the column containing rsid')
    model_gwas_file_parser.add_argument('--header', type=str, help='header file')
    model_gwas_file_parser.add_argument('-o', '--output', required=True, help='output file')
    model_gwas_file_parser.add_argument('--print', default=False, action='store_true', help='print output to stdout')

    # model-biobankuk
    model_biobankuk_parser = subparsers.add_parser('model-biobankuk', description='model-biobankuk builds a model from a biobankuk results')
    model_biobankuk_parser.add_argument('--code', '--phenocode', type=str, required=True, help='phenocode of phenotype form Uk Biobank')
    model_biobankuk_parser.add_argument('--sex', '--pheno_sex', type=str, default="both_sexes", help='pheno_sex of phenotype form Uk Biobank')
    model_biobankuk_parser.add_argument('--coding', type=str, default="", help='additional coding of phenotype form Uk Biobank')
    model_biobankuk_parser.add_argument('--output-directory', type=str, default='.', help='output directory')
    model_biobankuk_parser.add_argument('--index-file', type=str, help='path to Index file from PAN UKBiobank. It can be downloaded using gbe-get')
    model_biobankuk_parser.add_argument('--variant-metrics-file', type=str, help='path to annotation file. It can be downloaded from https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz')
    model_biobankuk_parser.add_argument('--index-url', type=str, default='https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz', help='url of index file for PAN UKBiobank.')
    model_biobankuk_parser.add_argument('--variant-metrics-url', type=str, default='https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz', help='url for variant summary metrics')
    model_biobankuk_parser.add_argument('--pvalue-threshold', type=float, default=1e-08, help='significance cut-off threshold. e.g. 1e-08')
    model_biobankuk_parser.add_argument('--clump-r2', type=float, default=0.25, help='clumping r2 threshold')
    model_biobankuk_parser.add_argument('--clump-kb', type=float, default=1000, help='clumping kb threshold')
    model_biobankuk_parser.add_argument('--population', type=str, default='EUR', help='population: meta, AFR, AMR, CSA, EUR, EAS, EUR, MID')
    model_biobankuk_parser.add_argument('--clumping-vcf', type=str, default='eur.phase3.biobank.set.vcf.gz', help='')
    model_biobankuk_parser.add_argument('--source-ref-vcf', type=str, default='dbsnp155.grch37.norm.vcf.gz', help='')
    model_biobankuk_parser.add_argument('--target-ref-vcf', type=str, default='dbsnp155.grch38.norm.vcf.gz', help='')
    model_biobankuk_parser.add_argument('--gene-positions', type=str, default='ensembl-genes.104.tsv', help='table with ensembl genes')
    model_biobankuk_parser.add_argument('--ignore-warnings', type=bool, default='False', help='')
    model_biobankuk_parser.add_argument('--test', default=False, action='store_true', help='print output to stdout')

    ### utils ###
    # vcf-index
    vcf_index_parser = subparsers.add_parser('vcf-index', description='vcf-index creates index for vcf file')
    vcf_index_parser.add_argument('-i', '--vcf', required=True, help='path to vcf.gz')

    # vcf-stat-baf
    vcf_stat_baf_parser = subparsers.add_parser('vcf-stat-baf', description='vcf-stat creates index for vcf file')
    vcf_stat_baf_parser.add_argument('-i', '--vcf', required=True, help='path to vcf.gz')
    vcf_stat_baf_parser.add_argument('-o', '--output-directory', type=str, default='.', help='output directory (default: .)')

    # vcf-stat-zygosity
    vcf_stat_zygosity_parser = subparsers.add_parser('vcf-stat-zygosity', description='vcf-stat creates index for vcf zygosity')
    vcf_stat_zygosity_parser.add_argument('-i', '--vcf', required=True, help='path to vcf.gz')
    vcf_stat_zygosity_parser.add_argument('-o', '--output', type=str, default='./zygosity.json', help='output directory (default: .)')
    
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
        tool = getattr(tools, parsed_args.tool.replace('-',''))
        tool.run(parsed_args)
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
