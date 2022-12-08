from polygenic.data.gwas import Gwas
from calendar import c
import sys
import os

import polygenic.tools.utils as utils
import polygenic.data.csv_accessor as csv_accessor

#    model_gwas_file_parser = subparsers.add_parser('model-gwas-file', description='model-gwas-file builds a model from a gwas results')
#    model_gwas_file_parser.add_argument('-i', '--gwas-file', required=True, help='gwas results file')
#    model_gwas_file_parser.add_argument('-c', '--chromosome-column-name', type=str, default='CHROM', help='name of the column containing chromosome')
#    model_gwas_file_parser.add_argument('-s', '--position-column-name', type=str, default='POS', help='name of the column containing position')
#    model_gwas_file_parser.add_argument('-r', '--ref-allele-column-name', type=str, default='REF', help='name of the column containing reference allele')
#    model_gwas_file_parser.add_argument('-a', '--alt-allele-column-name', type=str, default='ALT', help='name of the column containing alternate allele')
#    model_gwas_file_parser.add_argument('-e', '--effect-allele-column-name', type=str, default='EFFECT', help='name of the column containing effect allele')
#    model_gwas_file_parser.add_argument('-p', '--pvalue-column-name', type=str, default='PVALUE', help='name of the column containing p-values')
#    model_gwas_file_parser.add_argument('-b', '--beta-column-name', type=str, default='BETA', help='name of the column containing beta')
#    model_gwas_file_parser.add_argument('-d', '--rsid-column-name', type=str, default='RSID', help='name of the column containing rsid')
#    model_gwas_file_parser.add_argument('-h', '--header', type=str, help='header file')
#    model_gwas_file_parser.add_argument('-o', '--output', required=True, help='output file')
#    model_gwas_file_parser.add_argument('--print', default=False, action='store_true', help='print output to stdout')


def run(args):
    gwas = Gwas()
    gwas.load_gwas_from_csv(args.gwas_file, vars(args))
    gwas.clump()
    gwas.validate()
    gwas.plot_manhattan()

    # dbsnp37 = vcf_accessor.VcfAccessor("/home/marpiech/data/vcf/dbsnp155.grch37.norm.vcf.gz")
    # dbsnp38 = vcf_accessor.VcfAccessor("/home/marpiech/data/vcf/dbsnp155.grch38.norm.vcf.gz")

    # data = csv.get_data()

    # print(str(data.iloc[0:2, :]))
