import argparse
import logging

import sys
import os
import urllib.request

import configparser
import tabix

import yaml 

from polygenic.tools import pgscompute
from polygenic.tools import modelbiobankuk
from polygenic.tools import modelpgscat

# utils
# simulate
from polygenic.data.vcf_accessor import VcfAccessor
from polygenic.data.vcf_accessor import DataNotPresentError
from polygenic.model.utils import is_valid_path
from polygenic.model.utils import download
from polygenic.model.utils import read_header
from polygenic.model.utils import read_table
from polygenic.error.polygenic_exception import PolygenicException
# clumping
import subprocess
import re
# simlating
import random
import statistics
# saving
import json

logger = logging.getLogger('polygenicmaker')
config = configparser.ConfigParser()
config.read(os.path.dirname(__file__) + "/../polygenic/polygenic.cfg")

#######################
### vcf-index #########
#######################


def vcf_index(args):
    parser = argparse.ArgumentParser(
        description='vcf index indexes vcf file for query')  # todo dodać opis
    parser.add_argument('--vcf', type=str, required=True,
                        help='path to vcf file')
    parsed_args = parser.parse_args(args)
    VcfAccessor(parsed_args.vcf)
    return

#####################################################################################################
###                                                                                               ###
###                                   Global Biobank Engine                                       ###
###                                                                                               ###
#####################################################################################################

#######################
### gbe-index #########
#######################

def gbe_index(args):
    print("GBEINDEX")
    parser = argparse.ArgumentParser(
        description='polygenicmaker gbe-index downloads index of gwas results from pan.ukbb study')  # todo dodać opis
    parser.add_argument('--url', type=str, default='https://biobankengine.stanford.edu/static/degas-risk/degas_n_977_traits.tsv',
                        help='alternative url location for index')
    parser.add_argument('--output-directory', type=str, default='',
                        help='output directory')
    parsed_args = parser.parse_args(args)
    output_path = os.path.abspath(os.path.expanduser(
        parsed_args.output_directory)) + "/gbe_phenotype_manifest.tsv"
    download(parsed_args.url, output_path)
    return

###############
### gbe-get ###
###############

def gbe_get(parsed_args):
    url = "https://biobankengine.stanford.edu/static/PRS_map/" + parsed_args.code + ".tsv"
    output_directory = os.path.abspath(os.path.expanduser(parsed_args.output_directory))
    output_file_name = os.path.splitext(os.path.basename(url))[0]
    output_path = output_directory + "/" + output_file_name
    download(url=url, output_path=output_path,
             force=parsed_args.force, progress=True)
    return output_path

#######################
### gbe-model #########
#######################

def gbe_model(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker biobankuk-build-model constructs polygenic score model based on p value data')  # todo dodać opis
    parser.add_argument('-c','--code', required=True, type=str, help='path to PRS file from gbe. It can be downloaded using gbe-get')
    parser.add_argument('-o', '--output-directory', type=str, default='output', help='output directory')
    parser.add_argument('--gbe-index', type=str, default='gbe-index.1.3.1.tsv', help='gbe-index')
    parser.add_argument('--source-ref-vcf', type=str, default='dbsnp155.grch37.norm.vcf.gz', help='source reference vcf (hg19)')
    parser.add_argument('--target-ref-vcf', type=str, default='dbsnp155.grch38.norm.vcf.gz', help='source reference vcf (hg38)')
    parser.add_argument('--af-vcf', type=str, default='gnomad.3.1.vcf.gz', help='path to allele frequency vcf.')
    parser.add_argument('--af-field', type=str, default='AF_nfe', help='name of the INFO field with ALT allele frequency')
    parser.add_argument('--gene-positions', type=str, default='ensembl-genes.104.tsv', help='table with ensembl genes')
    parser.add_argument('-i', '--iterations', type=float, default=1000, help='simulation iterations for mean and sd')
    parser.add_argument('-f', '--force', action='store_true', help='overwrite downloaded file')
    parsed_args = parser.parse_args(args)

    if not is_valid_path(parsed_args.output_directory, is_directory=True): return
    path = gbe_get(parsed_args)

    if not is_valid_path(path): return
    data = read_table(path)

    if not is_valid_path(parsed_args.gbe_index): return
    gbe_index = read_table(parsed_args.gbe_index)
    info = [line for line in gbe_index if parsed_args.code in line['Trait name']][0]

    if not is_valid_path(parsed_args.gene_positions): return
    gene_positions = read_table(parsed_args.gene_positions)

    if not is_valid_path(parsed_args.af_vcf, possible_url = True): return
    af_vcf = VcfAccessor(parsed_args.af_vcf)

    if not is_valid_path(parsed_args.source_ref_vcf, possible_url = True): return
    source_vcf = VcfAccessor(parsed_args.source_ref_vcf)

    if not is_valid_path(parsed_args.target_ref_vcf, possible_url = True): return
    target_vcf = VcfAccessor(parsed_args.target_ref_vcf)

    data = [line for line in data if "rs" in line['ID']]
    for line in data: line.update({"rsid": line['ID']})
    data = [validate(line, validation_source = target_vcf) for line in data]

    for line in data:
        if not "rsid" in line:
            print(line)

    data = [line for line in data if "rs" in line['ID']]

    data = [add_annotation(
        line, 
        annotation_name = "af", 
        annotation_source = af_vcf, 
        annotation_source_field = parsed_args.af_field,
        default_value = 0.001) for line in data]

    symbols = [add_gene_symbol(
        line, 
        target_vcf, 
        gene_positions) for line in data]
    symbols = [symbol for symbol in symbols if symbol]

    description = dict()
    parameters = simulate_parameters(data)
    description["pmids"] = [33095761]
    description["info"] = info
    description["parameters"] = parameters
    description["genes"] = symbols

    model_path = parsed_args.output_directory + "/" + parsed_args.code + ".yml"
    write_model(data, description, model_path)
    return




def main(args=sys.argv[1:]):
    try:
        if args[0] == 'pgs-compute':
            pgscompute.main(args[1:])
        elif args[0] == 'model-biobankuk':
            modelbiobankuk.main(args[1:])
        elif args[0] == 'model-pgscat':
            modelpgscat.main(args[1:])
        elif args[0] == 'vcf-index':
            vcf_index(args[1:])
        else:
            print('ERROR: Please select proper tool name"')
            print("""
            Program: polygenicmaker (downloads gwas data, clumps and build polygenic scores)
            Contact: Marcin Piechota <piechota@intelliseq.com>
            Usage:   pgstk <command> [options]

            Command:
            pgs-compute             computes pgs score for vcf file
            model-biobankuk         prepare polygenic score model based on gwas results from biobankuk
            model-pgscat            prepare polygenic score model based on gwas results from PGS Catalogue
            vcf-index               prepare rsidx for vcf

            """)
    except PolygenicException as e:
        print("Analysis failed")
        print("ERROR: " + str(e))
    except RuntimeError as e:
        print("ERROR: " + str(e))

if __name__ == '__main__':
    main(sys.argv[1:])
