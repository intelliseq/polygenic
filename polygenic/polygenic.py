import argparse
import configparser
import glob
import importlib
import json
import logging
import os
import sys

from typing import Dict
from typing import List
from typing import Union
from typing import Iterable

MODULE_PATH = os.path.abspath(__file__).rsplit(os.path.sep, 4)[0]
sys.path.insert(0, MODULE_PATH)

from polygenic.version import __version__

from polygenic.data.data_accessor import DataAccessor
from polygenic.data.vcf_accessor import VcfAccessor
from polygenic.core.model import Model, SeqqlOperator

logger = logging.getLogger('polygenic')

def expand_path(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path)) if path else ''

def setup_logger(path):
    log_directory = os.path.dirname(os.path.abspath(os.path.expanduser(path)))
    if log_directory:
        try:
            os.makedirs(log_directory)
        except OSError:
            pass
    logger.setLevel(logging.DEBUG)
    logging_file_handler = logging.FileHandler(path)
    logging_file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logging_file_handler.setFormatter(formatter)
    logger.addHandler(logging_file_handler)

def main(args = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='')  # todo dodać opis
    parser.add_argument('-i', '--vcf', required=True, help='vcf.gz file with genotypes')
    parser.add_argument('-m', '--model', action='append', help="path to model (can be specified multiple times)")
    parser.add_argument('--parameters', type=str, help="parameters json (to be used in formula models)")
    parser.add_argument('-s', '--sample-name', type=str, help='sample name in vcf.gz to calculate')
    parser.add_argument('-o', '--output-directory', type=str, default='', help='output directory')
    parser.add_argument('-n', '--output-name-appendix', type=str, default='', help='appendix for output file names')
    parser.add_argument('-l', '--log-file', type=str, default='/var/logs/polygenic.log', help='path to log file')
    parser.add_argument('--af', type=str, default='', help='vcf file containing allele freq data')
    parser.add_argument('--af-field', type=str, default='AF',help='name of the INFO field to be used as allele frequency')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    parsed_args = parser.parse_args(args)
    setup_logger(parsed_args.log_file)

    out_dir = expand_path(parsed_args.output_directory)

    models = {}
    for model in parsed_args.model:
        model_path = expand_path(model)
        model_name = os.path.basename(model_path)
        model_info = {"path": model_path, "name": model_name}
        models[model_info["path"]] = model_info
    
    if not model_info:
        raise RuntimeError("No models loaded. Exiting.")

    vcf_accessor = VcfAccessor(expand_path(parsed_args.vcf))

    if parsed_args.af == "":
        allele_accessor = None
    else:    
        allele_accessor = VcfAccessor(expand_path(parsed_args.af))
    sample_names = vcf_accessor.get_sample_names()

    if "sample_name" in parsed_args and not parsed_args.sample_name is None:
        sample_names = [parsed_args.sample_name]

    parameters = {}
    if "parameters" in parsed_args and not parsed_args.parameters is None:
        with open(parsed_args.parameters) as parameters_json:
            parameters = json.load(parameters_json)

    appendix = parsed_args.output_name_appendix
    if appendix != "": 
        appendix = "-" + appendix

    for sample_name in sample_names:
        results_representations = {}
        for model_path, model_desc in models.items():
            if ".yml" in model_path:
                data_accessor = DataAccessor(
                    genotypes = vcf_accessor,
                    imputed_genotypes = vcf_accessor,
                    allele_frequencies =  allele_accessor,
                    sample_name = sample_name,
                    af_field_name = "AF_nfe",
                    parameters = parameters)
                model = SeqqlOperator.fromYaml(model_path)
                model.compute(data_accessor)
            with open(os.path.join(out_dir, f'{sample_name}-{model_desc["name"]}{appendix}.sample.json'), 'w') as f:
                json.dump(model.refine_results(), f, indent=2)
                #json.dump(results_representations, f, indent=4)



if __name__ == '__main__':
    main(sys.argv[1:])
