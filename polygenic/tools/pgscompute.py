import argparse
import sys
import os

from json import load as json_load
from json import dump as json_dump
from polygenic.version import __version__ as version
from polygenic.tools.utils import error_exit
from polygenic.tools.utils import setup_logger
from polygenic.tools.utils import expand_path
from datetime import datetime

from polygenic.data.data_accessor import DataAccessor
from polygenic.data.vcf_accessor import VcfAccessor
from polygenic.model.model import Model, SeqqlOperator
from polygenic.error.polygenic_exception import PolygenicException

def parse_args(args):
    parser = argparse.ArgumentParser(description='pgs-compute computes polygenic scores for genotyped sample in vcf format')
    parser.add_argument('-i', '--vcf', required=True, help='vcf.gz file with genotypes')
    parser.add_argument('-m', '--model', nargs='+', help="path to .yml model (can be specified multiple times with space as separator)")
    parser.add_argument('-p', '--parameters', type=str, help="parameters json (to be used in formula models)")
    parser.add_argument('-s', '--sample-name', type=str, help='sample name in vcf.gz to calculate')
    parser.add_argument('-o', '--output-directory', type=str, default='.', help='output directory')
    parser.add_argument('-n', '--output-name-appendix', type=str, help='appendix for output file names')
    parser.add_argument('-l', '--log-file', type=str, help='path to log file')
    parser.add_argument('--af', type=str, help='vcf file containing allele freq data')
    parser.add_argument('--af-field', type=str, default='AF',help='name of the INFO field to be used as allele frequency')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    parser.add_argument('--print', default=False, action='store_true', help='Print output to stdout')
    parsed_args = parser.parse_args(args)
    return parsed_args

def run(args):

    ### get models
    models = {}
    for model in args.model:
        model_path = expand_path(model)
        model_name = os.path.basename(model_path)
        model_info = {"path": model_path, "name": model_name}
        models[model_info["path"]] = model_info
    
    if not model_info:
        raise PolygenicError("No models were defined. Exiting.")

    ### get sample data
    vcf_accessor = VcfAccessor(expand_path(args.vcf))
    sample_names = vcf_accessor.get_sample_names()
    if "sample_name" in args and not args.sample_name is None:
        sample_names = [args.sample_name]

    ### get reference data 
    af_accessor = VcfAccessor(expand_path(args.af)) if args.af else None

    ### get parameters
    parameters = {}
    if "parameters" in args and not args.parameters is None:
        with open(args.parameters) as parameters_json:
            parameters = json_load(parameters_json)

    for sample_name in sample_names:
        results_representations = {}
        for model_path, model_desc in models.items():
            if ".yml" in model_path:
                data_accessor = DataAccessor(
                    genotypes = vcf_accessor,
                    allele_frequencies =  af_accessor,
                    sample_name = sample_name,
                    model_name = model_desc["name"],
                    af_field_name = args.af_field,
                    parameters = parameters)
                model = SeqqlOperator.fromYaml(model_path)
                model.compute(data_accessor)
                model.result["description"]["sample_name"] = sample_name
                model.result["description"]["model_name"] = model_desc["name"].replace('.yml', '')
            # output file name 
            appendix = "-" + args.output_name_appendix if args.output_name_appendix else ""
            output_path = os.path.join(expand_path(args.output_directory), f'{sample_name}-{model_desc["name"]}{appendix}-result.json')
            
            with open(output_path, 'w') as f:
                json_dump(model.refine_results(), f, indent=2)
            if args.print:
                json_dump(model.refine_results(), sys.stdout, indent=2)


def main(args = sys.argv[1:]):

    args = parse_args(args) 
    setup_logger(args.log_file) if args.log_file else setup_logger(args.output_directory + "/pgstk.log")

    try:
        run(args)
    except PolygenicException as e:
        error_exit(e)

if __name__ == '__main__':
    main(sys.argv[1:])
