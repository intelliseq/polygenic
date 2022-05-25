import sys
import os

from json import load as json_load
from json import dump as json_dump
from polygenic.version import __version__ as version
from polygenic.tools.utils import expand_path

from polygenic.data.data_accessor import DataAccessor
from polygenic.data.vcf_accessor import VcfAccessor
from polygenic.model.model import Model, SeqqlOperator
from polygenic.error.polygenic_exception import PolygenicException

def run(args):

    ### get models
    models = {}
    for model in args.model:
        model_path = expand_path(model)
        model_name = os.path.basename(model_path)
        model_info = {"path": model_path, "name": model_name}
        models[model_info["path"]] = model_info
    
    if not model_info:
        raise PolygenicException("No models were defined. Exiting.")

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
                if "description" not in model.result:
                    model.result["description"] = {}
                model.result["description"]["sample_name"] = sample_name
                model.result["description"]["model_name"] = model_desc["name"].replace('.yml', '')
            # output file name 
            appendix = "-" + args.output_name_appendix if args.output_name_appendix else ""
            output_path = os.path.join(expand_path(args.output_directory), f'{sample_name}-{model_desc["name"]}{appendix}-result.json').replace('.yml', '')
            if args.merge_outputs:
                results_representations[model_desc["name"].replace('.yml', '')] = model.refine_results()
            else:
                with open(output_path, 'w') as f:
                    json_dump(model.refine_results(), f, indent=2)
                if args.print:
                    json_dump(model.refine_results(), sys.stdout, indent=2)
        if args.merge_outputs:
            output_path = os.path.join(expand_path(args.output_directory), f'{sample_name}{appendix}-result.json')
            with open(output_path, 'w') as f:
                json_dump(results_representations, f, indent=2)
            if args.print:
                json_dump(results_representations, sys.stdout, indent=2)