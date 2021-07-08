import logging
import os
import configparser
import sys
import argparse
import glob
import importlib
import json

from typing import Dict
from typing import List
from typing import Union
from typing import Iterable

MODULE_PATH = os.path.abspath(__file__).rsplit(os.path.sep, 4)[0]
sys.path.insert(0, MODULE_PATH)

from polygenic.version import __version__
from polygenic.lib.data_access.data_accessor import VcfAccessor

from polygenic.lib.output import create_res_representation_for_model
from polygenic.lib.data_access.allele_frequency_accessor import AlleleFrequencyAccessor
from polygenic.lib.data_access.data_accessor import VcfAccessor
from polygenic.seqql.score import PolygenicRiskScore
from polygenic.seqql.score import Data
from polygenic.lib.data_access.dto import ModelDescriptionInfo


logger = logging.getLogger('polygenic')

class Sequery(object):

    def __init__(self, testConfigFile = os.path.dirname(__file__) + "/test_module_config.json", configFile = os.path.dirname(__file__) + "/wdltest.cfg", index = -1):
        self.logger = logging.getLogger(__name__)
        self.logger.info('Reading config')
        config = self.getConfig(configFile)
        testConfig = TestConfiguration(testConfigFile).getConfiguration()

        # Preparing cromwell
        self.logger.info('Preparing cromwell handler')
        self.cromwell = CromwellHandler(config)
        self.testRunner = TestRunner(testConfig, self.cromwell, index = index)

    def getConfig(self, configFile = os.path.dirname(__file__) + "/wdltest.cfg"):
        config = configparser.ConfigParser()
        config.read(configFile)
        return config

    def run(self):
        exitCode = self.testRunner.run()
        self.cromwell.stop()
        return exitCode

    def localrun(self):
        exitCode = self.testRunner.run()
        self.cromwell.stop()
        return exitCode

    def stop(self):
        self.cromwell.stop()

def expand_path(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path)) if path else ''

def load_models_when_mapping_absent(directories:Iterable[str], population:str, logger) -> Dict[str, ModelDescriptionInfo]:
    model_infos = {}
    for directory in directories:
        logger.info(f"Discovering models in {directory}")
        for model_path_ in glob.glob(os.path.join(directory, '*_{}_model.py'.format(population))):
            model_path = expand_path(model_path_)
            model_fname = os.path.basename(model_path)
            description_path = model_path.replace('.py', '.json')
            description = [description_path] if os.path.exists(description_path) else []
            plot_data_path_str = model_path.replace('_traits/', '_data/').replace('_model.py', '_data.json')
            plot_data_path = plot_data_path_str if os.path.exists(plot_data_path_str) else None
            model_infos[model_path] = ModelDescriptionInfo(
                model_fname = model_fname,
                model_path = model_path,
                desc_paths = description,
                plot_data_path = plot_data_path
            )
    return model_infos

def process_model(model_description_info:ModelDescriptionInfo, vcf_accessor:VcfAccessor, allele_freq_accessor:VcfAccessor, pop:str, sample_name:str):
    package_str = os.path.dirname(model_description_info.model_path)
    if package_str not in sys.path:
        sys.path.insert(0, package_str)
    module = importlib.import_module(model_description_info.model_fname.split('.')[0])
    pop_short = pop.replace('AF', '').lstrip('_')
    #if pop_short != module.trait_was_prepared_for_population:
    #    raise ImproperPopulationForModelError(f"You requested data for population {pop_short} while the model was prepared for {module.trait_was_prepared_for_population}")
    model: PolygenicRiskScore = module.model
    data = Data(vcf_accessor, allele_freq_accessor, sample_name, pop, model)
    return data.compute_model()

def parse_model(path: str):
    model_path = expand_path(path)
    model_fname = os.path.basename(model_path)
    description_path = model_path.replace('.py', '.json')
    description = [description_path] if os.path.exists(description_path) else []
    #plot_data_path_str = model_path.replace('_traits/', '_data/').replace('_model.py', '_data.json')
    #plot_data_path = plot_data_path_str if os.path.exists(plot_data_path_str) else None
    model_info = ModelDescriptionInfo(
        model_fname = model_fname,
        model_path = model_path,
        desc_paths = description,
        plot_data_path = None
    )
    return model_info

def main(args = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='')  # todo dodać opis
    parser.add_argument('--vcf', nargs='+', required=True,
                        help='Vcf file with genotypes')
    parser.add_argument('--log_file', type=str, default='polygenic.log')
    parser.add_argument('--out_dir', type=str, default="", help='Directory for result jsons.')
    parser.add_argument('--population', type=str, default='eas',
                        choices=['', 'nfe', 'eas', 'afr', 'amr', 'asj', 'fin', 'oth'],
                        help='''Population code:
        empty - use average allele frequency in all population,
        'nfe' - Non-Finnish European ancestry,
        'eas' - East Asian ancestry,
        'afr' - African-American/African ancestry,
        'amr' - Latino ancestry,
        'asj' - Ashkenazi Jewish ancestry, 
        'fin' - Finnish ancestry,
        'oth' - Other ancestry''')
    parser.add_argument('--traits_dirs', nargs='+', type=str, default=[],
                        help='Directories containing models from the inside of this repo')
    # parser.add_argument('--actionable_dir', type=str, default='',
    #                     help='Directory containing "actionable" models')  # default='vitalleo_actionable'
    parser.add_argument('--models_path', type=str, default='', help="Path to a directory containing models and corresponding_descriptions")
    parser.add_argument('--model', action='append', help="Path to model")
    parser.add_argument('--mapping_json', type=str, default='', help="A file containing mapping between models and descriptions")
    parser.add_argument('--af', type=str, default='',
                        help="A file containing allele freq data")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    parsed_args = parser.parse_args(args)

    ### setup logging ###
    log_directory = os.path.dirname(os.path.abspath(os.path.expanduser(parsed_args.log_file)))
    if log_directory:
        try:
            os.makedirs(log_directory)
        except OSError:
            pass
    logger.setLevel(logging.DEBUG)
    logging_file_handler = logging.FileHandler(parsed_args.log_file)
    logging_file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logging_file_handler.setFormatter(formatter)
    logger.addHandler(logging_file_handler)

    ###
    out_dir = expand_path(parsed_args.out_dir)

    population = 'AF' if not parsed_args.population else 'AF_' + parsed_args.population

    models_info = {}
    for model in parsed_args.model:
        model_info = parse_model(model)
        models_info[model_info.model_path] = model_info
    directory = str(os.path.join(parsed_args.models_path, '*_{}_model.py'.format(parsed_args.population)))
    for model in glob.glob(directory):
        model_info = parse_model(model)
        models_info[model_info.model_path] = model_info
    
    if not models_info:
        raise RuntimeError("No models loaded. Exiting.")

    vcf_accessor = VcfAccessor(expand_path(parsed_args.vcf[0]))

    if parsed_args.af == "":
        allele_accessor = None
    else:    
        allele_accessor = VcfAccessor(expand_path(parsed_args.af))#AlleleFrequencyAccessor(allele_freq_json_path=allele_freq_path)
    sample_names = vcf_accessor.sample_names
    
    for sample_name in sample_names:
         results_representations = {}
         for model_path, model_desc_info in models_info.items():
             print("POPULATION " + str(population))
             res = process_model(model_desc_info, vcf_accessor, allele_accessor, population, sample_name)
             results_representations[model_path] = create_res_representation_for_model(res, model_desc_info, parsed_args.population)
         with open(os.path.join(out_dir, f'{sample_name}.sample.json'), 'w') as f:
             json.dump(results_representations, f, indent=4)


if __name__ == '__main__':
    main(sys.argv[1:])
