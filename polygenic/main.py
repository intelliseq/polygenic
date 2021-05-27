#!/usr/bin/python3

import os
import sys
import glob
import json
import logging
import argparse
import importlib
from typing import Dict
from typing import List
from typing import Union
from typing import Iterable

MODULE_PATH = os.path.abspath(__file__).rsplit(os.path.sep, 4)[0]
sys.path.insert(0, MODULE_PATH)

from polygenic.version import version
from polygenic.lib.output import create_res_representation_for_model
from polygenic.lib.data_access.allele_frequency_accessor import AlleleFrequencyAccessor
from polygenic.lib.data_access.data_accessor import VcfAccessor
from polygenic.lib.model.seqql import PolygenicRiskScore
from polygenic.lib.model.seqql import Data
from polygenic.lib.data_access.dto import ModelDescriptionInfo


class ImproperPopulationForModelError(RuntimeError):
    pass


def expand_path(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path)) if path else ''


def create_model_description_info(models_paths:Dict[str, str], description_paths:Dict[str, str], mapping:Dict) -> ModelDescriptionInfo:
    model_fname:str = mapping['model']
    descriptions_fnames:List[str] = mapping['descriptions']
    return ModelDescriptionInfo(
        model_fname=model_fname,
        desc_paths={fname:description_paths[fname] for fname in descriptions_fnames},
        model_path=models_paths[model_fname]
    )


def process_model(model_description_info:ModelDescriptionInfo, initial_accessor:VcfAccessor, imputed_accessor:VcfAccessor, allele_freq_accessor:AlleleFrequencyAccessor, pop:str, sample_name:str):
    package_str = os.path.dirname(model_description_info.model_path)
    if package_str not in sys.path:
        sys.path.insert(0, package_str)
    module = importlib.import_module(model_description_info.model_fname.split('.')[0])
    pop_short = pop.replace('AF', '').lstrip('_')
    if pop_short != module.trait_was_prepared_for_population:
        raise ImproperPopulationForModelError(f"You requested data for population {pop_short} while the model was prepared for {module.trait_was_prepared_for_population}")
    model: PolygenicRiskScore = module.model
    data = Data(initial_accessor, imputed_accessor, allele_freq_accessor, sample_name, pop, model)
    return data.compute_model(), module.trait_was_prepared_for_population


def logical_xor(a, b) -> bool:
    return bool(a) != bool(b)


def load_models_when_mapping_present(directory:str, mappings_between_models_descriptions_and_plot_data:List[Dict[str, Union[str,List[str]]]], logger) -> Dict[str, ModelDescriptionInfo]:
    validate_data_completness(directory, mappings_between_models_descriptions_and_plot_data)
    ret = {}
    for mapping in mappings_between_models_descriptions_and_plot_data:
        try:
            plot_data_file_name = mapping['plot_data']
        except KeyError:
            logger.info(f'No plot data in {mapping}')
            plot_data_file_name = ''
        try:
            descriptions_filenames = mapping['descriptions']
        except KeyError:
            logger.info(f'No descriptions in {mapping}')
            descriptions_filenames = []
        plot_data_path = expand_path(os.path.join(directory, plot_data_file_name)) if plot_data_file_name else None
        descriptions_paths = [expand_path(os.path.join(directory, description)) for description in descriptions_filenames] if descriptions_filenames else []
        model_path = expand_path(os.path.join(directory, mapping['model']))
        ret[model_path] = ModelDescriptionInfo(
            model_fname=mapping['model'],
            model_path = model_path,
            desc_paths=descriptions_paths,
            plot_data_path=plot_data_path
        )
    return ret


def validate_data_completness(directory:str, mappings_between_models_descriptions_and_plot_data:List[Dict[str, Union[str,List[str]]]]):
    files_in_directory = [x for x in os.listdir(directory) if os.path.isfile(os.path.join(directory, x))]
    files_in_mappings_not_flattened = [item for mapping in mappings_between_models_descriptions_and_plot_data for item in mapping.values()]
    files_in_mappings = []
    for item in files_in_mappings_not_flattened:
        if isinstance(item, str):
            files_in_mappings.append(item)
        elif isinstance(item, list):
            files_in_mappings.extend(item)
        else:
            raise RuntimeError(f'Data not supported: {item}')
    files_missing_in_mapping = set(files_in_directory) - set(files_in_mappings)
    files_missing_in_data = set(files_in_mappings) - set(files_in_directory)
    if files_missing_in_mapping:
        raise RuntimeError(f"{', '.join(files_missing_in_mapping)} found in data, but not in mapping. Exiting.")
    if files_missing_in_data:
        raise RuntimeError(f"{', '.join(files_missing_in_data)} found in mapping, but not in data. Exiting.")


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


def abspath_to_model_dir_from_repo(directory:str) -> str:
    return expand_path(os.path.join(MODULE_PATH, 'src', 'main', 'resources', directory))


logger = logging.getLogger('description_language')


def main(args = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='')  # todo dodać opis
    parser.add_argument('--imputed_files', nargs='+', required=True,
                        help='Compressed vcfs')
    parser.add_argument('--log_file', type=str, default='/app_data/log.log')
    parser.add_argument('--curated_initial_vcfs', nargs='+',
                        help='Curated initial experimental file, before phasing, compressed and splitted by chromosome. If provided, will override imputed data.')
    parser.add_argument('--out_dir', type=str, default="", help='Directory for result jsons.')
    parser.add_argument('--population', type=str, default='nfe',
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
    parser.add_argument('--models_and_descriptions_path', type=str, default='', help="Path to a directory containing models and corresponding_descriptions")
    parser.add_argument('--mapping_json', type=str, default='', help="A file containing mapping between models and descriptions")
    parser.add_argument('--allele_freq_json', type=str, default='',
                        help="A file containing allele freq data")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(version))

    parsed_args = parser.parse_args(args)

    log_directory = os.path.dirname(os.path.abspath(os.path.expanduser(parsed_args.log_file)))
    if log_directory:
        try:
            os.makedirs(log_directory)
        except OSError:
            pass

    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(parsed_args.log_file)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    out_dir = expand_path(parsed_args.out_dir)

    population = 'AF' if not parsed_args.population else 'AF_' + parsed_args.population
    allele_freq_path = parsed_args.allele_freq_json or os.path.join(MODULE_PATH, 'src', 'main', 'resources', 'allele_frequencies', 'gnomad.json') #todo sprawdzić

    # mapping file should be available only when we supply this application with models from the outside
    if parsed_args.mapping_json and not parsed_args.models_and_descriptions_path:
        raise RuntimeError('You supplied mapping file, but no models')

    # coupling model files with descriptions and plot data
    directories = [abspath_to_model_dir_from_repo(trait_dir) for trait_dir in parsed_args.traits_dirs] or []
    model_files_info = {}
    if parsed_args.mapping_json:
        with open(expand_path(parsed_args.mapping_json)) as f:
            models_descriptions_mappings = json.load(f)
        model_files_info.update(
            load_models_when_mapping_present(parsed_args.models_and_descriptions_path, models_descriptions_mappings,
                                             logger))
    elif parsed_args.models_and_descriptions_path:
        directories.append(parsed_args.models_and_descriptions_path)
    model_files_info.update(load_models_when_mapping_absent(directories, parsed_args.population, logger))
    if not model_files_info:
        raise RuntimeError("No models loaded. Exiting.")

    # create data accessors
    initial_accessor = VcfAccessor([expand_path(path) for path in parsed_args.curated_initial_vcfs])
    imputed_accessor = VcfAccessor([expand_path(path) for path in parsed_args.imputed_files])
    allele_accessor = AlleleFrequencyAccessor(allele_freq_json_path=allele_freq_path)
    sample_names = initial_accessor.sample_names

    for sample_name in sample_names:
        results_representations = {}
        for model_path, model_desc_info in model_files_info.items():
            res, pop_in_model = process_model(model_desc_info, initial_accessor, imputed_accessor, allele_accessor, population, sample_name)
            results_representations[model_path] = create_res_representation_for_model(res, model_desc_info, pop_in_model)
        with open(os.path.join(out_dir, f'{sample_name}.sample.json'), 'w') as f:
            json.dump(results_representations, f, indent=4)


if __name__ == '__main__':
    main(sys.argv[1:])
