import os
import sys

module_path = os.path.abspath(__file__).rsplit(os.path.sep, 4)[0]
sys.path.insert(0, module_path)

import json
import argparse
import importlib
import subprocess
from typing import List
from typing import Dict
from typing import Union
from polygenic.lib.allele_frequency_sqlite_dao import read_multiple_rs
from polygenic.allele_frequency_gnomad import snp_data_to_dict

__version__ = '0.0.1'
model_pattern_end = 'model.py'


def discover_models(dir_to_be_checked: str, model_file_name_end: str = model_pattern_end) -> List[str]:
    model_import_strings = []
    for root, _, files in os.walk(dir_to_be_checked):
        for f in files:
            if f.endswith(model_file_name_end):
                model_import_strings.append(transform_path_to_import_string(os.path.join(root, f)))
    return model_import_strings


def transform_path_to_import_string(path: str) -> str:
    splitted = path.split(os.path.sep)
    first_part_of_import_string_index = splitted.index('src')
    return '.'.join(splitted[first_part_of_import_string_index:])[:-3]


def get_allele_freqs_using_gnomad_docker(rsids: List[str], docker_image_with_tag) -> Dict[
    str, Dict[str, Union[float, bool]]]:
    p = subprocess.run('docker run --rm {} --rsids {}'.format(docker_image_with_tag, ','.join(rsids)),
                       stdout=subprocess.PIPE, check=True, shell='bash')
    return json.loads(p.stdout)


if __name__ == '__main__':
    main_dir = os.path.abspath(__file__).rsplit('/', 2)[0]
    default_out = os.path.join(main_dir, 'resources', 'allele_frequencies', 'gnomad.json')
    default_model_dir = os.path.join(main_dir, 'resources')

    parser = argparse.ArgumentParser(description=('Creates allele frequency file for  rsids from given models ',
                                                  'using given gnomad docker image. ',
                                                  'If the allele frequency file exists and contains appropriate rsids,',
                                                  'this script does nothing (unless forced)'))
    parser.add_argument('--gnomad_image', type=str, default='')
    parser.add_argument('--output', type=str, default=default_out)
    parser.add_argument('--model_directory', type=str, default=default_model_dir,
                        help='Model file names should follow the pattern: "*{}"'.format(model_pattern_end))
    parser.add_argument('--force', action='store_true', help='Force creating file containing allele frequencies')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))

    args = parser.parse_args()

    rsids = set()
    for import_string in discover_models(args.model_directory):
        module = importlib.import_module(import_string)
        model = module.model
        rsids.update(model.get_rsids_used_by_model())

    previous_content = {}
    if not args.force:
        try:
            with open(args.output) as f:
                previous_content = json.load(f)
        except (IOError, json.decoder.JSONDecodeError):
            pass
    previous_rsids = previous_content.keys()

    if set(rsids) == set(previous_rsids):
        print('All rsids are present in current json file. Exiting.',
              'Run again with --force flag to force overwrite the file.')
        sys.exit(0)

    rsids_to_be_checked = list(set(rsids) - set(previous_rsids))
    print('Getting allele frequencies for {} new rsids'.format(len(rsids_to_be_checked)))
    if args.gnomad_image:
        new_allele_frequencies = get_allele_freqs_using_gnomad_docker(rsids_to_be_checked, args.gnomad_image)
    else:
        data_list = read_multiple_rs(list(rsids_to_be_checked))
        new_allele_frequencies = {data.rsid: snp_data_to_dict(data) for data in data_list}

    if len(rsids_to_be_checked) != len(new_allele_frequencies):
        found_rsids = set(new_allele_frequencies.keys())
        raise RuntimeError("Data for {} not found. Exiting".format(
            ', '.join(rsid for rsid in rsids_to_be_checked if rsid not in found_rsids)))
    previous_content.update(new_allele_frequencies)
    with open(args.output, 'w') as fout:
        json.dump(previous_content, fout, sort_keys=True, indent=4)
    print('Wrote allele frequencies for {} rsids'.format(len(previous_content)))
