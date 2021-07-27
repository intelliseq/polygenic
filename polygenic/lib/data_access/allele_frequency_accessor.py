import os
import json
import logging
from typing import Dict
from typing import Union
from functools import partial
from polygenic.lib.allele_frequency_sqlite_dao import read_multiple_rs

logger = logging.getLogger('description_language.' + __name__)


class DataNotUsefulInModels(RuntimeError):
    pass


class AlleleFrequencyAccessor(object):

    def __init__(self, allele_freq_json_path: str = '', allele_freq_database_path: str = ''):
        assert not (allele_freq_json_path and allele_freq_database_path)
        if allele_freq_json_path:
            with open(allele_freq_json_path) as f:
                content = json.load(f)
            self.reading_method = partial(get_allele_freq_from_json, allele_freq_dict=content)
        if allele_freq_database_path:
            assert os.path.isfile(allele_freq_database_path)
            if not allele_freq_database_path.endswith('allele_frequency_db_gnomad_2.1.1.sqlite'):
                raise RuntimeError("implemented only for gnomad 2.1.1")
            self.reading_method = partial(get_allele_freq_from_db, path=allele_freq_database_path)

    def __call__(self, rsid: str, population: str):
        return self.reading_method(rsid, population)


def get_allele_freq_from_db(rsid: str, population_name: str, path: str = ''):
    snp_data = read_multiple_rs([rsid], path=path)[0]
    ref_allele = snp_data.ref
    alt_allele = snp_data.alt
    alt_allele_freq = snp_data.__dict__[population_name]
    if not snp_data.is_unique:
        logger.info(
            f'{rsid} is multiallelic but only two alleles are provided. Only {ref_allele} and {alt_allele} were considered')
    return {alt_allele: alt_allele_freq, ref_allele: 1 - alt_allele_freq}


def get_allele_freq_from_json(rsid: str, population_name: str, allele_freq_dict:Dict[str, Dict[str, Union[str, float, bool]]]):
    snp_data = allele_freq_dict[rsid]
    ref_allele = snp_data['ref']
    alt_allele = snp_data['alt']
    alt_allele_freq = snp_data[population_name]
    if not snp_data['is_biallelic']:
        logger.info(
            f'{rsid} is multiallelic but only two alleles are provided. Only {ref_allele} and {alt_allele} were considered')
    return {alt_allele: alt_allele_freq, ref_allele: 1 - alt_allele_freq}


if __name__ == '__main__':
    afa = AlleleFrequencyAccessor(
        allele_freq_database_path='/home/wojtek/PycharmProjects/mobigen/src/main/resources/databases/allele_frequency_db_gnomad_2.1.1.sqlite')
    print(afa('rs925966', 'AF_nfe'))
    afa2 = AlleleFrequencyAccessor(allele_freq_json_path='/home/wojtek/PycharmProjects/mobigen/src/main/resources/allele_frequencies/gnomad.json')
    print(afa2('rs925966', 'AF_asj'))
