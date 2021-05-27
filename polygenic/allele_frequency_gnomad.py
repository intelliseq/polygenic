import os
import sys

module_path = os.path.abspath(__file__).rsplit(os.path.sep, 4)[0]
sys.path.insert(0, module_path)

import argparse
import json
from typing import Dict
from typing import Any
from polygenic.lib.allele_frequency_sqlite_dao import get_column_names
from polygenic.lib.allele_frequency_sqlite_dao import read_multiple_rs
from polygenic.lib.gnomad_parser import SNPdata

__version__ = '1.0.0'

population_dictionary = {'AF': 'Alternate allele frequency in samples',
                         'AF_nfe_seu': 'Alternate allele frequency in samples of Southern European ancestry',
                         'controls_AF_afr_male': 'Alternate allele frequency in male samples of African-American/African ancestry in the controls subset',
                         'non_topmed_AF_amr': 'Alternate allele frequency in samples of Latino ancestry in the non_topmed subset',
                         'AF_raw': 'Alternate allele frequency in samples, before removing low-confidence genotypes',
                         'AF_fin_female': 'Alternate allele frequency in female samples of Finnish ancestry',
                         'non_neuro_AF_asj_female': 'Alternate allele frequency in female samples of Ashkenazi Jewish ancestry in the non_neuro subset',
                         'non_neuro_AF_afr_male': 'Alternate allele frequency in male samples of African-American/African ancestry in the non_neuro subset',
                         'AF_afr_male': 'Alternate allele frequency in male samples of African-American/African ancestry',
                         'AF_afr': 'Alternate allele frequency in samples of African-American/African ancestry',
                         'non_neuro_AF_afr_female': 'Alternate allele frequency in female samples of African-American/African ancestry in the non_neuro subset',
                         'non_topmed_AF_amr_female': 'Alternate allele frequency in female samples of Latino ancestry in the non_topmed subset',
                         'non_topmed_AF_oth_female': 'Alternate allele frequency in female samples of Other ancestry in the non_topmed subset',
                         'AF_eas_female': 'Alternate allele frequency in female samples of East Asian ancestry',
                         'AF_afr_female': 'Alternate allele frequency in female samples of African-American/African ancestry',
                         'non_neuro_AF_female': 'Alternate allele frequency in female samples in the non_neuro subset',
                         'controls_AF_afr': 'Alternate allele frequency in samples of African-American/African ancestry in the controls subset',
                         'AF_nfe_onf': 'Alternate allele frequency in samples of Other Non-Finnish European ancestry',
                         'controls_AF_fin_male': 'Alternate allele frequency in male samples of Finnish ancestry in the controls subset',
                         'non_neuro_AF_nfe_nwe': 'Alternate allele frequency in samples of North-Western European ancestry in the non_neuro subset',
                         'AF_fin_male': 'Alternate allele frequency in male samples of Finnish ancestry',
                         'AF_nfe_female': 'Alternate allele frequency in female samples of Non-Finnish European ancestry',
                         'AF_amr': 'Alternate allele frequency in samples of Latino ancestry',
                         'non_topmed_AF_nfe_male': 'Alternate allele frequency in male samples of Non-Finnish European ancestry in the non_topmed subset',
                         'AF_eas': 'Alternate allele frequency in samples of East Asian ancestry',
                         'non_neuro_AF_nfe_female': 'Alternate allele frequency in female samples of Non-Finnish European ancestry in the non_neuro subset',
                         'non_neuro_AF_afr': 'Alternate allele frequency in samples of African-American/African ancestry in the non_neuro subset',
                         'controls_AF_raw': 'Alternate allele frequency in samples in the controls subset, before removing low-confidence genotypes',
                         'controls_AF_male': 'Alternate allele frequency in male samples in the controls subset',
                         'non_topmed_AF_male': 'Alternate allele frequency in male samples in the non_topmed subset',
                         'controls_AF_nfe_female': 'Alternate allele frequency in female samples of Non-Finnish European ancestry in the controls subset',
                         'non_neuro_AF_amr': 'Alternate allele frequency in samples of Latino ancestry in the non_neuro subset',
                         'non_neuro_AF_eas_female': 'Alternate allele frequency in female samples of East Asian ancestry in the non_neuro subset',
                         'AF_asj_male': 'Alternate allele frequency in male samples of Ashkenazi Jewish ancestry',
                         'controls_AF_nfe_male': 'Alternate allele frequency in male samples of Non-Finnish European ancestry in the controls subset',
                         'non_neuro_AF_fin': 'Alternate allele frequency in samples of Finnish ancestry in the non_neuro subset',
                         'AF_oth_female': 'Alternate allele frequency in female samples of Other ancestry',
                         'controls_AF_nfe': 'Alternate allele frequency in samples of Non-Finnish European ancestry in the controls subset',
                         'controls_AF_oth_female': 'Alternate allele frequency in female samples of Other ancestry in the controls subset',
                         'controls_AF_asj': 'Alternate allele frequency in samples of Ashkenazi Jewish ancestry in the controls subset',
                         'non_neuro_AF_amr_male': 'Alternate allele frequency in male samples of Latino ancestry in the non_neuro subset',
                         'controls_AF_nfe_nwe': 'Alternate allele frequency in samples of North-Western European ancestry in the controls subset',
                         'AF_nfe_nwe': 'Alternate allele frequency in samples of North-Western European ancestry',
                         'controls_AF_nfe_seu': 'Alternate allele frequency in samples of Southern European ancestry in the controls subset',
                         'non_neuro_AF_amr_female': 'Alternate allele frequency in female samples of Latino ancestry in the non_neuro subset',
                         'non_neuro_AF_nfe_onf': 'Alternate allele frequency in samples of Other Non-Finnish European ancestry in the non_neuro subset',
                         'non_topmed_AF_eas_male': 'Alternate allele frequency in male samples of East Asian ancestry in the non_topmed subset',
                         'controls_AF_amr_female': 'Alternate allele frequency in female samples of Latino ancestry in the controls subset',
                         'non_neuro_AF_fin_male': 'Alternate allele frequency in male samples of Finnish ancestry in the non_neuro subset',
                         'AF_female': 'Alternate allele frequency in female samples',
                         'non_neuro_AF_oth_male': 'Alternate allele frequency in male samples of Other ancestry in the non_neuro subset',
                         'non_topmed_AF_nfe_est': 'Alternate allele frequency in samples of Estonian ancestry in the non_topmed subset',
                         'non_topmed_AF_nfe_nwe': 'Alternate allele frequency in samples of North-Western European ancestry in the non_topmed subset',
                         'non_topmed_AF_amr_male': 'Alternate allele frequency in male samples of Latino ancestry in the non_topmed subset',
                         'non_topmed_AF_nfe_onf': 'Alternate allele frequency in samples of Other Non-Finnish European ancestry in the non_topmed subset',
                         'controls_AF_eas_male': 'Alternate allele frequency in male samples of East Asian ancestry in the controls subset',
                         'controls_AF_oth_male': 'Alternate allele frequency in male samples of Other ancestry in the controls subset',
                         'non_topmed_AF': 'Alternate allele frequency in samples in the non_topmed subset',
                         'controls_AF_fin': 'Alternate allele frequency in samples of Finnish ancestry in the controls subset',
                         'non_neuro_AF_nfe': 'Alternate allele frequency in samples of Non-Finnish European ancestry in the non_neuro subset',
                         'non_neuro_AF_fin_female': 'Alternate allele frequency in female samples of Finnish ancestry in the non_neuro subset',
                         'non_topmed_AF_nfe_seu': 'Alternate allele frequency in samples of Southern European ancestry in the non_topmed subset',
                         'controls_AF_eas_female': 'Alternate allele frequency in female samples of East Asian ancestry in the controls subset',
                         'non_topmed_AF_asj': 'Alternate allele frequency in samples of Ashkenazi Jewish ancestry in the non_topmed subset',
                         'controls_AF_nfe_onf': 'Alternate allele frequency in samples of Other Non-Finnish European ancestry in the controls subset',
                         'non_neuro_AF': 'Alternate allele frequency in samples in the non_neuro subset',
                         'non_topmed_AF_nfe': 'Alternate allele frequency in samples of Non-Finnish European ancestry in the non_topmed subset',
                         'non_topmed_AF_raw': 'Alternate allele frequency in samples in the non_topmed subset, before removing low-confidence genotypes',
                         'non_neuro_AF_nfe_est': 'Alternate allele frequency in samples of Estonian ancestry in the non_neuro subset',
                         'non_topmed_AF_oth_male': 'Alternate allele frequency in male samples of Other ancestry in the non_topmed subset',
                         'AF_nfe_est': 'Alternate allele frequency in samples of Estonian ancestry',
                         'non_topmed_AF_afr_male': 'Alternate allele frequency in male samples of African-American/African ancestry in the non_topmed subset',
                         'AF_eas_male': 'Alternate allele frequency in male samples of East Asian ancestry',
                         'controls_AF_eas': 'Alternate allele frequency in samples of East Asian ancestry in the controls subset',
                         'non_neuro_AF_eas_male': 'Alternate allele frequency in male samples of East Asian ancestry in the non_neuro subset',
                         'non_neuro_AF_asj_male': 'Alternate allele frequency in male samples of Ashkenazi Jewish ancestry in the non_neuro subset',
                         'controls_AF_oth': 'Alternate allele frequency in samples of Other ancestry in the controls subset',
                         'AF_nfe': 'Alternate allele frequency in samples of Non-Finnish European ancestry',
                         'non_topmed_AF_female': 'Alternate allele frequency in female samples in the non_topmed subset',
                         'non_neuro_AF_asj': 'Alternate allele frequency in samples of Ashkenazi Jewish ancestry in the non_neuro subset',
                         'non_topmed_AF_eas_female': 'Alternate allele frequency in female samples of East Asian ancestry in the non_topmed subset',
                         'non_neuro_AF_raw': 'Alternate allele frequency in samples in the non_neuro subset, before removing low-confidence genotypes',
                         'non_topmed_AF_eas': 'Alternate allele frequency in samples of East Asian ancestry in the non_topmed subset',
                         'non_topmed_AF_fin_male': 'Alternate allele frequency in male samples of Finnish ancestry in the non_topmed subset',
                         'AF_fin': 'Alternate allele frequency in samples of Finnish ancestry',
                         'AF_nfe_male': 'Alternate allele frequency in male samples of Non-Finnish European ancestry',
                         'controls_AF_amr_male': 'Alternate allele frequency in male samples of Latino ancestry in the controls subset',
                         'controls_AF_afr_female': 'Alternate allele frequency in female samples of African-American/African ancestry in the controls subset',
                         'controls_AF_amr': 'Alternate allele frequency in samples of Latino ancestry in the controls subset',
                         'AF_asj_female': 'Alternate allele frequency in female samples of Ashkenazi Jewish ancestry',
                         'non_neuro_AF_eas': 'Alternate allele frequency in samples of East Asian ancestry in the non_neuro subset',
                         'non_neuro_AF_male': 'Alternate allele frequency in male samples in the non_neuro subset',
                         'AF_asj': 'Alternate allele frequency in samples of Ashkenazi Jewish ancestry',
                         'controls_AF_nfe_est': 'Alternate allele frequency in samples of Estonian ancestry in the controls subset',
                         'non_topmed_AF_asj_female': 'Alternate allele frequency in female samples of Ashkenazi Jewish ancestry in the non_topmed subset',
                         'non_topmed_AF_oth': 'Alternate allele frequency in samples of Other ancestry in the non_topmed subset',
                         'non_topmed_AF_fin_female': 'Alternate allele frequency in female samples of Finnish ancestry in the non_topmed subset',
                         'AF_oth': 'Alternate allele frequency in samples of Other ancestry',
                         'non_neuro_AF_nfe_male': 'Alternate allele frequency in male samples of Non-Finnish European ancestry in the non_neuro subset',
                         'controls_AF_female': 'Alternate allele frequency in female samples in the controls subset',
                         'non_topmed_AF_fin': 'Alternate allele frequency in samples of Finnish ancestry in the non_topmed subset',
                         'non_topmed_AF_nfe_female': 'Alternate allele frequency in female samples of Non-Finnish European ancestry in the non_topmed subset',
                         'controls_AF_asj_male': 'Alternate allele frequency in male samples of Ashkenazi Jewish ancestry in the controls subset',
                         'non_topmed_AF_asj_male': 'Alternate allele frequency in male samples of Ashkenazi Jewish ancestry in the non_topmed subset',
                         'non_neuro_AF_oth': 'Alternate allele frequency in samples of Other ancestry in the non_neuro subset',
                         'AF_male': 'Alternate allele frequency in male samples',
                         'controls_AF_fin_female': 'Alternate allele frequency in female samples of Finnish ancestry in the controls subset',
                         'controls_AF_asj_female': 'Alternate allele frequency in female samples of Ashkenazi Jewish ancestry in the controls subset',
                         'AF_amr_male': 'Alternate allele frequency in male samples of Latino ancestry',
                         'AF_amr_female': 'Alternate allele frequency in female samples of Latino ancestry',
                         'AF_oth_male': 'Alternate allele frequency in male samples of Other ancestry',
                         'non_neuro_AF_nfe_seu': 'Alternate allele frequency in samples of Southern European ancestry in the non_neuro subset',
                         'non_topmed_AF_afr_female': 'Alternate allele frequency in female samples of African-American/African ancestry in the non_topmed subset',
                         'non_topmed_AF_afr': 'Alternate allele frequency in samples of African-American/African ancestry in the non_topmed subset',
                         'controls_AF': 'Alternate allele frequency in samples in the controls subset',
                         'non_neuro_AF_oth_female': 'Alternate allele frequency in female samples of Other ancestry in the non_neuro subset',
                         'controls_AF_popmax': 'Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry) in the controls subset',
                         'AF_popmax': 'Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry)',
                         'non_neuro_AF_popmax': 'Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry) in the non_neuro subset',
                         'non_topmed_AF_popmax': 'Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry) in the non_topmed subset'}


def get_available_populations() -> Dict[str, str]:
    column_names = get_column_names()
    return {name: population_dictionary[name] for name in column_names if name.startswith('AF')}


def snp_data_to_dict(snp_data: SNPdata) -> Dict[str, Any]:
    d = {field: value for field, value in zip(SNPdata._fields, snp_data) if not field=='rsid'}
    d['is_biallelic'] = bool(d['is_unique'])
    del d['is_unique']
    return d


if __name__ == '__main__':
    #parser = argparse.ArgumentParser(description='rsid-indexed allele frequencies from gnomAD 3.0')
    parser = argparse.ArgumentParser(description='rsid-indexed allele frequencies from gnomAD 2.1.1')
    parser.add_argument('--populations', action='store_true', help='List available populations')
    parser.add_argument('--rsids', type=str, help='Coma-separated rsids')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))
    args = parser.parse_args()

    if (args.populations and args.rsids) or (not args.rsids and not args.populations):
        raise RuntimeError('Please decide whether you want to receive population list or allele frequencies')
    elif args.populations:
        print(get_available_populations())
    else:
        rsid_list = args.rsids.split(',')
        if not all(rsid.startswith('rs') for rsid in rsid_list):
            raise RuntimeError('Did you provide valid rsids?')
        data_list = read_multiple_rs(rsid_list)
        print(json.dumps({data.rsid: snp_data_to_dict(data) for data in data_list}))
