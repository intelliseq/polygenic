"""
high level support for gwas files
"""
import logging
import sys
# import os

import pandas as pd
import numpy as np

from tqdm import tqdm
from collections import OrderedDict

import polygenic.data.csv_accessor as csv_accessor
# from polygenic.error.polygenic_exception import PolygenicException



logger = logging.getLogger('polygenic.data.' + __name__)

class GwasRow(object):
    """
    class for manipulating gwas records
    """

    def __init__(self, values: dict):
        super().__init__()
        self.__values = values

    def get(self, key: str):
        """
        return the value of a key
        """
        return self.__values[key]

class GwasData(object):

    """
    class for manipulating gwas data
    """

    def __init__(self):
        super().__init__()
        self.__gwas = {}
        self.__count = 0


    def load_gwas_data_from_csv(self, gwas_file_path: str, column_names: dict):
        logging.debug("Reading csv file: %s", gwas_file_path)
        csv = csv_accessor.CsvAccessor(gwas_file_path)
        logging.debug("Standardizing column names")
        csv.standardize_column_names(column_names)
        logging.debug("Organizing data")

        # organizing is implemente by array iteration in numpy for efficiency
        with tqdm(total = csv.get_data().shape[0], file = sys.stdout) as pbar:
            chromosome_vector = csv.get_data()['chromosome']
            position_vector = csv.get_data()['position']
            rsid_vector = csv.get_data()['rsid']
            ref_vector = csv.get_data()['ref']
            alt_vector = csv.get_data()['alt']
            effect_vector = csv.get_data()['effect']
            pvalue_vector = csv.get_data()['pvalue']
            beta_vector = csv.get_data()['beta']
            
            chromosomes = list(np.unique(chromosome_vector))
            for chromosome in chromosomes:
                if chromosome not in self.__gwas:
                    self.__gwas.update({chromosome: OrderedDict()})
            for index in list(range(csv.get_data().shape[0])):
                values = {'chromosome': chromosome_vector[index],
                            'position': position_vector[index],
                            'rsid': rsid_vector[index],
                            'ref': ref_vector[index],
                            'alt': alt_vector[index],
                            'effect': effect_vector[index],
                            'pvalue': pvalue_vector[index],
                            'beta': beta_vector[index]}
                self.__gwas[chromosome_vector[index]].update({position_vector[index]: GwasRow(values)})
                pbar.update(1)
#                self.__count += 1
#                self.__gwas[self.__count] = GwasRow(row)
            # for index, row in csv.get_data().iterrows():
            #     self.__count += 1
            #     #if row.get("chromosome") not in self.__gwas:
            #     #    self.__gwas.update({row.get("chromosome"): {}})
            #     #self.__gwas.get(row.get("chromosome")).update({row.get("position"): GwasRow(row)})
            #     pbar.set_description('Organizing gwas records (estimated): %d' % index)
            #     pbar.update(1)



