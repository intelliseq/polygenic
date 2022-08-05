"""
high level support for csv files
"""
import logging
import os

import pandas as pd
import numpy as np
from polygenic.error.polygenic_exception import PolygenicException

logger = logging.getLogger('polygenic.data.' + __name__)

class CsvAccessor(object):
    """
    class for reading csv files
    """
    def __init__(self, csv_path: str):
        super().__init__()
        self.__path = csv_path
        self.__delimiter = '\t'
        self.__rsid_column_index = None
        self.__chrom_column_index = None
        self.__pos_column_index = None
        self.__ref_column_index = None
        self.__alt_column_index = None
        self.__effect_column_index = None
        self.__pvalue_column_index = None
        self.__beta_column_index = None
        if not os.path.exists(self.__path):
            raise PolygenicException(f"Can not access {self.__path}")
        self.__data = self.read_data()

    def __find_index_of_column_by_name(self, name: str, equals_instead_of_contains: bool = True):
        """
        return the index of a column
        """
        for column_index, column_name in enumerate(self.__data.columns):
            if equals_instead_of_contains:
                if name.lower() == column_name.lower():
                    return column_index
            else:
                if name.lower() in column_name.lower():
                    return column_index
        return None

    def get_rsid_column_index(self, rsid_column_name: str = 'rsid'):
        """
        return the index of the rsid column
        """
        if self.__rsid_column_index is None:
            self.__rsid_column_index = self.__find_index_of_column_by_name(rsid_column_name)
        if self.__rsid_column_index is None:
            self.__rsid_column_index = self.__find_index_of_column_by_name('rsid')
        if self.__rsid_column_index is None:
            self.__rsid_column_index = self.__find_index_of_column_by_name('id')
        return self.__rsid_column_index

    def get_chrom_column_index(self, chrom_column_name: str = 'chromosome'):
        """
        return the index of the chrom column
        """
        if self.__chrom_column_index is None:
            self.__chrom_column_index = self.__find_index_of_column_by_name(chrom_column_name)
        if self.__chrom_column_index is None:
            self.__chrom_column_index = self.__find_index_of_column_by_name('chromosome')
        if self.__chrom_column_index is None:
            self.__chrom_column_index = self.__find_index_of_column_by_name('chrom')
        if self.__chrom_column_index is None:
            self.__chrom_column_index = self.__find_index_of_column_by_name('chr')
        return self.__chrom_column_index

    def get_pos_column_index(self, pos_column_name: str = 'position'):
        """
        return the index of the pos column
        """
        if self.__pos_column_index is None:
            self.__pos_column_index = self.__find_index_of_column_by_name(pos_column_name)
        if self.__pos_column_index is None:
            self.__pos_column_index = self.__find_index_of_column_by_name('position')
        if self.__pos_column_index is None:
            self.__pos_column_index = self.__find_index_of_column_by_name('pos')
        if self.__pos_column_index is None:
            self.__pos_column_index = self.__find_index_of_column_by_name('bp')
        return self.__pos_column_index

    def get_ref_column_index(self, ref_column_name: str = 'ref'):
        """
        return the index of the ref column
        """
        if self.__ref_column_index is None:
            self.__ref_column_index = self.__find_index_of_column_by_name(ref_column_name)
        return self.__ref_column_index

    def get_alt_column_index(self, alt_column_name: str = 'alt'):
        """
        return the index of the alt column
        """
        if self.__alt_column_index is None:
            self.__alt_column_index = self.__find_index_of_column_by_name(alt_column_name)
        return self.__alt_column_index
    
    def get_effect_column_index(self, effect_column_name: str = 'effect'):
        """
        return the index of the effect column
        """
        if self.__effect_column_index is None:
            self.__effect_column_index = self.__find_index_of_column_by_name(effect_column_name)
        return self.__effect_column_index
    
    def get_pvalue_column_index(self, pvalue_column_name: str = 'pvalue'):
        """
        return the index of the pvalue column
        """
        if self.__pvalue_column_index is None:
            self.__pvalue_column_index = self.__find_index_of_column_by_name(pvalue_column_name)
        return self.__pvalue_column_index

    def get_beta_column_index(self, beta_column_name: str = 'beta'):
        """
        return the index of the beta column
        """
        if self.__beta_column_index is None:
            self.__beta_column_index = self.__find_index_of_column_by_name(beta_column_name)
        return self.__beta_column_index

    def get_symbol_for_rsid(self, rsid):
        """
        return the symbol for a rsid
        """
        data = self.__data
        data = data.loc[data["rsid"] == rsid]
        if len(data.index) == 0:
            return None
        return data['symbol'].head(1).iloc[0]

    def get_column_names(self):
        """
        return the column names of the csv file
        """
        return self.__data.columns

    def get_data(self):
        """
        return the dataframe
        """
        return self.__data

    def read_data(self):
        """
        read the csv file and return a dataframe
        """
        return pd.read_csv(filepath_or_buffer = self.__path, sep = self.__delimiter)

    def get_symbol_for_genomic_position(self, chrom, pos):
        """
        return the symbol for a genomic position
        """
        data = self.__data
        data = data.loc[data["chromosome"] == str(chrom)]
        if len(data.index) == 0:
            return None
        data = data.assign(pos_start = abs(data["start"] - np.int64(pos)),
                           pos_end = abs(data["end"] - np.int64(pos)))
        data = data.assign(position = data[["pos_start", "pos_end"]].min(axis = 1))
        return data.sort_values(by=['pos_start'])['symbol'].head(1).iloc[0]
