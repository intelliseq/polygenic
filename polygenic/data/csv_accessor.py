import logging
import os 
import math

import pandas as pd
import numpy as np
from polygenic.error.polygenic_exception import PolygenicException

logger = logging.getLogger('polygenic.data.' + __name__)

class CsvAccessor(object):
    """
    high level support for csv files
    """
    def __init__(self, csv_path: str):
        super().__init__()
        self.__path = csv_path
        self.__delimiter = '\t'
        self.__rsid_column_index = None
        if not os.path.exists(self.__path):
            raise PolygenicException('Can not access {path}'.format(path = self.__path))
        self.__data = self.read_data()

    def __find_index_of_column_by_name(self, name):
        """
        return the index of a column
        """
        for column_index, column_name in enumerate(self.__data.columns):
            if name.lower() in column_name.lower():
                return column_index
        return None

    def get_rsid_column_index(self, rsid_column_name: str = 'rsid'):
        """
        return the index of the rsid column
        """
        if self.__rsid_column_index is None:
            self.__rsid_column_index = self.__find_index_of_column_by_name(rsid_column_name)
        return self.__rsid_column_index

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