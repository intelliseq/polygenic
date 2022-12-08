"""
high level support for gwas files
"""
import logging
import sys
import math
import csv as csvlib

from tqdm import tqdm

import pandas as pd
import numpy as np
from scipy.signal import argrelextrema

from plotnine import ggplot
from plotnine import geom_point, geom_hline
from plotnine import aes, theme, theme_void, theme_minimal, scale_fill_manual, scale_color_manual, scale_size_continuous
from plotnine import element_rect, element_line
from plotnine import *

from polygenic.data.csv_accessor import CsvAccessor
from polygenic.tools.data.chromsizes import Chromsizes
from polygenic.tools.data.colors import Colors as colors
from polygenic.error.polygenic_exception import PolygenicException

logger = logging.getLogger('polygenic.data.' + __name__)

class Gwas(object):

    """
    class for manipulating gwas data
    """

    COLUMN_NAMES = ['chromosome',
                    'position',
                    'gnomadid',
                    'rsid',
                    'ref',
                    'alt',
                    'effect',
                    'pvalue',
                    'beta',
                    'af']

    def __init__(self):
        super().__init__()
        self.__data = {}
        self.__size = 0
        self.__iteridx = 0
        self.__iterlist = 0

    def __get_size(self):
        if self.__size > 0:
            return self.__size
        for(positions) in self.__data.values():
            self.__size += len(positions)
        return self.__size

    def __as_list(self):
        data = []
        for chromosome, positions in self.__data.items():
            for position, values in positions.items():
                data.append(values)
        return data

    def __iter__(self):
        self.__iteridx = 0
        self.__iterlist = self.__as_list()
        return self.__iterlist[self.__iteridx]

    def __next__(self):
        self.__iteridx += 1
        if self.__iteridx > self.__get_size():
            raise StopIteration
        return self.__iterlist[self.__iteridx]

    def __get_column_from_csv(self, csv, column_name):
        if column_name in csv.get_column_names():
            return csv.get_data()[column_name]
        else:
            # return vector filled with NaNs
            return np.full(csv.get_data().shape[0], np.nan)

    def load(self, input_path: str, column_mappings: dict = {}):
        logging.debug("Reading csv file: %s", input_path)
        csv = CsvAccessor(input_path)
        logging.debug("Standardizing column names")
        csv.standardize_column_names(column_mappings)
        logging.debug("Organizing data")

        # organizing is implemented by array iteration in numpy for efficiency
        with tqdm(total = csv.get_data().shape[0], file = sys.stdout, leave=True) as pbar:
            chromosome_vector = self.__get_column_from_csv(csv, 'chromosome')
            position_vector = self.__get_column_from_csv(csv, 'position')
            rsid_vector = self.__get_column_from_csv(csv, 'rsid')
            ref_vector = self.__get_column_from_csv(csv, 'ref')
            alt_vector = self.__get_column_from_csv(csv, 'alt')
            effect_vector = self.__get_column_from_csv(csv, 'effect')
            pvalue_vector = self.__get_column_from_csv(csv, 'pvalue')
            beta_vector = self.__get_column_from_csv(csv, 'beta')
            af_vector = self.__get_column_from_csv(csv, 'af')
            
            chromosomes = list(np.unique(chromosome_vector))
            for chromosome in chromosomes:
                if chromosome not in self.__data:
                    self.__data.update({chromosome: {}})
            for index in list(range(csv.get_data().shape[0])):
                values = {'chromosome': chromosome_vector[index],
                            'position': position_vector[index],
                            'gnomadid': "-".join((str(chromosome_vector[index]), str(position_vector[index]), str(ref_vector[index]), str(alt_vector[index]))),
                            'rsid': rsid_vector[index],
                            'ref': ref_vector[index],
                            'alt': alt_vector[index],
                            'effect': effect_vector[index],
                            'pvalue': pvalue_vector[index],
                            'beta': beta_vector[index],
                            'af': af_vector[index]}
                print(index)
                print(chromosome_vector[index])
                print(position_vector[index])
                print(values)            
                self.__data[chromosome_vector[index]].update({int(position_vector[index]): values})
                pbar.set_description('Organizing records')
                pbar.update(1)
            for chromosome, positions in self.__data.items():
                self.__data[chromosome] = dict(sorted(positions.items()))
    
    def save(self, output_path: str):
        logging.debug("Writing csv file: %s", output_path)
        with open(output_path, 'w', encoding = 'UTF8') as file:
            writer = csvlib.writer(file, delimiter = '\t', quotechar = '', quoting = csvlib.QUOTE_NONE)
            writer.writerow(self.COLUMN_NAMES)
            for positions in self.__data.values():
                for position in positions.values():
                    writer.writerow(position.values())

    def __get_size(self):
        size = 0
        for(positions) in self.__data.values():
            size += len(positions)
        return size

    # def clump(self):
    #     if self.__clumped_data:
    #         return self.__clumped_data
    #     clumped_data = []
    #     print(self.__get_size())
    #     with tqdm(total = self.__get_size(), file = sys.stdout, leave=False) as pbar:
    #         for chromosome, positions in self.__data.items():
    #             p_sequence = np.array([])
    #             pos_sequence = np.array([])
    #             for position in positions.values():
    #                 pbar.set_description('Clumping records')
    #                 pbar.update(1)
    #                 p_sequence = np.append(p_sequence, -math.log(position.get("pvalue"),10))
    #                 pos_sequence = np.append(pos_sequence, position.get("position"))
    #             extrema = argrelextrema(p_sequence, np.greater_equal, order=int(len(positions) / 10))[0]
    #             for extremum in extrema:
    #                 if p_sequence[extremum] > 5:
    #                     clumped_data.append(positions[int(pos_sequence[extremum])])
    #     self.__clumped_data = clumped_data
    #     return clumped_data

    # def validate(self):
    #     """
    #     validate gwas data
    #     """
    #     grch37_ref_vcf = 'polygenic/tests/resources/largefiles/dbsnp155.grch37.norm.vcf.gz'
    #     grch38_ref_vcf = 'polygenic/tests/resources/largefiles/dbsnp155.grch38.norm.vcf.gz'
        

    # def __get_filtered_data(self, pvalue_threshold: float = 0.05):
    #     """
    #     return a filtered gwas data
    #     """
    #     filtered_data = []
        
    #     with tqdm(total = self.__get_size(), file = sys.stdout, leave=False) as pbar:
    #         for chromosome, positions in self.__data.items():
    #             for position, gwas_record in positions.items():
    #                 if gwas_record.get("pvalue") < pvalue_threshold:
    #                     filtered_data.append(gwas_record)
    #         pbar.set_description('Filtering records')
    #         pbar.update(1)            
    #     return filtered_data

    # def __get_filtered_data_for_manhattan_plot(self, evry_nth: int = 2, pvalue_threshold: float = 0.05):
    #     """
    #     return a filtered gwas data for manhattan plot
    #     """
    #     filtered_data = []
        
    #     with tqdm(total = self.__get_size(), file = sys.stdout, leave=False) as pbar:
    #         for chromosome, positions in self.__data.items():
    #             for position, record in positions.items():
    #                 if record.get("pvalue") < pvalue_threshold and int(position) % evry_nth == 0:
    #                     filtered_data.append(record)
    #         pbar.set_description('Filtering records')
    #         pbar.update(1)
        
    #     return filtered_data

    # def plot_manhattan(self):
    #     """
    #     plot manhattan plot
    #     """

    #     threshold = 5
        
    #     data = self.__get_filtered_data_for_manhattan_plot()
    #     clumped_data = self.clump()
    #     data_length = len(data)
    #     clumped_data_length = len(clumped_data)
    #     #data.append(clumped_data)
    #     data = pd.DataFrame(data)
    #     data['set'] = 'regular'
    #     clumped_data = pd.DataFrame(clumped_data)
    #     clumped_data['set'] = 'selected'
    #     data = pd.concat([data, clumped_data], ignore_index=True, sort=False)
    #     chromsizes = Chromsizes().chromsizes["GRCh37"]
    #     cumulative_chromsizes = Chromsizes().chromsizes["GRCh37" + "_cumulative"]
    #     # get summary length of all chromosomes
    #     summary_length = sum(chromsizes.values())
    #     print(data.shape)
    #     print(str(data.iloc[0:3]))
    #     print(clumped_data.shape)
    #     print(str(clumped_data.iloc[0:3]))
    #     min_beta = min(abs(clumped_data['beta']))
    #     max_beta = max(abs(clumped_data['beta']))
    #     # add cumulative position column to data
    #     data['cumulative_position'] = data.apply(lambda row: int(row["position"]) + cumulative_chromsizes[str(int(row["chromosome"]))], axis=1)
    #     # # add log10 pvalue column to data
    #     data['log10_pvalue'] = data.apply(lambda row: -math.log10(row["pvalue"]), axis=1)
    #     # add different color to every second chromosome
    #     data['color'] = data.apply(lambda row: '0' if row["chromosome"] in ['X', 'Y'] or int(row["chromosome"]) % 2 == 1 else '1', axis=1)
    #     # data['color'] = data.apply(lambda row: row['color'] if row['log10_pvalue'] < 8 else '2', axis=1)
    #     # data['color'] = data.apply(lambda row: row['color'] if row['log10_pvalue'] < 20 else '3', axis=1)
    #     data['color'] = data.apply(lambda row: row['color'] if row['set'] != 'selected' else '5', axis=1)
    #     min_beta_size = 0.1
    #     max_beta_size = 5
    #     diff_beta_size = max_beta_size - min_beta_size
    #     data['size'] = min_beta_size
    #     data['size'] = data.apply(lambda row: row['size'] if row['set'] != 'selected' else (min_beta + diff_beta_size * ((abs(row['beta']) - min_beta) / max_beta)), axis=1)
    
    
    #     color_dict = {'0': colors.grey, 
    #                 '1': colors.grey_dark, 
    #                 '2': colors.teal, 
    #                 '3': colors.teal_dark, 
    #                 '4': colors.teal_darker,
    #                 '5': colors.pumpkin_light}

    #     plot = (
    #         ggplot(data, aes('cumulative_position', 'log10_pvalue', color = 'color', size = 'size')) + 
    #         geom_point() +
    #         geom_hline(yintercept = 8, color = colors.grey, size = 0.5, linetype = 'dashed') +
    #         geom_hline(yintercept = 5, color = colors.pumpkin, size = 0.5, linetype = 'dashed') +
    #         scale_color_manual(values=color_dict) +
    #         scale_size_continuous(range=(min_beta_size,max_beta_size)) + 
    #         theme(
    #             legend_position = "none",
    #             line = element_line(color = colors.grey, size = 1),
    #             rect = element_rect(fill = colors.grey_lighter),
    #             panel_grid_major = element_blank(),
    #             panel_grid_minor = element_blank(),
    #             panel_border = element_blank(),
    #             panel_background = element_blank(),
    #             axis_line = element_blank(),
    #             axis_title_x = element_blank(),
    #             axis_title_y = element_blank(),
    #             axis_text_x = element_blank(),
    #             axis_text_y = element_blank(),
    #             axis_ticks = element_line(size = 1),
    #             axis_ticks_length = 5,
    #             axis_ticks_length_minor = 2,
    #             axis_ticks_major_x = element_blank(),
    #             axis_ticks_major_y = element_line(color = colors.grey_light),
    #             axis_ticks_minor_y = element_line(color = colors.grey_light),
    #         )
    #     )

    #     plot.save("/tmp/plot.png", height=6, width=8)



