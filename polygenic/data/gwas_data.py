"""
high level support for gwas files
"""
import logging
import sys
import math

from collections import OrderedDict
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
        csv = CsvAccessor(gwas_file_path)
        logging.debug("Standardizing column names")
        csv.standardize_column_names(column_names)
        logging.debug("Organizing data")

        # organizing is implemente by array iteration in numpy for efficiency
        with tqdm(total = csv.get_data().shape[0], file = sys.stdout, leave=True) as pbar:
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
                self.__gwas[chromosome_vector[index]].update({int(position_vector[index]): GwasRow(values)})
                pbar.set_description('Organizing records')
                pbar.update(1)
    
    def __get_size(self):
        size = 0
        for(positions) in self.__gwas.values():
            size += len(positions)
        return size

    def clump(self):
        print("clumping")
        clumped_data = []
        print(self.__get_size())
        with tqdm(total = self.__get_size(), file = sys.stdout, leave=False) as pbar:
            for chromosome, positions in self.__gwas.items():
                print("chromosome: " + str(chromosome))
                p_sequence = np.array([])
                pos_sequence = np.array([])
                for position in positions.values():
                    pbar.set_description('Clumping records')
                    pbar.update(1)
                    p_sequence = np.append(p_sequence, -math.log(position.get("pvalue"),10))
                    pos_sequence = np.append(pos_sequence, position.get("position"))
                # print("hello")
                # print(str(p_sequence[0:10]))
                extrema = argrelextrema(p_sequence, np.greater_equal, order=int(len(positions) / 10))[0]
                for extremum in extrema:
                    print("extremum: " + str(extremum))
                    print("position: " + str(pos_sequence[extremum]))
                    print(str(positions[pos_sequence[extremum]]))


        return clumped_data

    def __get_filtered_data(self, pvalue_threshold: float = 0.05):
        """
        return a filtered gwas data
        """
        filtered_data = []
        
        with tqdm(total = self.__get_size(), file = sys.stdout, leave=False) as pbar:
            for chromosome, positions in self.__gwas.items():
                for position, gwas_record in positions.items():
                    if gwas_record.get("pvalue") < pvalue_threshold:
                        filtered_data.append(gwas_record)
            pbar.set_description('Filtering records')
            pbar.update(1)            
        return filtered_data

    def __get_filtered_data_for_manhattan_plot(self, evry_nth: int = 2, pvalue_threshold: float = 0.05):
        """
        return a filtered gwas data for manhattan plot
        """
        filtered_data = []
        
        with tqdm(total = self.__get_size(), file = sys.stdout, leave=False) as pbar:
            for chromosome, positions in self.__gwas.items():
                for position, gwas_record in positions.items():
                    if gwas_record.get("pvalue") < pvalue_threshold and int(position) % evry_nth == 0:
                        filtered_data.append({
                            "pvalue": gwas_record.get("pvalue"),
                            "chromosome": gwas_record.get("chromosome"),
                            "position": gwas_record.get("position")
                        })
            pbar.set_description('Filtering records')
            pbar.update(1)
        filtered_data = pd.DataFrame(filtered_data)
        return filtered_data

    def plot_manhattan(self):
        """
        plot manhattan plot
        """
        
        data = self.__get_filtered_data_for_manhattan_plot()
        chromsizes = Chromsizes().chromsizes["GRCh37"]
        cumulative_chromsizes = Chromsizes().chromsizes["GRCh37" + "_cumulative"]
        # get summary length of all chromosomes
        summary_length = sum(chromsizes.values())

        print(str(data.iloc[0:3]))

        # add cumulative position column to data
        data['cumulative_position'] = data.apply(lambda row: int(row["position"]) + cumulative_chromsizes[str(int(row["chromosome"]))], axis=1)
        # # add log10 pvalue column to data
        data['log10_pvalue'] = data.apply(lambda row: -math.log10(row["pvalue"]), axis=1)
        # add different color to every second chromosome
        data['color'] = data.apply(lambda row: '0' if row["chromosome"] in ['X', 'Y'] or int(row["chromosome"]) % 2 == 1 else '1', axis=1)
        data['color'] = data.apply(lambda row: row['color'] if row['log10_pvalue'] < 8 else '2', axis=1)
        data['color'] = data.apply(lambda row: row['color'] if row['log10_pvalue'] < 20 else '3', axis=1)
        data['size'] = data.apply(lambda row: 0.5 if row['log10_pvalue'] < 8 else 1.5 if row['log10_pvalue'] > 18 else (row['log10_pvalue'] - 8) / 10 + 0.5, axis=1)
    
    
        color_dict = {'0': colors.grey, 
                    '1': colors.grey_dark, 
                    '2': colors.teal, 
                    '3': colors.teal_dark, 
                    '4': colors.teal_darker}

        plot = (
            ggplot(data, aes('cumulative_position', 'log10_pvalue', color = 'color', size = 'size')) + 
            geom_point() +
            geom_hline(yintercept = 8, color = colors.grey, size = 0.5, linetype = 'dashed') +
            scale_color_manual(values=color_dict) +
            scale_size_continuous(range=(0.5,2)) + 
            theme(
                legend_position = "none",
                line = element_line(color = colors.grey, size = 1),
                rect = element_rect(fill = colors.grey_lighter),
                panel_grid_major = element_blank(),
                panel_grid_minor = element_blank(),
                panel_border = element_blank(),
                panel_background = element_blank(),
                axis_line = element_blank(),
                axis_title_x = element_blank(),
                axis_title_y = element_blank(),
                axis_text_x = element_blank(),
                axis_text_y = element_blank(),
                axis_ticks = element_line(size = 1),
                axis_ticks_length = 5,
                axis_ticks_length_minor = 2,
                axis_ticks_major_x = element_blank(),
                axis_ticks_major_y = element_line(color = colors.grey_light),
                axis_ticks_minor_y = element_line(color = colors.grey_light),
            )
        )

        plot.save("/tmp/plot.png", height=6, width=8)



