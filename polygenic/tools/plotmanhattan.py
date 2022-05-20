import logging
import math
import csv
import pandas as pd

from polygenic.tools.data.chromsizes import Chromsizes
from polygenic.tools.data.colors import Colors as colors
from polygenic.error.polygenic_exception import PolygenicException

from plotnine import ggplot
from plotnine import geom_point, geom_hline
from plotnine import aes, theme, theme_void, theme_minimal, scale_fill_manual, scale_color_manual, scale_size_continuous
from plotnine import element_rect, element_line
from plotnine import *

def run(args):
    # read csv file 
    with open(args.tsv, 'r') as file:
        line = file.readline()
        # pass through header
        while line[0] == '#':
            line = file.readline()
        # read column names
        column_names = next(csv.reader([line.rstrip('\r\n')], delimiter = args.delimiter))
        pvalue_column_index = column_names.index(args.pvalue_column)
        chromosome_column_index = column_names.index(args.chromosome_column)
        position_column_index = column_names.index(args.position_column)
        # create empty dataframe with pvalue, chromosome and position columns
        data = pd.DataFrame(columns = [args.pvalue_column, args.chromosome_column, args.position_column])
        # row index
        row_index = 0
        while True:
            line = next(csv.reader([file.readline().rstrip('\r\n')], delimiter = args.delimiter))
            if len(line) < 2:
                break
            if not len(column_names) == len(line):
                logging.error("Line and first row have different lengths")
                raise PolygenicException("Line and header have different lengths. LineL {line}".format(line = str(line)))
            
            data.append({
                "position": line[position_column_index],
                "chromosome": line[chromosome_column_index],
                "pvalue": line[pvalue_column_index]
            }, ignore_index = True)

    # extract data
    data = tsv.get_data()
    # get chromosome sizes
    chromsizes = Chromsizes().chromsizes[args.genome_version]
    cumulative_chromsizes = Chromsizes().chromsizes[args.genome_version + "_cumulative"]
    # get summary length of all chromosomes
    summary_length = sum(chromsizes.values())
    # add cumulative position column to data
    data['cumulative_position'] = data.apply(lambda row: int(row[args.position_column]) + cumulative_chromsizes[row[args.chromosome_column]], axis=1)
    # add log10 pvalue column to data
    data['log10_pvalue'] = data.apply(lambda row: -math.log10(row[args.pvalue_column]), axis=1)
    # add different color to every second chromosome
    data['color'] = data.apply(lambda row: '0' if row[args.chromosome_column] in ['X', 'Y'] or int(row[args.chromosome_column]) % 2 == 1 else '1', axis=1)
    data['color'] = data.apply(lambda row: row['color'] if row['log10_pvalue'] < 8 else '2', axis=1)
    data['color'] = data.apply(lambda row: row['color'] if row['log10_pvalue'] < 20 else '3', axis=1)
    data['size'] = data.apply(lambda row: 0.5 if row['log10_pvalue'] < 8 else 1.5 if row['log10_pvalue'] > 18 else (row['log10_pvalue'] - 8) / 10 + 0.5, axis=1)
    logging.info(summary_length)
    
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

    plot = plot.save(args.output, height=6, width=8)

    return 0
