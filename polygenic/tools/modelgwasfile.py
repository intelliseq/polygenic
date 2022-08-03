import sys
import os

import polygenic.tools.utils as utils
import polygenic.data.csv_accessor as csv_accessor

def run(args):

    csv = csv_accessor.CsvAccessor(args.gwas_file)
    print(str(csv.get_column_names()))