from calendar import c
import sys
import os

import polygenic.tools.utils as utils
import polygenic.data.csv_accessor as csv_accessor

def run(args):

    csv = csv_accessor.CsvAccessor(args.gwas_file)
    rsid_idx = csv.get_rsid_column_index(rsid_column_name=args.rsid_column)
    chr_idx = csv.get_

    print(str(csv.get_column_names()))
    print(str())
