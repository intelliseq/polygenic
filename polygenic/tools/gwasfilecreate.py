from polygenic.data.gwas import Gwas
from calendar import c
import sys
import os

import polygenic.tools.utils as utils
import polygenic.data.csv_accessor as csv_accessor

def run(args):
    gwas = Gwas()
    gwas.load(args.input, vars(args))
    gwas.save(args.output)