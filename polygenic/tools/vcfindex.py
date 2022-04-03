import argparse
import sys
import logging

from polygenic.version import __version__ as version
from polygenic.tools import utils

from polygenic.data.vcf_accessor import VcfAccessor
from polygenic.error.polygenic_exception import PolygenicException

def parse_args(args):
    parser = argparse.ArgumentParser(description='vcf-index creates index for vcf file')
    parser.add_argument('-i', '--vcf', required=True, help='path to vcf.gz')
    parser.add_argument('-l', '--log-file', type=str, default='$HOME/.pgstk/log/pgstk.log', help='path to log file')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version)
    parser.add_argument('--print', default=False, action='store_true', help='print output to stdout')
    parsed_args = parser.parse_args(args)
    return parsed_args

def run(args):
    logging.getLogger().debug("vcf-index")
    VcfAccessor(args.vcf)
    return

def main(args = sys.argv[1:]):

    args = parse_args(args) 
    try:
        run(args)
    except PolygenicException as e:
        utils.error_exit(e)

if __name__ == '__main__':
    main(sys.argv[1:])
