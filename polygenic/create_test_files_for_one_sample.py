import gzip
import os
import argparse

__version__ = '0.0.1'


def get_first_sample_from_vcf_and_write_it_to_new_vcf(old_vcf_path: str, new_vcf_path: str, compressed=False,
                                                      separator='\t'):
    print("Processing {}".format(old_vcf_path))
    if compressed:
        f_in = gzip.open(old_vcf_path, 'rt')
        f_out = gzip.open(new_vcf_path, 'wt')
    else:
        f_in = open(old_vcf_path)
        f_out = open(new_vcf_path, 'wt')
    try:
        for line in f_in:
            if line.startswith('##'):
                f_out.write(line)
            else:
                f_out.write(separator.join(line.split(separator)[:10]))
                f_out.write('\n')
    finally:
        f_in.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Prepares input files for one (first) sample based on input files with many samples')
    parser.add_argument('--imputed_files_in', nargs='+', required=True,
                        help='Compressed impoted vcfs with many samples')
    parser.add_argument('--curated_initial_vcf_in', type=str, required=True,
                        help='Curated initial experimental file, for many samples')
    parser.add_argument('--imputed_directory_out', required=True,
                        help='Compressed impoted vcfs with one sample')
    parser.add_argument('--curated_initial_vcf_out', type=str, required=True,
                        help='Curated initial experimental file, for one sample')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))

    args = parser.parse_args()

    get_first_sample_from_vcf_and_write_it_to_new_vcf(args.curated_initial_vcf_in, args.curated_initial_vcf_out)

    if not os.path.isfile(args.imputed_directory_out):
        raise Exception('{} should be an empty directory'.format(args.imputed_directory_out))

    if not os.path.exists(args.imputed_directory_out):
        os.mkdir(args.imputed_directory_out)

    for f_path in args.imputed_files_in:
        f_name = os.path.basename(f_path)
        get_first_sample_from_vcf_and_write_it_to_new_vcf(f_path, os.path.join(args.imputed_directory_out, f_name),
                                                          compressed=True)
