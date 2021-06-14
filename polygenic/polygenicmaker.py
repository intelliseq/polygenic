import argparse
import logging

import sys
import os
import urllib.request

# biobankuk-index
import gzip
import io

# biobankuk-get
import progressbar
import os.path

logger = logging.getLogger('polygenicmaker')

def biobankuk_index(args):
    parser = argparse.ArgumentParser(description='polygenicmaker biobankuk-index downloads index of gwas results from pan.ukbb study')  # todo dodać opis
    parser.add_argument('--url', type=str, default='https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz', help='alternative url location for index')
    parser.add_argument('--output', type=str, default='', help='output directory')
    parsed_args = parser.parse_args(args)
    output_path = os.path.abspath(os.path.expanduser(parsed_args.output)) + "/phenotype_manifest.tsv"
    if os.path.isfile(output_path):
        print("Index exists: " + output_path)
        return
    print("Downloading from " + parsed_args.url)
    response = urllib.request.urlopen(parsed_args.url)
    decompressed_file = gzip.GzipFile(fileobj=response)
    with open(output_path, 'w') as outfile:
        outfile.write(str(decompressed_file.read(), 'utf-8'))
    return

def biobankuk_get(args):
    parser = argparse.ArgumentParser(description='polygenicmaker biobankuk-get downloads specific gwas result from pan.ukbb study')  # todo dodać opis
    parser.add_argument('--index', type=str, default='phenotype_manifest.tsv', help='path to phenotype_manifest.tsv index file. Can be downloaded using polygenicmaker biobankuk-index command')
    parser.add_argument('--phenocode', type=str, required=True, help='biobankUK phenotype code. Example 30600')
    parser.add_argument('--output', type=str, default='', help='output directory')
    parser.add_argument('--force', action='store_true', help='overwrite downloaded file')
    parsed_args = parser.parse_args(args)
    # checking index file for download url
    with open(parsed_args.index, 'r') as indexfile:
        firstline = indexfile.readline()
        phenocode_colnumber = firstline.split('\t').index("phenocode")
        aws_link_colnumber = firstline.split('\t').index("aws_link")
        while True:
            line = indexfile.readline()
            if not line:
                break
            if line.split('\t')[phenocode_colnumber] != parsed_args.phenocode:
                continue
            url = line.split('\t')[aws_link_colnumber]
            break
    # downloading
    if not url is None:
        logger.info("Downloading from " + url)
        output_directory = os.path.abspath(os.path.expanduser(parsed_args.output))
        output_file_name = os.path.splitext(os.path.basename(url))[0]
        output_path = output_directory + "/" + output_file_name
        print(parsed_args.force)
        if os.path.isfile(output_path) and parsed_args.force is False:
            print("File is laready downloaded")
            return
        logger.info("Saving to " + output_path)
        response = urllib.request.urlopen(url)
        file_size = 3.5 * int(response.getheader('Content-length'))
        decompressed_file = gzip.GzipFile(fileobj=response)
        if file_size is None:
            file_size = 7078686639
        else:
            bar = progressbar.ProgressBar(max_value = file_size).start()
            downloaded = 0
            with open(output_path, 'w') as outfile:
                while (bytes := decompressed_file.read(1024)):
                    outfile.write(str(bytes, 'utf-8'))
                    downloaded = downloaded + 1024
                    bar.update(min(downloaded, file_size))
            bar.update(file_size)
            bar.finish()
    return

def main(args = sys.argv[1:]):
    try:
        if args[0] == 'biobankuk-index':
            biobankuk_index(args[1:])
        elif args[0] == 'biobankuk-get':
            biobankuk_get(args[1:])
        else:
            raise Exception()
    except Exception as e:
        print("ERROR " + str(e))
        print("""
        Program: polygenicmaker (downloads gwas data, clumps and build polygenic scores)
        Contact: Marcin Piechota <piechota.marcin@gmail.com>

        Usage:   polygenicmaker <command> [options]
\
        Command: 
        biobankuk-index     downloads pan biobankuk index of gwas results

        """)

if __name__ == '__main__':
    main(sys.argv[1:])
