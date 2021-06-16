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

# bioabnkuk-build-model
## clumping
import subprocess
import re
import random
import statistics

logger = logging.getLogger('polygenicmaker')

#######################
### biobankuk-index ###
#######################

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

#######################
### biobankuk-get ###
#######################

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

#############################
### biobankuk-build-model ###
#############################

def biobankuk_build_model(args):
    parser = argparse.ArgumentParser(description='polygenicmaker biobankuk-build-model constructs polygenic score model based on p value data')  # todo dodać opis
    parser.add_argument('--data', type=str, required=True, help='path to biomarkers file from biobank uk. It can be downloaded using biobankuk-get')
    parser.add_argument('--output', type=str, default='', help='output directory')
    parser.add_argument('--anno', type=str, required=True, help='path to annotation file. It can be downloaded with biobank-get-anno')
    parser.add_argument('--pop', type=str, default='meta', help='population: meta, AFR, AMR, CSA, EAS, EUR, MID')
    parser.add_argument('--threshold', type=float, default=1e-08, help='population: meta, AFR, AMR, CSA, EAS, EUR, MID')
    parser.add_argument('--iterations', type=float, default=1000, help='simulation iterations for mean and sd')
    parsed_args = parser.parse_args(args)
    if not os.path.isdir(parsed_args.output):
        print("ERROR: " + parsed_args.output + " does not exists or is not directory")
        return
    if not os.path.isfile(parsed_args.data):
        print("ERROR: " + parsed_args.data + " does not exists or is not a file")
        return
    if not os.path.isfile(parsed_args.anno):
        print("ERROR: " + parsed_args.anno + " does not exists or is not a file")
        return
    #filter_pval(parsed_args)
    #clump(parsed_args)
    simulation_results = simulate(parsed_args)
    description = {
        'mean': simulation_results.mean,
        'sd': simulation_results.sd
    }
    save_model(parsed_args, description)
    return

def filter_pval(args):
    output_path = args.output + "/" + os.path.basename(args.data) + ".filtered"
    with open(args.data, 'r') as data, open(args.anno, 'r') as anno, open(output_path, 'w') as output:
        data_header = data.readline().rstrip().split('\t')
        anno_header = anno.readline().rstrip().split('\t')
        output.write('\t'.join(data_header + anno_header) + "\n")
        while True:
            try:
                data_line = data.readline().rstrip().split('\t')
                anno_line = anno.readline().rstrip().split('\t')
                if float(data_line[data_header.index('pval_' + args.pop)].replace('NA','1',1)) <= args.threshold:
                    output.write('\t'.join(data_line + anno_line) + "\n")
            except:
                break
    return

def clump(args):
    filtered_path = args.output + "/" + os.path.basename(args.data) + ".filtered"
    subprocess.call("plink" + 
        " --clump " + filtered_path + 
        " --clump-p1 " + str(args.threshold) +
        " --clump-r2 0.25 " + 
        " --clump-kb 1000 " + 
        " --clump-snp-field rsid " +
        " --clump-field pval_" + args.pop +
        " --vcf results/eur.phase3.biobank.set.vcf.gz " +
        " --allow-extra-chr", 
        shell=True)
    clumped_rsids = []
    with open("plink.clumped", 'r') as plink_file:
        while(line := plink_file.readline()):
            if ' rs' in line:
                line = re.sub(' +', '\t', line).rstrip().split('\t')
                clumped_rsids.append(line[3])
    try:
        os.remove("plink.clumped")
        os.remove("plink.log")
        os.remove("plink.nosex")
    except:
        pass
    filtered_path = args.output + "/" + os.path.basename(args.data) + ".filtered"
    clumped_path = args.output + "/" + os.path.basename(args.data) + ".clumped"
    with open(filtered_path, 'r') as filtered_file, open(clumped_path, 'w') as clumped_file:
        filtered_header = filtered_file.readline().rstrip().split('\t')
        clumped_file.write('\t'.join(filtered_header) + "\n")
        while True:
            try:
                filtered_line = filtered_file.readline().rstrip().split('\t')
                if filtered_line[filtered_header.index('rsid')] in clumped_rsids:
                    clumped_file.write('\t'.join(filtered_line) + "\n")
            except:
                break
    return

def simulate(args):
    clumped_path = args.output + "/" + os.path.basename(args.data) + ".clumped"
    random.seed(0)
    simulation_data = []
    with open(clumped_path, 'r') as clumped_file:
        clumped_header = clumped_file.readline().rstrip().split('\t')
        clumped_line = clumped_header
        while True:
            clumped_line = clumped_file.readline().rstrip().split('\t')
            if len(clumped_line) < 2:
                break
            rsid = clumped_line[clumped_header.index('rsid')]
            af = float(clumped_line[clumped_header.index('af_' + args.pop)])
            beta = float(clumped_line[clumped_header.index('beta_' + args.pop)])
            simulation_data.append({'rsid': rsid, 'af': af, 'beta': beta})

    randomized_beta_list = []
    for _ in range(args.iterations):
        randomized_beta_list.append(sum(map(lambda snp: randomize_beta(snp['beta'], snp['af']), simulation_data)))
    return {'mean': statistics.mean(randomized_beta_list), 'sd': statistics.stdev(randomized_beta_list)}

def randomize_beta(beta: float, af: float):
    first_allele_beta = beta if random.uniform(0, 1) < af else 0
    second_allele_beta = beta if random.uniform(0, 1) < af else 0
    return first_allele_beta + second_allele_beta

def save_model(args, description):
    return

def main(args = sys.argv[1:]):
    try:
        if args[0] == 'biobankuk-index':
            biobankuk_index(args[1:])
        elif args[0] == 'biobankuk-get':
            biobankuk_get(args[1:])
        elif args[0] == 'biobankuk-build-model':
            biobankuk_build_model(args[1:])
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
        biobankuk-index         downloads pan biobankuk index of gwas results
        biobankuk-get           downloads gwas results for given phenocode
        biobankuk-build-model   build polygenic score based on gwas results

        """)

if __name__ == '__main__':
    main(sys.argv[1:])
