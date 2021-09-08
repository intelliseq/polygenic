import argparse
import logging

import sys
import os
import urllib.request

import configparser
import tabix

# utils
# simulate
from polygenic.lib.data_access.vcf_accessor import VcfAccessor
from polygenic.lib.data_access.vcf_accessor import DataNotPresentError
from polygenic.lib.utils import is_valid_path
from polygenic.lib.utils import download
from polygenic.lib.utils import read_header
from polygenic.lib.utils import read_table

# clumping
import subprocess
import re
# simlating
import random
import statistics
# saving
import json

logger = logging.getLogger('polygenicmaker')
config = configparser.ConfigParser()
config.read(os.path.dirname(__file__) + "/../polygenic/polygenic.cfg")

####################
###   add af     ###
####################

def add_annotation(
      dict, 
      vcf: VcfAccessor, 
      column_name: str = 'rsid', 
      new_name: str = None, 
      backup_column: str = None,
      default = ""):
    print(str(dict))
    if new_name is None:
        new_name = column_name
    if backup_column is None:
        backup_column = column_name
    annotation_dict = None
    if annotation_dict is None and "CHR" in dict:
        annotation_dict = vcf.get_record_by_position(dict['CHR'], dict['POS'])
    if annotation_dict is None and "CHROM" in dict:
        annotation_dict = vcf.get_record_by_position(dict['CHROM'], dict['POS'])
    if annotation_dict is None and "chr_name" in dict:
        annotation_dict = vcf.get_record_by_position(dict['chr_name'], dict['chr_position'])
    if annotation_dict is None and "ID" in dict:
        annotation_dict = vcf.get_record_by_rsid(dict['ID'])
    if annotation_dict is None and "rsid" in dict:
        annotation_dict = vcf.get_record_by_rsid(dict['rsid'])
    if not annotation_dict is None:
        annotation_dict = annotation_dict.get_as_dict()
    if not annotation_dict is None and column_name in annotation_dict:
        print(str(annotation_dict[column_name]))
        dict[new_name] = annotation_dict[column_name]
    elif backup_column in dict:
        dict[new_name] = dict[backup_column]
    else:
        dict[new_name] = default
    return dict

def add_af(line, af_accessor: VcfAccessor, population: str = 'nfe', rsid_column_name: str = 'rsid'):
    try:
        af = af_accessor.get_af_by_pop(line['ID'], 'AF_' + population)
        line['af'] = af[line['ALT']]
    except DataNotPresentError:
        line['af'] = 0
    except Exception:
        # TODO convert alternate alleles
        af = af_accessor.get_af_by_pop(line['ID'], 'AF_' + population)
        line['af'] = af[list(af.keys())[1]]
    return line

####################
###   add rsid   ###
####################


def add_rsid(line, tabix_source: open, population: str = 'nfe', rsid_column_name: str = 'ID'):
    if "rs" in line[rsid_column_name]:
        line['rsid'] = line[rsid_column_name]
    else:
        try:
            records = tabix_source.query(
                line['CHROM'], int(line['POS']) - 1, int(line['POS']))
            for record in records:
                line['rsid'] = record[2]
                print(str(line))
                break
        except Exception:
            line['rsid'] = line[rsid_column_name]

    return line

####################
###   simulate   ###
####################


def simulate_parameters(data, iterations: int = 1000, coeff_column_name: str = 'BETA'):
    random.seed(0)

    randomized_beta_list = []
    for _ in range(iterations):
        randomized_beta_list.append(sum(map(lambda snp: randomize_beta(
            float(snp[coeff_column_name]), snp['af']), data)))
    minsum = sum(map(lambda snp: min(float(snp[coeff_column_name]), 0), data))
    maxsum = sum(map(lambda snp: max(float(snp[coeff_column_name]), 0), data))
    return {
        'mean': statistics.mean(randomized_beta_list),
        'sd': statistics.stdev(randomized_beta_list),
        'min': minsum,
        'max': maxsum
    }

###################
### write model ###
###################


def write_model(data, description, destination):

    with open(destination, 'w') as model_file:
        model_file.write(
            "from polygenic.seqql.score import PolygenicRiskScore\n")
        model_file.write("from polygenic.seqql.score import ModelData\n")
        model_file.write(
            "from polygenic.seqql.category import QuantitativeCategory\n")
        model_file.write("\n")
        model_file.write("model = PolygenicRiskScore(\n")
        model_file.write("categories=[\n")
        model_file.write("QuantitativeCategory(from_=" + str(description['min']) + ", to=" + str(
            description['mean'] - 1.645 * description['sd']) + ", category_name='Reduced'),\n")
        model_file.write("QuantitativeCategory(from_=" + str(description['mean'] - 1.645 * description['sd']) + ", to=" + str(
            description['mean'] + 1.645 * description['sd']) + ", category_name='Average'),\n")
        model_file.write("QuantitativeCategory(from_=" + str(description['mean'] + 1.645 * description['sd']) + ", to=" + str(
            description['max']) + ", category_name='Increased')\n")
        model_file.write("],\n")
        model_file.write("snips_and_coefficients={\n")
        snps = []
        for snp in data:
            snps.append("'" + snp['rsid'] + "': ModelData(effect_allele='" +
                        snp['ALT'] + "', coeff_value=" + snp['BETA'] + ")")
        model_file.write(",\n".join(snps))
        model_file.write("},\n")
        model_file.write("model_type='beta'\n")
        model_file.write(")\n")
        model_file.write("description = " + json.dumps(description, indent=4))

    return

#####################################################################################################
###                                                                                               ###
###                                   Utils                                                       ###
###                                                                                               ###
#####################################################################################################

#######################
### vcf-index #########
#######################


def vcf_index(args):
    parser = argparse.ArgumentParser(
        description='vcf index indexes vcf file for query')  # todo dodać opis
    parser.add_argument('--vcf', type=str, required=True,
                        help='path to vcf file')
    parsed_args = parser.parse_args(args)
    VcfAccessor(parsed_args.vcf)
    return


#####################################################################################################
###                                                                                               ###
###                                   Polygenic Score Catalogue                                   ###
###                                                                                               ###
#####################################################################################################

#######################
### pgs-index #########
#######################

def pgs_index(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker pgs-index downloads index of gwas results from Polgenic Score Catalogue')  # todo dodać opis
    parser.add_argument('--url', type=str, default='http://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/pgs_all_metadata_scores.csv',
                        help='alternative url location for index')
    parser.add_argument('--output', type=str, default='',
                        help='output directory')
    parsed_args = parser.parse_args(args)
    output_path = os.path.abspath(os.path.expanduser(
        parsed_args.output)) + "/pgs_manifest.tsv"
    download(parsed_args.url, output_path, force=True)
    return

#######################
### pgs-get ###########
#######################


def pgs_get(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker pgs-get downloads specific gwas result from polygenic score catalogue')  # todo dodać opis
    parser.add_argument('-c', '--code', type=str, required=False,
                        help='PGS score code. Example: PGS000814')
    parser.add_argument('-o', '--output-path', type=str,
                        default='', help='output directory')
    parser.add_argument('-f', '--force', action='store_true',
                        help='overwrite downloaded file')
    parsed_args = parser.parse_args(args)
    url = "http://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/" + \
        parsed_args.code + "/ScoringFiles/" + parsed_args.code + ".txt.gz"
    download(url=url, output_path=parsed_args.output_path,
             force=parsed_args.force, progress=True)
    return

#######################
### pgs-prepare #######
#######################


def pgs_prepare_model(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker pgs-prepare-model constructs polygenic score model')  # todo dodać opis
    parser.add_argument('-i', '--input', type=str,
                        required=True, help='path to PRS file from PGS.')
    parser.add_argument('-o', '--output-path', type=str,
                        required=True, help='path to output file.')
    parser.add_argument('--origin-reference-vcf', type=str,
                        required=True, help='path to rsid vcf.')
    parser.add_argument('--model-reference-vcf', type=str,
                        required=True, help='path to rsid vcf.')
    parser.add_argument('--af', type=str, required=True,
                        help='path to allele frequency vcf.')
    parser.add_argument('--pop', type=str, default='nfe',
                        help='population: meta, afr, amr, csa, eas, eur, mid')
    parser.add_argument('--iterations', type=float, default=1000,
                        help='simulation iterations for mean and sd')
    parsed_args = parser.parse_args(args)
    if not is_valid_path(parsed_args.input):
        return
    if not is_valid_path(parsed_args.origin_reference_vcf):
        return
    if not is_valid_path(parsed_args.model_reference_vcf):
        return
    if not is_valid_path(parsed_args.af):
        return
    data = read_table(parsed_args.input)
    af_vcf = VcfAccessor(parsed_args.af)
    origin_ref_vcf = VcfAccessor(parsed_args.origin_reference_vcf)
    model_ref_vcf = VcfAccessor(parsed_args.model_reference_vcf)

    clean_data = []
    for pgs_entry in data:
        origin_record = origin_ref_vcf.get_record_by_position(
            pgs_entry['chr_name'], pgs_entry['chr_position'])
        model_record = model_ref_vcf.get_record_by_rsid(origin_record.get_id())
        af_records = af_vcf.get_records_by_rsid(origin_record.get_id())
        af_record = None
        for record in af_records:
            if record.get_alt()[0] == pgs_entry['effect_allele']:
                af_record = record
                break
        if model_record is None and not af_record is None:
            model_record = af_record
        if origin_record is None or model_record is None:
            logger.warning("Variant {chromosome}:{position} is not present in reference.".format(
                chromosome=pgs_entry['chr_name'], position=pgs_entry['chr_position']))
            continue
        if not pgs_entry['reference_allele'] == origin_record.get_ref():
            logger.warning("Variant {chromosome}:{position} has mismatch nucleotide in reference.".format(
                chromosome=pgs_entry['chr_name'], position=pgs_entry['chr_position']))
            continue
        if not origin_record.get_ref() == model_record.get_ref():
            logger.warning("Variant {chromosome}:{position} has mismatch nucleotide between references (grch37 vs grch38).".format(
                chromosome=pgs_entry['chr_name'], position=pgs_entry['chr_position']))
            continue
        pgs_entry['rsid'] = origin_record.get_id()
        if af_record is None:
            pgs_entry['af'] = 0.001
        else:
            pgs_entry['af'] = af_record.get_af_by_pop(
                'AF_' + parsed_args.pop)[pgs_entry['effect_allele']]
        pgs_entry['ALT'] = pgs_entry['effect_allele']
        pgs_entry['BETA'] = pgs_entry['effect_weight']
        clean_data.append(pgs_entry)

    description = simulate_parameters(clean_data)
    description.update(read_header(parsed_args.input))
    write_model(clean_data, description, parsed_args.output_path)

    return

#####################################################################################################
###                                                                                               ###
###                                   Global Biobank Engine                                       ###
###                                                                                               ###
#####################################################################################################

#######################
### gbe-index #########
#######################


def gbe_index(args):
    print("GBEINDEX")
    parser = argparse.ArgumentParser(
        description='polygenicmaker gbe-index downloads index of gwas results from pan.ukbb study')  # todo dodać opis
    parser.add_argument('--url', type=str, default='https://biobankengine.stanford.edu/static/degas-risk/degas_n_977_traits.tsv',
                        help='alternative url location for index')
    parser.add_argument('--output', type=str, default='',
                        help='output directory')
    parsed_args = parser.parse_args(args)
    output_path = os.path.abspath(os.path.expanduser(
        parsed_args.output)) + "/gbe_phenotype_manifest.tsv"
    download(parsed_args.url, output_path)
    return

#######################
### gbe-get ###########
#######################


def gbe_get(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker biobankuk-get downloads specific gwas result from pan.ukbb study')  # todo dodać opis
    parser.add_argument('--code', type=str, required=False,
                        help='GBE phenotype code. Example: BIN1210')
    parser.add_argument('--output', type=str, default='',
                        help='output directory')
    parser.add_argument('--force', action='store_true',
                        help='overwrite downloaded file')
    parsed_args = parser.parse_args(args)
    url = "https://biobankengine.stanford.edu/static/PRS_map/" + parsed_args.code + ".tsv"
    output_directory = os.path.abspath(os.path.expanduser(parsed_args.output))
    output_file_name = os.path.splitext(os.path.basename(url))[0]
    output_path = output_directory + "/" + output_file_name
    download(url=url, output_path=output_path,
             force=parsed_args.force, progress=True)
    return

#######################
### gbe-prepare #######
#######################


def gbe_prepare_model(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker biobankuk-build-model constructs polygenic score model based on p value data')  # todo dodać opis
    parser.add_argument('--data', type=str, required=True,
                        help='path to PRS file from gbe. It can be downloaded using gbe-get')
    parser.add_argument('--output', type=str, default='',
                        help='output directory')
    parser.add_argument('--af', type=str, required=True,
                        help='path to allele frequency vcf. It can be downloaded with biobank-get-anno')
    parser.add_argument('--pop', type=str, default='nfe',
                        help='population: meta, AFR, AMR, CSA, EAS, EUR, MID')
    parser.add_argument('--iterations', type=float, default=1000,
                        help='simulation iterations for mean and sd')
    parsed_args = parser.parse_args(args)
    if not is_valid_path(parsed_args.output, is_directory=True):
        return
    if not is_valid_path(parsed_args.data):
        return
    if not is_valid_path(parsed_args.af):
        return
    data = read_table(parsed_args.data)
    af = VcfAccessor(parsed_args.af)
    rsid_vcf = VcfAccessor(config['urls']['hg19-rsids'])
    data = [line for line in data if "rs" in line['ID']]
    #data = [add_rsid(line, tabix_source) for line in data]
    #print(data[1])
    #add_annotation(data[1], rsid_vcf, 'ID')
    data = [add_annotation(line, rsid_vcf, 'ID', new_name='rsid', backup_column='ID') for line in data]
    #data = [add_af(line, af, population=parsed_args.pop) for line in data]
    #description = simulate_parameters(data)
    #model_path = parsed_args.output + "/" + \
    #    os.path.basename(parsed_args.data).split('.')[0] + ".py"
    #write_model(data, description, model_path)
    return


#######################
### biobankuk-index ###
#######################

def biobankuk_index(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker biobankuk-index downloads index of gwas results from pan.ukbb study')  # todo dodać opis
    parser.add_argument('--url', type=str, default='https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz',
                        help='alternative url location for index')
    parser.add_argument('--output', type=str, default='',
                        help='output directory')
    parsed_args = parser.parse_args(args)
    output_path = os.path.abspath(os.path.expanduser(
        parsed_args.output)) + "/panukbb_phenotype_manifest.tsv"
    download(parsed_args.url, output_path)
    return

#######################
### biobankuk-get #####
#######################


def biobankuk_get(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker biobankuk-get downloads specific gwas result from pan.ukbb study')  # todo dodać opis
    parser.add_argument('--index', type=str, default='panukbb_phenotype_manifest.tsv',
                        help='path to phenotype_manifest.tsv index file. Can be downloaded using polygenicmaker biobankuk-index command')
    parser.add_argument('--phenocode', type=str, required=False,
                        help='biobankUK phenotype code. Example: 30600')
    parser.add_argument('--pheno_sex', type=str, default='both_sexes',
                        help='biobankUK pheno_sex code. Example: both_sexes')
    parser.add_argument('--coding', type=str, required=True,
                        help='biobankUK phenotype code. Example 30600')
    parser.add_argument('--modifier', type=str, required=True,
                        help='biobankUK phenotype code. Example 30600')
    parser.add_argument('--output', type=str, default='',
                        help='output directory')
    parser.add_argument('--force', action='store_true',
                        help='overwrite downloaded file')
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
        output_directory = os.path.abspath(
            os.path.expanduser(parsed_args.output))
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
            bar = progressbar.ProgressBar(max_value=file_size).start()
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


def biobankuk_prepare_model(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker biobankuk-build-model constructs polygenic score model based on p value data')  # todo dodać opis
    parser.add_argument('--data', type=str, required=True,
                        help='path to biomarkers file from biobank uk. It can be downloaded using biobankuk-get')
    parser.add_argument('--h2', type=str, help='')
    parser.add_argument('--output', type=str, default='',
                        help='output directory')
    parser.add_argument('--anno', type=str, required=True,
                        help='path to annotation file. It can be downloaded with biobank-get-anno')
    parser.add_argument('--pop', type=str, default='meta',
                        help='population: meta, AFR, AMR, CSA, EAS, EUR, MID')
    parser.add_argument('--threshold', type=float, default=1e-08,
                        help='population: meta, AFR, AMR, CSA, EAS, EUR, MID')
    parser.add_argument('--iterations', type=float, default=1000,
                        help='simulation iterations for mean and sd')
    parsed_args = parser.parse_args(args)
    if not os.path.isdir(parsed_args.output):
        print("ERROR: " + parsed_args.output +
              " does not exists or is not directory")
        return
    if not os.path.isfile(parsed_args.data):
        print("ERROR: " + parsed_args.data +
              " does not exists or is not a file")
        return
    if not os.path.isfile(parsed_args.anno):
        print("ERROR: " + parsed_args.anno +
              " does not exists or is not a file")
        return
    # filter_pval(parsed_args)
    # clump(parsed_args)
    simulation_results = simulate(parsed_args)
    description = {
        'mean': simulation_results['mean'],
        'sd': simulation_results['sd'],
        'min': simulation_results['min'],
        'max': simulation_results['max'],
        'population': parsed_args.pop
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
                if float(data_line[data_header.index('pval_' + args.pop)].replace('NA', '1', 1)) <= args.threshold:
                    output.write('\t'.join(data_line + anno_line) + "\n")
            except:
                break
    return


def clump(args):
    filtered_path = args.output + "/" + \
        os.path.basename(args.data) + ".filtered"
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
    filtered_path = args.output + "/" + \
        os.path.basename(args.data) + ".filtered"
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
            beta = float(
                clumped_line[clumped_header.index('beta_' + args.pop)])
            simulation_data.append({'rsid': rsid, 'af': af, 'beta': beta})

    randomized_beta_list = []
    for _ in range(args.iterations):
        randomized_beta_list.append(
            sum(map(lambda snp: randomize_beta(snp['beta'], snp['af']), simulation_data)))
    minsum = sum(map(lambda snp: min(snp['beta'], 0), simulation_data))
    maxsum = sum(map(lambda snp: max(snp['beta'], 0), simulation_data))
    return {
        'mean': statistics.mean(randomized_beta_list),
        'sd': statistics.stdev(randomized_beta_list),
        'min': minsum,
        'max': maxsum
    }


def randomize_beta(beta: float, af: float):
    first_allele_beta = beta if random.uniform(0, 1) < af else 0
    second_allele_beta = beta if random.uniform(0, 1) < af else 0
    return first_allele_beta + second_allele_beta


def save_model(args, description):
    model_path = args.output + "/" + \
        os.path.basename(args.data).split('.')[0] + ".py"
    with open(model_path, 'w') as model_file:
        model_file.write(
            "from polygenic.seqql.score import PolygenicRiskScore\n")
        model_file.write("from polygenic.seqql.score import ModelData\n")
        model_file.write(
            "from polygenic.seqql.category import QuantitativeCategory\n")
        model_file.write("\n")
        model_file.write("model = PolygenicRiskScore(\n")
        model_file.write("categories=[\n")
        model_file.write("QuantitativeCategory(from_=" + str(description['min']) + ", to=" + str(
            description['mean'] - 1.645 * description['sd']) + ", category_name='Reduced'),\n")
        model_file.write("QuantitativeCategory(from_=" + str(description['mean'] - 1.645 * description['sd']) + ", to=" + str(
            description['mean'] + 1.645 * description['sd']) + ", category_name='Average'),\n")
        model_file.write("QuantitativeCategory(from_=" + str(description['mean'] + 1.645 * description['sd']) + ", to=" + str(
            description['max']) + ", category_name='Increased')\n")
        model_file.write("],\n")
        model_file.write("snips_and_coefficients={\n")
        clumped_path = args.output + "/" + \
            os.path.basename(args.data) + ".clumped"
        snps = []
        with open(clumped_path, 'r') as clumped_file:
            clumped_header = clumped_file.readline().rstrip().split('\t')
            while True:
                clumped_line = clumped_file.readline().rstrip().split('\t')
                if len(clumped_line) < 2:
                    break
                rsid = clumped_line[clumped_header.index('rsid')]
                allele = clumped_line[clumped_header.index('alt')]
                beta = str(
                    float(clumped_line[clumped_header.index('beta_' + args.pop)]))
                snps.append("'" + rsid + "': ModelData(effect_allele='" +
                            allele + "', coeff_value=" + beta + ")")
        model_file.write(",\n".join(snps))
        model_file.write("},\n")
        model_file.write("model_type='beta'\n")
        model_file.write(")\n")
        model_file.write("description = " + json.dumps(description, indent=4))

    return


def main(args=sys.argv[1:]):
    try:
        if args[0] == 'biobankuk-index':
            biobankuk_index(args[1:])
        elif args[0] == 'biobankuk-get':
            biobankuk_get(args[1:])
        elif args[0] == 'biobankuk-prepare':
            biobankuk_prepare_model(args[1:])
        elif args[0] == 'gbe-index':
            gbe_index(args[1:])
        elif args[0] == 'gbe-get':
            gbe_get(args[1:])
        elif args[0] == 'gbe-prepare':
            gbe_prepare_model(args[1:])
        elif args[0] == 'pgs-index':
            pgs_index(args[1:])
        elif args[0] == 'pgs-get':
            pgs_get(args[1:])
        elif args[0] == 'pgs-prepare':
            pgs_prepare_model(args[1:])
        elif args[0] == 'vcf-index':
            vcf_index(args[1:])
        else:
            print('ERROR: Please select proper tool name"')
            print("""
            Program: polygenicmaker (downloads gwas data, clumps and build polygenic scores)
            Contact: Marcin Piechota <piechota@intelliseq.com>
            Usage:   polygenicmaker <command> [options]

            Command:
            biobankuk-index         downloads pan biobankuk index of gwas results
            biobankuk-get           downloads gwas results for given phenocode
            biobankuk-build-model   build polygenic score based on gwas results

            """)
    except RuntimeError as e:
        print("ERROR: " + str(e))

if __name__ == '__main__':
    main(sys.argv[1:])
