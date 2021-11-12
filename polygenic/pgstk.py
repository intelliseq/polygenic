import argparse
import logging

import sys
import os
import urllib.request

import configparser
import tabix

import yaml 

# utils
# simulate
from polygenic.data.vcf_accessor import VcfAccessor
from polygenic.data.vcf_accessor import DataNotPresentError
from polygenic.core.utils import is_valid_path
from polygenic.core.utils import download
from polygenic.core.utils import read_header
from polygenic.core.utils import read_table
from polygenic.core.trial import PolygenicException
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

def validate(
    validated_line: dict,
    validation_source: VcfAccessor,
    invert_field: str = None):
    record = validation_source.get_record_by_rsid(validated_line['rsid'])
    if record is None:
        print("WARNING: Failed validation for " + validated_line['rsid'] + ". SNP not present in validation vcf.")
        validated_line["status"] = "WARNING: snp not present"
        return validated_line
    if not (validated_line['REF'] == record.get_ref()): 
        if (validated_line['REF'] == record.get_alt()[0] and validated_line['ALT'] == record.get_ref()):
            ref = validated_line['REF']
            alt = validated_line['ALT']
            validated_line['REF'] = alt
            validated_line['ALT'] = ref
            if invert_field is not None:
                validated_line[invert_field] = - float(validated_line[invert_field])
            print("WARNING: " + "Failed validation for " + validated_line['rsid'] + ". REF and ALT do not match. " + record.get_ref() + "/" + str(record.get_alt()) + " succesful invert!")
            validated_line["status"] = "WARNING: ref alt inverted"
            return validated_line
        else:
            print("ERROR: " + "Failed validation for " + validated_line['rsid'] + ". REF and ALT do not match. " + record.get_ref() + "/" + str(record.get_alt()))
            validated_line["status"] = "WARNING: ref alt do not match"
            return validated_line
    validated_line["status"] = "SUCCESS"
    return validated_line

def add_annotation(
    annotated_line: dict,
    annotation_name: str,
    annotation_source: VcfAccessor,
    annotation_source_field: str,
    default_value,
    rsid_field: str = 'rsid'):

    rsid = annotated_line[rsid_field]
    record = annotation_source.get_record_by_rsid(annotated_line[rsid_field])
    annotated_line[annotation_name] = default_value
    if record is None:
        print("WARNING: Failed annotation for " + rsid + ". No record in source file.")
        return annotated_line
    if not (annotated_line['REF'] == record.get_ref()) or not (annotated_line['ALT'] in record.get_alt()): 
        print("WARNING: Failed annotation for " + rsid + ". REF and ALT do not match. " + record.get_ref() + "/" + str(record.get_alt()))
        return annotated_line
    if annotation_source_field == "id" or annotation_source_field == "rsid":    
        info = record.get_id()
    else:
        info = record.get_info_field(annotation_source_field)
    if info is None:
        print("WARNING: Failed annotation for " + rsid + ". No " + annotation_source_field + " in source.")
        return annotated_line    
    annotated_line[annotation_name] = info
    return annotated_line

def add_gene_symbol(line, reference, genes):
    record = reference.get_record_by_rsid(line["rsid"])
    if not record: return None
    symbol = [(gene["symbol"], abs(int(record.get_pos()) - int(gene["start"]))) for gene in genes if record.get_chrom() == ("chr" + gene["chromosome"]) and abs(int(record.get_pos()) - int(gene["start"])) < 500000]
    if not symbol: return None
    symbol = sorted(symbol, key=lambda tup: tup[1])[0][0]
    return symbol

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
###   simulate   ###
####################


def simulate_parameters(data, iterations: int = 1000, coeff_column_name: str = 'BETA'):
    random.seed(0)

    randomized_beta_list = []
    for _ in range(iterations):
        randomized_beta_list.append(sum(map(lambda snp: randomize_beta(
            float(snp[coeff_column_name]), float(snp['af'])), data)))
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

        categories = dict()
        borders = [
            description["parameters"]['mean'] - 1.645 * description["parameters"]['sd'],
            description["parameters"]['mean'] + 1.645 * description["parameters"]['sd']
        ]
        categories["reduced"] = {"from": description["parameters"]['min'], "to": borders[0]}
        categories["average"] = {"from": borders[0], "to": borders[1]}
        categories["increased"] = {"from": borders[1], "to": description["parameters"]['max']}
        
        variants = dict()
        for snp in data:
            variant = dict()
            variant["effect_allele"] = snp['ALT']
            variant["effect_size"] = snp['BETA']
            variants[snp['rsid']] = variant

        model = {"score_model": {"categories": categories, "variants": variants}}
        model_file.write(yaml.dump(model, indent=2))

        description = {"description": description}
        model_file.write(yaml.dump(description, indent=2, default_flow_style=False))

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
    if not is_valid_path(parsed_args.input): return
    if not is_valid_path(parsed_args.origin_reference_vcf): return
    if not is_valid_path(parsed_args.model_reference_vcf): return
    if not is_valid_path(parsed_args.af): return
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
    parser.add_argument('--output-directory', type=str, default='',
                        help='output directory')
    parsed_args = parser.parse_args(args)
    output_path = os.path.abspath(os.path.expanduser(
        parsed_args.output_directory)) + "/gbe_phenotype_manifest.tsv"
    download(parsed_args.url, output_path)
    return

###############
### gbe-get ###
###############

def gbe_get(parsed_args):
    url = "https://biobankengine.stanford.edu/static/PRS_map/" + parsed_args.code + ".tsv"
    output_directory = os.path.abspath(os.path.expanduser(parsed_args.output_directory))
    output_file_name = os.path.splitext(os.path.basename(url))[0]
    output_path = output_directory + "/" + output_file_name
    download(url=url, output_path=output_path,
             force=parsed_args.force, progress=True)
    return output_path

#######################
### gbe-model #########
#######################

def gbe_model(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker biobankuk-build-model constructs polygenic score model based on p value data')  # todo dodać opis
    parser.add_argument('-c','--code', required=True, type=str, help='path to PRS file from gbe. It can be downloaded using gbe-get')
    parser.add_argument('-o', '--output-directory', type=str, default='output', help='output directory')
    parser.add_argument('--gbe-index', type=str, default='gbe-index.1.3.1.tsv', help='gbe-index')
    parser.add_argument('--source-ref-vcf', type=str, default='dbsnp155.grch37.norm.vcf.gz', help='source reference vcf (hg19)')
    parser.add_argument('--target-ref-vcf', type=str, default='dbsnp155.grch38.norm.vcf.gz', help='source reference vcf (hg38)')
    parser.add_argument('--af-vcf', type=str, default='gnomad.3.1.vcf.gz', help='path to allele frequency vcf.')
    parser.add_argument('--af-field', type=str, default='AF_nfe', help='name of the INFO field with ALT allele frequency')
    parser.add_argument('--gene-positions', type=str, default='ensembl-genes.104.tsv', help='table with ensembl genes')
    parser.add_argument('-i', '--iterations', type=float, default=1000, help='simulation iterations for mean and sd')
    parser.add_argument('-f', '--force', action='store_true', help='overwrite downloaded file')
    parsed_args = parser.parse_args(args)

    if not is_valid_path(parsed_args.output_directory, is_directory=True): return
    path = gbe_get(parsed_args)

    if not is_valid_path(path): return
    data = read_table(path)

    if not is_valid_path(parsed_args.gbe_index): return
    gbe_index = read_table(parsed_args.gbe_index)
    info = [line for line in gbe_index if parsed_args.code in line['Trait name']][0]

    if not is_valid_path(parsed_args.gene_positions): return
    gene_positions = read_table(parsed_args.gene_positions)

    if not is_valid_path(parsed_args.af_vcf, possible_url = True): return
    af_vcf = VcfAccessor(parsed_args.af_vcf)

    if not is_valid_path(parsed_args.source_ref_vcf, possible_url = True): return
    source_vcf = VcfAccessor(parsed_args.source_ref_vcf)

    if not is_valid_path(parsed_args.target_ref_vcf, possible_url = True): return
    target_vcf = VcfAccessor(parsed_args.target_ref_vcf)

    data = [line for line in data if "rs" in line['ID']]
    for line in data: line.update({"rsid": line['ID']})
    data = [validate(line, validation_source = target_vcf) for line in data]

    for line in data:
        if not "rsid" in line:
            print(line)

    data = [line for line in data if "rs" in line['ID']]

    data = [add_annotation(
        line, 
        annotation_name = "af", 
        annotation_source = af_vcf, 
        annotation_source_field = parsed_args.af_field,
        default_value = 0.001) for line in data]

    symbols = [add_gene_symbol(
        line, 
        target_vcf, 
        gene_positions) for line in data]
    symbols = [symbol for symbol in symbols if symbol]

    description = dict()
    parameters = simulate_parameters(data)
    description["pmids"] = [33095761]
    description["info"] = info
    description["parameters"] = parameters
    description["genes"] = symbols

    model_path = parsed_args.output_directory + "/" + parsed_args.code + ".yml"
    write_model(data, description, model_path)
    return


#######################
### biobankuk-index ###
#######################

def biobankuk_index(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker biobankuk-index downloads index of gwas results from pan.ukbb study')  # todo dodać opis
    parser.add_argument('--url', type=str, default='https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz',
                        help='alternative url location for index')
    parser.add_argument('--variant-metrics', type=str, default='https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz',
                        help='alternative url location for variant metrics')
    parser.add_argument('--output-directory', type=str, default='',
                        help='output directory')
    parsed_args = parser.parse_args(args)
    output_path = os.path.abspath(os.path.expanduser(
        parsed_args.output_directory)) + "/panukbb_phenotype_manifest.tsv"
    download(parsed_args.url, output_path)
    output_path = os.path.abspath(os.path.expanduser(
        parsed_args.output_directory)) + "/full_variant_qc_metrics.txt"
    download(parsed_args.variant_metrics, output_path)
    return

#######################
### biobankuk-get #####
#######################

def biobankuk_get(parsed_args):
    # checking index file for download url
    with open(parsed_args.index, 'r') as indexfile:
        firstline = indexfile.readline()
        phenocode_colnumber = firstline.split('\t').index("phenocode")
        pheno_sex_colnumber = firstline.split('\t').index("pheno_sex")
        coding_colnumber = firstline.split('\t').index("coding")
        aws_link_colnumber = firstline.split('\t').index("aws_link")
        while True:
            line = indexfile.readline()
            if not line:
                break
            if line.split('\t')[phenocode_colnumber] != parsed_args.code:
                continue
            if line.split('\t')[pheno_sex_colnumber] != parsed_args.sex:
                continue
            if line.split('\t')[coding_colnumber] != parsed_args.coding:
                continue
            url = line.split('\t')[aws_link_colnumber]
            break
    # downloading
    if not url is None:
        output_directory = os.path.abspath(os.path.expanduser(parsed_args.output_directory))
        output_file_name = os.path.splitext(os.path.basename(url))[0]
        output_path = output_directory + "/" + output_file_name
        logger.info("Downloading from " + url + " to " + output_path)
        download(url=url, output_path=output_path,
             force=False, progress=True)
        return output_path
    return None

#############################
### biobankuk-model #########
#############################

def biobankuk_model(args):
    parser = argparse.ArgumentParser(
        description='polygenicmaker biobankuk-model prepares polygenic score model based on p value data')
    parser.add_argument('--code', '--phenocode', type=str, required=True,
                        help='phenocode of phenotype form Uk Biobank')
    parser.add_argument('--sex', '--pheno_sex', type=str, default="both_sexes",
                        help='pheno_sex of phenotype form Uk Biobank')
    parser.add_argument('--coding', type=str, default="",
                        help='additional coding of phenotype form Uk Biobank')
    parser.add_argument('--index', type=str, required=True,
                        help='path to Index file from PAN UKBiobank. It can be downloaded using gbe-get')
    parser.add_argument('--output-directory', type=str, default='',
                        help='output directory')
    parser.add_argument('--variant-metrics', type=str, required=True,
                        help='path to annotation file. It can be downloaded from https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz')
    parser.add_argument('--threshold', type=float, default=1e-08,
                        help='significance cut-off threshold')
    parser.add_argument('--population', type=str, default='EUR',
                        help='population: meta, AFR, AMR, CSA, EAS, EUR, MID')
    parser.add_argument('--clumping-vcf', type=str, default='/tmp/polygenic/results/eur.phase3.biobank.set.vcf.gz',
                        help='')
    parser.add_argument('--source-ref-vcf', type=str, default='ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz',
                        help='')
    parser.add_argument('--target-ref-vcf', type=str, default='/tmp/marpiech/kenobi/resources/GRCh38.dbSNP155.chr.norm.rsidonly.vcf.gz',
                        help='')

    parsed_args = parser.parse_args(args)
    if not is_valid_path(parsed_args.output_directory, is_directory=True): return
    path = biobankuk_get(parsed_args)
    if not is_valid_path(path): return
    if not is_valid_path(parsed_args.variant_metrics): return
    if not is_valid_path(parsed_args.target_ref_vcf): return
    if not is_valid_path(parsed_args.target_ref_vcf, possible_url = True): return
    source_vcf = VcfAccessor(parsed_args.source_ref_vcf)
    target_vcf = VcfAccessor(parsed_args.target_ref_vcf)

    filter_pval(path, parsed_args)
    clump(path, parsed_args)

    data = read_table(path + ".clumped")
    for line in data: line.update({"rsid": line['chr'] + ":" + line['pos'] + "_" + line['ref'] + "_" + line['alt']})
    for line in data: line.update({"REF": line['ref'], "ALT": line['alt']})
    for line in data: line.update({"BETA": line["beta_" + parsed_args.population]})
    for line in data: line.update({"af": line["af_" + parsed_args.population]})

    data = [add_annotation(
        line, 
        annotation_name = "rsid", 
        annotation_source = source_vcf, 
        annotation_source_field = "rsid",
        default_value = "rs0") for line in data]

    data = [validate(line, validation_source = target_vcf) for line in data]

    description = simulate_parameters(data)

    model_path = path + ".yml"
    write_model(data, description, model_path)

    return


def filter_pval(path, parsed_args):
    output_path = path + ".filtered"
    with open(path, 'r') as data, open(parsed_args.variant_metrics, 'r') as anno, open(output_path, 'w') as output:
        data_header = data.readline().rstrip().split('\t')
        anno_header = anno.readline().rstrip().split('\t')
        output.write('\t'.join(data_header + anno_header) + "\n")
        while True:
            try:
                data_line = data.readline().rstrip().split('\t')
                anno_line = anno.readline().rstrip().split('\t')
                if float(data_line[data_header.index('pval_' + parsed_args.population)].replace('NA', '1', 1)) <= parsed_args.threshold:
                    output.write('\t'.join(data_line + anno_line) + "\n")
            except:
                break
    return


def clump(path, parsed_args):

    filtered_path = path + ".filtered"
    clumped_path = path + ".clumped"

    subprocess.call("plink" +
                    " --clump " + filtered_path +
                    " --clump-p1 " + str(parsed_args.threshold) +
                    " --clump-r2 0.25 " +
                    " --clump-kb 1000 " +
                    " --clump-snp-field rsid " +
                    " --clump-field pval_" + parsed_args.population +
                    " --vcf " + parsed_args.clumping_vcf + " " +
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
        elif args[0] == 'biobankuk-model':
            biobankuk_model(args[1:])
        elif args[0] == 'gbe-index':
            gbe_index(args[1:])
        elif args[0] == 'gbe-model':
            gbe_model(args[1:])
        elif args[0] == 'pgs-index':
            pgs_index(args[1:])
        elif args[0] == 'pgs-model':
            pgs_model(args[1:])
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
            biobankuk-model         prepare polygenic score model based on gwas results from biobankuk
            gbe-index               downloads Global Biobank Engine index
            gbe-model               prepare polygenic score model from GBE data
            pgs-index               downloads Polygenic Risk Score index
            pgs-model               prepare polygenic score model from PGS data
            vcf-index               prepare rsidx for vcf

            """)
    except polygenic.core.trial.PolygenicException as e:
        print("Analysis failed")
        print("ERROR: " + str(e))
    except RuntimeError as e:
        print("ERROR: " + str(e))

if __name__ == '__main__':
    main(sys.argv[1:])
