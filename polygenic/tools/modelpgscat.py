import argparse
import sys
import os
import re
import csv
from tqdm import tqdm

import polygenic.tools.utils as utils
from polygenic.data.vcf_accessor import VcfAccessor
from polygenic.error.polygenic_exception import PolygenicException

def parse_args(args):
    parser = argparse.ArgumentParser(description='pgstk model-biobankuk prepares polygenic score model based on p value data')
    parser.add_argument('--code', '--phenocode', type=str, required=True, help='phenocode of phenotype form Uk Biobank')
    parser.add_argument('--output-directory', type=str, default='.', help='output directory')
    parser.add_argument('--index-file', type=str, help='path to index file')
    parser.add_argument('--index-url', type=str, default='http://ftp.ebi.ac.uk/pub/databases/spot/pgs/metadata/pgs_all_metadata_scores.csv', help='url of index file for PAN UKBiobank.')
    parser.add_argument('--af', type=str, help='vcf file containing allele freq data', default='gnomad.3.1.vcf.gz')
    parser.add_argument('--af-field', type=str, default='AF',help='name of the INFO field to be used as allele frequency')
    parser.add_argument('--origin-genome-build', type=str)
    parser.add_argument('--source-ref-vcf', type=str, default='dbsnp155.grch37.norm.vcf.gz', help='')
    parser.add_argument('--target-ref-vcf', type=str, default='dbsnp155.grch38.norm.vcf.gz', help='')
    parser.add_argument('--gene-positions', type=str, default='ensembl-genes.104.tsv', help='table with ensembl genes')
    parser.add_argument('--ignore-warnings', type=bool, default='False', help='')
    parser.add_argument('-l', '--log-file', type=str, help='path to log file')
    parsed_args = parser.parse_args(args)
    parsed_args.index_file = parsed_args.index_file if parsed_args.index_file else parsed_args.output_directory + "/pgs_all_metdata_scores.csv"
    return parsed_args

def get_index(args):
    utils.download(args.index_url, os.path.abspath(args.index_file))
    return

def get_data(args):
    with open(args.index_file, 'r') as indexfile:
        firstline = next(csv.reader([indexfile.readline()]))
        phenocode_colnumber = firstline.index("Polygenic Score (PGS) ID")
        aws_link_colnumber = firstline.index("FTP link")
        while True:
            line = indexfile.readline()
            if not line:
                break
            line = next(csv.reader([line]))
            if line[phenocode_colnumber] != args.code:
                continue
            url = line[aws_link_colnumber]
            break
    if not url is None:
        output_directory = os.path.abspath(os.path.expanduser(args.output_directory))
        output_file_name = os.path.splitext(os.path.basename(url))[0]
        output_path = output_directory + "/" + output_file_name
        utils.download(url=url, output_path=output_path, force=False, progress=True)
        args.gwas_file = output_path
        return output_path
    return None

def get_info(args):
    with open(args.index_file, 'r') as indexfile:
        firstline = indexfile.readline()
        colnames = next(csv.reader([firstline]))
        phenocode_colnumber = colnames.index("Polygenic Score (PGS) ID")
        while True:
            line = indexfile.readline()
            if not line:
                break
            splitted_line = next(csv.reader([line]))
            if splitted_line[phenocode_colnumber] != args.code:
                continue
            output = dict()
            for i in range(len(colnames) - 1):
                output[colnames[i]] = splitted_line[i]
            header = utils.read_header(args.gwas_file)
            for key in header:
                output[key] = header[key]
            return output
        
    return dict()

def define_origin_genome_build(build: str):
    if build == "NR":
        raise PolygenicException("Origin genome build is not defined")
    if "19" in build:
        return "GRCh37"
    if "37" in build:
        return "GRCh37"
    if "36" in build:
        raise PolygenicException("Origin genome build GRCh36 is not supported")
    if "38" in build:
        return "GRCh38"

def validate_paths(args):
    if not utils.is_valid_path(args.output_directory, is_directory=True): raise PolygenicException("Output directory does not exists")
    if not utils.is_valid_path(args.gwas_file): raise PolygenicException("PRS file does not exists")
    if args.source_ref_vcf and not utils.is_valid_path(args.source_ref_vcf): raise PolygenicException("GRCh37 reference vcf does not exists")
    if not utils.is_valid_path(args.target_ref_vcf): PolygenicException("GRCh38 reference vcf does not exists")

def read_variants(args):
    data = utils.read_table(args.gwas_file)
    if "rsID" in data[0]:
        for line in data: line.update({"rsid": line["rsID"]}) 
    if "chr_name" in data[0]:
        for line in data: line.update({"chr": line["chr_name"]})
    if "chr_position" in data[0] and data[0]["chr_position"] != "":
        for line in data: line.update({"pos": line["chr_position"]})
    if "reference_allele" in data[0]:
        for line in data: line.update({"ref": line["reference_allele"]})
    if "other_allele" in data[0]:
        for line in data: line.update({"ref": line["other_allele"]})
    if "effect_allele" in data[0]:
        for line in data: line.update({"alt": line["effect_allele"]})
    if "effect_weight" in data[0]:
        for line in data: line.update({"beta": float(line["effect_weight"])})
    if "weight_type" in data[0]:
        for line in data: line.update({"weight_type": line["weight_type"]})
    if "allelefrequency_effect" in data[0]:
        for line in data: line.update({"af": float(line["allelefrequency_effect"])})
    if "chr" in data[0] and "pos" in data[0] and "ref" in data[0] and "alt" in data[0]:
        for line in data: line.update({"gnomadid": line['chr'] + ":" + line['pos'] + "_" + line['ref'] + "_" + line['alt']})
    return data

def run(args):

    args.logger.info("Checking for index")
    get_index(args) # download index
    gwas_file = get_data(args) # download gwas results
    validate_paths(args) # check if vcf files are correct

    ### check genome build
    description = dict()
    description["info"] = get_info(args)
    if args.origin_genome_build:
        origin_genome_build = define_origin_genome_build(args.origin_genome_build)
    else:
        origin_genome_build = define_origin_genome_build(description["info"]["Original Genome Build"])
    
    args.logger.info("Reading variants")
    data = read_variants(args) # read filtered variants
    origin_beta = utils.sum_beta(data)
    if origin_genome_build == "GRCh37":
        args.logger.info("Validating variants with GRCh37")
        data = utils.validate_with_source(data, args.source_ref_vcf, ignore_warnings = True) # validate if variants are present in hg19
        args.logger.info("Succesfully validated " + str(len(data)) + " variants")
        args.logger.info("Validating variants with GRCh38")
        data = utils.validate_with_source(data, args.target_ref_vcf, ignore_warnings = args.ignore_warnings, use_gnomadid = False) # validate if variants are present in hg38
    else:
        args.logger.info("Validating variants with GRCh38")
        data = utils.validate_with_source(data, args.target_ref_vcf, ignore_warnings = args.ignore_warnings) # validate if variants are present in hg38
    final_beta = utils.sum_beta(data)
    args.logger.info("Succesfully validated " + str(len(data)) + " variants")
    args.logger.info("Annotating with symbols")
    data = utils.annotate_with_symbols(data, args.gene_positions)
    genes = utils.get_gene_symbols(data)
    if "af" not in data[0]:
        args.logger.info("Annotating with af")
        data = utils.annotate_with_af(data, args.af, af_field = args.af_field, default_af = 0)
    
    description["arguments"] = utils.args_to_dict(args)
    description["parameters"] = utils.simulate_parameters(data)
    description["parameters"]["number of variants"] = len(data)
    description["parameters"]["origin_beta"] = origin_beta
    description["parameters"]["final_beta"] = final_beta
    description["pmid"] = description["info"]['Publication (PMID)']
    description["genes"] = {"symbols": genes}
    name = re.sub("[^0-9a-zA-Z]+", "_", description["info"]["Reported Trait"].lower().replace(" ", "_")) # trait name
    filename = "-".join(["pgscat", args.code, name, args.af_field]) + ".yml"
    model_path = "/".join([args.output_directory, filename]) # output path
    utils.write_model(data, description, model_path, included_fields_list = ['ref', 'gnomadid', 'af']) # writing model
    return

def main(args = sys.argv[1:]):

    args = parse_args(args) 
    args.logger = utils.setup_logger(args.log_file) if args.log_file else utils.setup_logger(args.output_directory + "/pgstk.log")

    try:
        run(args)
    except PolygenicException as e:
        utils.error_exit(e)

if __name__ == '__main__':
    main(sys.argv[1:])