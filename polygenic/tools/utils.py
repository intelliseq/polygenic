import os
import logging
import sys
import subprocess
import urllib
import re
import random
import statistics
import yaml
import gzip
import io
import os.path
import progressbar
import urllib.request
import numpy
from datetime import datetime
from polygenic.version import __version__ as version

from polygenic.data.vcf_accessor import VcfAccessor

def error_print(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def error_exit(e):
    time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    error_print(" -> polygenic " + version + " failed at " + time)
    error_print(" -> with command: pgstk " + (" ").join(sys.argv))
    error_print(" -> with message: ")
    error_print(str(e))
    exit(1)

def expand_path(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path)) if path else ''

def setup_logger(path):
    logger = logging.getLogger('pgstk')

    log_directory = os.path.dirname(os.path.abspath(os.path.expanduser(path)))
    if log_directory:
        try:
            os.makedirs(log_directory)
        except OSError:
            pass
    logger.setLevel(logging.DEBUG)
    logging_file_handler = logging.FileHandler(path)
    logging_file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logging_file_handler.setFormatter(formatter)
    logger.addHandler(logging_file_handler)

    return logger

### model tools
def download(url: str, output_path: str, force: bool=False, progress: bool=False):
    """Downloads file from url

    Keyword arguments:
    url -- url to file
    output_path -- path to output file
    force -- flag whether to overwrite downloaded file
    progress -- flag whether to present progress
    """
    logger = logging.getLogger('utils')

    if os.path.isfile(output_path) and not force:
        logger.warning("File already exists: " + output_path)
        return output_path
    logger.info("Downloading from " + url)
    response = urllib.request.urlopen(url)
    file_size = int(response.getheader('Content-length'))
    if file_size is None:
        progress = False
    if ".gz" in url or ".bgz" in url:
        subprocess.call("wget " + url + " -O " + output_path + ".gz",
                    shell=True)
        subprocess.call("gzip -d " + output_path + ".gz",
                    shell=True)
        return output_path
    else:
        response_data = response
    if progress: bar = progressbar.ProgressBar(max_value = file_size).start()
    downloaded = 0
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    bytebuffer = b''
    while (bytes := response_data.read(1024)):
            bytebuffer = bytebuffer + bytes
            downloaded = downloaded + 1024
            if not file_size == progressbar.UnknownLength: downloaded = min(downloaded, file_size)
            progress and bar.update(downloaded)
    with open(output_path, 'w') as outfile:
        outfile.write(str(bytebuffer, 'utf-8'))
    progress and bar.finish()
    return output_path

def is_valid_path(path: str, is_directory: bool = False, create_directory: bool = True, possible_url: bool = False):
    """Checks whether path is valid.

    Keyword arguments:
    path -- the path to file or directory
    is_directory -- flag if the targe is directory
    """
    if possible_url and "://" in path:
        return True
    if is_directory:
        if create_directory:
            try:
                os.makedirs(path, exist_ok=True)
            except:
                print("ERROR: Could not create " + path)
                return False
        if not os.path.isdir(path):
            print("ERROR: " + path + " does not exists or is not directory")
            return False
    else:
        if not os.path.isfile(path):
            print("ERROR: " + path + " does not exists or is not a file")
            return False
    return True

def clump(
    gwas_file, 
    reference,  
    clump_p1 = "1e-08", 
    clump_r2 = "0.25", 
    clump_kb = "1000",
    clump_snp_field = "rsid",
    clump_field = "pval_EUR"):

    clumped_path = gwas_file + ".clumped"

    subprocess.call("plink" +
                    " --clump " + gwas_file +
                    " --clump-p1 " + str(clump_p1) +
                    " --clump-r2 " + str(clump_r2) +
                    " --clump-kb " + str(clump_kb) +
                    " --clump-snp-field " + str(clump_snp_field) +
                    " --clump-field " + str(clump_field) +
                    " --vcf " + str(reference) + " " +
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
    
    with open(gwas_file, 'r') as filtered_file, open(clumped_path, 'w') as clumped_file:
        filtered_header = filtered_file.readline().rstrip().split('\t')
        clumped_file.write('\t'.join(filtered_header) + "\n")
        while True:
            try:
                filtered_line = filtered_file.readline().rstrip().split('\t')
                if filtered_line[filtered_header.index('rsid')] in clumped_rsids:
                    clumped_file.write('\t'.join(filtered_line) + "\n")
            except:
                break
    return clumped_path

def read_header(file_path: str):
    """Reads header into dictionary. First row is treated as keys for dictionary.

    Keyword arguments:
    path -- the path to .tsv file
    """
    header = {}
    with open(file_path, 'r') as file:
        while True:
            line = file.readline().rstrip()
            if line[0] == '#':
                if line[1] == ' ':
                    key,value = line[2:].split(' = ')
                    header[key] = value
            else:
                break
    return header


def read_table(file_path: str, delimiter: str = '\t'):
    """Reads table into dictionary. First row is treated as keys for dictionary.

    Keyword arguments:
    path -- the path to .tsv file
    """
    logger = logging.getLogger('utils')


    table = []
    with open(file_path, 'r') as file:
        line = file.readline()
        while line[0] == '#':
            line = file.readline()
        header = line.rstrip('\r\n').split(delimiter)
        while True:
            line = file.readline().rstrip('\r\n').split(delimiter)
            if len(line) < 2:
                break
            if not len(header) == len(line):
                logger.error("Line and header have different leangths")
                raise RuntimeError("Line and header have different leangths. LineL {line}".format(line = str(line)))
            line_dict = {}
            for header_element, line_element in zip(header, line):
                line_dict[header_element] = line_element
            table.append(line_dict)
    return table

def write_data(data: list, file_path: str, delimiter: str = '\t'):
    """Reads table into dictionary. First row is treated as keys for dictionary.

    Keyword arguments:
    path -- the path to .tsv file
    """

    keys = set()
    for line in data: keys.update(set(line.keys()))
    with open(file_path, 'w') as file:
        file.write(delimiter.join(keys) + os.linesep)
        for line in data:
            values = [str(line[key]) if key in line else "" for key in keys]
            file.write(delimiter.join(values) + os.linesep)
    return file_path

def validate(
    validated_line: dict,
    validation_source: VcfAccessor,
    invert_field: str = None,
    ignore_warnings: bool = False,
    strict: bool = True):
    record = None
    snpid = None
    if "id" in validated_line:
        record = validation_source.get_record_by_rsid(validated_line['id'])
        snpid = validated_line['id']
    if "rsid" in validated_line and record is None:
        record = validation_source.get_record_by_rsid(validated_line['rsid'])
        snpid = validated_line['rsid']
    if "gnomadid" in validated_line and record is None:
        record = validation_source.get_record_by_rsid(validated_line['gnomadid'])
        snpid = validated_line['gnomadid']
    if record is None:
        error_print("ERROR: Failed validation for " + snpid + ". SNP not present in validation vcf.")
        validated_line["status"] = "ERROR"
        return None if strict else validated_line
    if not (validated_line['ref'] == record.get_ref()): 
        if (validated_line['ref'] == record.get_alt()[0] and validated_line['alt'] == record.get_ref()):
            ref = validated_line['ref']
            alt = validated_line['alt']
            validated_line['ref'] = alt
            validated_line['alt'] = ref
            if invert_field is not None:
                validated_line[invert_field] = - float(validated_line[invert_field])
            error_print("WARNING: " + "Failed validation for " + validated_line['rsid'] + ". REF and ALT do not match. " + record.get_ref() + "/" + str(record.get_alt()) + " succesful invert!")
            validated_line["status"] = "WARNING"
            return validated_line if ignore_warnings else None
        else:
            error_print("ERROR: " + "Failed validation for " + validated_line['rsid'] + ". REF and ALT do not match. " + record.get_ref() + "/" + str(record.get_alt()))
            validated_line["status"] = "ERROR"
            return None if strict else validated_line
    validated_line["status"] = "SUCCESS"
    return validated_line


def validate_with_source(data, source_file, ignore_warnings = False):
    source_accessor = VcfAccessor(source_file)
    data = [validate(
        validated_line = line,
        validation_source = source_accessor,
        ignore_warnings = ignore_warnings) for line in data
    ]
    data = [line for line in data if line]
    return data

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
        'max': maxsum,
        'percentiles': get_percentiles(randomized_beta_list)
    }

def randomize_beta(beta: float, af: float):
    first_allele_beta = beta if random.uniform(0, 1) < af else 0
    second_allele_beta = beta if random.uniform(0, 1) < af else 0
    return first_allele_beta + second_allele_beta

def get_percentiles(value_list: list):
    value_array = numpy.array(value_list)
    percentiles = {}
    for i in range(101):
        percentiles[str(i)] = numpy.percentile(value_array, i)
    return percentiles

def write_model(
    data, 
    description, 
    destination, 
    id_field = 'rsid',
    effect_allele_field = 'alt',
    effect_size_field = 'beta',
    included_fields_list = []):

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
            variant["effect_allele"] = snp[effect_allele_field]
            variant["effect_size"] = snp[effect_size_field]
            for field in included_fields_list:
                variant[field] = snp[field]
            variants[snp[id_field]] = variant

        model = {"score_model": {"categories": categories, "variants": variants}}
        model_file.write(yaml.dump(model, indent=2))

        description = {"description": description}
        model_file.write(yaml.dump(description, indent=2, default_flow_style=False))

    return
