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
import csv
from datetime import datetime
from polygenic.version import __version__ as version

from polygenic.data.vcf_accessor import VcfAccessor
from polygenic.data.csv_accessor import CsvAccessor
from tqdm import tqdm

def expand_path(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path)) if path else ''

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
    
    # replace apostrophe in url
    url = url.replace("'", "%27")

    # convert all nonalphanumeric characters to underscores except slashes, dots and dashes in output path
    output_path = re.sub(r'[^a-zA-Z0-9/.-]', '_', output_path).lower()

    if os.path.isfile(output_path) and not force:
        logger.warning("File already exists: " + output_path)
        return output_path
    logger.info("Downloading from " + url)
    response = urllib.request.urlopen(url)
    file_size = int(response.getheader('Content-length'))
    if file_size is None:
        progress = False
    if ".gz" in url or ".bgz" in url:
        subprocess.call("wget '" + url + "' -O '" + output_path + ".gz'",
                    shell=True)
        subprocess.call("gzip -d '" + output_path + ".gz'",
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

### lasso helper function
def is_in_the_same_window(first_snp, second_snp, window_size = 1000):
    contig_field_name = 'chr'
    if 'chrom' in first_snp:
        contig_field_name = 'chrom'
    return first_snp[contig_field_name] == second_snp[contig_field_name] and \
           abs(int(first_snp['pos']) - int(second_snp['pos'])) <= window_size

### lasso helper function
def flatten(lst):
    return sum(lst, [])

### lasso helper function
def local_minima(lst):
    if len(lst) < 2:
        return list(range(len(lst)))
    
    lmin = True
    curr = list()
    local_minima = list()
    
    for i, x in enumerate(lst[:-1]):
        if lst[i + 1] == x and lmin:
            curr.append(i)
        if lst[i + 1] < x:
            lmin = True
            curr = list()
        if lst[i + 1] > x:
            if lmin:
                local_minima.append(curr + [i])
            lmin = False

    if lst[-1] < lst[-2] or (lmin and lst[-1] == lst[-2]):
        local_minima.append(curr + [len(lst) - 1])
    
    return flatten(local_minima)

### lasso helper function
def invert(a):
    return [-x for x in a]

### lasso helper function
def local_maxima(a):
    return local_minima(invert(a))

def lasso_clump(
    gwas_file, 
    clump_p = "1e-08", 
    clump_kb = "200",
    clump_snp_field = "rsid",
    clump_field = "pval_EUR",
    return_file = True):

    # define output path
    clumped_path = gwas_file + ".clumped2"

    # read gwas file into dictionary
    gwas_dict = read_table(gwas_file)

    # cluster snps into regions
    clumped_rsids = []
    region = []
    for snp in gwas_dict:
        if float(snp[clump_field]) > float(clump_p):
            continue
        if float(snp[clump_field]) == 0:
            continue
        if not region:
            region.append(snp)
            continue
        if is_in_the_same_window(region[-1], snp, window_size = int(clump_kb) * 1000):
            region.append(snp)
            continue
        pvalues = [float(x[clump_field]) for x in region]
        min_index = pvalues.index(min(pvalues))
        clumped_rsids.append(region[min_index][clump_snp_field])
        region = []

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

    # stdout, stderr = process.communicate()
    clumped_rsids = []

    # if plink.clumped does not exist return None
    if not os.path.isfile("plink.clumped"):
        return None

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

def sum_beta(data):
    sum = 0
    for line in data:
        if 'beta' in line:
            sum += abs(float(line['beta']))
    return sum

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
            line = next(csv.reader([file.readline().rstrip('\r\n')], delimiter = delimiter))
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
def invert_nucleotides(obj):
    if type(obj).__name__ == "str": return _invert_nucleotides_str(obj)
    if type(obj).__name__ == "list": return _invert_nucleotides_list(obj)
    return None

def _invert_nucleotides_list(nucleotides: list):
    return [_invert_nucleotides_str(nucleotide) if type(nucleotide).__name__ == "str" else None for nucleotide in nucleotides ]

def _invert_nucleotides_str(nucleotide: str):
    return {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }.get(nucleotide, None)

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

def get_record(
    line: dict,
    vcf_accessor: VcfAccessor,
    use_gnomadid: bool = True):
    record = None
    if "id" in line:
        record = vcf_accessor.get_record_by_rsid(line['id'])
        snpid = ine['id']
        return record, snpid
    if "rsid" in line and record is None:
        record = vcf_accessor.get_record_by_rsid(line['rsid'])
        snpid = line['rsid']
        return record, snpid
    if "gnomadid" in line and use_gnomadid and record is None:
        record = vcf_accessor.get_record_by_rsid(line['gnomadid'])
        snpid = line['gnomadid']    
        return record, snpid
    return None, None

def validate(
    validated_line: dict,
    validation_source: VcfAccessor,
    invert_field: str = "beta",
    ignore_warnings: bool = False,
    strict: bool = True,
    use_gnomadid: bool = True,
    verbose = False):
    record, snpid = get_record(validated_line, validation_source, use_gnomadid)
    if record is None:
        if verbose: error_print("ERROR: " + "Failed validation for " + validated_line['rsid'] + ". No id in source")
        validated_line["status"] = "ERROR"
        return None if strict else validated_line
    if not "ref" in validated_line:
        validated_line['ref'] = record.get_ref()
    if not (validated_line['ref'] == record.get_ref() and validated_line['alt'] in record.get_alt()): 
        message = "REF mismatch. Line: " + validated_line['ref'] + "/" + validated_line['alt'] + " . Reference: " + record.get_ref() + "/" + str(record.get_alt()) + "."
        if validated_line['ref'] == invert_nucleotides(validated_line['alt']):
            if verbose: error_print("ERROR: " + message + " REF and ALT are complementary. Not sure whether to swap (beta will be negative) or invert (beta will be unchanged) nucleotides.")
            validated_line["status"] = "ERROR"
            return None if strict else validated_line
        if (invert_nucleotides(validated_line['ref']) == record.get_ref() and invert_nucleotides(validated_line['alt']) in record.get_alt()):
            validated_line['ref'] = invert_nucleotides(validated_line['ref'])
            validated_line['alt'] = invert_nucleotides(validated_line['alt'])
            if verbose: error_print("WARNING: " + message + " REF and ALT are on '-' strand notation. Succesfull inversion.")
            validated_line["status"] = "WARNING. Minus strand notation"
            if "pos" not in validated_line: validated_line["pos"] = record.get_pos()
            if validated_line["pos"] == "": validated_line["pos"] = record.get_pos()
            validated_line["gnomadid"] = record.get_chrom().replace("chr", "") + ":" + record.get_pos() + "_" + record.get_ref() + "_" + validated_line['alt']
            return validated_line if ignore_warnings else None
        if (validated_line['ref'] in record.get_alt() and validated_line['alt'] == record.get_ref()):
            ref = validated_line['ref']
            alt = validated_line['alt']
            validated_line['ref'] = alt
            validated_line['alt'] = ref
            if invert_field is not None:
                validated_line[invert_field] = - float(validated_line[invert_field])
            if verbose: error_print("WARNING: " + message + " REF and ALT are swapped. Succesful swap.")
            validated_line["status"] = "WARNING. Inverted ref with alt"
            if "pos" not in validated_line: validated_line["pos"] = record.get_pos()
            if validated_line["pos"] == "": validated_line["pos"] = record.get_pos()
            validated_line["gnomadid"] = record.get_chrom().replace("chr", "") + ":" + record.get_pos() + "_" + record.get_ref() + "_" + validated_line['alt']
            return validated_line if ignore_warnings else None
        if (invert_nucleotides(validated_line['ref']) in record.get_alt() and invert_nucleotides(validated_line['alt']) == record.get_ref()):
            ref = invert_nucleotides(validated_line['ref'])
            alt = invert_nucleotides(validated_line['alt'])
            validated_line['ref'] = alt
            validated_line['alt'] = ref
            if invert_field is not None:
                validated_line[invert_field] = - float(validated_line[invert_field])
            if verbose: error_print("WARNING: " + message + " REF and ALT are on '-' strand notation and swapped. Succesfull swap and inversion.")
            validated_line["status"] = "WARNING. Inverted ref with alt"
            if "pos" not in validated_line: validated_line["pos"] = record.get_pos()
            if validated_line["pos"] == "": validated_line["pos"] = record.get_pos()
            validated_line["gnomadid"] = record.get_chrom().replace("chr", "") + ":" + record.get_pos() + "_" + record.get_ref() + "_" + validated_line['alt']
            return validated_line if ignore_warnings else None
        if verbose: error_print("ERROR: " + message + ". REF and ALT do not match.")
        validated_line["status"] = "ERROR"
        return None if strict else validated_line
    if record.get_id():
        validated_line["rsid"] = record.get_id()
    if "pos" not in validated_line: validated_line["pos"] = record.get_pos()
    if validated_line["pos"] == "": validated_line["pos"] = record.get_pos()
    validated_line["gnomadid"] = record.get_chrom().replace("chr", "") + ":" + record.get_pos() + "_" + record.get_ref() + "_" + validated_line['alt']
    validated_line["status"] = "SUCCESS"
    return validated_line

def args_to_dict(args):
    args_dict = vars(args)
    out_dict = dict()
    for arg in args_dict:
        if type(args_dict[arg]).__name__ in ["str", "bool", "int", "float"]:
            out_dict[arg] = args_dict[arg]
    return out_dict

def validate_with_source(data, source_path, ignore_warnings = False, use_gnomadid = True):
    source_accessor = VcfAccessor(source_path)
    validated_data = list()
    for line in tqdm(data):
        validated_data.append(validate(
            validated_line = line,
            validation_source = source_accessor,
            ignore_warnings = ignore_warnings,
            use_gnomadid = use_gnomadid))
    validated_data = [line for line in validated_data if line]
    return validated_data

def add_af(line, source_accessor, af_field = "af", default_af = 0):
    record, snpid = get_record(line, source_accessor)
    if not record is None:
        line["af"] = float(record.get_info_field(af_field))
    else:
        line["af"] = float(default_af)
    return line

def annotate_with_af(data, source_path, af_field = "af", default_af = 0):
    source_accessor = VcfAccessor(source_path)
    annotated_data = list()
    for line in tqdm(data):
        annotated_data.append(add_af(line, source_accessor, af_field, default_af))
    return annotated_data

def annotate_with_symbols(data, source_path):
    source_accessor = CsvAccessor(source_path)
    annotated_data = list()
    for line in tqdm(data):
        if "symbol" in line:
            annotated_data.append(line)
            continue
        if "gnomadid" in line:
            chrom = line["gnomadid"].split(":")[0]
            pos = line["gnomadid"].split(":")[1].split("_")[0]
            line["symbol"] = source_accessor.get_symbol_for_genomic_position(chrom, pos)
            annotated_data.append(line)
            continue
        annotated_data.append(line)
    return annotated_data

def get_gene_symbols(data):
    genes = set()
    for line in data:
        if "symbol" in line:
            genes.add(line["symbol"])
    return list(genes)

def simulate_parameters(data, iterations: int = 1000, coeff_column_name: str = 'beta'):
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
        percentiles[str(i)] = str(numpy.percentile(value_list, i))
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
            variant["effect_size"] = float(snp[effect_size_field])
            if "symbol" in snp:
                variant["symbol"] = snp["symbol"]
            for field in included_fields_list:
                variant[field] = snp[field]
            variants[snp[id_field]] = variant

        model = {"score_model": {"categories": categories, "variants": variants}}
        model_file.write(yaml.dump(model, indent=2))

        description = {"description": description}
        model_file.write(yaml.dump(description, indent=2, default_flow_style=False))

    return
