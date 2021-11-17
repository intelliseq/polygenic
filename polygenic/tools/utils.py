import os
import logging
import sys
import subprocess
import urllib
import re

from polygenic.data.vcf_accessor import VcfAccessor

def error_print(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def error_exit(e):
    time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    error_print(" -> polygenic " + version + " failed at " + time)
    error_print(" -> with command: pgstk " + sys.argv.join(" "))
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

def clump(gwas_file, reference, clump_field = "pval_EUR", threshold = "1e-08"):

    filtered_path = gwas_file + ".filtered"
    clumped_path = gwas_file + ".clumped"

    subprocess.call("plink" +
                    " --clump " + filtered_path +
                    " --clump-p1 " + str(threshold) +
                    " --clump-r2 0.25 " +
                    " --clump-kb 1000 " +
                    " --clump-snp-field rsid " +
                    " --clump-field " + clump_field +
                    " --vcf " + reference + " " +
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
        header = line.rstrip().split(delimiter)
        while True:
            line = file.readline().rstrip().split(delimiter)
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


def validate_with_source(data, source_file):
    source_accessor = VcfAccessor(source_file)
    data = [validate(
        validated_line = line,
        validation_source = source_accessor) for line in data
    ]
    # data = [add_annotation(
    #      line, 
    #      annotation_name = "rsid", 
    #      annotation_source = source_vcf, 
    #      annotation_source_field = "rsid",
    #      default_value = "rs0") for line in data]