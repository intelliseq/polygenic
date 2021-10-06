import gzip
import io
import logging
import os
import os.path
import progressbar
import sys
import urllib.request
import subprocess

def merge(old_dict, new_dict):
    if new_dict is None:
        return old_dict
    for item in new_dict:
        if item in old_dict:
            if type(new_dict[item]) is int or type(new_dict[item]) is float or type(new_dict[item]) is list:
                old_dict[item] += new_dict[item]
            if type(new_dict[item]) is dict:
                merge(old_dict[item], new_dict[item])
        else:
            old_dict[item] = new_dict[item]

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

def download(url: str, output_path: str, force: bool=False, progress: bool=False):
    """Downloads file from url

    Keyword arguments:
    url -- url to file
    output_path -- path to output file
    force -- flag whether to overwrite downloaded file
    progress -- flag whether to present progress
    """
    logger = logging.getLogger('utils')

    print(url)

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
    #elif ".bgz" in url:
    #    response_data = gzip.GzipFile(fileobj = response)
    #    file_size = progressbar.UnknownLength
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