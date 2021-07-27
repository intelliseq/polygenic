import logging
import re
import subprocess
from typing import List
from typing import Dict
from polygenic.lib.data_access.dto import SnpDataManySamples


SEARCH_CMD_PATTERN = 'bcftools view {} {}:{}'
GENOTYPE_STRING_RE = re.compile('\\sGT(:[A-Z]{2,3})*(.+)')

logger = logging.getLogger('description_language.' + __name__)


def extract_line_containing_rs_id(file_path: str, rsid_end: str, id_: str) -> str:
    """
    :param file_path:
    :param rsid_end: 3 last characters in the rsid
    :param id_: rsid without "rs"
    :return:
    """
    p = subprocess.Popen(SEARCH_CMD_PATTERN.format(file_path, rsid_end, id_), stdout=subprocess.PIPE, shell=True)
    output, _ = p.communicate()
    return output.decode("utf-8").strip().rsplit('\n', 1)[-1]

def get_genotypes(vcf_line: str, sample_names: List[str], genotype_string_re=GENOTYPE_STRING_RE) -> SnpDataManySamples:
    splitted = vcf_line.split(None)
    rsid = splitted[2]
    ref = splitted[3]
    alts = splitted[4].split(',')

    #print("SPL: " + str(splitted) + " rsid: " + str(rsid) + " ref: " + str(ref)) 
    #print("SPL-1: " + str(splitted[9:])) 
    #print("SEARCH: " + str(genotype_string_re.search(splitted[8:])))
    #print("GROUP: " + str(genotype_string_re.search(splitted[8:]).group(0)))

    #genotype_string_all_samples = genotype_string_re.search(splitted[-1]).group(2).strip()
    genotype_string_splitted = splitted[9:]#genotype_string_all_samples.split()
    assert len(genotype_string_splitted) == len(sample_names)
    genotype_strings = {name: genotype_string.split(':')[0] for name, genotype_string in
                        zip(sample_names, genotype_string_splitted)}
    return SnpDataManySamples(ref, alts, {sample_name: genotype_string_to_letters(ref, alts, genotype_string, rsid) for sample_name, genotype_string
            in genotype_strings.items()})


def genotype_string_to_letters(ref: str, alts: List[str], genotype_string: str, rsid: str) -> List[str]:
    # genotype_string - sth line 1|0
    if '|' in genotype_string:
        indices = genotype_string.split('|')
    else:
        indices = genotype_string.split('/')
    possible_genotypes = [ref] + alts
    genotype = []
    for index in indices:
        try:
            genotype.append(possible_genotypes[int(index)])
        except ValueError:
            logger.debug('Invalid genotype for {}: {}'.format(rsid, index))
    return genotype


if __name__ == '__main__':
    # print(get_sample_names('~/PycharmProjects/imputing/example_data/imputed_genotype'))
    # concatenate_vcfs('~/PycharmProjects/imputing/example_data/imputed_genotype', '/home/wojtek/logs/all.gz')
    # make_index('/home/wojtek/logs/all.gz')
    # print(os.listdir('/home/wojtek/logs/'))
    sample_names = get_sample_names('~/PycharmProjects/imputing/example_data/imputed_genotype')
    extracted = extract_line_containing_rs_id('/home/wojtek/logs/all.gz', 'rs568182971')
    print(extracted)
    print(get_genotypes(extracted, sample_names, 'rs568182971'))
