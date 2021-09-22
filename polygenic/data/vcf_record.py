import logging
from typing import List
from typing import Dict
from typing import Tuple
from polygenic.lib.polygenic_exception import PolygenicException

logger = logging.getLogger('description_language.' + __name__)

class VcfRecord(object):

    def __init__(self, vcf_line:str, sample_names:List[str] = []):
        super().__init__()
        self.__dict = self.__parse(vcf_line)
        self.__sample_names = sample_names

    def __parse(self, line):
        splitted_line = line.split("\t")

        parsed_line = {
             "CHROM": splitted_line[0],
             "POS": splitted_line[1],
             "ID": splitted_line[2],
             "REF": splitted_line[3],
             "ALT": splitted_line[4],
             "QUAL": splitted_line[5],
             "FILTER": splitted_line[6],
             "INFO": splitted_line[7]
        }

        if len(splitted_line) > 8:
             parsed_line["FORMAT"] = splitted_line[8]
             parsed_line["SAMPLES"] = [sample.rstrip() for sample in splitted_line[9:]]

        for field in parsed_line["INFO"].split(';'):
            if '=' in field:
                parsed_line[field.split('=')[0]] = field.split('=')[1]

        return parsed_line

    def get_genotype(self, sample_name) -> list:
        samples = self.__dict["SAMPLES"]
        if not self.__sample_names is None:
            idx = self.__sample_names.index(sample_name)
        else:
            idx = None
        if not samples is None and not idx is None:
            sample = samples[idx]
            sample = sample.split(":")[0]
            alleles = [sample[:1], sample[2:]]
            alleles = [self.recode_allele(allele) for allele in alleles]
            return alleles
        return None

    def recode_allele(self, allele) -> str:
        if int(allele) == 0:
            return self.get_ref()
        else:
            return self.get_alt()[int(allele) - 1]

    def get_as_dict(self) -> dict:
        return self.__dict

    def get_chrom(self) -> str:
        return self.__dict["CHROM"]

    def get_pos(self) -> str:
        return self.__dict["POS"]

    def get_id(self) -> str:
        return self.__dict["ID"]

    def get_alt(self) -> List[str]:
        return self.__dict["ALT"].split(",")
  
    def get_ref(self) -> str:
        return self.__dict["REF"]

    def get_info(self) -> str:
        return self.__dict["INFO"]

    def get_format(self) -> str:
        return self.__dict["FORMAT"]

    def is_imputed(self) -> bool:
        return self.get_format().find("DS") != -1

    def get_info_field(self, name) -> str:
        for field in self.get_info().split(";"):
            field_name = field.split("=")[0]
            if field_name == name:
                return field.split("=")[1]
        return None

    def get_af_by_pop(self, af_field_name) -> Dict[str, float]:
        if self.get_info_field(af_field_name) is None:
            raise RuntimeError("No {field} field in allele frequency vcf ".format(field = af_field_name))
        af = {}
        counter = 0
        sumfreq = 0
        for allele in self.get_alt():
            freq = float(self.get_info_field(af_field_name).split(",")[counter])
            sumfreq = sumfreq + freq
            af[allele] = freq
            counter = counter + 1
        af[self.get_ref()] = 1 - freq
        return af

    def get_genotype_by_af(self, af_field_name):
        allele_frequencies = self.get_af_by_pop(af_field_name)
        return list(max(self.__hardy_weinberg_diploid_frequencies(allele_frequencies).items(), key=lambda x: x[1])[0])

    def __hardy_weinberg_diploid_frequencies(self, allele_frequencies: Dict[str, float]) -> Dict[Tuple[str], float]:
        assert len(allele_frequencies) == 2
        allele0 = list(allele_frequencies.keys())[0]
        allele1 = list(allele_frequencies.keys())[1]
        return {
            (allele0, allele0): allele_frequencies[allele0] ** 2,
            (allele1, allele1): allele_frequencies[allele1] ** 2,
            (allele0, allele1): allele_frequencies[allele1] * allele_frequencies[allele0] * 2
        }