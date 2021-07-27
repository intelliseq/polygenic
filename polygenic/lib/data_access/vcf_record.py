import logging
from typing import List
from typing import Dict
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
             "chrom": splitted_line[0],
             "pos": splitted_line[1],
             "id": splitted_line[2],
             "ref": splitted_line[3],
             "alt": splitted_line[4],
             "qual": splitted_line[5],
             "filter": splitted_line[6],
             "info": splitted_line[7]
        }
        if len(splitted_line) > 8:
             parsed_line["format"]: splitted_line[8]
             parsed_line["samples"]: splitted_line[9:]
        return parsed_line

    def get_chrom(self) -> str:
        return self.__dict["chrom"]

    def get_pos(self) -> str:
        return self.__dict["pos"]

    def get_id(self) -> str:
        return self.__dict["id"]

    def get_alt(self) -> List[str]:
        return self.__dict["alt"].split(",")
  
    def get_ref(self) -> str:
        return self.__dict["ref"]

    def get_info(self) -> str:
        return self.__dict["info"]

    def get_format(self) -> str:
        return self.__dict["format"]

    def is_imputed(self) -> bool:
        return self.get_format().find("DS") != -1

    def get_info_field(self, name) -> str:
        for field in self.get_info().split(";"):
            field_name = field.split("=")[0]
            if field_name == name:
                return field.split("=")[1]

    def get_af_by_pop(self, population_name) -> Dict[str, float]:
        af = {}
        counter = 0
        sumfreq = 0
        for allele in self.get_alt():
            #print(population_name)
            #print(self.get_info())
            freq = float(self.get_info_field(population_name).split(",")[counter])
            sumfreq = sumfreq + freq
            af[allele] = freq
            counter = counter + 1
        af[self.get_ref()] = 1 - freq
        return af