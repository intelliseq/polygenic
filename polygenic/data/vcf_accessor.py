import logging
import pathlib
from typing import Dict
from typing import List
from typing import Union
from polygenic.data import mobigen_utils
from polygenic.data.snp_data import SnpData
from polygenic.data.vcf_record import VcfRecord

# rsidx
import os
from gzip import open as gzopen
import polygenic.rsidx as rsidx
import sqlite3
import tabix

class VcfAccessor(object):
    def __init__(self, vcf_path:str):
        super().__init__()
        self.__path = vcf_path
        self.__is_remote = True if "://" in self.__path else False
        if not self.__is_remote and not os.path.exists(self.__path):
            raise RuntimeError('Can not access {path}'.format(path = self.__path))
        if not self.__is_remote and not os.path.exists(self.__path + '.tbi'):
            raise RuntimeError('Can not access tabix index for {path}'.format(path = self.__path))
        self.__tabix = tabix.open(self.__path)
        if not self.__is_remote and not os.path.exists(self.__path + '.idx.db'):
            logging.debug('Creating index for {path}'.format(path = self.__path))
            self.create_rsidx_index(self.__path)
        if not self.__is_remote:
            self.__sample_names = self.get_sample_names()
        else:
            self.__sample_names = []
        #self.__rsidx_conn = sqlite3.connect(self.__path + '.idx.db')
        self.__data: Dict[str, Dict[str:SnpData]] = {}  # dictionary rsid:{sample_name:ModelSnpData}

    ### create rsidx index
    @staticmethod
    def create_rsidx_index(path):
        with sqlite3.connect(path + '.idx.db') as vcf_index, gzopen(path, 'rt') as vcf_file:
            rsidx.index.index(vcf_index, vcf_file)

    def get_record_by_gnomadid(self, gnomadid) -> VcfRecord:
        splitted_gnomadid = gnomadid.split("-")
        chromosome = splitted_gnomadid[0]
        position = splitted_gnomadid[1]
        ref = splitted_gnomadid[2]
        alt = splitted_gnomadid[3]
        records = self.get_records_by_position(chromosome, position)
        if records:
            for record in records:
                if record.get_ref() == ref and alt in record.get_alt():
                    return record
        return None

    def get_record_by_position(self, chromosome, position) -> VcfRecord:
        records = self.get_records_by_position(chromosome, position)
        if records:
            return records[0]
        return None

    def get_records_by_position(self, chromosome, position) -> List:
        records = []
        try:
            records = self.__tabix.query(chromosome, int(position) - 1, int(position))
        except:
            try:
                if "chr" in chromosome:
                    records = self.__tabix.query(chromosome.replace("chr",""), int(position) - 1, int(position))
                else:
                    records = self.__tabix.query("chr" + chromosome, int(position) - 1, int(position))
            except:
                pass
        vcf_records = []
        if records:
            for record in records:
                vcf_record = VcfRecord("\t".join(record), self.__sample_names)
                if vcf_record.get_pos() == position:
                    vcf_records.append(vcf_record)
        return vcf_records

    def get_records_by_rsid(self, rsid) -> List:
        if ":" in rsid:
            position = (rsid.split("_")[0]).split(':')
            return self.get_records_by_position(position[0], position[1])
        if rsid.count("-") == 3:
            record = self.get_record_by_gnomadid(rsid)
            if record:
                return [record]
        vcf_records = []
        try:
            with sqlite3.connect(self.__path + '.idx.db') as dbconn:
                if "rs" in rsid:
                    for line in rsidx.search.search([rsid], dbconn, self.__path):
                        vcf_records.append(VcfRecord(line, self.__sample_names))
        except KeyError:
            pass
        return vcf_records

    def get_record_by_rsid(self, rsid, ref = "", alt = "", allow_invert = False) -> VcfRecord:
        records = self.get_records_by_rsid(rsid)
        if records:
            if "rs" in rsid:
                for record in records:
                    if record.get_id() == rsid:
                        return record
            if ":" in rsid and "_" in rsid:
                ref = rsid.split("_")[1]
                alt = rsid.split("_")[2]
                for record in records:
                    if ref == record.get_ref() and alt in record.get_alt():
                        return record
                    if allow_invert and alt == record.get_ref() and ref in record.get_alt():
                        return record
            if rsid.count("-") == 3:
                if records:
                    return records[0]                      
        return None

    def __get_record_for_rsid(self, rsid) -> VcfRecord:
        return VcfRecord(self.__get_vcf_line_for_rsid(rsid), self.__sample_names)

    def get_sample_names(self) -> List[str]:
        logging.debug('Getting sample names')
        sample_names_for_all_files = []
        with gzopen(self.__path) as vcf_file:
            for line in vcf_file:
                line = line.decode("utf-8")
                if line.find("#CHROM") != -1:
                    break
        if line.find('FORMAT') == -1:
            return None
        samples_string = line.split('FORMAT')[1].strip()
        sample_names_for_all_files.append(samples_string.split())
        assert all(sample_names == sample_names_for_all_files[-1] for sample_names in sample_names_for_all_files[:-1])
        return sample_names_for_all_files[-1]

    def __get_data_for_given_rsid(self, rsid, imputed:bool = False) -> Dict[str, SnpData]:
        line = self.__get_vcf_line_for_rsid(rsid)
        if not line:
            logging.debug(f'No line for rsid {rsid} found')
            raise DataNotPresentError
        if VcfRecord(line).is_imputed() == imputed:
            data = mobigen_utils.get_genotypes(line, self.__sample_names)
            self.__data[rsid] = {sample_name: SnpData(data.ref, data.alts, genotype) for sample_name, genotype in data.genotypes.items()}
        else:
            raise DataNotPresentError
        return self.__data[rsid]


    def get_af_by_pop(self, rsid:str, population_name:str) -> Dict[str, float]:
        return self.__get_record_for_rsid(rsid).get_af_by_pop(population_name)


    def get_data_for_sample(self, sample_name:str, rsid:str, imputed:bool = False) -> SnpData:
        try:
            return self.__data[rsid][sample_name]
        except KeyError:
            try:
                return self.__get_data_for_given_rsid(rsid, imputed)[sample_name]
            except DataNotPresentError:
                return None

    def __get_vcf_line_for_rsid(self, rsid:str) -> Union[None, str]:
        try:
            with sqlite3.connect(self.__path + '.idx.db') as dbconn:
                for line in rsidx.search.search([rsid], dbconn, self.__path):
                    return line
        except KeyError:
            print("Record " + str(rsid) + " not found")
            raise DataNotPresentError
        raise DataNotPresentError

    def get_allele_freq_from_db(rsid: str, population_name: str):
        record = self.__get_record_for_rsid(rsid)
        ref_allele = record.get_ref()
        alt_allele = record.get_alt()
        alt_allele_freq = record.get_af_by_pop(population_name)
        if not len(alt_allele) == 1:
            logging.info(
                f'{rsid} is multiallelic but only two alleles are provided. Only {ref_allele} and {alt_allele} were considered')
        return {alt_allele: alt_allele_freq, ref_allele: 1 - alt_allele_freq}

class DataNotPresentError(RuntimeError):
    pass

def path_to_fname_stem(path:str) -> str:
    return pathlib.PurePath(path).name.split('.')[0]