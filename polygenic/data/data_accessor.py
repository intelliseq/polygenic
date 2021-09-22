import logging
from polygenic.data.vcf_accessor import VcfAccessor
from polygenic.data.vcf_record import VcfRecord

logger = logging.getLogger('description_language.' + __name__)

class DataAccessor(object):
    def __init__(self, 
        genotypes: VcfAccessor,
        imputed_genotypes: VcfAccessor,
        allele_frequencies: VcfAccessor,
        sample_name: str,
        af_field_name: str = "AF_nfe",
        parameters = {}):
        self.__genotypes = genotypes
        self.__imputed_genotypes = imputed_genotypes
        self.__allele_frequencies = allele_frequencies
        self.__sample_name = sample_name
        self.__af_field_name = af_field_name
        self.__parameters = parameters

    def get_parameters(self) -> dict:
        return(self.__parameters)

    def get_genotype_by_rsid(self, rsid) -> VcfRecord:
        genotype = {"rsid": rsid}
        if not self.__genotypes is None:
            record = self.__genotypes.get_record_by_rsid(rsid)
            if not record is None:
                genotype["genotype"] = record.get_genotype(self.__sample_name)
                genotype["source"] = "genotyping"
                genotype["ref"] = record.get_ref()
                return genotype
        if record is None and not self.__imputed_genotypes is None:
            record = self.__imputed_genotypes.get_record_by_rsid(rsid)
            if not record is None:
                genotype["genotype"] = record.get_genotype(self.__sample_name)
                genotype["source"] = "imputing"
                genotype["ref"] = record.get_ref()
                return genotype
        if record is None and not self.__allele_frequencies is None:
            record = self.__allele_frequencies.get_record_by_rsid(rsid)
            if not record is None:
                genotype["genotype"] = record.get_genotype_by_af(self.__af_field_name)
                genotype["source"] = "af"
                genotype["ref"] = record.get_ref()
                return genotype
        genotype["genotype"] = [None, None]
        genotype["source"] = "missing"
        genotype["ref"] = None
        return genotype

            

