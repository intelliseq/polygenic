import os
from unittest import TestCase

# VcfAccessor
from polygenic.data.vcf_accessor import VcfAccessor

# CsvAccessor
from polygenic.data.csv_accessor import CsvAccessor

# DataAccessor
from polygenic.data.data_accessor import DataAccessor

# Utils
import polygenic.tools.utils as utils

class VcfAccessorTest(TestCase):

    def testGetRecordByPosition(self):
        vcf = VcfAccessor("polygenic/tests/resources/vcf/test-af.vcf.gz")
        record = vcf.get_record_by_position("chr1", "10147")
        self.assertEqual('C', record.get_ref())
    
    def testGetRecordByIdPositon(self):
        vcf = VcfAccessor("polygenic/tests/resources/vcf/test-af.vcf.gz")
        record = vcf.get_record_by_rsid("chr1:10147_C_A")
        self.assertEqual('C', record.get_ref())

    def testMultiallelicByPosition(self):
        vcf = VcfAccessor("polygenic/tests/resources/vcf/test-dbsnp.vcf.gz")
        record = vcf.get_record_by_rsid("1:164507787_A_AC")
        self.assertEqual('A', record.get_ref())

class DataAccessorTest(TestCase):

    def testGetGenotypeByRsid(self):
        vcf = VcfAccessor("polygenic/tests/resources/vcf/test-vcf-general.vcf.gz")
        data_accessor = DataAccessor(vcf)
        genotype = data_accessor.get_genotype_by_rsid("rs201694901")
        self.assertEqual(True, genotype["phased"])
        self.assertEqual('T', genotype["genotype"][0])
        genotype = data_accessor.get_genotype_by_rsid("rs1292226269")
        self.assertEqual(False, genotype["phased"])
        self.assertEqual('A', genotype["genotype"][0])
        genotype = data_accessor.get_genotype_by_rsid("rs1570391830")
        self.assertEqual(None, genotype["phased"])
        self.assertEqual(None, genotype["genotype"][0])

class CsvAccessorTest(TestCase):

    def testColumnNames(self):
        csv = CsvAccessor("polygenic/tests/resources/csv/test.csv")
        column_names = csv.get_column_names()
        self.assertTrue("ensembl_id" in column_names)
        gene = csv.get_symbol_for_genomic_position("1", 100133140)
        self.assertEqual("TRMT13", gene)

class UtilsTest(TestCase):

    def testDownloadGzip(self):
        url = "https://github.com/marpiech/scalable-genomics/raw/main/test/resources/archive/gzip.gz"
        path = "/tmp/gzip"
        return_path = utils.download(url, path)
        self.assertEqual(return_path, path)
        url = "https://github.com/marpiech/scalable-genomics/raw/main/test/resources/archive/gzip.bgz"
        path = "/tmp/gzip"
        return_path = utils.download(url, path)
        self.assertEqual(return_path, path)
    
    def testWriteData(self):
        test_directory = "/tmp/polygenic/test"
        os.makedirs(test_directory, exist_ok = True)
        path = test_directory + "/test.tsv"
        data = [{"a": 1, "b": 2, "d": 6}, {"a": 3, "c":"leleu", "d": 6}]
        utils.write_data(data, path)
        data = utils.read_table(path)
        self.assertEqual(data[0]["d"], "6")

    def testInvertNucleotide(self):
        nucleotide = 'A'
        self.assertEqual('T', utils.invert_nucleotides(nucleotide))
        nucleotides = ['A','G']
        self.assertEqual(['T', 'C'], utils.invert_nucleotides(nucleotides))

    def testLassoClump(self):
        test_directory = "/tmp/polygenic/test"
        test_file = "polygenic/tests/resources/tsv/snps.validated"
        utils.lasso_clump(test_file)
        