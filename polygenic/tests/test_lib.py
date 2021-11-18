import os
from unittest import TestCase

# VcfAccessor
from polygenic.data.vcf_accessor import VcfAccessor

# Utils
from polygenic.tools.utils import download
from polygenic.tools.utils import read_table
from polygenic.tools.utils import write_data

class VcfAccessorTest(TestCase):

    def testGetRecordByPosition(self):
        vcf = VcfAccessor("polygenic/tests/resources/vcf/test.af.vcf.gz")
        record = vcf.get_record_by_position("chr1", "10147")
        self.assertEqual('C', record.get_ref())
    
    def testGetRecordByIdPositon(self):
        vcf = VcfAccessor("polygenic/tests/resources/vcf/test.af.vcf.gz")
        record = vcf.get_record_by_rsid("chr1:10147_C_A")
        self.assertEqual('C', record.get_ref())

    def testMultiallelicByPosition(self):
        vcf = VcfAccessor("polygenic/tests/resources/vcf/test.dbsnp.vcf.gz")
        record = vcf.get_record_by_rsid("1:164507787_A_AC")
        self.assertEqual('A', record.get_ref())

class UtilsTest(TestCase):

    def testDownloadGzip(self):
        url = "https://github.com/marpiech/scalable-genomics/raw/main/test/resources/archive/gzip.gz"
        path = "/tmp/gzip"
        return_path = download(url, path)
        self.assertEqual(return_path, path)
        url = "https://github.com/marpiech/scalable-genomics/raw/main/test/resources/archive/gzip.bgz"
        path = "/tmp/gzip"
        return_path = download(url, path)
        self.assertEqual(return_path, path)
    
    def testWriteData(self):
        test_directory = "/tmp/polygenic/test"
        os.makedirs(test_directory, exist_ok = True)
        path = test_directory + "/test.tsv"
        data = [{"a": 1, "b": 2, "d": 6}, {"a": 3, "c":"leleu", "d": 6}]
        write_data(data, path)
        data = read_table(path)
        self.assertEqual(data[0]["d"], "6")
