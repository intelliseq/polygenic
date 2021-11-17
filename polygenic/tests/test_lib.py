from unittest import TestCase

# VcfAccessor
from polygenic.data.vcf_accessor import VcfAccessor

# Utils
from polygenic.tools.utils import download

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
