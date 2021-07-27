from unittest import TestCase
from polygenic.lib.data_access.vcf_accessor import VcfAccessor

class VcfAccessorTest(TestCase):

    def testGetRecordByPosition(self):
        vcf = VcfAccessor("test/resources/vcf/af.vcf.gz")
        record = vcf.get_record_by_position("chr22", "38936618")
        self.assertEqual('C', record.get_ref())