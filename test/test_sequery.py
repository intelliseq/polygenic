from unittest import TestCase
from polygenic import polygenic

class TestSequery(TestCase):
    def testCommandLine(self):
        polygenic.main(["--vcf", "/home/marpiech/data/clustered_204800980122_R01C02.vcf.gz", "--log_file", "/dev/null", "--model", "/home/marpiech/data/breast_cancer_eas_model.py", "--population", "eas", "--out_dir", "/tmp/polygenic"])
        self.assertEqual('1', '1')

#class TestSequeryShort(TestCase):
#    def testCommandLine(self):
#        polygenic.main(["--version"])
#        self.assertEqual('1', '1')

if __name__ == "__main__":
    unittest.main()
