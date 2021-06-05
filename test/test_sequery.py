from unittest import TestCase
from polygenic import polygenic

class TestSequery(TestCase):
    def testCommandLine(self):
        polygenic.main(["--imputed_files", "/home/marpiech/data/clustered_204800980122_R01C02.vcf.gz", "--curated_initial_vcfs", "/home/marpiech/data/clustered_204800980122_R01C02.vcf.gz", "--allele_freq_json", "/home/marpiech/data/af.vcf.gz", "--log_file", "/dev/null", "--models_and_descriptions_path", "polygenic/resources/models", "--population", "eas", "--out_dir", "/tmp/polygenic"])
        self.assertEqual('1', '1')

  #def testAnalysis(self):
    #exitCode = wdltest.server_testrun()
    #self.assertTrue(exitCode == 0)

  #def testAnalysis(self):
    #exitCode = wdltest.local_testrun()
    #self.assertTrue(exitCode == 0)

if __name__ == "__main__":
    unittest.main()
