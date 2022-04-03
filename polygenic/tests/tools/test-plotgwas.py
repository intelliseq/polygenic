from unittest import TestCase
from polygenic import pgstk

class PlotGwasTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(PlotGwasTest, self).__init__(*args, **kwargs)
        self.output_directory = "/tmp/polygenic/test"
        path(self.output_directory).mkdir(parents=True, exist_ok=True)

    def testPlotGwas(self):
        pgstk.main([
            "model-gbe",
            "--code", "2395",
            "--sex", "both_sexes",
            "--coding", "4", 
            "--output-directory", self.output_directory,
            "--pvalue-threshold", "1e-08",
            "--clumping-vcf", "polygenic/tests/resources/largefiles/eur.phase3.biobank.set.vcf.gz",
            "--source-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch37.norm.vcf.gz",
            "--target-ref-vcf", "polygenic/tests/resources/largefiles/dbsnp155.grch38.norm.vcf.gz",
            "--gene-positions", "polygenic/tests/resources/largefiles/ensembl-genes.104.tsv"
        ])

        result_path = self.output_directory + "/biobankuk-2395-both_sexes-4-hair_balding_pattern-EUR-1e-08.yml"

        with open(result_path, 'r') as output:
            data = output.read()
            header = list(filter(lambda line: "score_model:" in line, data.split('\n')))
            self.assertEqual(1, len(header))