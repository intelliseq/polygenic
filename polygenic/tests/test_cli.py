from unittest import TestCase
from polygenic import polygenic
from polygenic import polygenicmaker
from polygenic import vcfstat

from pathlib import Path as path

import tabix
import os
import configparser

class PolygenicTest(TestCase):

    def test(self):
        self.assertEqual(1, 1)

    def testPolygenicYml(self):
        polygenic.main([
            #'--vcf', '/home/marpiech/data/gtcgenome/clusterd.id.vcf.gz',
            #'--vcf', '/home/marpiech/data/mygenome-bgi/bgi-genotype.nochr.vcf.gz',
            '--vcf', '/tmp/marpiech/kenobi/mygenome/bgi/bgi-genotype.nochr.rsid.vcf.gz',
            '--log-file', '/dev/null',
            #'--model', 'test/resources/model/gc_prodia.yml',  
            #'--model', 'test/resources/model/hirisplex.yml',
            '--model', 'polygenic/tests/resources/model/breast_prodia.yml',
            '--output-directory', '/tmp/polygenic',
            '--output-name-appendix', 'yml'])
            #"--vcf", "test/resources/vcf/my.vcf.gz",
            #"--sample", "yfyfy", 
            #"--log-file", "/dev/null",
            #"--model", "test/resources/model/variantset.yml",  
            #"--model", "test/resources/model/brest_prodia.yml",  
            #"--model", "test/resources/model/polygenicscore.yml",  
            #"--output-directory", "/tmp/polygenic",
            #"--output-name-appendix", "yml",
            #"--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')

    def testPolygenicMat(self):
        polygenic.main([
            "--vcf", "/tmp/marpiech/kenobi/polygenic/illu_merged-imputed.vcf.gz",
            "--model", "/tmp/marpiech/kenobi/polygenic/described-model.yml",
            "--output-directory", "/tmp/polygenic"])
        self.assertEqual('1', '1')

    def testPolygenicGc(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz",
            "--sample", "yfyfy", 
            "--log-file", "/dev/null",
            #"--model", "test/resources/model/variantset.yml",  
            #"--model", "test/resources/model/gc_prodia.yml",
            "--model", "test/resources/model/breast_prodia.yml",
            #"--model", "test/resources/model/polygenicscore.yml",  
            "--output-directory", "/tmp/polygenic",
            "--output-name-appendix", "yml",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')

    def testPolygenicCoreWithAf(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz",
            "--sample", "yfyfy", 
            "--log-file", "/dev/null",
            "--model", "test/resources/model/scaled_eas_model.py", 
            "--population", "eas", 
            "--output-directory", "/tmp/polygenic",
            "--output-name-appendix", "bambala",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')
    
    def testPolygenicForBiobankModel(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz", 
            "--log-file", "/dev/null",
            "--model", "test/resources/model/biomarkers-30600-both_sexes-irnt.tsv.py", 
            "--population", "eas", 
            "--output-directory", "/tmp/",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')

    def testPolygenicForGbeModel(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz", 
            "--log-file", "/dev/null",
            "--model", "results/model/BIN1210.py", 
            "--population", "nfe", 
            "--output-directory", "/tmp/",
            "--af", "test/resources/vcf/af.vcf.gz"])
        self.assertEqual('1', '1')

    def testPolygenicParameters(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz", 
            "--model", "test/resources/model/test-params.yml",
            "--parameters", "test/resources/json/test-params.json",
            "--output-directory", "/tmp/polygenic/",
            "--log-file", "/dev/null"])

    def testPolygenicGeneSymbol(self):
        polygenic.main([
            "--vcf", "test/resources/vcf/my.vcf.gz", 
            "--model", "test/resources/model/test-gene-symbol.yml",
            "--output-directory", "/tmp/polygenic/",
            "--log-file", "/dev/null"])

class PolygenicMakerTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(PolygenicMakerTest, self).__init__(*args, **kwargs)
        self.output_directory = "/tmp/polygenic/test"
        path(self.output_directory).mkdir(parents=True, exist_ok=True)

    def testTabix(self):
        config = configparser.ConfigParser()
        config.read(os.path.dirname(__file__) + "/../polygenic/polygenic.cfg")
        url = config['urls']['hg19-rsids']
        tb = tabix.open(url)
        print(type(tb).__name__)
        records = tb.query("16", 1650945, 1650946)
        for record in records:
            print(record[:3])

    def testGbeIndex(self):
         polygenicmaker.main([
             "gbe-index",
             "--output", "/tmp/polygenic/output"
         ])

    def testGbeModel(self):
        
        result_path = self.output_directory + '/HC710.yml'
        os.remove(result_path) if os.path.exists(result_path) else None

        polygenicmaker.main([
            "gbe-model",
            "--code", "HC710",
            "--gbe-index", "/tmp/marpiech/kenobi/polygenic/gbe-index.1.3.1.tsv",
            "--gene-positions", "/tmp/marpiech/kenobi/polygenic/ensembl-genes.104.tsv",
            "--source-ref-vcf", "/tmp/marpiech/kenobi/polygenic/dbsnp138.37.norm.vcf.gz",
            "--target-ref-vcf", "/tmp/marpiech/kenobi/polygenic/dbsnp138.38.norm.vcf.gz",
            "--af-vcf", "/tmp/marpiech/kenobi/polygenic/gnomad.3.1.vcf.gz",
            "--af-field", "AF_nfe",
            "--output-directory", self.output_directory
        ])
        with open(result_path, 'r') as output:
            data = output.read()
            header = list(filter(lambda line: "score_model:" in line, data.split('\n')))
            self.assertEqual(1, len(header))
            print(str(data))

    def testINI78(self):
        
        result_path = self.output_directory + '/INI78.yml'
        os.remove(result_path) if os.path.exists(result_path) else None

        polygenicmaker.main([
            "gbe-model",
            "--code", "INI78",
            "--gbe-index", "/tmp/marpiech/kenobi/polygenic/gbe-index.1.3.1.tsv",
            "--gene-positions", "/tmp/marpiech/kenobi/polygenic/ensembl-genes.104.tsv",
            "--source-ref-vcf", "/tmp/marpiech/kenobi/polygenic/dbsnp138.37.norm.vcf.gz",
            "--target-ref-vcf", "/tmp/marpiech/kenobi/polygenic/dbsnp138.38.norm.vcf.gz",
            "--af-vcf", "/tmp/marpiech/kenobi/polygenic/gnomad.3.1.vcf.gz",
            "--af-field", "AF_nfe",
            "--output-directory", self.output_directory
        ])
        
        with open(result_path, 'r') as output:
            data = output.read()
            header = list(filter(lambda line: "score_model:" in line, data.split('\n')))
            self.assertEqual(1, len(header))
            print(str(data))

    def testBiobankukIndex(self):
        polygenicmaker.main([
            "biobankuk-index",
            "--output-directory", "/tmp/polygenic/results"
        ])

    def testBiobankukModel(self):
        polygenicmaker.main([
            "biobankuk-model",
            "--code", "2395",
            "--sex", "both_sexes",
            "--coding", "4", 
            "--index", "/tmp/polygenic/results/panukbb_phenotype_manifest.tsv",
            "--output-directory", "/tmp/polygenic/results/model",
            "--variant-metrics", "/tmp/polygenic/results/full_variant_qc_metrics.txt",
            "--threshold", "1e-08",
            "--source-ref-vcf", "/tmp/polygenic/results/ALL.2of4intersection.20100804.genotypes.vcf.gz",
            "--target-ref-vcf", "/tmp/marpiech/kenobi/resources/GRCh38.dbSNP155.chr.norm.rsidonly.vcf.gz"
        ])

    def testPgsIndex(self):
        polygenicmaker.main([
            "pgs-index",
            "--output", "/tmp/polygenic"
        ])

    def testPgsGet(self):
        polygenicmaker.main([
            "pgs-get",
            "--code", "PGS000004",
            "--output-path", "/tmp/polygenic/PGS000004.txt"
        ])

    def testPgsPrepare(self):
        polygenicmaker.main([
            "pgs-prepare",
            "--input", "/tmp/polygenic/PGS000004.txt",
            "--output-path", "/tmp/polygenic/PGS000004.py",
            "--af", "/home/marpiech/data/af.vcf.gz",
            "--origin-reference-vcf", "/tmp/dbsnp/grch37/00-common_all.vcf.gz",
            "--model-reference-vcf", "/tmp/dbsnp/grch38/00-common_all.vcf.gz"
        ])

    def testVcfIndex(self):
        polygenicmaker.main([
            "vcf-index",
            "--vcf", "/tmp/dbsnp/grch38/00-common_all.vcf.gz"
        ])

class VcfstatTest(TestCase):

    def testBaf(self):
        vcfstat.main([
            "baf",
            "--vcf", "/home/marpiech/data/clustered_204800980122_R01C02.vcf.gz",
            "--output-directory", "/tmp/baf"
        ])
        self.assertEqual('1', '1')

    def testZygosity(self):
        vcfstat.main([
            "zygosity",
            "--vcf", "/home/marpiech/data/clustered_204800980122_R01C02.vcf.gz",
            "--output-file", "/tmp/baf/stats.json"
        ])
        self.assertEqual('1', '1')

class Debug(TestCase):

    def testDebug(self):
        polygenic.main([
            '--vcf', '/tmp/marpiech/tmp/B102_genotyped-by-vcf.vcf.gz',
            '--model', '/data/projects/BIOINFO-221/gbe-INI78-bone-density.yml',
            '--output-directory', '/tmp/polygenic'])