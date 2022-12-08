import os
from unittest import TestCase

from polygenic.data.polars_frame import PolarsFrame

class PolarsFrameTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(PolarsFrameTest, self).__init__(*args, **kwargs)
        os.makedirs("/tmp/polygenic/tests", exist_ok=True)

    def test_polarsframe(self):
        pf = PolarsFrame("polygenic/tests/resources/tsv/biobankuk-volume-10klines.tsv")
        self.assertTrue("chr" in pf.get_cols(), "Column 'chr' not found in the dataframe")
        pf.convert_to_gwas({"af": "af_meta_hq", "beta": "beta_meta_hq", "se": "se_meta_hq", "pvalue": "pval_meta_hq"})
        self.assertTrue("chromosome" in pf.get_cols(), "Column 'chromosome' not found in the dataframe " + str(pf.get_cols()) + " after gwas conversion")
        df = pf.get_dataframe()
        df.to_pandas().to_csv("/tmp/polygenic/tests/biobankuk-volume-10klines-polarsframe.tsv", sep = '\t', index=False)

    def test_largepolarsframe(self):
        pf = PolarsFrame("/home/marpiech/data/dat/continuous-IBil-both_sexes-irnt.tsv.gz")
        self.assertTrue("chr" in pf.get_cols(), "Column 'chr' not found in the dataframe")
    