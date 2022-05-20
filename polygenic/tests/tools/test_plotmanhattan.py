from unittest import TestCase
from polygenic import pgstk
from pathlib import Path as path

class PlotManhattanTest(TestCase):

    def __init__(self, *args, **kwargs):
        super(PlotManhattanTest, self).__init__(*args, **kwargs)
        self.output_directory = "/tmp/polygenic/test"
        path(self.output_directory).mkdir(parents=True, exist_ok=True)

    def test_plot_manhattan(self):
        pgstk.main([
            "--log-level", "INFO",
            "--log-stdout",
            "plot-manhattan",
            "--tsv", "/home/marpiech/data/tsv/continuous-IBil-both_sexes-irnt.tsv.gz",
            #"--tsv", "polygenic/tests/resources/tsv/biobankuk-test.tsv.gz",
            "--output", f"{self.output_directory}/biobankuk-test.png"
        ])