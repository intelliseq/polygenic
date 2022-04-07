from unittest import TestCase
from polygenic import pgstk
from io import StringIO as io
from unittest.mock import patch

class PgstkTest(TestCase):

    def testHelp(self):

        ### test help
        with patch('sys.stdout', new=io()) as helpOutput:
            with self.assertRaises(SystemExit):
                pgstk.main([
                    "--help"
                ])
            self.assertIn("pgstk", helpOutput.getvalue())