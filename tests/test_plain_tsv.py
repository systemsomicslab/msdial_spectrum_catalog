import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.parsers import read_tsv


class PlainTsvTests(unittest.TestCase):
    def test_long_spectrum_cells_and_literal_quotes_preserve_rows(self):
        # MS-DIAL's plain TSV exporter does not CSV-quote its fields. A long
        # spectrum cell must not prevent ingest even when only other columns
        # are needed. A literal opening quote must not consume subsequent rows.
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "Invented.tsv"
            long_cell = "100:1 " * 40000
            for preamble, prefix in (("", None), ("metadata\tignored\n", "Peak ID\t")):
                with self.subTest(prefix=prefix):
                    path.write_text(preamble + "Peak ID\tName\tSpectrum\tEmpty\n"
                                    '1\t"Invented opening quote\t' + long_cell + "\t\n"
                                    '2\t"Invented closed quote"\t50:1\t\n', encoding="utf-8")
                    rows = list(read_tsv(path, prefix))
                    self.assertEqual(len(rows), 2)
                    self.assertEqual(rows[0][1]["Spectrum"], long_cell)
                    self.assertEqual(rows[0][1]["Name"], '"Invented opening quote')
                    self.assertEqual(rows[1][1]["Name"], '"Invented closed quote"')
                    self.assertEqual(rows[1][1]["Empty"], "")
