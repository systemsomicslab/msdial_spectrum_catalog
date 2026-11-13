import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.annotator_settings import read_annotator_settings
from msdial_spectrum_catalog.msemblator_export import annotator_tiers


class AnnotatorSettingsTests(unittest.TestCase):
    def test_projection_and_export_tiers_preserve_tsv_semantics(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            path = output / "msp_annotator_settings.tsv"
            path.write_text(
                'msp_file_path\tlibrary_kind\tannotator_id\tevidence_tier\n'
                'X:\\synthetic\\ignored.msp\tacquired\t"strict"\t"MSP high quality"\n'
                'ignored\t\trelaxed\tlow\n'
                'ignored\tpredicted\tcalculated\thigh\n'
                'ignored\t\tunknown\tother\n', encoding="utf-8")
            rows = read_annotator_settings(path)
            self.assertEqual(set(rows[0]), {"annotator_id", "evidence_tier", "library_kind"})
            self.assertEqual(rows[0]["annotator_id"], "strict")
            self.assertEqual(annotator_tiers(output), {
                "strict": "VS21 strict", "relaxed": "VS21 relaxed",
                "calculated": "GLACIER in-silico"})
            self.assertEqual(read_annotator_settings(output / "absent.tsv"), [])
