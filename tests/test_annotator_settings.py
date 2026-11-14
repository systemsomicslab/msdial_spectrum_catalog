import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

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
            self.assertEqual(set(rows[0]), {"annotator_id", "evidence_tier", "library_kind", "tier"})
            self.assertEqual(rows[0]["annotator_id"], "strict")
            self.assertEqual(annotator_tiers(output), {
                "strict": "VS21 strict", "relaxed": "VS21 relaxed",
                "calculated": "GLACIER in-silico"})
            self.assertEqual(read_annotator_settings(output / "absent.tsv"), [])

    def test_export_never_hashes_or_opens_settings_outside_shared_reader(self):
        from test_msemblator_export import MsemblatorExportTest
        from msdial_spectrum_catalog import msemblator_export
        fixture = MsemblatorExportTest()
        fixture.setUp()
        self.addCleanup(fixture.tearDown)
        settings = fixture.output / "msp_annotator_settings.tsv"
        original = Path.open
        reader_active = False
        shared = msemblator_export.read_annotator_settings

        def read(path):
            nonlocal reader_active
            reader_active = True
            try:
                return shared(path)
            finally:
                reader_active = False

        def checked_open(path, *args, **kwargs):
            if path == settings:
                self.assertTrue(reader_active, "settings opened outside shared reader")
            return original(path, *args, **kwargs)

        with patch.object(msemblator_export, "read_annotator_settings", side_effect=read), \
                patch.object(Path, "open", checked_open):
            summary, _, _ = fixture._export("unresolved")
        self.assertNotIn("msp_annotator_settings.tsv", summary["sources"])
