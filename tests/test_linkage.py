from contextlib import closing
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.linkage import export_sample_attributes, ingest_unit, read_unit_metadata
from ws1_fixture import write_unit


class LinkageTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.unit = write_unit(self.root)
        self.db = self.root / "catalog.sqlite"

    def test_unit_sample_attributes_and_tool_run(self):
        report = ingest_unit(self.db, self.unit, granularity="representative")
        self.assertTrue(report.valid, report.errors)
        self.assertEqual(report.metadata_matched, 2)
        with closing(sqlite3.connect(self.db)) as c, c:
            self.assertEqual(c.execute("SELECT metadata_matched_by FROM sample ORDER BY sample_name").fetchall(),
                             [("stem",), ("analysis_csv_basename",)])
            self.assertEqual(c.execute("SELECT repository_sample_id, class_id FROM sample ORDER BY sample_name").fetchall(),
                             [("sample_a", "Case"), ("sample_b.d", "Control")])
            self.assertEqual(c.execute("SELECT COUNT(*) FROM sample_attribute").fetchone()[0], 2)
            self.assertEqual(c.execute("SELECT COUNT(*) FROM run_metadata_source").fetchone()[0], 4)
            self.assertEqual(c.execute("SELECT COUNT(*) FROM artifact WHERE relative_path LIKE '%manifest%' OR relative_path LIKE '%metadata_reviewed%' OR relative_path LIKE '%annotator_settings%'").fetchone()[0], 0)
            self.assertEqual(c.execute("SELECT tool_name FROM annotation_tool_run").fetchone()[0], "msdial.console")
            self.assertEqual(c.execute("SELECT COUNT(*) FROM annotation_tool_run_library").fetchone()[0], 2)
            self.assertEqual(c.execute("SELECT target_omics_applied FROM analysis_unit").fetchone()[0], "Lipidomics")

    def test_identifier_disagreement_writes_nothing(self):
        path = self.unit / "campaign-record.json"
        data = json.loads(path.read_text())
        data["unit_key"] = "different"
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "disagree"):
            ingest_unit(self.db, self.unit)
        self.assertFalse(self.db.exists())

    def test_manual_disagreement_writes_nothing(self):
        with self.assertRaises(ValueError):
            ingest_unit(self.db, self.unit, expected={"accession": "OTHER"})
        self.assertFalse(self.db.exists())

    def test_whitelisted_path_refused_but_ignored_fields_not_stored(self):
        path = self.unit / "output/run-manifest.json"
        data = json.loads(path.read_text())
        data["annotation_pipeline_profile"] = "Q:/not-allowed"
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "path-like"):
            ingest_unit(self.db, self.unit)
        self.assertFalse(self.db.exists())

    def test_handoff_attributes_and_sources(self):
        path = self.root / "handoff.json"
        path.write_text(json.dumps(dict(schema="msdial-repository-reanalysis-handoff.v1",
            unit_attributes={"species": "Invented"}, technical_settings={"collision": [10, 20]}, ignored="X:/ignored")))
        ingest_unit(self.db, self.unit, handoff=path, granularity="representative")
        with closing(sqlite3.connect(self.db)) as c, c:
            self.assertEqual(c.execute("SELECT COUNT(*) FROM analysis_unit_attribute WHERE source LIKE 'handoff.%'").fetchone()[0], 2)
            self.assertEqual(c.execute("SELECT relative_label FROM run_metadata_source WHERE role='campaign_handoff'").fetchone()[0], "handoff.json")

    def test_metadata_relink_preserves_run_and_export(self):
        first = ingest_unit(self.db, self.unit, granularity="representative")
        path = self.unit / "output/SYN001_sample_metadata_reviewed.tsv"
        path.write_text(path.read_text(encoding="utf-8-sig").replace("Invented species", "Another invented species"), encoding="utf-8-sig")
        second = ingest_unit(self.db, self.unit, granularity="representative")
        self.assertEqual(first.run_id, second.run_id)
        out = self.root / "attributes.tsv"
        self.assertEqual(export_sample_attributes(self.db, out, unit="unit-a"), 2)
        self.assertIn("Another invented species", out.read_text())

    def test_multiple_metadata_files_are_error(self):
        (self.unit / "output/OTHER_sample_metadata_reviewed.tsv").touch()
        with self.assertRaises(ValueError):
            ingest_unit(self.db, self.unit)
        self.assertFalse(self.db.exists())

    def test_missing_settings_and_metadata_warn(self):
        (self.unit / "output/msp_annotator_settings.tsv").unlink()
        (self.unit / "output/SYN001_sample_metadata_reviewed.tsv").unlink()
        report = ingest_unit(self.db, self.unit, granularity="representative")
        self.assertIn("annotator settings missing", report.warnings)
        self.assertIn("reviewed sample metadata missing", report.warnings)
