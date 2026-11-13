from contextlib import closing
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from msdial_spectrum_catalog.catalog_build import ingest_devset, privacy_scan, size_probe
from msdial_spectrum_catalog.devset import canonical_digest
from msdial_spectrum_catalog.storage import initialize
from ws1_fixture import write_unit


class CatalogBuildTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def selection(self, unit):
        content = dict(campaigns=[], exclusions=[], units=[dict(unit_key=unit.name, workspace=str(unit), campaign_id="synthetic-campaign", handoff_path=None)])
        selection = dict(schema="msdial-annotation-devset.v1", **content, content_sha256=canonical_digest(content))
        path = self.root / "devset.json"
        path.write_text(json.dumps(selection))
        return path

    def test_private_name_refusal_writes_nothing(self):
        unit = write_unit(self.root, private=True)
        selection = self.selection(unit)
        target = self.root / "never-created/catalog.sqlite"
        with self.assertRaisesRegex(ValueError, "suffix"):
            ingest_devset(target, selection)
        self.assertFalse(target.parent.exists())

    def test_build_records_reports_pin_and_validation(self):
        unit = write_unit(self.root)
        target = self.root / "build/catalog.sqlite"
        with patch("msdial_spectrum_catalog.catalog_build.git_commit", return_value="synthetic"):
            result = ingest_devset(target, self.selection(unit))
        self.assertEqual(result["failures"], 0, result["units"])
        self.assertTrue(result["units"][0]["valid"])
        self.assertEqual(result["table_counts"]["analysis_run"], 1)
        run = json.loads((target.parent / "run.json").read_text())
        self.assertEqual(run["git_commit"], "synthetic")
        self.assertIn("pin_worktree", run)
        self.assertTrue(run["local_only"])

    def test_one_unit_failure_is_isolated(self):
        unit = write_unit(self.root)
        selection = self.selection(unit)
        data = json.loads(selection.read_text())
        second = write_unit(self.root, key="unit-b")
        data["units"].append(dict(unit_key=second.name, workspace=str(second), campaign_id="synthetic-campaign"))
        data["content_sha256"] = canonical_digest({key: data[key] for key in ("campaigns", "exclusions", "units")})
        selection.write_text(json.dumps(data))
        from msdial_spectrum_catalog.linkage import ingest_unit
        def run(database, unit_dir, **kwargs):
            if Path(unit_dir).name == "unit-a":
                raise RuntimeError("SENTINEL_PRIVATE_VALUE")
            return ingest_unit(database, unit_dir, **kwargs)
        with patch("msdial_spectrum_catalog.catalog_build.git_commit", return_value="synthetic"), patch("msdial_spectrum_catalog.catalog_build.ingest_unit", side_effect=run):
            result = ingest_devset(self.root / "build/catalog.sqlite", selection)
        self.assertEqual(result["failures"], 1)
        self.assertEqual(result["table_counts"]["analysis_run"], 1)
        self.assertNotIn("SENTINEL_PRIVATE_VALUE", repr(result))

    def test_probe_four_variants_counts_vacuum_and_cleanup(self):
        unit = write_unit(self.root)
        out = self.root / "probe/size-report.json"
        with patch("msdial_spectrum_catalog.catalog_build.git_commit", return_value="synthetic"):
            result = size_probe([unit], out)
        variants = result["units"][0]["variants"]
        self.assertEqual(len(variants), 4)
        self.assertTrue(all(v["valid"] for v in variants))
        self.assertTrue(all(v["bytes"] > 0 and v["row_counts"]["analysis_run"] == 1 for v in variants))
        self.assertEqual(variants[2]["row_counts"]["alignment_member"], 0)
        self.assertEqual(variants[3]["row_counts"]["mztab_record"], 0)
        self.assertEqual(list(out.parent.glob("*.sqlite*")), [])

    def test_privacy_scan_all_text_columns_and_exemptions(self):
        db = self.root / "catalog.sqlite"
        initialize(db)
        with closing(sqlite3.connect(db)) as c, c:
            c.execute("INSERT INTO study(study_id,repository,accession,title) VALUES ('s','synthetic','a',?)", ("X:\\private\\Invented.msp",))
        scan = privacy_scan(db)
        self.assertEqual(scan["hits"], 1)
        self.assertEqual(scan["columns"], {"study.title": 1})
        self.assertNotIn("private", repr(scan))
        with closing(sqlite3.connect(db)) as c, c:
            c.execute("UPDATE study SET title='Invented'")
            c.execute("INSERT INTO analysis_unit(analysis_unit_id,study_id,external_unit_id) VALUES ('u','s','u')")
            c.execute("INSERT INTO analysis_run(run_id,analysis_unit_id,run_fingerprint,output_directory) VALUES ('r','u','fp',?)", ("X:\\exempt\\Invented.msp",))
        self.assertTrue(privacy_scan(db)["valid"])

    def test_analysis_output_refused(self):
        with self.assertRaisesRegex(ValueError, "analysis root"):
            ingest_devset(self.root / "MSDIAL_Public_Reanalysis/analysis/catalog.sqlite", "missing.json")
