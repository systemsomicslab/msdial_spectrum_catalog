from contextlib import closing
import sqlite3
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.linkage import ingest_unit, match_annotator_library, read_annotator_settings
from ws1_fixture import write_unit


class RunAnnotatorTests(unittest.TestCase):
    def test_allowlisted_reader_tiers_sentinel_and_unresolved(self):
        with tempfile.TemporaryDirectory() as tmp:
            unit = write_unit(tmp)
            rows = read_annotator_settings(unit / "output/msp_annotator_settings.tsv")
            self.assertEqual([r["tier"] for r in rows], ["strict", "relaxed", "in_silico", "unknown"])
            self.assertEqual(rows[1]["library_kind"], "acquired")
            self.assertEqual(rows[3]["library_kind"], "")
            self.assertNotIn("SENTINEL_DO_NOT_READ", repr(rows))
            db = Path(tmp) / "catalog.sqlite"
            report = ingest_unit(db, unit, granularity="representative")
            self.assertEqual(report.library_identity_unresolved, 1)
            self.assertNotIn("SENTINEL_DO_NOT_READ", repr(report))
            with closing(sqlite3.connect(db)) as c, c:
                self.assertEqual(c.execute("SELECT COUNT(*) FROM run_annotator").fetchone()[0], 4)
                self.assertEqual(c.execute("SELECT COUNT(*) FROM run_annotator WHERE library_sha256 IS NULL").fetchone()[0], 1)
                self.assertNotIn("SENTINEL_DO_NOT_READ", "\n".join(c.iterdump()))

    def test_name_role_matching_and_ambiguity(self):
        lib = dict(name="Invented.msp", role="msp", sha256="a" * 64, library_kind="acquired")
        self.assertEqual(match_annotator_library("MSP:Invented_high", [lib]), lib)
        self.assertIsNone(match_annotator_library("MSP:Invented_high", [lib, lib]))
        self.assertIsNone(match_annotator_library("msp_high_quality", [lib, lib]))
        self.assertIsNone(match_annotator_library("LbmDB:Invented", [lib]))
