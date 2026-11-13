from contextlib import closing
import json
import sqlite3
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from msdial_spectrum_catalog import storage
from msdial_spectrum_catalog.annotation import LibraryIdentity, ToolRunInput, record_tool_run
from msdial_spectrum_catalog.schema import MIGRATIONS, SCHEMA_SQL


class ToolRunIdentityTests(unittest.TestCase):
    def test_reject_library_paths_and_invalid_hashes(self):
        for name, digest in (("D:\\x\\y.msp", "a" * 64), ("/a/b.msp", "a" * 64),
                             ("D:relative.msp", "a" * 64), ("invented.msp", "g" * 64),
                             ("invented.msp", "A" * 64)):
            with self.subTest(name=name), self.assertRaises(ValueError):
                LibraryIdentity(name, digest)

    def test_libraries_kernel_scope_roundtrip_and_order_independence(self):
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "catalog.sqlite"
            storage.initialize(db)
            with storage.transaction(db) as c:
                c.execute("INSERT INTO study(study_id, repository, accession) VALUES ('s','synthetic','a')")
                c.execute("INSERT INTO analysis_unit(analysis_unit_id,study_id,external_unit_id) VALUES ('u','s','u')")
                for run in ("r1", "r2"):
                    c.execute("INSERT INTO analysis_run(run_id,analysis_unit_id,run_fingerprint,output_directory) VALUES (?,'u','fingerprint','.')", (run,))
            libraries = (LibraryIdentity("Invented-a.msp", "a" * 64, role="msp", distribution="public"),
                         LibraryIdentity("Invented-b.msp", "b" * 64, role="msp", library_kind="predicted"))
            one = record_tool_run(db, ToolRunInput("invented", libraries=libraries, kernel_identity="invented.cosine?tol=0.01", scope_run_ids=("r1", "r2")))
            two = record_tool_run(db, ToolRunInput("invented", libraries=libraries[::-1], kernel_identity="invented.cosine?tol=0.01", scope_run_ids=("r2", "r1")), run_id="r1")
            self.assertEqual(one, two)
            with closing(storage.connect(db)) as c:
                row = c.execute("SELECT run_id, kernel_identity FROM annotation_tool_run").fetchone()
                self.assertEqual(tuple(row), (None, "invented.cosine?tol=0.01"))
                self.assertEqual(c.execute("SELECT COUNT(*) FROM annotation_tool_run_library").fetchone()[0], 2)
                self.assertEqual(c.execute("SELECT COUNT(*) FROM annotation_tool_run_scope").fetchone()[0], 2)
