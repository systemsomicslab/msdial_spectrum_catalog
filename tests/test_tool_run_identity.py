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

    def test_main_v10_upgrade_preserves_all_rows_and_legacy_tool_id(self):
        # Read the committed public baseline, never another checkout's working files.
        schema_source = subprocess.check_output(["git", "show", "36f1fce:src/msdial_spectrum_catalog/schema.py"], text=True)
        annotation_source = subprocess.check_output(["git", "show", "36f1fce:src/msdial_spectrum_catalog/annotation.py"], text=True)
        schema_ns = {"__name__": "baseline_schema"}
        # Dataclasses resolve their defining module; use an ordinary existing module name.
        schema_ns["__name__"] = __name__
        exec(compile(schema_source, "baseline_schema", "exec"), schema_ns)
        annotation_ns = {"__name__": __name__, "__package__": "msdial_spectrum_catalog"}
        exec(compile(annotation_source, "baseline_annotation", "exec"), annotation_ns)
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "catalog.sqlite"
            with patch.object(storage, "SCHEMA_VERSION", 10), patch.object(storage, "SCHEMA_SQL", schema_ns["SCHEMA_SQL"]), patch.object(storage, "MIGRATIONS", schema_ns["MIGRATIONS"]), patch.object(storage, "AddColumn", schema_ns["AddColumn"]), patch.object(storage, "RunSql", schema_ns["RunSql"]):
                storage.initialize(db)
                with storage.transaction(db) as c:
                    c.execute("INSERT INTO study(study_id, repository, accession) VALUES ('s','synthetic','a')")
                    c.execute("INSERT INTO analysis_unit(analysis_unit_id,study_id,external_unit_id) VALUES ('u','s','u')")
                    c.execute("INSERT INTO analysis_run(run_id,analysis_unit_id,run_fingerprint,output_directory) VALUES ('r','u','fp','.')")
                    c.execute("INSERT INTO sample(sample_id,analysis_unit_id,sample_name) VALUES ('sample','u','invented')")
                old = annotation_ns["record_tool_run"](db, annotation_ns["ToolRunInput"]("invented", tool_version="1", parameters={"key": 4}), run_id="r")
                with closing(storage.connect(db)) as c:
                    before = {row[0]: ([col[1] for col in c.execute('PRAGMA table_info("' + row[0] + '")')],
                                       [tuple(r) for r in c.execute('SELECT * FROM "' + row[0] + '"')])
                              for row in c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name != 'catalog_meta'").fetchall()}
            storage.initialize(db)
            with closing(storage.connect(db)) as c:
                self.assertEqual(storage.read_schema_version(c), 13)
                for table, (columns, rows) in before.items():
                    actual = [tuple(r) for r in c.execute('SELECT ' + ','.join('"' + col + '"' for col in columns) + ' FROM "' + table + '"')]
                    for row in rows:
                        self.assertIn(row, actual, table)
                self.assertEqual(c.execute("SELECT ingest_granularity FROM analysis_run").fetchone()[0], "full")
                self.assertEqual(c.execute("PRAGMA foreign_key_check").fetchall(), [])
            new = record_tool_run(db, ToolRunInput("invented", tool_version="1", parameters={"key": 4}), run_id="r")
            self.assertEqual(old, new)
