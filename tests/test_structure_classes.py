from contextlib import closing
import sqlite3
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.catalog_build import privacy_scan
from msdial_spectrum_catalog.structure_classes import load_structure_classes


class StructureClassesTests(unittest.TestCase):
    def test_three_column_lookup_identity_no_path_and_idempotency(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Invented.icd"
            path.write_bytes((Path(__file__).parent / "fixtures/ws1/class_lookup.tsv").read_bytes())
            db = Path(tmp) / "catalog.sqlite"
            one = load_structure_classes(db, path)
            two = load_structure_classes(db, path)
            self.assertEqual((one["assignments_added"], two["assignments_added"]), (2, 0))
            self.assertEqual(two["terms_added"], 0)
            self.assertNotIn(str(path), repr(two))
            with closing(sqlite3.connect(db)) as c, c:
                self.assertEqual(c.execute("SELECT COUNT(*) FROM annotation_tool_run").fetchone()[0], 1)
                self.assertEqual(c.execute("SELECT library_name, sha256 FROM annotation_tool_run_library").fetchone(), ("Invented.icd", one["sha256"]))
                self.assertEqual(c.execute("SELECT DISTINCT source_kind FROM structure_class_assignment").fetchone()[0], "lookup")
                self.assertNotIn(str(path), "\n".join(c.iterdump()))
            self.assertTrue(privacy_scan(db)["valid"])

    def test_bad_layout_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Invented.icd"
            path.write_text("wrong\n")
            db = Path(tmp) / "catalog.sqlite"
            with self.assertRaises(ValueError):
                load_structure_classes(db, path)
            self.assertFalse(db.exists())

    def test_emitted_lineage_can_supply_missing_lookup_depth(self):
        from msdial_spectrum_catalog.classes import ClassTermInput, upsert_lineage
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "Invented.icd"
            path.write_bytes((Path(__file__).parent / "fixtures/ws1/class_lookup.tsv").read_bytes())
            db = Path(tmp) / "catalog.sqlite"
            load_structure_classes(db, path)
            upsert_lineage(db, [ClassTermInput("classyfire.chemont", "unknown", "class", "Invented class", "CHEMONTID:9999999")])
            load_structure_classes(db, path)
            with closing(sqlite3.connect(db)) as c:
                self.assertEqual(c.execute("SELECT level FROM class_term").fetchone()[0], "class")
                self.assertEqual(c.execute("SELECT COUNT(*) FROM structure_class_assignment").fetchone()[0], 2)
