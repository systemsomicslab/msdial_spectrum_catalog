from contextlib import closing
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog.annotation import ToolRunInput, record_tool_run
from msdial_spectrum_catalog.classes import (ClassPredictionInput, ClassTermInput,
    list_class_predictions, record_class_predictions, upsert_lineage)
from msdial_spectrum_catalog.linkage import ingest_unit
from msdial_spectrum_catalog.storage import connect
from ws1_fixture import write_unit


class ClassesTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.db = Path(self.tmp.name) / "catalog.sqlite"
        ingest_unit(self.db, write_unit(self.tmp.name), granularity="representative")
        self.tool = record_tool_run(self.db, ToolRunInput("invented.classifier"))
        with closing(connect(self.db)) as c:
            self.subject = c.execute("SELECT spectrum_id FROM spectrum LIMIT 1").fetchone()[0]
        self.lineage = [ClassTermInput("classyfire.chemont", "synthetic-1", level, "Invented " + level,
                                      "CHEMONTID:" + str(i).zfill(7))
                        for i, level in enumerate(("kingdom", "superclass", "class", "subclass", "level5"))]

    def test_lineage_predictions_posteriors_and_idempotency(self):
        prediction = ClassPredictionInput("spectrum", self.subject, self.lineage,
                                         {t.level: 0.91 for t in self.lineage}, formula_context="C2H6 [M+H]+",
                                         formula_rank=1, source_row_ref="synthetic.tsv:sha256:2")
        self.assertEqual(record_class_predictions(self.db, self.tool, [prediction]), 5)
        record_class_predictions(self.db, self.tool, [prediction])
        rows = list_class_predictions(self.db, self.subject)
        self.assertEqual(len(rows), 5)
        self.assertEqual([r["level"] for r in rows if r["is_most_specific"]], ["level5"])
        self.assertTrue(all(r["posterior"] == 0.91 for r in rows))
        with closing(connect(self.db)) as c:
            self.assertEqual(c.execute("SELECT COUNT(*) FROM class_term_parent").fetchone()[0], 4)

    def test_npclassifier_multiple_parent_edges(self):
        def term(level, name):
            return ClassTermInput("npclassifier", "unknown", level, name)
        for parent in ("Invented superclass A", "Invented superclass B"):
            upsert_lineage(self.db, [term("pathway", "Invented pathway"), term("superclass", parent), term("class", "Invented class")])
        with closing(connect(self.db)) as c:
            self.assertEqual(c.execute("SELECT COUNT(*) FROM class_term_parent WHERE term_id IN (SELECT term_id FROM class_term WHERE level='class')").fetchone()[0], 2)

    def test_invalid_levels_and_order(self):
        for lineage in (self.lineage[::-1], [ClassTermInput("npclassifier", "unknown", "invalid", "Invented")]):
            with self.assertRaises(ValueError):
                upsert_lineage(self.db, lineage)

    def test_invalid_posterior_and_subject_and_tool(self):
        for value in (1.2, -0.1, float("nan")):
            with self.assertRaises(ValueError):
                record_class_predictions(self.db, self.tool, [ClassPredictionInput("spectrum", self.subject, self.lineage, {"class": value})])
        with self.assertRaises(ValueError):
            record_class_predictions(self.db, "missing", [])
        with self.assertRaises(ValueError):
            record_class_predictions(self.db, self.tool, [ClassPredictionInput("spectrum", "missing", self.lineage, {})])
