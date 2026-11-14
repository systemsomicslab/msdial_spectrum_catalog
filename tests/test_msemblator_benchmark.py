"""Tests for the benchmark's helpers and its leak and citation self-check. No RDKit is needed."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "msemblator_benchmark.py"
_spec = importlib.util.spec_from_file_location("msemblator_benchmark", SCRIPT)
benchmark = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(benchmark)

CITATION = "Library-pos.msp (sha256 0123abcd)"
CLEAN = f"# Benchmark\n\nReference library: {CITATION}.\n\n| metric | rate |\n|---|---|\n| hit@1 | 0.5 |\n"


class SelfCheckTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def _check(self, text: str, identities=("Caffeine",)):
        path = self.base / "benchmark_summary.md"
        path.write_text(text, encoding="utf-8")
        return benchmark.self_check([path], list(identities), CITATION)

    def test_clean_report_passes(self):
        hits, uncited, _ = self._check(CLEAN + "| Caffeine | 1 |\n")
        self.assertEqual(hits, 0)
        self.assertEqual(uncited, [])

    def test_planted_absolute_library_path_is_a_hit(self):
        hits, _, detail = self._check(CLEAN + "library: X:\\private\\Library-pos.msp\n")
        self.assertGreater(hits, 0)
        self.assertGreater(detail["benchmark_summary.md"]["absolute_or_library_path"], 0)

    def test_planted_msp_peak_list_is_a_hit(self):
        hits, _, detail = self._check(CLEAN + "Num Peaks: 2\n100.0123\t999\n120.5\t10\n")
        self.assertGreater(detail["benchmark_summary.md"]["msp_record_text"], 0)
        self.assertGreater(hits, 0)

    def test_planted_account_detail_is_a_hit(self):
        hits, _, detail = self._check(CLEAN + "run by someone@example.org\n")
        self.assertGreater(detail["benchmark_summary.md"]["account_detail"], 0)
        self.assertGreater(hits, 0)

    def test_identity_without_citation_is_flagged(self):
        _, uncited, _ = self._check("| Caffeine | rank 1 |\n")
        self.assertEqual(uncited, ["benchmark_summary.md"])

    def test_numbers_in_a_table_are_not_peak_lists(self):
        hits, _, _ = self._check(CLEAN + "| structure hit@5 | 6 | 0.5 | 0.19-0.81 |\n")
        self.assertEqual(hits, 0)


class HelperTest(unittest.TestCase):
    def test_wilson(self):
        interval = benchmark.wilson(5, 10)
        self.assertEqual(interval["rate"], 0.5)
        self.assertAlmostEqual(interval["low"], 0.2366, places=3)
        self.assertAlmostEqual(interval["high"], 0.7634, places=3)
        self.assertIsNone(benchmark.wilson(0, 0)["rate"])

    def test_hill(self):
        self.assertEqual(benchmark.hill("H2O"), "H2O")
        self.assertEqual(benchmark.hill("O4C2H2"), "C2H2O4")
        self.assertEqual(benchmark.hill("C6H12O6"), "C6H12O6")
        self.assertEqual(benchmark.hill("ClNa"), "ClNa")
        self.assertIsNone(benchmark.hill(""))
        self.assertIsNone(benchmark.hill("C6H12O6+"))

    def test_reference_parent_fallback(self):
        full = {"AAAAAAAAAAAAAA-BBBBBBBBBB-N": ("Flavones", "CHEMONTID:1")}
        blocks = {"AAAAAAAAAAAAAA": {("Flavones", "CHEMONTID:1")},
                  "CCCCCCCCCCCCCC": {("A", "1"), ("B", "2")}}
        self.assertEqual(benchmark.reference_parent("AAAAAAAAAAAAAA-BBBBBBBBBB-N", full, blocks)[2], "full_inchikey")
        self.assertEqual(benchmark.reference_parent("AAAAAAAAAAAAAA-XXXXXXXXXX-N", full, blocks)[2],
                         "first_block_unique")
        self.assertIsNone(benchmark.reference_parent("CCCCCCCCCCCCCC-XXXXXXXXXX-N", full, blocks))

    def test_mz_bins(self):
        self.assertEqual(benchmark.mz_bin(299.9), "<300")
        self.assertEqual(benchmark.mz_bin(300.0), "300-650")
        self.assertEqual(benchmark.mz_bin(651.0), ">650")


if __name__ == "__main__":
    unittest.main()
