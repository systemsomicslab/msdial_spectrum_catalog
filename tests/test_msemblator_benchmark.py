"""Tests for the benchmark's helpers and its leak and citation self-check. No RDKit is needed."""

from __future__ import annotations

import importlib.util
import math
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

    def test_fallback_text_is_not_a_citation(self):
        path = self.base / "benchmark_summary.md"
        path.write_text("Reference library: reference library not recorded.\n| Caffeine | 1 |\n", encoding="utf-8")
        _, uncited, _ = benchmark.self_check([path], ["Caffeine"], None)
        self.assertEqual(uncited, ["benchmark_summary.md"])
        path.write_text("| hit@1 | 0.5 |\n", encoding="utf-8")  # no identity named: nothing to cite
        self.assertEqual(benchmark.self_check([path], ["Caffeine"], None)[1], [])

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

    def test_chunk_exclusion(self):
        self.assertEqual(benchmark.chunk_exclusion({"exit_code": 0, "failed_steps": []}), "")
        self.assertIn("formula:msbuddy", benchmark.chunk_exclusion({"exit_code": 1, "failed_steps": ["formula:msbuddy"]}))
        self.assertTrue(benchmark.chunk_exclusion({"exit_code": 2, "failed_steps": None}).startswith("driver exit 2"))
        self.assertTrue(benchmark.chunk_exclusion({"exit_code": 0, "failed_steps": ["x"]}).startswith("degraded"))


def row(name, **values):
    base = {"name": name, "unit_key": "u1", "ion_mode": "Positive", "instrument_family": "QTOF", "mz_bin": "<300",
            "excluded": "", "formula_evaluable": True, "structure_evaluable": True, "reference_in_db": True,
            "formula_returned": True, "structure_returned": True, "structure_step": "sent",
            "rank_formula": 1, "rank_structure": 1, "sirius_formula_rank": 1, "sirius_structure_rank": 1,
            "top1_structure_margin": 0.5, "reference_direct_parent": "", "canopus_c1": False,
            "canopus_in_all_classifications": False, **{f"canopus_c2_p{p}": False for p in benchmark.CANOPUS_P}}
    return base | values


class SummaryTest(unittest.TestCase):
    CONTEXT = {"plan": {"run_id": "r"}, "databases": [], "canopus_columns": {}, "icd": {}, "not_run": 0,
               "excluded_chunks": [{"chunk_id": "u1/benchmark/c001", "exit_code": 1,
                                    "failed_steps": ["formula:msbuddy"], "spectra": 1, "reason": "driver exit 1"}]}

    def test_excluded_chunks_and_unevaluable_references_stay_out_of_the_rates(self):
        rows = [row("hit"),
                row("miss", rank_formula=math.inf, rank_structure=math.inf),
                row("failed chunk", excluded="driver exit 1", rank_formula=math.inf, rank_structure=math.inf),
                row("charged reference", formula_evaluable=False, rank_formula=math.inf),
                row("no inchikey", structure_evaluable=False, rank_structure=math.inf)]
        summary = benchmark.summarise(rows, self.CONTEXT, "b", {"text": "lib", "libraries": []})
        structure = summary["structure"]["denominators"]["N_sent"]["hit@1"]
        formula = summary["formula"]["denominators"]["N_sent"]["hit@1"]
        self.assertEqual((structure["hits"], structure["n"]), (2, 3))
        self.assertEqual((formula["hits"], formula["n"]), (2, 3))
        counts = summary["counts"]
        self.assertEqual(counts["benchmark_spectra_evaluated"], 4)
        self.assertEqual(counts["benchmark_spectra_in_excluded_chunks"], 1)
        self.assertEqual(counts["formula_reference_unevaluable"], 1)
        self.assertEqual(counts["structure_reference_unevaluable"], 1)
        self.assertEqual(summary["excluded_chunks"][0]["chunk_id"], "u1/benchmark/c001")
        self.assertIn("u1/benchmark/c001", benchmark.markdown(summary, None))
        self.assertEqual(summary["counts"]["icd_coverage"]["n"], 3)  # the keyless reference is left out
        self.assertEqual(summary["unit_macro_average"]["structure_hit@1"], {"average": 0.6667, "units": 1})
        cited = [r | {"reference_name": r["name"], "reference_inchikey": "", "canopus_most_specific": ""}
                 for r in rows]
        table = benchmark.markdown(summary, cited)
        self.assertIn("| excluded |", table)
        self.assertIn("| driver exit 1 |", table)

    def test_unit_counts_follow_each_metric(self):
        rows = [row("a"), row("b", unit_key="u2", structure_evaluable=False, rank_structure=math.inf)]
        summary = benchmark.summarise(rows, {**self.CONTEXT, "excluded_chunks": []}, "b",
                                      {"text": "lib", "libraries": []})
        self.assertEqual(summary["unit_macro_average"]["structure_hit@1"]["units"], 1)
        self.assertEqual(summary["unit_macro_average"]["formula_hit@1"]["units"], 2)


if __name__ == "__main__":
    unittest.main()
