"""Tests for the MS-Emblator batch: chunking, the hash chain, resume, retry, the lock and the stop rule."""

from __future__ import annotations

import csv
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

from msdial_spectrum_catalog import msemblator_batch as batch
from msdial_spectrum_catalog.msemblator_driver import sha256

STEM = "AlignResult-2026"
SETTINGS = (
    "annotator_id\tmsp_file_path\tevidence_tier\tlibrary_kind\n"
    "msp_high_quality\tX:\\private\\lib.msp\tMSP high quality\tacquired\n"
    "msp_low_quality\tX:\\private\\lib.msp\tMSP low quality candidate\tacquired\n"
)

# A driver stand-in: it writes a record whose input_sha256 is the chunk's, and exits with the code named
# in FAKE_EXITS (a comma list consumed one call at a time through a counter file).
FAKE_DRIVER = r'''
import hashlib, json, os, sys
from pathlib import Path
args = sys.argv[1:]
inp = Path(args[args.index("--input") + 1])
out = Path(args[args.index("--output") + 1])
counter = Path(os.environ["FAKE_COUNTER"])
calls = int(counter.read_text()) if counter.exists() else 0
counter.write_text(str(calls + 1))
codes = [int(c) for c in os.environ.get("FAKE_EXITS", "0").split(",")]
code = codes[min(calls, len(codes) - 1)]
if out.exists() and any(out.iterdir()):
    sys.exit(3)
out.mkdir(parents=True, exist_ok=True)
record = {"exit_code": code, "input_sha256": hashlib.sha256(inp.read_bytes()).hexdigest(),
          "started_at": "2026-11-13T00:00:%02d+00:00" % calls, "finished_at": "2026-11-13T00:00:%02d.5+00:00" % calls,
          "duration_seconds": 0.5, "failed_steps": [], "login_ok": os.environ.get("FAKE_LOGIN") != "bad",
          "service_error_lines": 0, "stage_seconds": {"structure.sirius": 2.0}, "name_map": {"sent_to_structure": 2},
          "error": None if code == 0 else "fake failure"}
(out / "msemblator_run.json").write_text(json.dumps(record))
sys.exit(code)
'''


def write_unit(output: Path, strict: int, relaxed: int) -> None:
    """A unit output with `strict` VS21-strict and `relaxed` VS21-relaxed representative spectra."""
    output.mkdir(parents=True)
    (output / "msp_annotator_settings.tsv").write_text(SETTINGS, encoding="utf-8")
    rows, records = [], []
    # alignment IDs written out of order, to prove the chunks are in alignment-ID order
    spots = list(range(strict + relaxed, 0, -1))
    for spot in spots:
        annotator = "msp_high_quality" if spot <= strict else "msp_low_quality"
        rows.append([str(spot), "true", annotator, "430", "ReferenceSpectrum", f"name{spot}", "C1", "KEY", ""])
        records.append(f"NAME: Unknown\nPRECURSORMZ: {100 + spot}.1\nPRECURSORTYPE: [M+H]+\nIONMODE: Positive\n"
                       f"COMMENT: |PEAKID={spot}|ISOTOPE=M+0|\nNum Peaks: 3\n50.0\t1\n60.0\t2\n70.0\t3\n\n")
    with (output / f"{STEM}.mdcandidate.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["alignment_master_id", "is_representative", "annotator_id", "annotation_tag_vs1",
                         "evidence_source", "name", "formula", "inchikey", "ontology"])
        writer.writerows(rows)
    (output / f"{STEM}.mdmsp").write_text("".join(records), encoding="utf-8")


class BatchTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        write_unit(self.base / "units" / "small" / "output", strict=2, relaxed=3)
        write_unit(self.base / "units" / "large" / "output", strict=1, relaxed=7)
        self.units = [{"unit_key": "large", "output_dir": str(self.base / "units" / "large" / "output")},
                      {"unit_key": "small", "output_dir": str(self.base / "units" / "small" / "output")}]
        self.root = self.base / "annotation" / "msemblator"
        self.run_dir = self.root / "run1"
        self.run_dir.mkdir(parents=True)
        self.plan = batch.make_plan(self.units, self.run_dir, "run1", chunk_size=3, chunk_size_origin="test")
        self.msemblator = self.base / "msemblator"
        (self.msemblator / "config").mkdir(parents=True)
        (self.msemblator / "config" / "msemblator_parameter_file.yaml").write_text("x\n", encoding="utf-8")
        self.fake = self.base / "fake_driver.py"
        self.fake.write_text(FAKE_DRIVER, encoding="utf-8")
        self.counter = self.base / "calls.txt"

    def tearDown(self):
        self.temporary.cleanup()

    def _run(self, exits="0", **kwargs):
        os.environ["FAKE_COUNTER"] = str(self.counter)
        os.environ["FAKE_EXITS"] = exits
        try:
            return batch.run_batch(self.run_dir / "plan.json", self.msemblator, python=Path(sys.executable),
                                   driver=self.fake, hash_libraries=False, **kwargs)
        finally:
            for name in ("FAKE_COUNTER", "FAKE_EXITS", "FAKE_LOGIN"):
                os.environ.pop(name, None)

    def _calls(self):
        return int(self.counter.read_text()) if self.counter.exists() else 0

    def test_plan_counts_order_and_hash_chain(self):
        plan = self.plan
        self.assertEqual(plan["totals"]["records_by_selection"], {"benchmark": 3, "unresolved": 10})
        # benchmark first; within a selection the larger unit last
        self.assertEqual([c["chunk_id"] for c in plan["chunks"]],
                         ["large/benchmark/c000", "small/benchmark/c000", "small/unresolved/c000",
                          "large/unresolved/c000", "large/unresolved/c001", "large/unresolved/c002"])
        for entry in plan["exports"]:
            selection = json.loads((self.run_dir / entry["selection_json"]).read_text(encoding="utf-8"))
            self.assertEqual(entry["source_msp_sha256"], selection["msp_sha256"])
            self.assertEqual(entry["source_msp_sha256"], sha256(self.run_dir / entry["msp"]))
        for chunk in plan["chunks"]:
            self.assertEqual(chunk["chunk_sha256"], sha256(self.run_dir / chunk["chunk"]))
        self.assertTrue(any(e["manifest"].endswith("_benchmark_manifest.local-vs21.tsv") for e in plan["exports"]))
        names = []
        for chunk in plan["chunks"]:
            text = (self.run_dir / chunk["chunk"]).read_text(encoding="utf-8")
            names += [line[6:] for line in text.splitlines() if line.startswith("NAME: ")]
        self.assertEqual(len(names), len(set(names)))
        first = (self.run_dir / "input/large/large_unresolved_c000.msp").read_text(encoding="utf-8")
        self.assertEqual([line[6:] for line in first.splitlines() if line.startswith("NAME: ")],
                         ["large_2", "large_3", "large_4"])

    def test_chunking_is_deterministic(self):
        other = self.base / "annotation" / "msemblator" / "run2"
        other.mkdir()
        again = batch.make_plan(self.units, other, "run2", chunk_size=3, chunk_size_origin="test")
        self.assertEqual([c["chunk_sha256"] for c in again["chunks"]],
                         [c["chunk_sha256"] for c in self.plan["chunks"]])

    def test_a_plan_is_written_once(self):
        with self.assertRaises(SystemExit):
            batch.make_plan(self.units, self.run_dir, "run1", chunk_size=3, chunk_size_origin="test")

    def test_run_chains_records_to_chunks_and_resume_skips_done(self):
        self.assertEqual(self._run(), 0)
        self.assertEqual(self._calls(), 6)
        for chunk in self.plan["chunks"]:
            record = json.loads((self.run_dir / chunk["output_dir"] / "msemblator_run.json").read_text())
            self.assertEqual(record["input_sha256"], chunk["chunk_sha256"])
        self.assertEqual(self._run(), 0)
        self.assertEqual(self._calls(), 6, "a resumed run must skip exit-0 chunks")
        summary = batch.status(self.run_dir / "plan.json")
        self.assertEqual(summary["chunks"]["done"], 6)
        self.assertEqual(summary["overlapping_chunks"], [])
        run_json = json.loads((self.run_dir / "run.json").read_text(encoding="utf-8"))
        self.assertEqual(len(run_json["invocations"]), 2)
        self.assertIn("catalog_commit", run_json["invocations"][0])

    def test_failed_chunk_is_set_aside_and_retried_once(self):
        self.assertEqual(self._run(exits="1,0", limit=1), 0)
        output = self.run_dir / self.plan["chunks"][0]["output_dir"]
        self.assertEqual(self._calls(), 2)
        self.assertTrue(output.with_name(output.name + ".failed-1").is_dir())
        self.assertEqual(json.loads((output / "msemblator_run.json").read_text())["exit_code"], 0)

    def test_a_chunk_that_fails_twice_is_left_failed(self):
        self._run(exits="2", limit=1)
        self.assertEqual(self._calls(), 2)
        summary = batch.status(self.run_dir / "plan.json")
        self.assertEqual(summary["chunks"]["failed"], 1)
        self.assertEqual(summary["failures"][0]["reason"], "fake failure")

    def test_held_lock_refuses_a_second_runner(self):
        lock = self.root / batch.LOCK_NAME
        acquired, _ = batch.acquire_lock(lock, "other-run")
        self.assertTrue(acquired)
        try:
            self.assertEqual(self._run(), 3)
            self.assertEqual(self._calls(), 0)
        finally:
            batch.release_lock(lock)
        self.assertFalse(lock.exists())

    def test_stale_lock_of_a_dead_process_is_taken_over(self):
        lock = self.root / batch.LOCK_NAME
        lock.write_text(json.dumps({"pid": 999999999, "host": batch.socket.gethostname(), "run_id": "old"}))
        self.assertEqual(self._run(limit=1), 0)
        self.assertFalse(lock.exists())

    def test_consecutive_login_failures_stop_the_batch(self):
        os.environ["FAKE_LOGIN"] = "bad"
        self.assertEqual(self._run(exits="2", max_attempts=1), 1)
        self.assertEqual(self._calls(), 3)
        run_json = json.loads((self.run_dir / "run.json").read_text(encoding="utf-8"))
        self.assertIn("consecutive", run_json["invocations"][-1]["outcome"])

    def test_refused_driver_stops_the_batch(self):
        self.assertEqual(self._run(exits="3"), 1)
        self.assertEqual(self._calls(), 1)

    def test_refuses_the_analysis_root(self):
        analysis = self.base / "analysis"
        (analysis / batch.CAMPAIGN_MARKER).mkdir(parents=True)
        target = analysis / "mb_post" / "run"
        target.mkdir(parents=True)
        with self.assertRaises(SystemExit):
            batch.make_plan(self.units, target, "run", chunk_size=3, chunk_size_origin="test")
        self.assertFalse((target / "plan.json").exists())

    def test_backfill_rankers_marks_the_record(self):
        self._run(limit=1)
        output = self.run_dir / self.plan["chunks"][0]["output_dir"]
        (output / "msemblator.log").write_text(
            "2026-11-13T00:00:00.000000+00:00 Structure ranking model: catboostranker_final_all.pkl\n",
            encoding="utf-8")
        (self.msemblator / "models" / "structure").mkdir(parents=True)
        (self.msemblator / "models" / "structure" / "catboostranker_final_all.pkl").write_bytes(b"m")
        model = self.msemblator / "models" / "structure" / "catboostranker_final_all.pkl"
        wrong = self.base / "wrong.txt"
        wrong.write_text("0" * 64 + " *models/structure/catboostranker_final_all.pkl\n", encoding="utf-8")
        with self.assertRaises(SystemExit):
            batch.backfill_rankers(self.run_dir / "plan.json", self.msemblator, wrong)
        self.assertNotIn("rankers", json.loads((output / "msemblator_run.json").read_text(encoding="utf-8")))
        listed = self.base / "installed-models-sha256.txt"
        listed.write_text(f"{sha256(model)} *models/structure/catboostranker_final_all.pkl\n", encoding="utf-8")
        changed = batch.backfill_rankers(self.run_dir / "plan.json", self.msemblator, listed)
        self.assertEqual(changed, [self.plan["chunks"][0]["chunk_id"]])
        record = json.loads((output / "msemblator_run.json").read_text(encoding="utf-8"))
        self.assertEqual(record["rankers"]["structure"]["ranker"], "catboost")
        self.assertEqual(record["rankers_backfilled"]["from"], "msemblator.log")
        self.assertEqual(record["rankers_backfilled"]["verified_against_sha256"], sha256(listed))
        self.assertTrue(record["rankers_backfilled"]["match"])
        self.assertEqual(batch.backfill_rankers(self.run_dir / "plan.json", self.msemblator, listed), [])

    def test_throughput_and_chunk_size(self):
        self._run()
        result = batch.throughput(self.run_dir / "plan.json", {"devset": 94849})
        self.assertEqual(result["spectra"], 13)
        self.assertEqual(result["structure_sirius_seconds_per_sent_spectrum"], 1.0)
        self.assertEqual(result["chunk_size"], 300)
        self.assertEqual(batch.chunk_size_from(10.0), 90)
        self.assertEqual(batch.chunk_size_from(100.0), 20)
        self.assertIn("devset", result["extrapolation"])
        self.assertTrue((self.run_dir / "throughput.json").is_file())


if __name__ == "__main__":
    unittest.main()
