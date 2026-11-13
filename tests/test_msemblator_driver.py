"""Tests for the MS-Emblator driver, against a stub MS-Emblator checkout: no SIRIUS is needed."""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

from msdial_spectrum_catalog import msemblator_driver as driver

DRIVER = Path(driver.__file__)

# A stand-in for MS-Emblator's cli.main(): it renames NAMEs the way MS-Emblator does, drops the last record
# from the structure step (no formula), writes the four CSVs, a CANOPUS summary and MetFrag files, and
# prints the stage lines the real workflows print. Environment variables switch on the failure paths.
STUB_CLI = '''
import os, sys, time, logging
from msemblator.paths import WORK_DIR, METFRAG_DIR
from msemblator.runners.sirius_cmd import sirius_login


def _read(path):
    blocks = [b for b in open(path, encoding="utf-8").read().split("\\n\\n") if b.strip()]
    return [dict(line.split(": ", 1) for line in b.strip().splitlines() if ": " in line) for b in blocks]


def main():
    args = sys.argv[1:]
    input_path = args[args.index("--input") + 1]
    output = args[args.index("--output") + 1]
    if os.environ.get("STUB_RAISE"):
        raise ValueError("stub exception")
    records = _read(input_path)
    save = WORK_DIR / "save_folder" / "formula_fixed_msp"
    save.mkdir(parents=True, exist_ok=True)
    with open(save / "id_change.msp", "w", encoding="utf-8") as handle:
        for i, record in enumerate(records, start=1):
            comment = record.get("COMMENT", "").rstrip("|")
            handle.write(f"NAME: {i}\\nFORMULA: \\nCOMMENT: {comment}|ORIGNAME={record['NAME']}|\\nNum Peaks: 0\\n\\n")
    with open(save / "formula_fixed.msp", "w", encoding="utf-8") as handle:
        for i in range(1, len(records)):
            handle.write(f"NAME: {i}\\nFORMULA: C6H12O6\\n\\n")
    print("No SIRIUS credentials given; using the login stored by `sirius login`." if False else "")
    if os.environ.get("STUB_LOGIN"):
        try:
            sirius_login("x", "<stored-sirius-login>", "<stored-sirius-login>")
        except RuntimeError as error:
            print(f"\\nError: {error}")
            return 1
    print("Running both formula and structure elucidation...")
    print("Running formula elucidation")
    for tool in ("SIRIUS", "MS-FINDER", "msbuddy"):
        print(f"{tool} processing start")
        time.sleep(0.01)
        print(f"{tool} processing complete")
    print("Processing completed in 0.5 seconds.")
    print("Logged in as: someone@example.org")
    print('{"sid": "abc123", "subscriberId": "x"}')
    print("Running structure elucidation")
    print("SIRIUS processing start")
    logging.basicConfig(level=logging.INFO)
    logging.info("SIRIUS processing time: 2.50 seconds")
    print("SIRIUS processing complete")
    logging.info("MetFrag processing time: 1.25 seconds")
    logging.info("MSFinder processing time: 0.75 seconds")
    logging.info("Summary generation time: 0.10 seconds")
    METFRAG_DIR.mkdir(parents=True, exist_ok=True)
    (METFRAG_DIR / "example_paramater.txt").write_text("FragmentPeakMatchAbsoluteMassDeviation = 0.01\\n")
    (METFRAG_DIR / "1.csv").write_text("Score\\n1\\n")
    sirius_out = WORK_DIR / "structure" / "sirius_output"
    sirius_out.mkdir(parents=True, exist_ok=True)
    (sirius_out / "structure_identifications_top-100.tsv").write_text("a\\tb\\n1\\t2\\n")
    (sirius_out / "sirius_project.sirius").write_bytes(b"project")
    if not os.environ.get("STUB_NO_CANOPUS"):
        (sirius_out / "canopus_formula_summary.tsv").write_text("mappingFeatureId\\n1\\n")
    for name in ("formula_summary.csv", "formula_score.csv", "structure_score.csv", "structure_summary.csv"):
        with open(os.path.join(output, name), "w", encoding="utf-8") as handle:
            handle.write("filename,rank\\n")
            for record in records:
                handle.write(f"{record['NAME']},1\\n")
    if os.environ.get("STUB_FAIL_METFRAG"):
        print("\\nError: structure elucidation completed with failed steps: MetFrag")
        return 1
    return 0
'''

STUB_PATHS = '''
from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
WORK_DIR = PROJECT_ROOT / "work"
METFRAG_DIR = PROJECT_ROOT / "external" / "metfrag"
FORMULA_MODEL_FILES = ("pipeline_sirius_score.pkl",)
STRUCTURE_MODEL_FILES = ("pipeline_csi_fingerid.pkl",)
RANKER_MODEL_FILES = ("xgbranker_final_all.pkl",)
'''

STUB_LOGIN_MODULE = "def sirius_login(sirius_directory, username, password):\n    raise AssertionError('not patched')\n"


def make_stub_root(base: Path) -> Path:
    root = base / "msemblator"
    package = root / "src" / "msemblator"
    for folder in ("cli", "runners", "workflows"):
        (package / folder).mkdir(parents=True)
        (package / folder / "__init__.py").write_text("", encoding="utf-8")
    (package / "__init__.py").write_text("", encoding="utf-8")
    (package / "paths.py").write_text(STUB_PATHS, encoding="utf-8")
    (package / "cli" / "msemblator.py").write_text(STUB_CLI, encoding="utf-8")
    for module in ("runners/sirius_cmd.py", "runners/sirius_struc_cmd.py", "workflows/struc_main.py"):
        (package / module).write_text(STUB_LOGIN_MODULE, encoding="utf-8")
    (root / "config").mkdir()
    (root / "config" / "msemblator_parameter_file.yaml").write_text("MS1: qtof\n", encoding="utf-8")
    (root / "external" / "sirius" / "app").mkdir(parents=True)
    (root / "external" / "sirius" / "app" / "sirius_rest_service-6.3.2-boot.jar").write_bytes(b"")
    (root / "external" / "metfrag").mkdir(parents=True)
    (root / "external" / "metfrag" / "MetFragCommandLine-2.5.0.jar").write_bytes(b"jar")
    (root / "models" / "formula").mkdir(parents=True)
    (root / "models" / "formula" / "pipeline_sirius_score.pkl").write_bytes(b"model")
    (root / "library").mkdir()
    (root / "library" / "metfrag_StructureDB.txt").write_text("InChIKey1\n", encoding="utf-8")
    return root


def write_input(base: Path, names: list[str]) -> tuple[Path, Path]:
    msp = base / "unit1_unresolved_c000.msp"
    with msp.open("w", encoding="utf-8", newline="\n") as handle:
        for index, name in enumerate(names, start=1):
            handle.write(f"NAME: {name}\nPRECURSORMZ: 100.0\nCOMMENT: |UNIT=unit1|ALIGNMENTID={index}|\n"
                         f"Num Peaks: 3\n50.0\t1\n60.0\t2\n70.0\t3\n\n")
    manifest = base / "unit1_unresolved_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["name", "unit_id", "alignment_id", "tier", "reference_name", "reference_inchikey"])
        for index, name in enumerate(names, start=1):
            writer.writerow([name, "unit1", index, "unassigned", "", ""])
    return msp, manifest


class DriverTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.root = make_stub_root(self.base)
        self.names = ["unit1_1", "unit1_2", "unit1_3"]
        self.msp, self.manifest = write_input(self.base, self.names)
        self.output = self.base / "out" / "c000"

    def tearDown(self):
        self.temporary.cleanup()

    def _run(self, **env):
        environment = dict(os.environ, PYTHONIOENCODING="utf-8", **env)
        command = [sys.executable, "-B", str(DRIVER), "--msemblator-root", str(self.root), "--input", str(self.msp),
                   "--output", str(self.output), "--mode", "2", "--manifest", str(self.manifest),
                   "--model-assets", "zenodo_record=123,md5=abc"]
        return subprocess.run(command, capture_output=True, text=True, encoding="utf-8", env=environment,
                              timeout=120)

    def _record(self):
        return json.loads((self.output / "msemblator_run.json").read_text(encoding="utf-8"))

    def test_clean_run_keeps_work_metfrag_name_map_and_record(self):
        result = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)
        record = self._record()
        self.assertEqual(record["exit_code"], 0)
        self.assertEqual(record["sirius_version"], "6.3.2")
        self.assertEqual(record["failed_steps"], [])
        self.assertTrue(record["name_map_ok"])
        self.assertEqual(record["name_map"]["records"], 3)
        self.assertEqual(record["name_map"]["skipped_no_formula"], 1)
        self.assertEqual(record["canopus_summary_files"], ["canopus_formula_summary.tsv"])
        self.assertIn("work/structure/sirius_output/sirius_project.sirius", record["harvest"])
        self.assertIn("metfrag/1.csv", record["harvest"])
        self.assertEqual(set(record["outputs"]), set(driver.CSV_OUTPUTS))
        self.assertEqual(record["input_sha256"], driver.sha256(self.msp))
        self.assertEqual(record["stage_seconds"]["structure.sirius"], 2.5)
        self.assertEqual(record["stage_seconds"]["structure.metfrag"], 1.25)
        self.assertIn("formula.sirius", record["stage_seconds"])
        self.assertIn("formula.total", record["stage_seconds"])
        self.assertEqual(record["environment"]["model_assets"], {"zenodo_record": "123", "md5": "abc"})
        self.assertIn("third_party_data", record["environment"])
        self.assertEqual(record["environment"]["model_sha256"]["formula/pipeline_sirius_score.pkl"],
                         driver.sha256(self.root / "models" / "formula" / "pipeline_sirius_score.pkl"))
        # work/ moved out of the checkout; MetFrag copied, the checkout's own files left behind
        self.assertFalse((self.root / "work").exists())
        self.assertTrue((self.output / "work" / "save_folder" / "formula_fixed_msp" / "id_change.msp").is_file())
        self.assertTrue((self.output / "metfrag" / "example_paramater.txt").is_file())
        self.assertFalse((self.output / "metfrag" / "MetFragCommandLine-2.5.0.jar").exists())
        with (self.output / "name_map.tsv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([r["name"] for r in rows], self.names)
        self.assertEqual([r["msemblator_id"] for r in rows], ["1", "2", "3"])
        self.assertEqual([r["alignment_id"] for r in rows], ["1", "2", "3"])
        self.assertEqual(rows[-1]["structure_step"], "skipped_no_formula")
        self.assertNotIn("reference_name", rows[0])

    def test_no_account_detail_or_placeholder_is_written(self):
        self._run()
        for name in ("msemblator_run.json", "msemblator.log", "name_map.tsv"):
            text = (self.output / name).read_text(encoding="utf-8")
            self.assertNotIn("someone@example.org", text, name)
            self.assertNotIn("abc123", text, name)
            self.assertNotIn(driver.PLACEHOLDER, text, name)
        self.assertIn("<account>", (self.output / "msemblator.log").read_text(encoding="utf-8"))

    def test_refuses_a_busy_work_folder(self):
        (self.root / "work" / "formula").mkdir(parents=True)
        (self.root / "work" / "formula" / "left.txt").write_text("x", encoding="utf-8")
        result = self._run()
        self.assertEqual(result.returncode, 3)
        self.assertFalse(self.output.exists())

    def test_refuses_a_non_empty_output(self):
        self.output.mkdir(parents=True)
        (self.output / "old.txt").write_text("x", encoding="utf-8")
        result = self._run()
        self.assertEqual(result.returncode, 3)
        self.assertFalse((self.output / "msemblator_run.json").exists())

    def test_failed_steps_are_a_top_level_list(self):
        result = self._run(STUB_FAIL_METFRAG="1")
        self.assertEqual(result.returncode, 1)
        record = self._record()
        self.assertEqual(record["failed_steps"], ["structure:MetFrag"])
        self.assertEqual(record["exit_code"], 1)

    def test_missing_canopus_is_exit_4(self):
        result = self._run(STUB_NO_CANOPUS="1")
        self.assertEqual(result.returncode, 4)
        self.assertEqual(self._record()["canopus_summary_files"], [])

    def test_exception_is_exit_2_and_work_still_moves(self):
        result = self._run(STUB_RAISE="1")
        self.assertEqual(result.returncode, 2)
        self.assertIn("stub exception", self._record()["error"])

    def test_missing_login_is_exit_2(self):
        result = self._run(STUB_LOGIN="1")
        self.assertEqual(result.returncode, 2)
        record = self._record()
        self.assertFalse(record["login_ok"])
        self.assertTrue(record["error"].startswith("login:"))

    def test_name_map_must_cover_the_manifest(self):
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            handle.write("name\tunit_id\talignment_id\nunit1_1\tunit1\t1\n")
        result = self._run()
        self.assertEqual(result.returncode, 2)
        record = self._record()
        self.assertFalse(record["name_map_ok"])
        self.assertTrue(record["error"].startswith("name_map:"))


class ParseLogTest(unittest.TestCase):
    def test_sections_and_failed_steps(self):
        lines = [
            "2026-11-13T00:00:00.000000+00:00 Running formula elucidation",
            "2026-11-13T00:00:01.000000+00:00 SIRIUS processing start",
            "2026-11-13T00:00:05.000000+00:00 SIRIUS processing complete",
            "2026-11-13T00:00:06.000000+00:00 Processing completed in 6.0 seconds.",
            "2026-11-13T00:00:07.000000+00:00 Running structure elucidation",
            "2026-11-13T00:00:08.000000+00:00 2026-11-13 09:00:08 - INFO - SIRIUS processing time: 30.00 seconds",
            "2026-11-13T00:00:09.000000+00:00 ERROR: HTTP 403 Forbidden from the web service",
            "2026-11-13T00:00:09.000000+00:00 401.2001\t1234",
            "2026-11-13T00:00:10.000000+00:00 Error: formula elucidation completed with failed steps: SIRIUS, msbuddy",
        ]
        parsed = driver.parse_log(lines)
        self.assertEqual(parsed["stage_seconds"]["formula.sirius"], 4.0)
        self.assertEqual(parsed["stage_seconds"]["formula.summary"], 1.0)
        self.assertEqual(parsed["stage_seconds"]["structure.sirius"], 30.0)
        self.assertEqual(parsed["failed_steps"], ["formula:SIRIUS", "formula:msbuddy"])
        self.assertEqual(parsed["service_error_lines"], 1)

    def test_third_party_data_records_msbuddy_files_with_their_source(self):
        import importlib

        with tempfile.TemporaryDirectory() as folder:
            package = Path(folder) / "msbuddy"
            (package / "data").mkdir(parents=True)
            (package / "__init__.py").write_text("", encoding="utf-8")
            (package / "load.py").write_text(
                "URL = 'https://github.com/example/msbuddy/releases/download/v1/formula_db_v1.joblib'\n",
                encoding="utf-8")
            (package / "data" / "formula_db_v1.joblib").write_bytes(b"db")
            sys.path.insert(0, folder)
            importlib.invalidate_caches()
            sys.modules.pop("msbuddy", None)
            try:
                data = driver.third_party_data()
            finally:
                sys.path.remove(folder)
                sys.modules.pop("msbuddy", None)
        files = data["msbuddy"]["files"]
        self.assertEqual(files[0]["file"], "formula_db_v1.joblib")
        self.assertEqual(files[0]["bytes"], 2)
        self.assertEqual(files[0]["source_url"],
                         "https://github.com/example/msbuddy/releases/download/v1/formula_db_v1.joblib")
        self.assertEqual(len(files[0]["sha256"]), 64)
        self.assertIn("fetched_at_file_mtime", files[0])

    def test_parse_assets(self):
        self.assertEqual(driver.parse_assets("zenodo_record=1, md5=ab"), {"zenodo_record": "1", "md5": "ab"})
        self.assertEqual(driver.parse_assets(None), {})

    def test_scrub(self):
        text = driver.scrub('Logged in as: a.b@c-d.ac.jp "sid": "xyz" refresh_token=abc')
        self.assertNotIn("a.b@c-d.ac.jp", text)
        self.assertNotIn("xyz", text)
        self.assertNotIn("abc", text)


if __name__ == "__main__":
    unittest.main()
