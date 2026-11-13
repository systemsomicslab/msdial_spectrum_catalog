"""Cross-language fixture drift and deterministic regeneration checks."""
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from msdial_spectrum_catalog import similarity_mn_tools as mt

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "similarity-fixtures"
TESTS = ROOT / "tests"


class FixtureChecks(unittest.TestCase):
    def test_fixture_drift(self):
        result = mt.check_fixtures(FIXTURES)
        self.assertEqual(result, {"ok": True, "files": 4, "cases": 56, "kernel_expectations": 312})

    def test_byte_identical_regeneration(self):
        with tempfile.TemporaryDirectory(prefix="mn-fixtures-", dir=TESTS) as tmp:
            target = Path(tmp)
            mt.write_fixtures(target)
            self.assertEqual(sorted(p.name for p in target.iterdir()), sorted(p.name for p in FIXTURES.iterdir()))
            for original in FIXTURES.iterdir():
                self.assertEqual(original.read_bytes(), (target / original.name).read_bytes(), original.name)

    def test_detect_corrupt_hash_and_rehashed_wrong_score(self):
        with tempfile.TemporaryDirectory(prefix="mn-drift-", dir=TESTS) as tmp:
            target = Path(tmp)
            for original in FIXTURES.iterdir():
                shutil.copyfile(original, target / original.name)
            path = target / "kernels.json"
            value = json.loads(path.read_bytes())
            value["cases"][0]["expected"][0]["score"] = .25
            path.write_bytes(mt.json_bytes(value))
            with self.assertRaisesRegex(ValueError, "manifest drift"):
                mt.check_fixtures(target)
            manifest_path = target / "manifest.json"
            manifest = json.loads(manifest_path.read_bytes())
            next(e for e in manifest["files"] if e["path"] == path.name)["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest_path.write_bytes(mt.json_bytes(manifest))
            with self.assertRaisesRegex(ValueError, "score drift"):
                mt.check_fixtures(target)

    def test_crlf_is_reported_explicitly_for_every_fixture_file(self):
        with tempfile.TemporaryDirectory(prefix="mn-crlf-", dir=TESTS) as tmp:
            target = Path(tmp)
            for original in FIXTURES.iterdir():
                with self.subTest(file=original.name):
                    mt.write_fixtures(target)
                    path = target / original.name
                    path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
                    with self.assertRaisesRegex(ValueError, "CRLF line endings in " + original.name):
                        mt.check_fixtures(target)

    def test_degenerate_fixture_contracts(self):
        value = json.loads((FIXTURES / "kernels.json").read_bytes())
        case = next(c for c in value["cases"] if c["id"] == "moddot-zero-nan")
        guarded = [e for e in case["expected"]
                   if e["kernel_id"].startswith("gnps.") and "prep=mn&" in e["kernel_id"]]
        self.assertEqual(len(guarded), 2)
        for expected in guarded:
            self.assertEqual(expected["implementations"], ["python"])
            self.assertEqual(expected["score"], "NaN")
            if "assign=greedy" in expected["kernel_id"]:
                self.assertEqual(expected["matched_peak_count"], 0)
        entropy = [e for e in case["expected"] if e["kernel_id"].startswith("entropy.")]
        self.assertEqual(len(entropy), 2)
        for expected in entropy:
            self.assertIsNone(expected["score"])

    def test_json_encoding_and_precision_contract(self):
        def precision(spectrum):
            for mz in [spectrum["precursor_mz"]] + [p[0] for p in spectrum["peaks"]]:
                self.assertGreaterEqual(Decimal(str(mz)).normalize().as_tuple().exponent, -5)
            for _, intensity in spectrum["peaks"]:
                self.assertLessEqual(len(Decimal(str(intensity)).normalize().as_tuple().digits), 6)
        for path in FIXTURES.glob("*.json"):
            raw = path.read_bytes()
            self.assertNotIn(b"\r", raw)
            self.assertTrue(raw.endswith(b"\n"))
            value = json.loads(raw)
            self.assertEqual(mt.json_bytes(value), raw)
            for case in value.get("cases", []):
                for field in ("a", "b", "input"):
                    if isinstance(case.get(field), dict) and "peaks" in case[field]:
                        precision(case[field])

    def test_vs20_provenance_and_record_cap(self):
        value = json.loads((FIXTURES / "kernels-vs20.json").read_bytes())
        self.assertEqual(value["source"], {"filename": mt.VS20_NAME, "md5": mt.VS20_MD5,
                         "doi": "10.5281/zenodo.21904103", "license": "CC BY 4.0"})
        indices = [i for case in value["cases"] for i in case["source_record_indices"]]
        self.assertEqual(len(indices), 16)
        self.assertEqual(indices, sorted(set(indices)))
        self.assertEqual(len(value["cases"]), 8)

    def test_smoke_parser_and_cli_on_synthetic_msp_only(self):
        # Never use a unit-output record in tests. Even this temporary synthetic
        # MSP is created only below tests and removed by TemporaryDirectory.
        with tempfile.TemporaryDirectory(prefix="mn-smoke-", dir=TESTS) as tmp:
            path = Path(tmp) / "synthetic.msp"
            path.write_text("NAME: synthetic-discard-me\nPRECURSORMZ: 250\n"
                "COMMENT: IGNORE=discard|PEAKID=7|OTHER=discard\nINCHIKEY: discard\n"
                "SMILES: discard\nONTOLOGY: discard\nNum Peaks: 2\n100 3\n150 4\n\n"
                "NAME: synthetic-empty\nPRECURSORMZ: 300\nNum Peaks: 0\n\n"
                "NAME: synthetic-shift\nPRECURSORMZ: 264\nCOMMENT: PEAKID=8\n"
                "Num Peaks: 2\n100 3\n164 4\n", encoding="utf-8")
            parsed = list(mt.read_smoke_msp(path))
            self.assertEqual([p[1] for p in parsed], ["7", None, "8"])
            self.assertNotIn("discard", repr(parsed))
            args = [sys.executable, "-B", "-m", "msdial_spectrum_catalog.similarity_mn_tools",
                    "smoke", "--msp", str(path)]
            first = subprocess.run(args, capture_output=True, text=True, check=True, timeout=20)
            second = subprocess.run(args, capture_output=True, text=True, check=True, timeout=20)
            a, b = json.loads(first.stdout), json.loads(second.stdout)
            self.assertGreaterEqual(a.pop("runtime_s"), 0)
            b.pop("runtime_s")
            self.assertEqual(a, b)
            self.assertEqual((a["records"], a["records_with_peaks"], a["pair_count"]), (3, 2, 1))
            self.assertEqual(a["fraction_k1_ge_0_95_and_k3_lt_0_5"], 0)
            self.assertAlmostEqual(a["kernels"][0]["min"], 1, delta=1e-12)
            self.assertEqual(a["kernels"][0]["fraction_matched_1"], 0)
            self.assertIsNone(a["kernels"][-1]["fraction_matched_1"])
            empty = mt.smoke(path, max_dprec=1)
            self.assertEqual(empty["pair_count"], 0)
            self.assertTrue(all(k["fraction_matched_1"] is None for k in empty["kernels"]))
            self.assertNotIn("PEAKID", first.stdout)
            self.assertNotIn("discard", first.stdout)


if __name__ == "__main__":
    unittest.main()
