"""WS8 P3 C# parity vectors (tools/ws8_p3/cjson_vectors.py): stable output, a loader vector for every error code."""
import base64
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile
import unittest

from tools.ws8_p3 import cjson_vectors as cv
from tools.ws8_p3 import rule_format as rf

SCRIPT = Path(cv.__file__).resolve()
READER_CODES = ("bom", "invalid_utf8", "syntax", "top_level_not_object", "duplicate_key", "nonfinite",
                "nesting_depth", "number_range", "negative_zero", "number_form", "float_precision", "control_char",
                "non_ascii")


def unb64(s):
    return base64.b64decode(s.encode("ascii"), validate=True)


class VectorFileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = cv.render()
        cls.vec = json.loads(cls.data)

    def test_output_is_stable_across_runs_and_hash_seeds(self):
        self.assertEqual(cv.render(), self.data)
        with tempfile.TemporaryDirectory() as d:
            outs = []
            for seed in ("0", "4242"):
                out = Path(d) / seed / "cjson-vectors.json"
                env = dict(os.environ, PYTHONHASHSEED=seed)
                proc = subprocess.run([sys.executable, "-B", str(SCRIPT), "write", str(out)], env=env,
                                      capture_output=True, text=True, timeout=300)
                self.assertEqual(proc.returncode, 0, proc.stderr)
                self.assertIn(hashlib.sha256(self.data).hexdigest(), proc.stdout)
                outs.append(out.read_bytes())
            self.assertEqual(outs[0], self.data)
            self.assertEqual(outs[1], self.data)

    def test_file_is_strict_json_in_storage_form(self):
        obj = rf.load_strict(self.data)
        self.assertEqual(rf.storage_bytes(obj), self.data)
        self.assertEqual((obj["format"], obj["format_version"]), (cv.FORMAT, 1))
        self.assertEqual(obj["engine_semantics"], rf.ENGINE_SEMANTICS)
        self.assertEqual(obj["error_codes"], dict(sorted(rf.ERROR_CODES.items())))
        self.assertEqual(set(obj["sections"]) | {"format", "format_version", "generator", "oracle",
                                                 "engine_semantics", "sections", "error_codes"}, set(obj))

    def test_every_error_code_has_a_loader_vector(self):
        loader = self.vec["loader_errors"]
        codes = [c["code"] for c in loader["cases"] if c["code"] is not None]
        self.assertEqual(sorted(codes), sorted(rf.ERROR_CODES))
        self.assertEqual(len(codes), len(set(codes)))
        base = {p: unb64(b) for p, b in loader["base_files"].items()}
        self.assertEqual(rf.load_files(base, draft=False).errors, [])
        for case in loader["cases"]:
            with self.subTest(case=case["name"]):
                files = dict(base)
                for p, b in case["patch"].items():
                    if b is None:
                        files.pop(p)
                    else:
                        files[p] = unb64(b)
                errors = [[e.code, e.file, e.pointer] for e in rf.load_files(files, draft=False).errors]
                self.assertEqual(errors, case["errors"])
                self.assertEqual(errors, sorted(errors, key=lambda e: (e[1], e[2], e[0])))
                if case["code"] is not None:
                    self.assertIn(case["code"], {e[0] for e in errors})

    def test_strict_vectors(self):
        cases = self.vec["strict_json"]
        self.assertEqual(len({c["name"] for c in cases}), len(cases))
        self.assertTrue(set(READER_CODES) <= {c.get("code") for c in cases})
        for c in cases:
            with self.subTest(c["name"]):
                raw = unb64(c["input_b64"])
                if c["expect"] == "ok":
                    canon = rf.canonical_bytes(rf.load_strict(raw))
                    self.assertEqual(canon.decode("ascii"), c["canonical"])
                    self.assertEqual(hashlib.sha256(canon).hexdigest(), c["sha256"])
                else:
                    with self.assertRaises(rf.RuleFormatError) as cm:
                        rf.load_strict(raw)
                    self.assertEqual(cm.exception.code, c["code"])
        by = {c["name"]: c for c in cases}
        self.assertEqual(by["key_order_code_points"]["canonical"], '{"":7,"B":3,"_":4,"a":2,"a ":6,"a0":5,"b":1}')
        self.assertEqual(by["float_trailing_zeros"]["canonical"], '{"a":[0.05,5.0,10.0,1.5,100.0]}')
        self.assertEqual((by["depth_128"]["expect"], by["depth_129"]["code"]), ("ok", "nesting_depth"))
        self.assertEqual(by["float_negative_zero_places"]["code"], "negative_zero")
        self.assertEqual(by["float_negative_tiny"]["code"], "number_range")

    def test_float_lexemes(self):
        for c in self.vec["float_lexeme"]:
            with self.subTest(c["lexeme"]):
                x = float(c["lexeme"])
                self.assertEqual(struct.pack(">d", x).hex(), c["ieee_hex"])
                self.assertEqual(c["mantissa"] / 10 ** c["scale"], x)
                self.assertLess(abs(c["mantissa"]), 10 ** 15)
                self.assertLessEqual(c["scale"], 18)
                self.assertEqual(repr(x), c["canonical"])

    def test_number_text_and_bits(self):
        self.assertEqual(cv.fixed(-1e-9, 6), "0.000000")
        self.assertEqual(cv.fixed(-0.0, 6), "0.000000")
        self.assertEqual(cv.fixed(2.5, 0), "2")
        self.assertEqual(cv.fixed(0.0625, 3), "0.062")
        self.assertIsNone(cv.result_double(float("nan"), 6))
        self.assertIsNone(cv.result_double(float("-inf"), 6))
        self.assertEqual(cv.bits(cv.from_hex("fff8000000000000")), "7ff8000000000000")
        self.assertEqual(cv.bits(-0.0), "8000000000000000")
        for c in self.vec["fixed_decimal"]:
            self.assertEqual(cv.fixed(cv.from_hex(c["ieee_hex"]), c["d"]), c["text"])
            self.assertFalse(re.fullmatch(r"-[0.]+", c["text"]), c)
        for c in self.vec["ieee_bits"]:
            self.assertEqual(cv.bits(cv.from_hex(c["ieee_hex_in"])), c["bits"])
        nan_in = {c["ieee_hex_in"] for c in self.vec["ieee_bits"] if c["bits"] == "7ff8000000000000"}
        self.assertTrue({"7ff8000000000000", "fff8000000000000", "7ff0000000000001"} <= nan_in)

    def test_strings_and_instrument_class(self):
        for c in self.vec["string_escape"]:
            s = cv.from_utf16_hex(c["text_utf16_hex"])
            self.assertEqual(json.dumps(s, ensure_ascii=True), c["json"])
        by = {c["text_utf16_hex"]: c["json"] for c in self.vec["string_escape"]}
        self.assertEqual(by["d83dde00"], '"\\ud83d\\ude00"')
        self.assertEqual(by["007f"], '"\\u007f"')
        classes = {None if c["text_utf16_hex"] is None else cv.from_utf16_hex(c["text_utf16_hex"]): c["class"]
                   for c in self.vec["instrument_class"]}
        self.assertEqual(classes[None], "unknown")
        self.assertEqual(classes[" \t\r\n"], "unknown")
        self.assertEqual(classes["\x0b"], "other")
        self.assertEqual(classes["qtof"], "QTOF-type")
        self.assertEqual(classes["q-exactive"], "Orbitrap/FT")
        self.assertEqual(classes["toﬀ"], "other")  # ASCII-only upper case: no TOF

    def test_hashes(self):
        for c in self.vec["file_hashes"]:
            files = {p: unb64(b) for p, b in c["files"].items()}
            self.assertEqual({p: rf.file_sha256(b) for p, b in files.items()}, c["file_sha256"])
            self.assertEqual(rf.ruleset_sha256(files), c["ruleset_sha256"])
        for c in self.vec["rule_semantic"]:
            self.assertEqual(rf.rule_semantic_sha256(c["rule"], c["resolved"], c["engine_semantics"]), c["sha256"])
        by = {c["name"]: c["sha256"] for c in self.vec["rule_semantic"]}
        self.assertEqual(by["review_fields_ignored"], by["flavonoid.o_glycoside.neg"])
        self.assertNotEqual(by["engine_semantics_bumped"], by["flavonoid.o_glycoside.neg"])
        for c in self.vec["ruleset_load"]:
            res = rf.load_files({p: unb64(b) for p, b in c["files"].items()}, draft=False)
            self.assertEqual(res.errors, [])
            self.assertEqual((res.ruleset_sha256, res.file_hashes), (c["ruleset_sha256"], c["file_sha256"]))
            for rid, r in c["rules"].items():
                self.assertEqual(res.semantic_sha256(rid), r["rule_semantic_sha256"])
                self.assertEqual(rf.rule_resolved(res.rules[rid], res.tables), r["resolved"])

    def test_masses_against_hand_sums(self):
        m = self.vec["masses"]
        self.assertEqual(m["electron_pda"], 549_000_000)
        formulas = {c["formula"]: c["pda"] for c in m["formula_pda"]}
        self.assertEqual(formulas["C6H10O5"], 6 * 12 * 10 ** 12 + 10 * 1_007_825_032_000 + 5 * 15_994_914_620_000)
        self.assertEqual(formulas["H2O"], 2 * 1_007_825_032_000 + 15_994_914_620_000)
        ions = {(c["formula"], c["charge"]): c for c in m["ions"]}
        self.assertEqual(ions[("H", 1)]["pda"], 1_007_825_032_000 - 549_000_000)
        self.assertEqual(ions[("C7H3O4", -1)]["mz_text"], "151.003683")
        prec = {(c["neutral"], c["adduct"]): c for c in m["adduct_precursors"]}
        self.assertEqual(prec[("C15H10O5", "[M+H]+")]["mz_text"], "271.060099")
        self.assertEqual(prec[("C21H21O11", "[M]+")]["mz_text"], "449.107837")
        two = prec[("C15H18O9", "[2M-H]-")]
        self.assertEqual(two["precursor_pda"], 2 * two["neutral_pda"] - 1_007_825_032_000 + 549_000_000)
        for c in m["ions"] + m["adduct_precursors"]:
            pda = c["pda"] if "pda" in c else c["precursor_pda"]
            self.assertEqual(struct.pack(">d", pda / 1e12).hex(), c["mz_ieee_hex"])
            self.assertLess(pda, 2 ** 53)
        self.assertIn({"lexeme": "1.0000000000001", "code": "mass_precision"}, m["lexeme_pda"])
        self.assertEqual({c["code"] for c in m["formula_errors"]}, {"bad_formula", "formula_not_hill",
                                                                     "unknown_element"})

    def test_check_command_and_source_md(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "v" / "cjson-vectors.json"
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cv.main(["write", str(out), "--source-md"]), 0)
                self.assertEqual(cv.main(["check", str(out)]), 0)
            self.assertEqual(out.read_bytes(), self.data)
            md = (out.parent / "SOURCE.md").read_text(encoding="ascii")
            self.assertIn(hashlib.sha256(self.data).hexdigest(), md)
            self.assertIsNone(re.search(r"\b[A-Za-z]:[\\/]", md))
            self.assertNotIn("Users", md)
            out.write_bytes(self.data + b" ")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cv.main(["check", str(out)]), 1)


if __name__ == "__main__":
    unittest.main()
