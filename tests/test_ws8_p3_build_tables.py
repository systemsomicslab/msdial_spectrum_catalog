"""WS8 P3 shared-table builder (tools/ws8_p3/build_tables.py), on temporary directories and synthetic rows only."""
import contextlib
import copy
import io
from pathlib import Path
import tempfile
import unittest

from tools.ws8_p3 import build_tables as bt
from tools.ws8_p3 import rule_format as rf

SYN_CITE = [{"doi": "10.1000/synthetic.1", "locator": "Table 1"}]


def quiet(fn, *args):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = fn(*args)
    return code, out.getvalue(), err.getvalue()


def synthetic_entry(sid, formula, ions=None):
    return {"scaffold_id": sid, "name": sid, "formula": formula, "class": "synthetic", "ions": ions or {},
            "origin": "as_published", "source": SYN_CITE}


class CrossCheckTests(unittest.TestCase):
    def setUp(self):
        self.el = bt.load_elements()

    def test_shipped_rows_pass(self):
        bt.check(self.el)

    def test_synthetic_entries_against_synthetic_prints(self):
        entries = [synthetic_entry("syn_a", "C15H10O5", {"rda_13a": {"pos": ["C7H5O4"], "neg": []}}),
                   synthetic_entry("syn_b", "C16H12O6")]
        good = dict(entries=entries, y0={"syn_a": "C15H11O5"}, plus_c={"syn_b": "C17H13O6"},
                    nominal_pos={"syn_b": 301}, nominal_neg={"syn_a": 269}, printed={"C7H5O4": (1, 153.0182)})
        bt.check(self.el, **good)
        bad_cases = {
            "y0": dict(good, y0={"syn_a": "C15H13O5"}),
            "plus_c": dict(good, plus_c={"syn_b": "C17H15O6"}),
            "nominal_pos": dict(good, nominal_pos={"syn_b": 302}),
            "nominal_neg": dict(good, nominal_neg={"syn_a": 271}),
            "printed": dict(good, printed={"C7H5O4": (1, 153.0190)}),
            "unknown entry": dict(good, y0={"syn_z": "C15H11O5"}),
        }
        for label, kw in bad_cases.items():
            with self.subTest(label):
                with self.assertRaises(bt.CheckError):
                    bt.check(self.el, **kw)
        with self.assertRaises(rf.RuleFormatError):
            bt.check(self.el, entries=[synthetic_entry("syn_c", "C15H10O5", {"k": {"pos": ["H7C7O3"]}})],
                     y0={}, plus_c={}, nominal_pos={}, nominal_neg={}, printed={})


class WriteAndCompareTests(unittest.TestCase):
    def write_set(self, root):
        code, out, _ = quiet(bt.main, ["write", str(root)])
        self.assertEqual(code, 0, out)
        return out

    def test_write_produces_a_loadable_table_set_in_storage_form(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / "RuleSet"
            out = self.write_set(root)
            self.assertEqual(sorted(line.split("\t")[1] for line in out.splitlines()),
                             sorted(bt.built_files()))
            self.assertTrue(all(line.startswith("write\t") for line in out.splitlines()))
            for rel, obj in bt.built_files().items():
                self.assertEqual((root / rel).read_bytes(), rf.storage_bytes(obj), rel)
            res = rf.validate(root, draft=True)
            self.assertEqual(res.errors, [])
            self.assertEqual([u.name for u in res.tables.units][:4], ["Hex", "dHex", "Pen", "HexA"])
            self.assertIn(bt.SCAFFOLD_LIST, res.tables.scaffolds)
            self.assertEqual(len(res.tables.scaffolds[bt.SCAFFOLD_LIST].entries), len(bt.ENTRIES))
            self.assertEqual(res.rule_set_version, bt.DEFAULT_RULE_SET_VERSION)
            # the second run changes nothing
            out = self.write_set(root)
            self.assertTrue(all(line.startswith("unchanged\t") for line in out.splitlines()))
            code, out, _ = quiet(bt.main, ["compare", str(root)])
            self.assertEqual(code, 0)
            self.assertIn("OK built rows reproduced", out)
            self.assertNotIn("extra\t", out)

    def test_merged_rows_are_reported_and_never_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / "RuleSet"
            self.write_set(root)
            units = rf.load_strict((root / rf.UNITS_FILE).read_bytes())
            units["units"].append({"name": "synunit", "class": "acyl", "residue": "C2H2O", "origin": "as_published",
                                   "source": SYN_CITE})
            rf.write_storage(root / rf.UNITS_FILE, units)
            adducts = rf.load_strict((root / rf.ADDUCTS_FILE).read_bytes())
            adducts["adducts"].append({"adduct": "[M+H-2H2O]+", "aliases": [], "charge": 1, "multiplier": 1,
                                       "add": "H", "remove": "H4O2", "common": False, "origin": "as_published",
                                       "source": SYN_CITE})
            rf.write_storage(root / rf.ADDUCTS_FILE, adducts)
            extended = {p: (root / p).read_bytes() for p in bt.built_files()}
            code, out, _ = quiet(bt.main, ["compare", str(root)])
            self.assertEqual(code, 0, out)
            self.assertIn("extra\ttables/units.json\tsynunit", out)
            self.assertIn("extra\ttables/adducts.json\t[M+H-2H2O]+", out)
            self.assertIn("built=12\tequal=12\tdiffer=0\tmissing=0\textra=1\tbuilt_rows_lead=true", out)
            # write refuses and leaves every file as it was, including the one it could have created
            (root / rf.ELEMENTS_FILE).unlink()
            code, out, err = quiet(bt.main, ["write", str(root)])
            self.assertEqual(code, 1)
            self.assertIn("differs\ttables/units.json", out)
            self.assertIn("refused", err)
            self.assertFalse((root / rf.ELEMENTS_FILE).exists())
            for p, b in extended.items():
                if p != rf.ELEMENTS_FILE:
                    self.assertEqual((root / p).read_bytes(), b, p)

    def test_compare_reports_changed_missing_and_reordered_rows(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / "RuleSet"
            self.write_set(root)
            doc = rf.load_strict((root / rf.UNITS_FILE).read_bytes())
            doc["units"][4]["residue"] = "C3H2O4"  # malonyl changed
            del doc["units"][5]  # acetyl missing
            rf.write_storage(root / rf.UNITS_FILE, doc)
            lines, ok = bt.compare(root)
            self.assertFalse(ok)
            self.assertIn("differ\ttables/units.json\tmalonyl", lines)
            self.assertIn("missing\ttables/units.json\tacetyl", lines)
            # an extension row placed before the built rows changes the canonical unit order
            self.write_set(Path(d) / "Other")
            doc = rf.load_strict((Path(d) / "Other" / rf.UNITS_FILE).read_bytes())
            doc["units"].insert(0, {"name": "synunit", "class": "acyl", "residue": "C2H2O",
                                    "origin": "as_published", "source": SYN_CITE})
            rf.write_storage(Path(d) / "Other" / rf.UNITS_FILE, doc)
            lines, ok = bt.compare(Path(d) / "Other")
            self.assertFalse(ok)
            self.assertTrue(any("built_rows_lead=false" in line for line in lines))
            # a different version pin is a differing ruleset.json
            lines, ok = bt.compare(root, "0.2.0-draft")
            self.assertIn("file\truleset.json\tdiffers", lines)
            (root / rf.ELEMENTS_FILE).unlink()
            lines, ok = bt.compare(root)
            self.assertIn("missing_file\ttables/elements.json", lines)

    def test_built_rows_are_independent_of_callers(self):
        files = bt.built_files()
        before = rf.canonical_bytes(files)
        copy.deepcopy(files)[rf.UNITS_FILE]["units"].clear()
        self.assertEqual(rf.canonical_bytes(bt.built_files()), before)
        self.assertEqual(bt.built_files("0.0.1-x")[rf.RULESET_FILE]["rule_set_version"], "0.0.1-x")


if __name__ == "__main__":
    unittest.main()
