"""WS8 P3 MS-FINDER off/on measurement (tools/ws8_p3/p3_msfinder.py), on synthetic records, synthetic library
snippets and hand-built candidate lists only."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from tools.ws8_p2 import p2_diagnostics as p2
from tools.ws8_p3 import p3_msfinder as ms

REVIEWED = {"saponin.a.neg": dict(family="saponin", motif_id="m1", version=1),
            "glycoside.b.pos": dict(family="glycoside", motif_id="m2", version=1)}


def cand(block, score, total=None, onto=""):
    return dict(block=block, ontology=onto, score=score, total=score if total is None else total, term=0.0)


class SelectionTest(unittest.TestCase):
    def reps(self):
        return [dict(record_id=i, key=f"K{i}", mode="neg" if i % 2 else "pos", cluster_id=f"c{i // 2}",
                     instrument_type="" if i % 3 else "LC-ESI-QTOF") for i in range(40)]

    def test_selects_reviewed_outcomes_in_db_and_counts_the_rest(self):
        outcomes = {1: [("saponin.a.neg", True, False)], 2: [("glycoside.b.pos", False, True)],
                    3: [("other.rule.neg", True, False)], 4: [("saponin.a.neg", False, False)],
                    5: [("saponin.a.neg", True, False)]}
        db = {f"K{i}" for i in range(40) if i != 5}
        rows, counts = ms.select_records(self.reps(), outcomes, {1, 3, 5}, db, REVIEWED, seed=7, control_n=10)
        sel = [r for r in rows if r["set"] == "selected"]
        self.assertEqual([(r["record_id"], r["family"], r["outcome"]) for r in sel],
                         [(1, "saponin", "fired"), (2, "glycoside", "veto_blocked")])
        self.assertEqual(counts["reviewed_outcome_not_in_db"], 1)
        self.assertEqual(counts["in_db"], 39)
        self.assertEqual(counts["in_db_no_reviewed_outcome"], 37)
        self.assertEqual(counts["in_db_no_rule_fired"], 36)  # record 3 fired a non-reviewed rule
        self.assertEqual(counts["control"], 10)
        self.assertTrue(all("key" not in r for r in rows))
        self.assertTrue(all(r["record_id"] not in (1, 2, 5) for r in rows if r["set"] == "control"))

    def test_control_sample_is_seeded(self):
        a, _ = ms.select_records(self.reps(), {}, set(), {f"K{i}" for i in range(40)}, REVIEWED, seed=3, control_n=5)
        b, _ = ms.select_records(self.reps(), {}, set(), {f"K{i}" for i in range(40)}, REVIEWED, seed=3, control_n=5)
        c, _ = ms.select_records(self.reps(), {}, set(), {f"K{i}" for i in range(40)}, REVIEWED, seed=4, control_n=5)
        self.assertEqual(a, b)
        self.assertNotEqual([r["record_id"] for r in a], [r["record_id"] for r in c])

    def test_smoke_limit(self):
        outcomes = {i: [("saponin.a.neg", True, False)] for i in range(1, 40, 2)}
        rows, counts = ms.select_records(self.reps(), outcomes, set(), {f"K{i}" for i in range(40)}, REVIEWED,
                                         seed=1, limit=8)
        self.assertEqual(counts["selected"], 8)
        self.assertEqual(counts["control"], 2)

    def test_console_instrument_class(self):
        self.assertEqual(ms.console_class(""), "unknown")
        self.assertEqual(ms.console_class("LC-ESI-QTOF"), "QTOF-type")
        self.assertEqual(ms.console_class("Q Exactive"), "Orbitrap/FT")


class LibraryTest(unittest.TestCase):
    def test_short_keys_icd_and_obo_class_level(self):
        with tempfile.TemporaryDirectory() as d:
            esd, icd, obo = Path(d) / "s.esd", Path(d) / "c.icd", Path(d) / "c.obo"
            esd.write_text("Title\tInChIkey\tShort InChIKey\tSMILES\nx\tAAAA-BBBB-N\tAAAA\tC\ny\tCCCC-DDDD-N\tCCCC\tCC\n",
                           encoding="utf-8")
            icd.write_text("InChIKey\tDirect parent name\tDirect parent ID\nAAAA-BBBB-N\tLeaf\tCHEMONTID:0000004\n"
                           "AAAA-XXXX-N\tOther\tCHEMONTID:0000003\n", encoding="utf-8")
            obo.write_text("".join(f"[Term]\nid: CHEMONTID:000000{i}\nname: T{i}\nis_a: {p} ! x\n\n" for i, p in
                                   [(0, "CHEMONTID:9999999"), (1, "CHEMONTID:0000000"), (2, "CHEMONTID:0000001"),
                                    (3, "CHEMONTID:0000002"), (4, "CHEMONTID:0000003")]), encoding="utf-8")
            self.assertEqual(ms.read_short_keys(esd), {"AAAA", "CCCC"})
            by_block, name_to_id = ms.read_icd(icd)
            self.assertEqual(by_block["AAAA"], ("Leaf", "CHEMONTID:0000004"))  # first row per first block
            parents, names = ms.read_obo(obo)
            self.assertEqual(names["T2"], "CHEMONTID:0000002")
            self.assertEqual(ms.class_level("CHEMONTID:0000004", parents), "CHEMONTID:0000002")
            self.assertEqual(ms.class_level("CHEMONTID:0000002", parents), "CHEMONTID:0000002")
            self.assertIsNone(ms.class_level("CHEMONTID:0000001", parents))
            self.assertIsNone(ms.class_level(None, parents))


class MetricTest(unittest.TestCase):
    def test_rank_is_tie_pessimistic_and_order_free(self):
        off = [cand("T", 5), cand("W1", 4), cand("W2", 3), cand("T", 2)]
        self.assertEqual(ms.rank_of(off, "T"), 1)
        self.assertIsNone(ms.rank_of(off, "Z"))
        self.assertEqual(ms.rank_of([cand("W1", 5), cand("T", 5)], "T"), 2)
        self.assertEqual(ms.rank_of([cand("T", 5), cand("W1", 5)], "T"), 2)
        self.assertEqual([c["block"] for c in ms.canonical([cand("B", 1), cand("A", 1), cand("C", 2)])], ["C", "A", "B"])

    def test_boost_wrong(self):
        off = [cand("T", 5), cand("W1", 4), cand("W2", 3)]
        self.assertTrue(ms.boost_wrong(off, [cand("W1", 6, 4), cand("T", 5), cand("W2", 3)], "T"))
        self.assertTrue(ms.boost_wrong(off, [cand("W1", 5, 4), cand("T", 5), cand("W2", 3)], "T"))  # tie counts
        self.assertFalse(ms.boost_wrong([cand("W1", 6), cand("T", 5)], [cand("T", 7, 5), cand("W1", 6)], "T"))
        self.assertFalse(ms.boost_wrong([cand("W1", 6), cand("T", 5)], [cand("W1", 7, 6), cand("T", 5)], "T"))  # was above
        self.assertFalse(ms.boost_wrong([cand("W1", 6)], [cand("W1", 7)], "T"))  # truth absent

    def test_spectrum_row_flags(self):
        parents = {"C4": "C3", "C3": "C2", "C2": "C1", "C1": "C0"}
        rec = dict(record_id=1, set="selected", family="f", instrument_class="unknown", cluster_id="c1")
        off = [cand("W", 6, onto="Leaf"), cand("T", 5)]
        on = [cand("T", 6.5, 5), cand("W", 6, onto="Leaf")]
        row = ms.spectrum_row(rec, "T", "C3", on, off, parents, {"Leaf": "C4"}, {})
        self.assertEqual((row["rank"], row["top1"], row["top5"], row["boost_wrong"]), (1, True, True, False))
        self.assertTrue(row["rank_changed"] and row["order_changed"])
        self.assertFalse(row["base_score_mismatch"])
        self.assertFalse(row["class_correct"])  # the top-1 has no ontology
        off_row = ms.spectrum_row(rec, "T", "C3", off, None, parents, {"Leaf": "C4"}, {})
        self.assertTrue(off_row["class_correct"])
        self.assertNotIn("boost_wrong", off_row)
        tied = ms.spectrum_row(rec, "T", "C3", [cand("T", 5), cand("W", 6, onto="Leaf")], [cand("W", 6), cand("T", 5)],
                               parents, {}, {})
        self.assertFalse(tied["order_changed"])  # the written order does not matter
        self.assertTrue(ms.spectrum_row(rec, "T", None, [cand("T", 5, 4)], [cand("T", 5)], parents, {}, {})["base_score_mismatch"])

    def test_cell_record_and_cluster_level(self):
        rows = [dict(rank=1, top1=True, top5=True, found=True, class_defined=True, class_correct=True, cluster_id="a",
                     boost_wrong=False, rank_changed=False, order_changed=False),
                dict(rank=None, top1=False, top5=False, found=False, class_defined=False, class_correct=False,
                     cluster_id="a", boost_wrong=False, rank_changed=False, order_changed=False),
                dict(rank=7, top1=False, top5=False, found=True, class_defined=True, class_correct=False, cluster_id="b",
                     boost_wrong=True, rank_changed=True, order_changed=True)]
        c = ms.cell(rows, on=True)
        self.assertEqual((c["n"], c["top1"], c["not_found"], c["rank_le10"], c["median_rank_found"]), (3, 1, 1, 2, 4.0))
        self.assertEqual((c["class_n"], c["class_correct"], c["class_rate"]), (2, 1, 0.5))
        self.assertEqual(c["boost_wrong_rate"], round(1 / 3, 6))
        self.assertEqual(c["c_top1_rate"], 0.25)  # (1/2 + 0/1) / 2
        self.assertEqual(c["c_class_rate"], 0.5)
        self.assertIsNone(ms.cell(rows, on=False)["boost_wrong_rate"])

    def test_choose_weight(self):
        w = ms.WEIGHTS
        self.assertEqual(ms.choose_weight({"off": 3, "w0.5": 4, "w1.0": 5, "w2.0": 5, "w4.0": 2}, w), 1.0)
        self.assertIsNone(ms.choose_weight({"off": 5, "w0.5": 5, "w1.0": 4, "w2.0": 5, "w4.0": 1}, w))
        self.assertEqual(ms.choose_weight({"off": 0, "w0.5": 1, "w1.0": 1, "w2.0": 1, "w4.0": 1}, w), 0.5)

    def test_suppression_only_vs21_by_denominator(self):
        c = dict(n=4, n_clusters=3, top1=2, top1_rate=0.5, class_n=4, class_correct=1, class_rate=0.25,
                 c_top1_rate=0.3, c_class_n=3, c_class_rate=0.1)
        self.assertEqual(ms.suppress_cell(c, "vs20")["top1"], 2)
        s = ms.suppress_cell(c, "vs21")
        self.assertIsNone(s["top1"])
        self.assertIsNone(s["c_class_rate"])
        self.assertEqual((s["n"], s["class_n"]), (4, 4))
        big = ms.suppress_cell(dict(c, n=9, n_clusters=9, class_n=3, c_class_n=6), "vs21")
        self.assertEqual(big["top1"], 2)
        self.assertIsNone(big["class_rate"])
        self.assertEqual(big["c_class_rate"], 0.1)


class FilesTest(unittest.TestCase):
    def test_structure_result_parse(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "Structure result-1.txt"
            head = ("File path\tFile name\tTitle\tMS1 count\tMSMS count\tPrecursor mz\tPrecursor type\tStructure rank\t"
                    "Total score\tDatabases\tFormula\tOntology\tInChIKey\tSMILES\tRanking score\tBase rank\tDiagnostic rule term")
            body = ["x\trec5\tq\t0\t3\t1\t[M-H]-\tA\t9\tdb\tC1\tOnt\tAAAA-B-N\tC\t9\t2\t1",
                    "x\trec5\tq\t0\t3\t1\t[M-H]-\tB\t8\tdb\tC1\t\tBBBB-B-N\tC\t8\t1\t0",
                    "x\trec7.local-vs21\tq\t0\t3\t1\t[M-H]-\tC\t8\tdb\tC1\t\tCCCC-B-N\tC\t8\t1\t0"]
            p.write_text("\n".join([head] + body) + "\n", encoding="ascii")
            r = ms.parse_structure_result(p)
            self.assertEqual([c["block"] for c in r[5]], ["AAAA", "BBBB"])
            self.assertEqual((r[5][0]["term"], r[5][0]["score"], r[5][0]["total"]), (1.0, 9.0, 9.0))
            self.assertIn(7, r)

    def test_method_ini(self):
        self.assertNotIn("DiagnosticRule", ms.method_ini(None))
        on = ms.method_ini(2.0)
        self.assertIn("DiagnosticRuleEvaluation=TRUE\r\n", on)
        self.assertIn("DiagnosticRuleRankingWeight=2.0\r\n", on)
        self.assertIn("DiagnosticRuleRankingScope=Reviewed\r\n", on)
        self.assertIn("PubChemNeverUse=True\r\n", on)
        self.assertTrue(on.startswith(ms.method_ini(None)))
        self.assertEqual((ms.variant_name(None), ms.variant_name(0.5)), ("off", "w0.5"))

    def test_query_msp(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "rec3.msp"
            ms.write_query_msp(p, dict(record_id=3, precursor=301.5, adduct="[M-H]-", mode="neg", instrument_type="",
                                       collision_energy="20", mz=[100.0, 150.25], intensity=[5.0, 10.0]))
            text = p.read_text(encoding="utf-8")
            self.assertTrue(text.startswith("NAME: rec3\nPRECURSORMZ: 301.5\nPRECURSORTYPE: [M-H]-\nIONMODE: Negative\n"))
            self.assertNotIn("INSTRUMENTTYPE", text)
            self.assertIn("Num Peaks: 2\n100.0\t5.0\n150.25\t10.0\n", text)

    def test_validate_job_refuses_wrong_rule_set(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            (out / "Structure result-1.txt").write_text("h\n", encoding="ascii")
            self.assertIn("result", ms.validate_job(out, None, 1))
            m = dict(rule_set=dict(version="0.1.0", sha256="x", pinned=True), settings=dict(evaluation=True, ranking_weight=1.0,
                     ranking_scope="reviewed"), counts=dict(errors=0, refused=0, queries=1), build_identity="b",
                     engine_semantics="e")
            (out / "Diagnostic rule run-1.json").write_text(json.dumps(m), encoding="utf-8")
            with self.assertRaises(RuntimeError):
                ms.validate_job(out, None, 1)  # an off job must not evaluate rules
            with self.assertRaises(RuntimeError):
                ms.validate_job(out, 1.0, 1)
            m["rule_set"] = dict(version=ms.EXPECTED_RULESET_VERSION, sha256=ms.EXPECTED_RULESET_SHA256, pinned=True)
            for k in ("evaluated", "fired", "veto_blocked", "none", "skipped"):
                m["counts"][k] = 0
            (out / "Diagnostic rule run-1.json").write_text(json.dumps(m), encoding="utf-8")
            self.assertEqual(ms.validate_job(out, 1.0, 1)["build_identity"], "b")
            with self.assertRaises(RuntimeError):
                ms.validate_job(out, 2.0, 1)


class LoadWaitTest(unittest.TestCase):
    def test_wait_loaded(self):
        class Proc:
            def __init__(self, code=None):
                self.code, self.killed = code, False

            def poll(self):
                return self.code

            def kill(self):
                self.killed = True

        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / "console.log"
            log.write_text("Loading library files..\n" + ms.LOADED_MARK + "\n", encoding="utf-8")
            self.assertTrue(ms.wait_loaded(Proc(), log))
            log.write_text("Loading library files..\nEI fragment DB is missing.\n", encoding="utf-8")
            self.assertFalse(ms.wait_loaded(Proc(code=255), log))
            proc = Proc()
            with self.assertRaises(RuntimeError):
                ms.wait_loaded(proc, log, timeout=0, poll=0)
            self.assertTrue(proc.killed)


class GuardTest(unittest.TestCase):
    def test_run_id(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            ms.check_out(root / "p3ms-20261125-abcdef1", False, "abcdef1234", out_root=root)
            ms.check_out(root / "p3mssmoke-20261125-abcdef1", True, "abcdef1234", out_root=root)
            for name, smoke in (("p3ms-20261125-abcdef1", True), ("p3mssmoke-20261125-abcdef1", False),
                                ("p3cs-20261125-abcdef1", False), ("p3ms-20261125-1234567", False)):
                with self.assertRaises(ValueError):
                    ms.check_out(root / name, smoke, "abcdef1234", out_root=root)

    def test_vs21_guard_refusals(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            args = SimpleNamespace(vs21=False, selection_run=None, console=None, ruleset_dir=None)
            with self.assertRaises(PermissionError):
                ms.vs21_guard(args, out_root=root)
            args.vs21 = True
            with self.assertRaises(PermissionError):
                ms.vs21_guard(args, out_root=root)  # no selection run
            for name in ("p3mssmoke-20261125-abcdef1", "p3ms-20261125-abcdef1"):
                (root / name).mkdir()
            args.selection_run = root / "p3mssmoke-20261125-abcdef1"
            with self.assertRaises(PermissionError):
                ms.vs21_guard(args, out_root=root)  # a smoke never selects
            sel = root / "p3ms-20261125-abcdef1"
            args.selection_run = sel
            with self.assertRaises(PermissionError):
                ms.vs21_guard(args, out_root=root)  # no files
            res = sel / "selection_result.json"
            p2.write_json(res, dict(binding=True, smoke=False, chosen_weight=None))
            p2.write_json(sel / "run.json", dict(status="complete", parameters=dict(source="vs20"),
                                                 files={"selection_result.json": p2.sha256(res)}))
            with self.assertRaisesRegex(PermissionError, "skipped"):
                ms.vs21_guard(args, out_root=root)
            p2.write_json(sel / "run.json", dict(status="running", parameters=dict(source="vs20"),
                                                 files={"selection_result.json": p2.sha256(res)}))
            with self.assertRaises(PermissionError):
                ms.vs21_guard(args, out_root=root)

    def test_vs21_parse_refuses_without_selection(self):
        with contextlib.redirect_stderr(io.StringIO()):
            self._parse_refusals()
        _, a = ms.parse_args(["select", "--out", "x", "--console", "c", "--ruleset-dir", "r", "--chemont-obo", "o"])
        self.assertEqual(a.weights, ms.WEIGHTS)
        _, a = ms.parse_args(["select", "--out", "x", "--limit", "20", "--console", "c", "--ruleset-dir", "r",
                              "--chemont-obo", "o"])
        self.assertEqual(a.weights, ms.SMOKE_WEIGHTS)

    def _parse_refusals(self):
        with self.assertRaises(SystemExit):
            ms.parse_args(["select", "--out", "x", "--source", "vs21", "--vs21", "--console", "c", "--ruleset-dir", "r",
                           "--chemont-obo", "o"])
        with self.assertRaises(SystemExit):
            ms.parse_args(["select", "--out", "x", "--source", "vs20", "--vs21", "--console", "c", "--ruleset-dir", "r",
                           "--chemont-obo", "o"])
        with self.assertRaises(SystemExit):
            ms.parse_args(["select", "--out", "x", "--workers", "9", "--console", "c", "--ruleset-dir", "r",
                           "--chemont-obo", "o"])


if __name__ == "__main__":
    unittest.main()
