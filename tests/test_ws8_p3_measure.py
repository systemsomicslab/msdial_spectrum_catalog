"""WS8 P3 rule-level measurement (tools/ws8_p3/p3_measure.py), on synthetic spectra, hand-built counts and
public molecules only. The console integration test runs only when WS8_P3_CONSOLE and WS8_P3_RULESET are set."""
from collections import Counter
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import random
import re
import tempfile
from types import SimpleNamespace
import unittest

from tools.ws8_p2 import p2_diagnostics as p2
from tools.ws8_p3 import p3_measure as m

HAS_CHEM = importlib.util.find_spec("rdkit") is not None and importlib.util.find_spec("pyarrow") is not None
CONSOLE = os.environ.get("WS8_P3_CONSOLE")
RULESET = os.environ.get("WS8_P3_RULESET")
# The console's peak-line grammar (PeakBlockCheck.PairLine).
NUMBER = r"[0-9]*\.?[0-9]+(?:[eE]\+[0-9]+)?"
PAIR = re.compile(r"^[ \t]*(" + NUMBER + r")[ \t]+(" + NUMBER + r")[ \t]*$")
ADDUCTS = [dict(adduct="[M+H]+", aliases=["M+H"], charge=1), dict(adduct="[M+Na]+", aliases=["M+Na"], charge=1),
           dict(adduct="[M-H]-", aliases=["M-H"], charge=-1),
           dict(adduct="[M+HCOO]-", aliases=["[M+FA-H]-", "M+FA-H"], charge=-1)]
MOTIFS = [dict(motif_id="carbonyl", smarts="[CX3]=O"), dict(motif_id="phenol", smarts="c[OX2H]")]
CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
THEOBROMINE = "CN1C=NC2=C1C(=O)NC(=O)N2C"
THEOPHYLLINE = "CN1C2=C(C(=O)N(C1=O)C)NC=N2"
QUERCETIN = "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O"
PROTOPINE = "CN1CCC2=CC3=C(C=C2C(=O)CC4=C(C1)C5=C(C=C4)OCO5)OCO3"


def msp_record(smiles=None, mz="195.0877", ptype="[M+H]+", mode="Positive", peaks=((138.0662, 100), (110.0713, 40)),
               comment=None, extra=()):
    lines = ["NAME: synthetic"]
    if mz is not None:
        lines.append(f"PRECURSORMZ: {mz}")
    if ptype is not None:
        lines.append(f"PRECURSORTYPE: {ptype}")
    if mode is not None:
        lines.append(f"IONMODE: {mode}")
    lines += ["INSTRUMENTTYPE: LC-ESI-QTOF", "COLLISIONENERGY: 20"]
    if smiles is not None:
        lines.append(f"SMILES: {smiles}")
    if comment is not None:
        lines.append(f"COMMENT: {comment}")
    lines += list(extra)
    lines.append(f"Num Peaks: {len(peaks)}")
    lines += [f"{a}\t{b}" for a, b in peaks]
    return "\n".join(lines) + "\n\n"


def quiet(fn, *args, **kw):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return fn(*args, **kw)


class ReaderTests(unittest.TestCase):
    def test_reader_equals_p2_iter_msp(self):
        text = (msp_record(CAFFEINE, comment="MassBank EU") + "NAME: no meta peaks\n\n" +
                msp_record(QUERCETIN, mz="301.0354", ptype="[M-H]-", mode="N",
                           peaks=((151.0031, 100), (178.998, 20), ("bad", 3), (-1, 5), (12.0, 0)), comment="GNPS in silico") +
                "PRECURSORMZ: 100\nNum Peaks: 1\n50.1 1e2\n")
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "x.msp"
            path.write_text("﻿" + text, encoding="utf-8")
            mine = list(m.iter_msp_p3(path, m.load_provenance()))
            theirs = list(p2.iter_msp(path))
        self.assertEqual([x[:3] for x in mine], theirs)
        self.assertEqual([x[3] for x in mine], [dict(bucket="massbank", in_silico=False), dict(bucket="gnps", in_silico=True),
                                                dict(bucket="unknown", in_silico=False)])

    def test_provenance_buckets_and_unavailable(self):
        t = m.load_provenance()
        self.assertEqual(m.provenance_of({"comment": "ReSpect PS012345"}, t)["bucket"], "respect")
        self.assertEqual(m.provenance_of({"comment": "lab standard"}, t)["bucket"], "other")
        self.assertEqual(m.provenance_of({"comment": ""}, t)["bucket"], "unknown")

    def test_provenance_table_rejects_unlisted_field(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "p.tsv"
            p.write_text("order\tfield\tpattern\toutcome\tnote\n1\tname\tx\tgnps\t\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                m.load_provenance(p)

    def test_alias_map_and_charge(self):
        a = m.alias_map(ADDUCTS)
        self.assertEqual(a["[M+FA-H]-"], ("[M+HCOO]-", -1))
        self.assertEqual(a["M+H"], ("[M+H]+", 1))
        self.assertEqual(m.adduct_charge("[M+2H]2+", None), 2)
        self.assertEqual(m.adduct_charge("[M+Foo]+", None), 1)
        self.assertEqual(m.adduct_charge("[M-H]-", a["[M-H]-"]), 1)

    def test_msp_numbers_follow_console_grammar(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "a.msp"
            m.write_msp(path, [dict(record_id=7, precursor=1e-5 + 300, adduct="[M+H]+", mode="pos", instrument_type="",
                                    collision_energy="", mz=[100.123456789, 2e-5], intensity=[1e-05, 1e20])])
            lines = path.read_text(encoding="utf-8").splitlines()
        peaks = lines[lines.index("Num Peaks: 2") + 1:lines.index("Num Peaks: 2") + 3]
        self.assertTrue(all(PAIR.match(x) for x in peaks), peaks)
        self.assertEqual(lines[0], "NAME: rec7")
        self.assertIn("IONMODE: Positive", lines)

    def test_shard_blocks(self):
        self.assertEqual(m.shard_blocks(10, 3), [(0, 4), (4, 8), (8, 10)])
        self.assertEqual(m.shard_blocks(2, 8), [(0, 1), (1, 2)])
        self.assertEqual(m.shard_blocks(0, 8), [(0, 0)])

    def test_p2_parity_is_a_hard_stop(self):
        theirs = {"vs20:pos:input": 3, "vs20:pos:nominal_mass": 1, "vs21:pos:input": 9}
        self.assertEqual(m.p2_parity(Counter({"vs20:pos:input": 3, "vs20:pos:nominal_mass": 1}), theirs, "vs20")["result"], "match")
        with self.assertRaises(RuntimeError):
            m.p2_parity(Counter({"vs20:pos:input": 3}), theirs, "vs20")


@unittest.skipUnless(HAS_CHEM, "RDKit and pyarrow are required")
class PrepareTests(unittest.TestCase):
    def build(self, d):
        pos = (msp_record(CAFFEINE, comment="MassBank") +
               msp_record(CAFFEINE, peaks=((138.0662, 100), (110.0713, 40), (83.0604, 5))) +
               msp_record(CAFFEINE, mz="217.0696", ptype="[M+Na]+") +
               msp_record(CAFFEINE, peaks=((138.0, 100), (110.0, 40))) +
               msp_record(CAFFEINE, peaks=()) +
               msp_record(CAFFEINE, mz=None) +
               msp_record(CAFFEINE, mode="Negative") +
               msp_record(CAFFEINE, ptype=None) +
               msp_record("not a smiles") +
               msp_record(CAFFEINE, mz="98.0475", ptype="[M+2H]2+") +
               msp_record(CAFFEINE, comment="in silico prediction") +
               msp_record(CAFFEINE, ptype="[M+Foo]+") +
               msp_record(THEOBROMINE, mz="181.0720") + msp_record(THEOPHYLLINE, mz="181.0720"))
        neg = msp_record(QUERCETIN, mz="301.0354", ptype="[M+FA-H]-", mode="Negative", peaks=((151.0031, 100), (178.9980, 20)))
        (Path(d) / "pos.msp").write_text(pos, encoding="utf-8")
        (Path(d) / "neg.msp").write_text(neg, encoding="utf-8")
        return [("pos", Path(d) / "pos.msp", None), ("neg", Path(d) / "neg.msp", None)]

    def run_prepare(self, d, out):
        out.mkdir()
        return quiet(m.prepare_source, "vs20", self.build(d), out, None, 1, MOTIFS, m.alias_map(ADDUCTS),
                     m.load_provenance(), 2)

    def test_prepare_counts_units_clusters_and_msp(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "run"
            counts, extra, stats, keys, hashes = self.run_prepare(d, out)
            self.assertEqual(counts["vs20:pos:input"], 14)
            for k in ("nominal_mass", "no_peaks", "precursor", "mode_mismatch", "missing_precursor_type"):
                self.assertEqual(counts["vs20:pos:" + k], 1, k)
            self.assertEqual(counts["vs20:structure_rejected_records"], 1)
            self.assertEqual(extra["vs20:charge_ge2_records"], 1)
            self.assertEqual(extra["vs20:in_silico_records"], 1)
            self.assertEqual(extra["vs20:adduct_unmapped_records"], 1)
            # caffeine [M+H]+ (2 records, one unit), [M+Na]+, [M+Foo]+, theobromine, theophylline, quercetin [M+HCOO]-
            self.assertEqual(stats["representatives"], 6)
            self.assertEqual(stats["clusters"], 5)  # theobromine and theophylline share formula, adduct and motifs
            rows = list(p2.parquet_rows(p2.row_path(out, "prepared", "vs20")))
            caf = [r for r in rows if r["adduct"] == "[M+H]+" and r["formula"] == "C8H10N4O2"]
            self.assertEqual(len(caf), 1)
            self.assertEqual(len(caf[0]["mz"]), 3)
            self.assertEqual(caf[0]["spectrum_count"], 2)  # the in-silico record is excluded
            self.assertEqual([r["adduct"] for r in rows if r["mode"] == "neg"], ["[M+HCOO]-"])
            self.assertIn("phenol", next(r for r in rows if r["mode"] == "neg")["motifs"])
            iso = [r for r in rows if r["formula"] == "C7H8N4O2"]
            self.assertEqual(len({r["cluster_id"] for r in iso}), 1)
            self.assertEqual(len({r["key"] for r in iso}), 2)
            index = m.read_rows(out / "msp/vs20/representatives-index.tsv")
            self.assertEqual(sorted(int(r["record_id"]) for r in index), sorted(r["record_id"] for r in rows))
            for f in sorted((out / "msp/vs20").glob("*.msp")):
                text = f.read_text(encoding="utf-8")
                self.assertNotIn("SMILES", text.upper())
                self.assertNotIn("COMMENT", text.upper())

    def test_prepare_is_deterministic(self):
        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            self.run_prepare(d1, Path(d1) / "run")
            self.run_prepare(d2, Path(d2) / "run")
            for name in ("representatives-00.msp", "representatives-01.msp", "representatives-index.tsv"):
                self.assertEqual((Path(d1) / "run/msp/vs20" / name).read_bytes(), (Path(d2) / "run/msp/vs20" / name).read_bytes())
            a = list(p2.parquet_rows(p2.row_path(Path(d1) / "run", "prepared", "vs20")))
            b = list(p2.parquet_rows(p2.row_path(Path(d2) / "run", "prepared", "vs20")))
            self.assertEqual(a, b)

    def test_smoke_limit_samples_by_ordinal(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "run"
            out.mkdir()
            counts, *_ = quiet(m.prepare_source, "vs20", self.build(d), out, 0.5, 1, MOTIFS, m.alias_map(ADDUCTS),
                               m.load_provenance(), 1)
            self.assertEqual(counts["vs20:pos:smoke_unsampled"], 7)
            self.assertEqual(counts["vs20:pos:input"], 14)


def rep(rid, key, mode, mz, peaks, adduct="[M+H]+"):
    return dict(record_id=rid, key=key, mode=mode, precursor=mz, adduct=adduct, instrument_type="", collision_energy="",
                mz=[p[0] for p in peaks], intensity=[p[1] for p in peaks])


class ChimericTests(unittest.TestCase):
    def reps(self):
        return [rep(0, "A", "pos", 300.0, [(100.0, 50), (150.0, 100), (300.0, 999)]),
                rep(2, "B", "pos", 310.0, [(120.0, 10), (200.0, 40)]),
                rep(4, "A", "pos", 305.0, [(90.0, 10)]),
                rep(6, "C", "pos", 500.0, [(250.0, 10)]),
                rep(1, "D", "neg", 301.0, [(80.0, 10)])]

    def test_partner_window_mode_key_and_scaling(self):
        rows, stats = m.build_chimeric(self.reps(), 7)
        self.assertEqual(stats["no_partner"], 2)  # record 6 (nothing within 12.5 Da) and the lone negative record
        by = {(r["record_id"], r["level"]): r for r in rows}
        r = by[(0, 0.3)]
        self.assertEqual(r["partner_record_id"], 2)  # record 4 shares the first block
        # The partner's highest peak becomes 0.3 x the primary base peak (300.0 is the precursor, excluded).
        self.assertAlmostEqual(dict(zip(r["mz"], r["intensity"]))[200.0], 30.0)
        self.assertEqual(r["mz"], sorted(r["mz"]))
        self.assertEqual(by[(4, 0.1)]["partner_record_id"], 2)
        self.assertNotIn((1, 0.1), by)

    def test_chimeric_is_seeded(self):
        reps = [rep(i * 2, f"K{i}", "pos", 300 + i * 0.5, [(100.0 + i, 10)]) for i in range(40)]
        a, _ = m.build_chimeric(reps, 11)
        b, _ = m.build_chimeric(list(reversed(reps)), 11)
        c, _ = m.build_chimeric(reps, 12)
        key = lambda rows: [(r["record_id"], r["level"], r["partner_record_id"]) for r in sorted(rows, key=lambda r: (r["record_id"], r["level"]))]
        self.assertEqual(key(a), key(b))
        self.assertNotEqual(key(a), key(c))


def evidence_line(name, i, outcomes, status="evaluated", iclass="QTOF-type"):
    res = None if status != "evaluated" else dict(instrument_class=iclass, outcomes=[
        dict(rule_id=r, version=v, fired=f, n_hypotheses=n, veto_blocked=False, chosen=dict(adduct="[M+H]+") if f else None)
        for r, v, f, n in outcomes])
    return json.dumps(dict(collision_energy_bin="x", file_name=name, input={}, record_index=i, result=res, status=status)) + "\n"


class JoinTests(unittest.TestCase):
    def shard(self, d, lines_by_pass):
        for p, lines in lines_by_pass.items():
            (Path(d) / f"Diagrules evidence-{p}-20261123000000.jsonl").write_text("".join(lines), encoding="utf-8")

    def test_join_and_collect(self):
        with tempfile.TemporaryDirectory() as d:
            ids = [("vs20", 0), ("vs20", 2), ("vs20", 4)]
            self.shard(d, {"base": [evidence_line("a.msp", 0, [("r1", 1, True, 2)]), evidence_line("a.msp", 1, [("r1", 1, False, 1)]),
                                    evidence_line("a.msp", 2, [], status="skipped:profile")],
                           "mz_only_curated": [evidence_line("a.msp", i, []) for i in range(3)]})
            ev = m.read_evidence(d, "a.msp", ids, ("base", "mz_only_curated"))
            by_rule, notev, inst = m.collect([ev], {"r1": 1})
            self.assertEqual([(e[0], e[1]) for e in by_rule["base"]["r1"]], [(("vs20", 0), True), (("vs20", 2), False)])
            self.assertEqual(notev["base"], {("vs20", 4)})
            self.assertEqual(inst[("vs20", 0)], "QTOF-type")
            with self.assertRaises(RuntimeError):
                m.collect([ev], {"r1": 2})

    def test_join_refuses_missing_or_foreign_lines(self):
        with tempfile.TemporaryDirectory() as d:
            self.shard(d, {"base": [evidence_line("a.msp", 0, [])]})
            with self.assertRaises(RuntimeError):
                m.read_evidence(d, "a.msp", [("vs20", 0), ("vs20", 2)], ("base",))
            with self.assertRaises(RuntimeError):
                m.read_evidence(d, "b.msp", [("vs20", 0)], ("base",))
            with self.assertRaises(RuntimeError):
                m.read_evidence(d, "a.msp", [("vs20", 0)], ("floor_x2",))


def rule(rid="r1", motif="M", origin="as_published", axis=0.01, mode="adduct_ambiguous", competitors=("Q",), **kw):
    return dict(dict(rule_id=rid, version=1, ion_mode="pos", motif_id=motif, competitors=list(competitors), measure=True,
                     source_checked=True, primary_mode=mode, origin=origin,
                     axis={"adduct_ambiguous": dict(occupancy=axis, lo_da=100, hi_da=101, n_hypotheses=1),
                           "mz_only": dict(occupancy=axis, lo_da=100, hi_da=101, n_hypotheses=1)}), **kw)


def scenario(tp_clusters=160, fp_clusters=0, members=2, neg_clusters=6000, decoy=0, neg_cond=0, comp=0, chim_fp=0):
    """Hand-built truth and fires: each fired cluster has `members` records that all fire."""
    truth, entries, rid = {}, [], 0

    def add(cluster, motifs, fired, cond=True, inst="QTOF-type"):
        nonlocal rid
        uid = ("vs20", rid)
        truth[uid] = dict(mode="pos", motifs=frozenset(motifs), cluster=cluster, instrument=inst, energy="eV:20-<40",
                          bucket="massbank", source="vs20")
        if cond:
            entries.append((uid, fired, 1, False, "[M+H]+"))
        rid += 2
        return uid

    for c in range(tp_clusters):
        for _ in range(members):
            add(f"tp{c}", {"M"}, True)
    for c in range(fp_clusters):
        for _ in range(members):
            add(f"fp{c}", {"Q"} if c < comp else set(), True)
    for c in range(neg_clusters):
        add(f"n{c}", set(), False, cond=c < neg_cond)
    for c in range(30):
        add(f"pm{c}", {"M"}, False, cond=False)
    by_rule = {p: {"r1": list(entries)} for p in m.PASSES}
    fired = [e for e in entries if e[1]]
    for d in m.DECOY_PASSES:
        by_rule[d] = {"r1": [(e[0], i < decoy, 1, False, None) for i, e in enumerate(fired)]}
    chim = dict(by_rule={lv: {"base": {"r1": [((u[0], u[1], lv), True, 1, False, None) for u, *_ in fired]}}
                         for lv in m.CHIMERIC_LEVELS},
                universe={lv: [(u[0], u[1], lv) for u in truth] for lv in m.CHIMERIC_LEVELS})
    if chim_fp:
        negs = [u for u in truth if not truth[u]["motifs"]][:chim_fp]
        for lv in m.CHIMERIC_LEVELS:
            chim["by_rule"][lv]["base"]["r1"] += [((u[0], u[1], lv), True, 1, False, None) for u in negs]
    return truth, by_rule, chim


class MeasureTests(unittest.TestCase):
    def core(self, rules=None, sc=None, **kw):
        truth, by_rule, chim = sc or scenario()
        return m.measure_core(rules or [rule()], truth, by_rule, {}, ["vs20"], replicates=200, chimeric=chim, **kw)

    def gate_row(self, res, p="base"):
        return next(r for r in res["precision"] if r["pass_id"] == p and (r["source"], r["instrument"], r["energy"]) == ("pooled", "all", "all"))

    def test_record_and_cluster_counts(self):
        res = self.core(sc=scenario(tp_clusters=120, fp_clusters=3, members=2, neg_cond=10))
        r = self.gate_row(res)
        self.assertEqual((r["fires"], r["tp"], r["fp"]), (246, 240, 6))
        self.assertEqual((r["c_fires"], r["c_tp"], r["c_fp"]), (123, 120, 3))
        self.assertEqual(r["n_cond"], 256)
        self.assertEqual(r["n_positive"], 270)
        self.assertEqual((r["c_positive"], r["c_cond_positive"]), (150, 120))
        self.assertAlmostEqual(r["coverage"], 240 / 270)
        self.assertEqual((r["wilson_lo"], r["wilson_hi"]), p2.wilson(240, 246))
        self.assertEqual((r["c_wilson_lo"], r["c_wilson_hi"]), p2.wilson(120, 123))
        st = p2.statistics(120, 3, 30, r["c_universe"] - 150 - 3)
        self.assertAlmostEqual(r["c_prior_stress"], st["prior_stress_precision"])
        self.assertLessEqual(r["c_boot_lo"], 120 / 123)
        self.assertGreaterEqual(r["c_boot_hi"], 120 / 123)

    def test_strata_cells(self):
        res = self.core()
        base = [r for r in res["precision"] if r["pass_id"] == "base"]
        self.assertEqual({(r["instrument"], r["energy"]) for r in base},
                         {("all", "all"), ("QTOF-type", "all"), ("all", "eV:20-<40"), ("QTOF-type", "eV:20-<40")})
        other = [r for r in res["precision"] if r["pass_id"] == "floor_x2"]
        self.assertEqual(len(other), 1)

    def test_decoy_fdr_is_mean_over_shifts(self):
        res = self.core(sc=scenario(tp_clusters=120, members=1, decoy=6))
        r = self.gate_row(res)
        self.assertEqual(r["decoy_fires_plus_5"], 6)
        self.assertAlmostEqual(r["decoy_fdr"], 6 / 120)
        self.assertEqual(res["gates"]["r1"]["status"], "fail")
        self.assertIn("decoy_fdr", res["gates"]["r1"]["failed"])

    def test_gate_not_estimable_below_100_clusters(self):
        res = self.core(sc=scenario(tp_clusters=99, members=3))
        self.assertEqual(res["gates"]["r1"]["status"], "not_estimable")

    def test_gate_incomplete_when_measured_criteria_met(self):
        g = self.core()["gates"]["r1"]
        self.assertEqual(g["status"], "incomplete")
        self.assertEqual(g["failed"], [])
        self.assertEqual(set(g["pending"]), {"synthetic_insource", "heldout"})

    def test_gate_fails_on_wilson(self):
        g = self.core(sc=scenario(tp_clusters=120, fp_clusters=4))["gates"]["r1"]
        self.assertEqual(g["status"], "fail")
        self.assertIn("cluster_wilson_lo", g["failed"])

    def test_gate_non_conditional_by_axis_and_by_records(self):
        g = self.core(rules=[rule(axis=0.2)])["gates"]["r1"]
        self.assertTrue(g["non_conditional"])
        self.assertIn("occupancy", g["failed"])
        res = self.core(sc=scenario(neg_cond=400))
        occ = next(o for o in res["occupancy"] if o["precondition_mode"] == "adduct_ambiguous")
        self.assertAlmostEqual(occ["record_occupancy"], 400 / 6000)
        self.assertEqual(res["gates"]["r1"]["status"], "fail")

    def test_gate_prior_stress(self):
        # Few negatives in the universe: LR+ cannot be high enough at half the prevalence.
        g = self.core(sc=scenario(tp_clusters=120, neg_clusters=5))["gates"]["r1"]
        self.assertIn("prior_stress", g["failed"])

    def test_synthetic_chimeric_precision(self):
        res = self.core(sc=scenario(chim_fp=20))
        s = [r for r in res["synthetic"] if r["level"] == 0.3][0]
        self.assertEqual((s["fires"], s["fp"]), (340, 20))
        self.assertIn("synthetic_chimeric", res["gates"]["r1"]["failed"])

    def test_empty_chimeric_set_is_pending_not_passed(self):
        truth, by_rule, _ = scenario()
        chim = dict(by_rule={}, universe={})
        g = m.measure_core([rule()], truth, by_rule, {}, ["vs20"], replicates=10, chimeric=chim)["gates"]["r1"]
        self.assertIsNone(g["criteria"]["synthetic_chimeric"]["ok"])
        self.assertIn("synthetic_chimeric", g["pending"])

    def test_unevaluable_criteria_are_not_failures(self):
        truth, by_rule, chim = scenario(tp_clusters=0, members=1)
        g = m.measure_core([rule()], truth, by_rule, {}, ["vs20"], replicates=10, chimeric=chim)["gates"]["r1"]
        self.assertEqual(g["status"], "not_estimable")
        self.assertNotIn("cluster_wilson_lo", g["failed"])
        self.assertNotIn("decoy_fdr", g["failed"])
        self.assertIn("decoy_fdr", g["pending"])

    def test_competitor_fpr(self):
        res = self.core(sc=scenario(tp_clusters=120, fp_clusters=2, comp=2))
        c = res["competitors"][0]
        self.assertEqual((c["competitor_motif_id"], c["n_cond"], c["fires"], c["c_fires"]), ("Q", 4, 4, 2))

    def test_tuned_rules_and_not_measured(self):
        rules = [rule("r1", origin="tuned_vs20"), rule("r2", measure=False), rule("r3", source_checked=False), rule("r4")]
        res = self.core(rules=rules, refused=("r4",))
        self.assertIn("tuned_vs20", res["gates"]["r1"]["criteria"]["heldout"]["pending"])
        self.assertEqual([res["gates"][r]["status"] for r in ("r2", "r3", "r4")], ["not_measured"] * 3)
        self.assertEqual([res["gates"][r]["reason"] for r in ("r2", "r3", "r4")], ["measure_false", "source_unchecked", "refused"])
        self.assertFalse(any(r["rule_id"] != "r1" for r in res["precision"]))

    def test_mz_only_primary_uses_mz_only_pass(self):
        truth, by_rule, chim = scenario()
        by_rule["mz_only_curated"] = {"r1": by_rule["base"]["r1"][:50]}
        g = m.measure_core([rule(mode="mz_only")], truth, by_rule, {}, ["vs20"], replicates=10, chimeric=chim)["gates"]["r1"]
        self.assertEqual(g["gate_pass"], "mz_only_curated")
        self.assertEqual(g["status"], "not_estimable")

    def test_pooled_and_source_cells_with_two_sources(self):
        truth, by_rule, chim = scenario()
        res = m.measure_core([rule()], truth, by_rule, {}, ["vs20", "vs21"], replicates=10, chimeric=chim)
        self.assertEqual({r["source"] for r in res["precision"]}, {"pooled", "vs20", "vs21"})

    def test_not_evaluated_records_are_counted(self):
        truth, by_rule, chim = scenario()
        notev = {"mz_only_curated": {("vs20", 0), ("vs20", 2)}}
        res = m.measure_core([rule()], truth, by_rule, notev, ["vs20"], replicates=10, chimeric=chim)
        self.assertEqual(self.gate_row(res, "mz_only_curated")["n_not_evaluated"], 2)

    def test_measure_core_formula_free(self):
        """Circularity control 1: truth carries motifs and clusters only; no formula enters the fires or the gate."""
        truth, by_rule, chim = scenario()
        for v in truth.values():
            v["formula"] = "garbage"
        a = m.measure_core([rule()], truth, by_rule, {}, ["vs20"], replicates=20, chimeric=chim)
        for v in truth.values():
            v["formula"] = "C1"
        b = m.measure_core([rule()], truth, by_rule, {}, ["vs20"], replicates=20, chimeric=chim)
        self.assertEqual(a, b)


class StatisticsTests(unittest.TestCase):
    def test_bootstrap_deterministic_and_degenerate(self):
        t = Counter({(True, 1): 90, (False, 1): 10})
        a = m.cluster_bootstrap(t, 200, 300, random.Random("s"))
        b = m.cluster_bootstrap(t, 200, 300, random.Random("s"))
        self.assertEqual(a, b)
        self.assertLess(a["boot_lo"], 0.9)
        self.assertGreater(a["boot_hi"], 0.9)
        allpos = m.cluster_bootstrap(Counter({(True, 2): 50}), 50, 100, random.Random(1))
        self.assertEqual((allpos["boot_lo"], allpos["boot_hi"], allpos["c_boot_lo"]), (1.0, 1.0, 1.0))
        self.assertIsNone(m.cluster_bootstrap(Counter(), 10, 5, random.Random(1))["boot_lo"])

    def test_bootstrap_record_level_weights_units(self):
        r = m.cluster_bootstrap(Counter({(True, 9): 1, (False, 1): 1}), 2, 400, random.Random(3))
        self.assertLessEqual(r["boot_lo"], 0.9)
        self.assertEqual(r["boot_hi"], 1.0)

    def test_axis_occupancy(self):
        pda = 10 ** 12
        tol = 10 ** 10
        one = m.axis_occupancy([int(100.5 * pda)], tol)
        self.assertAlmostEqual(one["occupancy"], 0.02)
        self.assertEqual((one["lo_da"], one["hi_da"]), (100, 101))
        two = m.axis_occupancy([int(100.5 * pda), int(100.505 * pda), int(102.5 * pda)], tol)
        self.assertAlmostEqual(two["occupancy"], (0.025 + 0.02) / 3)
        edge = m.axis_occupancy([100 * pda], tol)
        self.assertAlmostEqual(edge["occupancy"], 0.01)
        self.assertIsNone(m.axis_occupancy([], tol)["occupancy"])

    def test_stress(self):
        self.assertAlmostEqual(m.stress(19.0, 0.05), 0.5)
        self.assertIsNone(m.stress(None, 0.1))


class OutputTests(unittest.TestCase):
    def outputs(self, d, local=False):
        truth, by_rule, chim = scenario(fp_clusters=2)
        rules = [rule(), rule("r2", measure=False)]
        res = m.measure_core(rules, truth, by_rule, {}, ["vs20"], replicates=50, chimeric=chim)
        doc = m.write_aggregates(d, res, "p3cs-20261123-abcdef1", rules, local)
        m.write_report(d, doc, dict(prepared=dict(vs20=dict(provenance_unavailable=True)), p3_exclusions={}), local)
        return doc

    def test_outputs_are_deterministic_and_pass_the_audit(self):
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            self.outputs(a)
            self.outputs(b)
            names = sorted(p.name for p in Path(a).iterdir())
            self.assertEqual(names, sorted(p.name for p in Path(b).iterdir()))
            self.assertIn("review_proposal-p3cs-20261123-abcdef1.json", names)
            for n in names:
                self.assertEqual((Path(a) / n).read_bytes(), (Path(b) / n).read_bytes(), n)
            self.assertEqual(m.audit_aggregates(a), "pass")
            doc = json.loads((Path(a) / "rules_precision.json").read_text(encoding="utf-8"))
            self.assertEqual(doc["label"], "ws8-p3-class-rules")
            self.assertEqual(doc["status_counts"], {"incomplete": 1, "not_measured": 1})
            report = (Path(a) / "report.md").read_text(encoding="utf-8")
            self.assertIn("provenance_unavailable", report)
            self.assertIn("in-source", report)
            prop = json.loads((Path(a) / "review_proposal-p3cs-20261123-abcdef1.json").read_text(encoding="utf-8"))
            self.assertEqual(prop["passed"], [])

    def test_audit_catches_forbidden_columns_and_keys(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "x.csv").write_text("rule_id,first_block\n", encoding="utf-8")
            with self.assertRaises(AssertionError):
                m.audit_aggregates(d)
            (Path(d) / "x.csv").write_text("rule_id\nBSYNRYMUTXBXSQ-UHFFFAOYSA-N\n", encoding="utf-8")
            with self.assertRaises(AssertionError):
                m.audit_aggregates(d)

    def test_vs21_outputs_named_local_and_suppressed(self):
        with tempfile.TemporaryDirectory() as d:
            self.outputs(d, local=True)
            names = sorted(p.name for p in Path(d).iterdir())
            self.assertTrue(all(".local-vs21." in n for n in names), names)
            rows = m.read_rows(Path(d) / "rules_precision.local-vs21.csv")
            gate = next(r for r in rows if (r["pass_id"], r["instrument"], r["energy"]) == ("base", "all", "all"))
            self.assertEqual(gate["c_fp"], "")  # 2 fired FP clusters: a cell below 5
            self.assertEqual(gate["lr_plus"], "")
            self.assertNotEqual(gate["c_precision"], "")  # its denominator is 162
            doc = json.loads((Path(d) / "rules_precision.local-vs21.json").read_text(encoding="utf-8"))
            self.assertNotIn("value", doc["rules"][0]["criteria"]["cluster_wilson_lo"])

    def test_suppress(self):
        row = dict(rule_id="r", version=1, fires=3, tp=3, fp=0, precision=1.0, n_positive=40, coverage=0.5, lr_plus=4.0)
        s = m.suppress(row)
        self.assertEqual((s["fires"], s["tp"], s["fp"], s["precision"], s["coverage"], s["lr_plus"]), (None, None, 0, None, 0.5, None))
        self.assertEqual(s["version"], 1)


class GuardTests(unittest.TestCase):
    def files(self, d):
        rs = Path(d) / "RuleSet"
        (rs / "rules").mkdir(parents=True)
        (rs / "rules/a.json").write_text("{}", encoding="utf-8")
        console = Path(d) / "console.exe"
        console.write_bytes(b"MZ")
        return rs, console

    def test_vs21_refused_without_flag_freeze_or_on_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            rs, console = self.files(d)
            with self.assertRaises(PermissionError):
                m.vs21_guard(SimpleNamespace(vs21=False, freeze=None, ruleset_dir=rs, console=console))
            with self.assertRaises(PermissionError):
                m.vs21_guard(SimpleNamespace(vs21=True, freeze=None, ruleset_dir=rs, console=console))
            freeze = Path(d) / "freeze.json"
            freeze.write_text(json.dumps(m.freeze_hashes(rs, console)), encoding="utf-8")
            args = SimpleNamespace(vs21=True, freeze=freeze, ruleset_dir=rs, console=console)
            self.assertTrue(m.vs21_guard(args))
            (rs / "rules/a.json").write_text('{"x": 1}', encoding="utf-8")
            with self.assertRaisesRegex(PermissionError, "ruleset_files"):
                m.vs21_guard(args)
            (rs / "rules/a.json").write_text("{}", encoding="utf-8")
            console.write_bytes(b"MZ2")
            with self.assertRaisesRegex(PermissionError, "console"):
                m.vs21_guard(args)

    def test_cli_refuses_vs21_before_touching_anything(self):
        with tempfile.TemporaryDirectory() as d:
            rs, _ = self.files(d)
            out = Path(d) / "p3cs-20261123-abcdef1"
            with self.assertRaises(SystemExit):
                quiet(m.main, ["prepare", "--sources", "vs20,vs21", "--ruleset-dir", str(rs), "--out", str(out)])
            self.assertFalse(out.exists())
            with self.assertRaises(SystemExit):
                quiet(m.main, ["prepare", "--sources", "vs21", "--vs21", "--ruleset-dir", str(rs), "--out", str(out)])

    def test_run_directory_checks(self):
        with tempfile.TemporaryDirectory() as root:
            ok = Path(root) / "p3cs-20261123-abcdef1"
            self.assertEqual(m.check_out(ok, None, "abcdef1234", root), ok.resolve())
            self.assertEqual(m.check_out(Path(root) / "p3cssmoke-20261123-abcdef1", 0.01, "abcdef1234", root).name,
                             "p3cssmoke-20261123-abcdef1")
            bad = [(Path(root) / "x" / ok.name, None, "abcdef1"), (Path(root) / "p3-20261123-abcdef1", None, "abcdef1"),
                   (Path(root) / "p3cssmoke-20261123-abcdef1", None, "abcdef1"), (ok, 0.01, "abcdef1"),
                   (ok, None, "1234567"), (Path(root) / "p3cs-2026112-abcdef1", None, "abcdef1")]
            for out, limit, commit in bad:
                with self.assertRaises(ValueError, msg=str(out)):
                    m.check_out(out, limit, commit, root)
            self.assertTrue(m.check_freeze_path(Path(root) / "freeze-20261123-abcdef1.json", root))
            with self.assertRaises(ValueError):
                m.check_freeze_path(Path(root) / "frozen.json", root)

    def test_run_state(self):
        args = SimpleNamespace(sources=["vs20"], workers=8, limit=None, seed=1, bootstrap=10)
        m.check_run_state({}, args)
        with self.assertRaises(ValueError):
            m.check_run_state(dict(status="complete"), args)
        with self.assertRaises(ValueError):
            m.check_run_state(dict(parameters=dict(sources=["vs20"], workers=4, limit=None, seed=1, bootstrap=10)), args)


@unittest.skipUnless(HAS_CHEM and CONSOLE and RULESET and Path(CONSOLE or "-").is_file(),
                     "needs RDKit, pyarrow, WS8_P3_CONSOLE and WS8_P3_RULESET")
class ConsoleIntegrationTests(unittest.TestCase):
    def test_prepare_fire_measure_on_synthetic_records(self):
        rule_id = "alkaloid_benzylisoquinoline.protopine.pos"
        with tempfile.TemporaryDirectory() as d:
            pos = (msp_record(PROTOPINE, mz="354.1336", peaks=((336.1230, 30), (206.0812, 100), (149.0597, 60), (188.0706, 10))) +
                   msp_record(CAFFEINE) + msp_record(THEOBROMINE, mz="181.0720") +
                   msp_record(QUERCETIN, mz="303.0499", peaks=((153.0182, 100), (229.0495, 20))))
            neg = msp_record(QUERCETIN, mz="301.0354", ptype="[M-H]-", mode="Negative", peaks=((151.0031, 100), (178.9980, 20)))
            (Path(d) / "pos.msp").write_text(pos, encoding="utf-8")
            (Path(d) / "neg.msp").write_text(neg, encoding="utf-8")
            out = Path(d) / "run"
            out.mkdir()
            motif_rows, adducts = m.read_ruleset_tables(RULESET)
            args = SimpleNamespace(sources=["vs20"], out=out, workers=2, limit=None, ruleset_dir=Path(RULESET),
                                   console=Path(CONSOLE), seed=m.SEED, bootstrap=20)
            run = {}
            quiet(m.prepare, args, run, inputs={"vs20": [("pos", Path(d) / "pos.msp", None), ("neg", Path(d) / "neg.msp", None)]},
                  p2_run=None)
            self.assertEqual(run["prepared"]["vs20"]["representatives"], 5)
            quiet(m.fire, args, run)
            self.assertTrue(all(v["counts"]["base"]["errors"] == 0 for v in run["fire"].values()))
            quiet(m.measure, args, run)
            fires = list(p2.parquet_rows(out / "fires.parquet"))
            hit = [f for f in fires if f["rule_id"] == rule_id and f["pass_id"] == "base" and f["set"] == "representatives"]
            self.assertEqual(len(hit), 1)
            self.assertTrue(hit[0]["fired"])
            doc = json.loads((out / "rules_precision.json").read_text(encoding="utf-8"))
            g = next(r for r in doc["rules"] if r["rule_id"] == rule_id)
            self.assertEqual(g["status"], "not_estimable")
            self.assertEqual(g["passes"]["base"]["fires"], 1)
            self.assertEqual(run["acceptance_artifact_audit"], "pass")


if __name__ == "__main__":
    unittest.main()
