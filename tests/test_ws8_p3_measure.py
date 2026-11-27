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


HAS_ARROW = importlib.util.find_spec("pyarrow") is not None
C26_SCHEMA = [("ingest_lib", "s"), ("normalized_smiles", "s"), ("inchikey", "s"), ("inchikey14", "s"),
              ("molecular_formula", "s"), ("ionization_mode", "s"), ("instrument_type", "s"), ("adduct", "s"),
              ("adduct_orig", "s"), ("precursor_mz", "f"), ("precursor_error_ppm", "f"), ("ms2_mzs", "l"),
              ("ms2_normalized_intensities", "l"), ("num_peaks", "i"), ("base_peak_intensity", "f"),
              ("collision_energy_ev", "l"), ("collision_energy_orig", "s"), ("collision_energy_orig_units", "s")]


def c26_row(smiles=CAFFEINE, mode="positive", adduct="[M+H]+", mz=195.0877, ppm=1.5, lib="gnps", ik14="AAAAAAAAAAAAAA",
            peaks=((138.0662, 1.0), (110.0713, 0.4)), ce=(20.0,), inst="LC-ESI-QTOF"):
    """A synthetic CASMI26-shaped row (public molecules, made-up keys)."""
    return dict(ingest_lib=lib, normalized_smiles=smiles, inchikey=None, inchikey14=ik14, molecular_formula=None,
                ionization_mode=mode, instrument_type=inst, adduct=adduct, adduct_orig=adduct, precursor_mz=mz,
                precursor_error_ppm=ppm, ms2_mzs=[a for a, _ in peaks], ms2_normalized_intensities=[b for _, b in peaks],
                num_peaks=len(peaks), base_peak_intensity=1.0, collision_energy_ev=None if ce is None else list(ce),
                collision_energy_orig=None if ce is None else ",".join(f"{x:g}" for x in ce) + " eV",
                collision_energy_orig_units="eV")


def write_c26(path, rows):
    import pyarrow as pa
    import pyarrow.parquet as pq
    types = dict(s=pa.string(), f=pa.float64(), i=pa.int64(), l=pa.list_(pa.float64()))
    schema = pa.schema([(name, types[t]) for name, t in C26_SCHEMA])
    pq.write_table(pa.Table.from_pydict({n: [r[n] for r in rows] for n, _ in C26_SCHEMA}, schema=schema), path,
                   row_group_size=3)
    return path


C26_ROWS = [c26_row(ce=(10.0, 20.0, 40.0)),                                                    # 0 kept, merged CE
            c26_row(QUERCETIN, "negative", "[M-H]-", 301.0354, peaks=((151.0031, 1.0), (178.998, 0.2))),  # 1 neg
            c26_row(lib="enveda-180"),                                                         # 2 excluded library
            c26_row(smiles=None),                                                              # 3 no structure
            c26_row(ik14=None),                                                                # 4 no first block
            c26_row(ppm=None),                                                                 # 5 ppm missing
            c26_row(ppm=12.0),                                                                 # 6 ppm above 10
            c26_row(ppm=-9.9, inst=None),                                                      # 7 kept, no instrument
            c26_row(mode=None),                                                                # 8 mode missing
            c26_row(peaks=((138.0, 1.0), (110.0, 0.4)))]                                       # 9 nominal mass


class NamingTests(unittest.TestCase):
    def test_one_suffix_rule(self):
        self.assertEqual([m.suffix(*x) for x in (("vs20",), ("vs21",), ("c26",), ("vs20", "vs21"), ("vs20", "vs21", "c26"))],
                         ["", ".local-vs21", ".local-vs21.local-c26", ".local-vs21", ".local-vs21.local-c26"])
        self.assertEqual(m.row_path(Path("o"), "prepared", "c26"), Path("o") / "prepared.local-vs21.local-c26.parquet")
        self.assertEqual(m.row_path(Path("o"), "prepared", "vs21"), p2.row_path(Path("o"), "prepared", "vs21"))
        self.assertEqual(m.row_path(Path("o"), "prepared", "vs20"), p2.row_path(Path("o"), "prepared", "vs20"))
        self.assertEqual([m.local_vs21(x) for x in ("vs20", "vs21", "c26")], [False, True, True])

    def test_c26_inputs(self):
        self.assertEqual(m.library_inputs("c26", SimpleNamespace(c26_train=Path("x.parquet"))),
                         [("pos", Path("x.parquet"), m.C26_SHA256), ("neg", Path("x.parquet"), m.C26_SHA256)])
        self.assertEqual(m.C26_TRAIN.parent.parent.parent, m.OUT_ROOT.parent / "external")
        self.assertEqual(m.library_inputs("c26")[0][1], m.C26_TRAIN)


@unittest.skipUnless(HAS_ARROW, "pyarrow is required")
class C26ReaderTests(unittest.TestCase):
    def test_rows_become_msp_tuples_of_one_mode(self):
        with tempfile.TemporaryDirectory() as d:
            path = write_c26(Path(d) / "train.parquet", C26_ROWS)
            pos = list(m.iter_c26(path, "pos", batch_size=4))
            neg = list(m.iter_c26(path, "neg", batch_size=4))
        self.assertEqual([x[3]["ordinal"] for x in pos], [0, 2, 3, 4, 5, 6, 7, 8, 9])
        self.assertEqual([x[3]["ordinal"] for x in neg], [1])
        meta, peaks, decimals, prov = pos[0]
        self.assertEqual(meta, dict(precursormz="195.0877", precursortype="[M+H]+", ionmode="positive",
                                    instrumenttype="LC-ESI-QTOF", collisionenergy="10,20,40 eV", smiles=CAFFEINE))
        self.assertEqual((peaks, decimals), ([(138.0662, 1.0), (110.0713, 0.4)], [4, 4]))
        self.assertEqual({k: prov[k] for k in ("bucket", "in_silico", "exclude", "merged_ce")},
                         dict(bucket="gnps", in_silico=False, exclude=None, merged_ce=True))
        self.assertEqual([x[3]["exclude"] for x in pos],
                         [None, "excluded_ingest_lib", "missing_structure", "missing_structure",
                          "precursor_error_ppm_missing", "precursor_error_ppm_above_limit", None, None, None])
        self.assertEqual(pos[6][0]["instrumenttype"], "")
        self.assertEqual(pos[7][0]["ionmode"], "missing")
        self.assertIsNone(p2.ion_mode(pos[7][0]["ionmode"], "pos"))  # P2's mode check counts it
        self.assertTrue(p2.nominal_mass(pos[8][1], pos[8][2]))
        self.assertFalse(p2.nominal_mass(pos[0][1], pos[0][2]))
        self.assertEqual(neg[0][0]["ionmode"], "negative")
        self.assertFalse(neg[0][3]["merged_ce"])

    def test_decimals_from_shortest_repr(self):
        self.assertEqual([m.decimals_of(x) for x in (138.0, 138.0662, 1e-05, 2.5e20)], [1, 4, 0, 1])


@unittest.skipUnless(HAS_CHEM, "RDKit and pyarrow are required")
class C26PrepareTests(unittest.TestCase):
    def ruleset(self, d):
        t = Path(d) / "RuleSet" / "tables"
        t.mkdir(parents=True)
        (t / "motifs.json").write_text(json.dumps(dict(motifs=MOTIFS)), encoding="utf-8")
        (t / "adducts.json").write_text(json.dumps(dict(adducts=ADDUCTS)), encoding="utf-8")
        return t.parent

    def test_prepare_three_sources(self):
        rows = [c26_row(CAFFEINE),                                                             # 0 VS20 overlap
                c26_row(THEOBROMINE, mz=181.0720),                                             # 1 VS21 overlap
                c26_row(THEOPHYLLINE, mz=181.0720, ce=(10.0, 30.0)),                           # 2 kept, merged CE
                c26_row(THEOPHYLLINE, mz=181.0720, ppm=-9.9, inst=None, lib="mona"),           # 3 kept, same unit
                c26_row(QUERCETIN, "negative", "[M-H]-", 301.0354, peaks=((151.0031, 1.0), (178.998, 0.2))),  # 4 neg
                c26_row(PROTOPINE, mz=354.1336, peaks=((206.0812, 1.0), (149.0597, 0.6))),   # 5 kept
                c26_row(PROTOPINE, lib="enveda-180"), c26_row(smiles=None), c26_row(ppm=None), c26_row(ppm=11.0),
                c26_row(mode=None), c26_row(peaks=((138.0, 1.0), (110.0, 0.4)))]
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "vs20.msp").write_text(msp_record(CAFFEINE), encoding="utf-8")
            (Path(d) / "vs21.msp").write_text(msp_record(THEOBROMINE, mz="181.0720"), encoding="utf-8")
            (Path(d) / "empty.msp").write_text("", encoding="utf-8")
            c26 = write_c26(Path(d) / "train.parquet", rows)
            out = Path(d) / "run"
            out.mkdir()
            inputs = {"vs20": [("pos", Path(d) / "vs20.msp", None), ("neg", Path(d) / "empty.msp", None)],
                      "vs21": [("pos", Path(d) / "vs21.msp", None), ("neg", Path(d) / "empty.msp", None)],
                      "c26": [("pos", c26, None), ("neg", c26, None)]}
            args = SimpleNamespace(sources=["vs20", "vs21", "c26"], out=out, limit=None, workers=1,
                                   ruleset_dir=self.ruleset(d), c26_train=None)
            guarded, saved = [], m.vs21_guard
            m.vs21_guard = lambda a: guarded.append(a) or True
            try:
                run = {}
                quiet(m.prepare, args, run, inputs=inputs, p2_run=None)
            finally:
                m.vs21_guard = saved
            self.assertEqual(len(guarded), 2)  # vs21 and c26 are held out
            ex, p3x = run["exclusions"], run["p3_exclusions"]
            self.assertEqual((ex["c26:vs20_overlap_records"], ex["c26:vs21_overlap_records"]), (1, 1))
            self.assertEqual(ex.get("vs21:vs20_overlap_records", 0), 0)
            self.assertEqual((ex["c26:pos:input"], ex["c26:neg:input"]), (11, 1))
            self.assertEqual((ex["c26:pos:mode_mismatch"], ex["c26:pos:nominal_mass"]), (1, 1))
            for k, v in (("excluded_ingest_lib", 1), ("missing_structure", 1), ("precursor_error_ppm_missing", 1),
                         ("precursor_error_ppm_above_limit", 1), ("merged_collision_energy_records", 1),
                         ("bucket:gnps", 3), ("bucket:mona", 1)):
                self.assertEqual(p3x["c26:" + k], v, k)
            self.assertEqual(run["p2_parity"]["c26"], dict(result="not_checked", reason="no P2 run for c26"))
            self.assertEqual(run["prepared"]["c26"]["representatives"], 3)
            self.assertNotIn("instrument_types", run["prepared"]["c26"])
            dump = json.dumps(run)
            for secret in (CAFFEINE, THEOPHYLLINE, "AAAAAAAAAAAAAA", "QTOF"):
                self.assertNotIn(secret, dump)
            reps = list(p2.parquet_rows(out / "prepared.local-vs21.local-c26.parquet"))
            self.assertEqual(sorted(r["record_id"] for r in reps), [2 * 2, 2 * 4 + 1, 2 * 5])
            theo = next(r for r in reps if r["record_id"] == 4)
            self.assertEqual((theo["energy"], theo["spectrum_count"], theo["source"]), ("merged", 2, "c26"))
            self.assertTrue((out / "prepared.local-vs21.parquet").is_file())
            self.assertTrue((out / "prepared.parquet").is_file())
            names = sorted(p.name for p in (out / "msp" / "c26").iterdir())
            self.assertTrue(names and all(".local-vs21.local-c26." in n for n in names), names)
            itypes = m.read_rows(out / "instrument_types.local-vs21.local-c26.tsv")
            self.assertEqual([(r["instrument_type"], r["records"], r["instrument_class"]) for r in itypes],
                             [("", "1", "other"), ("LC-ESI-QTOF", "3", "QTOF-type")])
            self.assertFalse(list(out.glob("p3raw*")))


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


INS_ADDUCTS = [dict(adduct="[M+H]+", aliases=["M+H"], charge=1, multiplier=1, add="H", remove=None),
               dict(adduct="[M-H]-", aliases=["M-H"], charge=-1, multiplier=1, add=None, remove="H"),
               dict(adduct="[M+NH4]+", aliases=["M+NH4"], charge=1, multiplier=1, add="H4N", remove=None),
               dict(adduct="[M+HCOO]-", aliases=["[M+FA-H]-"], charge=-1, multiplier=1, add="CHO2", remove=None),
               dict(adduct="[M+H-H2O]+", aliases=[], charge=1, multiplier=1, add="H", remove="H2O")]
INS_HEADER = "parent_motif\tparent_adduct\tchange\tproduct_adduct\tsource\tnote\n"
INS_ROWS = ("*\t[M+H]+\t-H2O\t[M+H-H2O]+\tsynthetic\tloss\n"
            "*\t[M+H]+\t-NH3\t[M+H-NH3]+\tsynthetic\tloss needing N\n"
            "*\t[M+NH4]+\t-NH3\t[M+H]+\tsynthetic\tadduct conversion\n"
            "*\t[M+HCOO]-\t-CH2O2\t[M-H]-\tsynthetic\tadduct conversion\n"
            "gly\t[M-H]-\t-C6H10O5\t[M-H-Hex]-\tsynthetic\tmotif-restricted\n")


def elements(d):
    t = Path(d) / "tables"
    t.mkdir(parents=True, exist_ok=True)
    doc = {"electron": 0.000549, "elements": {"C": 12.0, "H": 1.007825032, "N": 14.003074, "O": 15.99491462},
           "format": "msfinder-diagnostic-elements", "format_version": 1}
    (t / "elements.json").write_text(json.dumps(doc, sort_keys=True), encoding="utf-8")
    return m.load_elements(d)


def transforms(d, rows=INS_ROWS):
    path = Path(d) / "insource-p3.tsv"
    path.write_text(INS_HEADER + rows, encoding="utf-8")
    return m.load_insource(path, {"gly", "M"}, INS_ADDUCTS, elements(d))


def parent(rid, formula, adduct, mz, peaks, motifs=(), mode=None):
    return dict(record_id=rid, key=f"K{rid}", formula=formula, motifs=list(motifs), adduct=adduct, precursor=mz,
                mode=mode or ("pos" if adduct.endswith("+") else "neg"), instrument_type="QTOF", collision_energy="20",
                mz=[a for a, _ in peaks], intensity=[b for _, b in peaks])


class InsourceTests(unittest.TestCase):
    def reps(self):
        return [parent(0, "C8H10N4O2", "[M+H]+", 195.0877, [(110.07, 40), (138.0662, 100), (177.08, 7), (180.0, 5)]),
                parent(2, "C6H12O6", "[M+H]+", 181.0707, [(85.03, 30), (163.06, 100)]),
                parent(4, "C6H12O6", "[M+NH4]+", 198.0972, [(85.03, 30), (163.06, 100), (181.07, 50)]),
                parent(1, "C21H20O11", "[M-H]-", 447.0933, [(284.03, 100), (285.04, 80)], motifs=["gly"]),
                parent(3, "C15H10O6", "[M-H]-", 285.0405, [(133.03, 100)]),
                parent(5, "C6H12O6", "[M+HCOO]-", 225.0616, [(179.056, 100), (89.02, 20)])]

    def test_rows_load_and_reject_bad_rows(self):
        with tempfile.TemporaryDirectory() as d:
            t = transforms(d)
            self.assertEqual([x["row"] for x in t], [0, 1, 2, 3, 4])
            self.assertAlmostEqual(t[0]["change_pda"] / 1e12, 18.010564684)
            self.assertEqual((t[0]["product_adduct"], t[0]["product_mapped"]), ("[M+H-H2O]+", True))
            self.assertEqual((t[1]["product_adduct"], t[1]["product_mapped"]), ("[M+H-NH3]+", False))
            self.assertEqual(t[3]["parent_adduct"], "[M+HCOO]-")
            for bad in ("nope\t[M+H]+\t-H2O\tx\ts\tn\n", "*\t[M+H]+\tH2O\tx\ts\tn\n",
                        "*\t[M+H]+\t-H2O\t[M-H]-\ts\tn\n", "*\t[M+Q]+\t-H2O\tx\ts\tn\n",
                        "*\t[M+H]+\t-HOH\tx\ts\tn\n", "*\t[M+H]+\t-SO3\tx\ts\tn\n"):  # repeated symbol; S not in the table
                with self.assertRaises(ValueError, msg=bad):
                    transforms(d, bad)

    def test_loss_row_truncates_and_adds_the_precursor_peak(self):
        with tempfile.TemporaryDirectory() as d:
            out, stats = m.build_insource(self.reps(), transforms(d), 7)
        r = next(x for x in out if x["record_id"] == 0 and x["row"] == 0)
        new = 195.0877 - 18.010564684
        self.assertAlmostEqual(r["precursor"], new)
        self.assertEqual(r["adduct"], "[M+H-H2O]+")
        # 180.0 lies above new + 0.01 and goes; 177.08 lies within 0.01 of the new precursor and is replaced.
        self.assertEqual([round(x, 4) for x in r["mz"]], [110.07, 138.0662, round(new, 4)])
        self.assertEqual(r["intensity"][-1], 100)
        self.assertEqual(r["mz"], sorted(r["mz"]))

    def test_atom_containment_including_adduct_atoms(self):
        with tempfile.TemporaryDirectory() as d:
            out, stats = m.build_insource(self.reps(), transforms(d), 7)
        got = {(x["record_id"], x["row"]) for x in out}
        self.assertIn((0, 1), got)            # caffeine [M+H]+ holds NH3
        self.assertNotIn((2, 1), got)         # hexose [M+H]+ has no N
        self.assertEqual(stats["by_row"]["1"]["atoms_missing"], 1)
        self.assertIn((4, 2), got)            # the NH4 of [M+NH4]+ supplies the N of -NH3
        self.assertIn((5, 3), got)            # formate supplies C and O of -CH2O2
        self.assertTrue(m.ion_contains("C6H12O6", dict(multiplier=1, add={"H": 1}, remove={}, change={"H": 2, "O": 1})))
        self.assertFalse(m.ion_contains("CH4", dict(multiplier=1, add={}, remove={"H": 1}, change={"H": 2, "O": 1})))
        self.assertTrue(m.ion_contains("C20H18NO4+", dict(multiplier=1, add={}, remove={}, change={"C": 1, "O": 2})))

    def test_mode_adduct_motif_and_star_rows(self):
        with tempfile.TemporaryDirectory() as d:
            out, stats = m.build_insource(self.reps(), transforms(d), 7)
        got = {(x["record_id"], x["row"]) for x in out}
        self.assertEqual({k for k in got if k[1] == 0}, {(0, 0), (2, 0)})  # [M+H]+ only, not [M+NH4]+ or negatives
        self.assertEqual({k for k in got if k[1] == 4}, {(1, 4)})          # the gly-restricted row skips record 3
        self.assertEqual(next(x for x in out if x["row"] == 3)["adduct"], "[M-H]-")
        self.assertEqual(stats["by_row"]["1"]["adduct_unmapped"], 1)
        self.assertEqual(stats["records"], len(out))

    def test_no_peaks_left(self):
        with tempfile.TemporaryDirectory() as d:
            reps = [parent(0, "C8H10N4O2", "[M+H]+", 195.0877, [(190.0, 100)])]
            out, stats = m.build_insource(reps, transforms(d), 7)
        self.assertEqual(stats["by_row"]["0"]["no_peaks_left"], 1)
        self.assertNotIn((0, 0), {(x["record_id"], x["row"]) for x in out})

    def test_deterministic_and_capped_by_seed(self):
        with tempfile.TemporaryDirectory() as d:
            t = transforms(d)
            a, _ = m.build_insource(self.reps(), t, 7)
            b, _ = m.build_insource(list(reversed(self.reps())), t, 7)
            self.assertEqual(a, b)
            c1, s1 = m.build_insource(self.reps(), t, 7, cap=3)
            c2, _ = m.build_insource(self.reps(), t, 7, cap=3)
            self.assertEqual((len(c1), s1["capped"], s1["before_cap"]), (3, True, len(a)))
            self.assertEqual(c1, c2)
            with tempfile.TemporaryDirectory() as o1, tempfile.TemporaryDirectory() as o2:
                m.write_insource(o1, "vs20", a, 2)
                m.write_insource(o2, "vs20", b, 2)
                for name in ("insource-00.msp", "insource-01.msp", "insource-index.tsv"):
                    self.assertEqual((Path(o1) / "msp/vs20" / name).read_bytes(), (Path(o2) / "msp/vs20" / name).read_bytes())
                peaks = [x for x in (Path(o1) / "msp/vs20/insource-00.msp").read_text(encoding="utf-8").splitlines()
                         if x and x[0].isdigit()]
                self.assertTrue(all(PAIR.match(x) for x in peaks))

    def test_truth_stays_the_parents(self):
        """A recast fires for the rule's motif: a TP only when the parent has the motif; otherwise an FP."""
        truth, by_rule, chim = scenario()
        rows = [dict(row=0, mode="pos", change_text="-H2O", product_adduct="[M+H-H2O]+"),
                dict(row=1, mode="neg", change_text="-CO2", product_adduct="[M-H-CO2]-")]
        pos_m = [u for u in sorted(truth) if "M" in truth[u]["motifs"]][:10]
        neg_m = [u for u in sorted(truth) if not truth[u]["motifs"]][:3]
        ents = [((u[0], u[1], 0), True, 1, False, None) for u in pos_m + neg_m]
        ins = dict(by_rule={0: {"base": {"r1": ents}}}, universe={0: [e[0] for e in ents] + [("vs20", 99999 * 2, 0)]},
                   notev={0: {"base": {("vs20", 99999 * 2, 0)}}}, rows=rows)
        truth[("vs20", 99999 * 2)] = dict(truth[neg_m[0]])
        res = m.measure_core([rule()], truth, by_rule, {}, ["vs20"], replicates=10, chimeric=chim, insource=ins)
        r0 = next(r for r in res["synthetic"] if r["set"] == "insource" and r["transform"] == 0)
        self.assertEqual((r0["fires"], r0["tp"], r0["fp"], r0["c_fp"], r0["n_not_evaluated"], r0["n_universe"]),
                         (13, 10, 3, 3, 1, 14))
        self.assertEqual((r0["change"], r0["product_adduct"]), ("-H2O", "[M+H-H2O]+"))
        self.assertFalse(any(r["set"] == "insource" and r["transform"] == 1 for r in res["synthetic"]))  # other mode
        g = res["gates"]["r1"]
        self.assertIn("synthetic_insource", g["failed"])
        self.assertEqual(g["criteria"]["synthetic_insource"]["unit"], "cluster")
        ins["by_rule"] = {0: {"base": {"r1": [(e[0], False, 1, False, None) for e in ents]}}}
        g = m.measure_core([rule()], truth, by_rule, {}, ["vs20"], replicates=10, chimeric=chim, insource=ins)["gates"]["r1"]
        self.assertTrue(g["criteria"]["synthetic_insource"]["ok"])  # no fires: passes, no longer pending
        self.assertNotIn("synthetic_insource", g["pending"])

    @unittest.skipUnless(RULESET, "WS8_P3_RULESET is not set")
    def test_shipped_table_loads_against_the_rule_set(self):
        motif_rows, adducts = m.read_ruleset_tables(RULESET)
        t = m.load_insource(m.INSOURCE_TABLE, {r["motif_id"] for r in motif_rows}, adducts, m.load_elements(RULESET))
        self.assertEqual(len(t), len(m.read_rows(m.INSOURCE_TABLE)))


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
        res = self.core(sc=scenario(tp_clusters=160, members=1, decoy=9))
        r = self.gate_row(res)
        self.assertEqual((r["decoy_fires_plus_5"], r["c_decoy_fires_plus_5"]), (9, 9))
        self.assertAlmostEqual(r["decoy_fdr"], 9 / 160)
        self.assertEqual(res["gates"]["r1"]["status"], "fail")
        self.assertIn("decoy_fdr", res["gates"]["r1"]["failed"])

    def test_decoy_fdr_exactly_at_the_bound_meets_it(self):
        """Regression: 6/120 per shift is exactly 0.05, which the bound allows; a float mean must not fail it."""
        self.assertTrue(m.decoy_fdr_within([6, 6, 6], 120, 0.05))
        self.assertTrue(m.decoy_fdr_within([5, 6, 7], 120, 0.05))
        self.assertFalse(m.decoy_fdr_within([6, 6, 7], 120, 0.05))
        self.assertTrue(m.decoy_fdr_within([0, 0, 0], 1, 0.05))

    def test_gate_counts_decoy_fires_in_clusters(self):
        """Regression: the gate reads the cluster-level decoy FDR (PR 908); the record-level value is report-only."""
        res = self.core(sc=scenario(tp_clusters=160, members=4, decoy=36))  # 36 records = 9 clusters
        r = self.gate_row(res)
        self.assertAlmostEqual(r["decoy_fdr"], 36 / 640)
        self.assertAlmostEqual(r["c_decoy_fdr"], 9 / 160)
        crit = res["gates"]["r1"]["criteria"]["decoy_fdr"]
        self.assertEqual((crit["value"], crit["unit"], crit["record_level"]), (r["c_decoy_fdr"], "cluster", r["decoy_fdr"]))
        res = self.core(sc=scenario(tp_clusters=160, members=4, decoy=28))  # 28 records = 7 clusters: 0.0438
        self.assertGreater(self.gate_row(res)["decoy_fdr"], 0.0)
        self.assertNotIn("decoy_fdr", res["gates"]["r1"]["failed"])

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

    def test_gate_counts_chimeric_precision_in_clusters(self):
        """Regression: 20 false chimeric fires from one cluster fail record-level precision (report only) but are one
        cluster against 160 at cluster level, which the gate reads."""
        truth, by_rule, chim = scenario()
        negs = [u for u in sorted(truth) if not truth[u]["motifs"]][:20]
        for u in negs:
            truth[u]["cluster"] = "n_shared"
        for lv in m.CHIMERIC_LEVELS:
            chim["by_rule"][lv]["base"]["r1"] += [((u[0], u[1], lv), True, 1, False, None) for u in negs]
        res = m.measure_core([rule()], truth, by_rule, {}, ["vs20"], replicates=10, chimeric=chim)
        srow = [r for r in res["synthetic"] if r["level"] == 0.5][0]
        self.assertLess(srow["precision"], 0.95)
        self.assertEqual((srow["c_fires"], srow["c_fp"]), (161, 1))
        g = res["gates"]["r1"]
        self.assertNotIn("synthetic_chimeric", g["failed"])
        self.assertEqual(g["criteria"]["synthetic_chimeric"]["unit"], "cluster")

    def test_chimeric_not_evaluated_per_level(self):
        """Regression: chimeric records a pass did not evaluate are carried into rules_synthetic.csv."""
        truth, by_rule, chim = scenario()
        chim["notev"] = {0.3: {"base": {("vs20", 0, 0.3), ("vs20", 2, 0.3)}}}
        res = m.measure_core([rule()], truth, by_rule, {}, ["vs20"], replicates=10, chimeric=chim)
        self.assertEqual({r["level"]: r["n_not_evaluated"] for r in res["synthetic"]}, {0.1: 0, 0.3: 2, 0.5: 0})
        self.assertIn("n_not_evaluated", m.SYNTHETIC_FIELDS)

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


def no_fire_insource(truth):
    """An in-source set with recast records in the positive mode and no fires, so the criterion is measured."""
    uids = [(u[0], u[1], 0) for u in sorted(truth)[:20]]
    return dict(by_rule={0: {"base": {"r1": []}}}, universe={0: uids}, notev={},
                rows=[dict(row=0, mode="pos", change_text="-H2O", product_adduct="[M+H-H2O]+")])


class TierTests(unittest.TestCase):
    """Tier 1 (strict) and tier 2 (provisional, registered 2026-11-24). With VS21 in the sources the held-out
    criterion of an as_published rule is met (the gate is then on VS20 and VS21 pooled), so `pass` is reachable."""

    def gate(self, sc, sources=("vs20", "vs21"), rules=None):
        truth, by_rule, chim = sc
        res = m.measure_core(rules or [rule()], truth, by_rule, {}, list(sources), replicates=10, chimeric=chim,
                             insource=no_fire_insource(truth))
        return res, res["gates"]["r1"]

    def test_rule_passing_tier_1(self):
        res, g = self.gate(scenario())
        self.assertEqual((g["status"], g["status_provisional"], g["label"]), ("pass", "pass", "pass_strict"))
        self.assertEqual((g["failed"], g["pending"], g["provisional"]["failed"]), ([], [], []))
        with tempfile.TemporaryDirectory() as d:
            m.write_aggregates(d, res, "p3cs-20261124-abcdef1", [rule()], "")
            prop = json.loads((Path(d) / "review_proposal-p3cs-20261124-abcdef1.json").read_text(encoding="utf-8"))
            self.assertEqual([x["rule_id"] for x in prop["passed_strict"]], ["r1"])
            self.assertEqual(prop["passed_strict"][0]["criteria"]["n_clusters"]["threshold"], 100)
            self.assertEqual(prop["passed_provisional_only"], [])
            row = next(r for r in m.read_rows(Path(d) / "rules_precision.csv") if r["gate_row"] == "True")
            self.assertEqual((row["status_strict"], row["status_provisional"], row["label"]), ("pass", "pass", "pass_strict"))

    def test_rule_passing_tier_2_only(self):
        res, g = self.gate(scenario(tp_clusters=40, members=1))  # n = 40 clusters, 0 FP
        self.assertEqual((g["status"], g["status_provisional"], g["label"]),
                         ("not_estimable", "pass", "pass_provisional_only"))
        self.assertGreaterEqual(g["provisional"]["criteria"]["cluster_wilson_lo"]["value"], 0.90)
        with tempfile.TemporaryDirectory() as d:
            m.write_aggregates(d, res, "p3cs-20261124-abcdef1", [rule()], "")
            prop = json.loads((Path(d) / "review_proposal-p3cs-20261124-abcdef1.json").read_text(encoding="utf-8"))
            self.assertEqual(prop["passed_strict"], [])
            listed = prop["passed_provisional_only"][0]
            self.assertEqual((listed["rule_id"], listed["tier"]), ("r1", "provisional"))
            self.assertEqual(listed["criteria"]["n_clusters"]["threshold"], 30)
            self.assertEqual(listed["criteria"]["cluster_wilson_lo"]["threshold"], 0.90)

    def test_rule_failing_both_tiers(self):
        res, g = self.gate(scenario(tp_clusters=160, fp_clusters=20))
        self.assertEqual((g["status"], g["status_provisional"], g["label"]), ("fail", "fail", "fail"))
        self.assertIn("cluster_wilson_lo", g["provisional"]["failed"])

    def test_prior_stress_is_report_only_in_tier_2(self):
        # Few negatives: prior stress falls below 0.90, which fails tier 1 but not tier 2.
        res, g = self.gate(scenario(tp_clusters=160, neg_clusters=5))
        self.assertIn("prior_stress", g["failed"])
        p = g["provisional"]["criteria"]["prior_stress"]
        self.assertTrue(p["report_only"])
        self.assertFalse(p["ok"])
        self.assertNotIn("prior_stress", g["provisional"]["failed"])
        self.assertEqual((g["status_provisional"], g["label"]), ("pass", "pass_provisional_only"))

    def test_incomplete_provisional_in_a_vs20_only_run(self):
        res, g = self.gate(scenario(tp_clusters=40, members=1), sources=("vs20",))
        self.assertEqual((g["status"], g["status_provisional"], g["label"]),
                         ("not_estimable", "incomplete", "incomplete_provisional"))
        self.assertEqual(g["provisional"]["pending"], ["heldout"])

    def test_tuned_rules(self):
        res, g = self.gate(scenario(), rules=[rule(origin="tuned_vs20")])
        self.assertEqual(g["label"], "incomplete")  # no VS21 fires: VS21-alone criterion pending in both tiers
        self.assertIn("heldout", g["provisional"]["pending"])
        res, g = self.gate(scenario(), rules=[rule(origin="tuned_pooled")])
        self.assertEqual((g["status"], g["status_provisional"]), ("incomplete", "incomplete"))

    def test_c26_is_a_third_source_of_the_per_source_criterion(self):
        truth, by_rule, chim = scenario(tp_clusters=160, fp_clusters=6, members=1)
        for rec in truth.values():
            c = rec["cluster"]
            if c.startswith("fp") or c.startswith("tp") and int(c[2:]) < 30:
                rec.update(source="c26", bucket="gnps")
            elif c.startswith("tp") and int(c[2:]) < 40:
                rec["source"] = "vs21"
        res, g = self.gate((truth, by_rule, chim), sources=("vs20", "vs21", "c26"))
        ps = g["criteria"]["per_source_precision"]["value"]
        self.assertEqual(set(ps), {"vs20", "vs21", "c26"})
        self.assertEqual((ps["c26"]["c_fires"], ps["c26"]["c_precision"]), (36, 30 / 36))
        self.assertIn("per_source_precision", g["failed"])
        self.assertIn("per_source_precision", g["provisional"]["failed"])
        self.assertIn("C26", g["criteria"]["heldout"]["note"])
        self.assertIn("c26:gnps", g["by_provenance"])
        self.assertEqual({r["source"] for r in res["precision"]}, {"pooled", "vs20", "vs21", "c26"})


def c26_split(sc, n_c26):
    """The first n_c26 TP clusters move to C26; everything else stays VS20."""
    truth, by_rule, chim = sc
    for rec in truth.values():
        c = rec["cluster"]
        if c.startswith("tp") and int(c[2:]) < n_c26:
            rec.update(source="c26", bucket="gnps")
    return truth, by_rule, chim


class ExploratoryTierTests(unittest.TestCase):
    """Tier 3 (exploratory, registered 2026-11-27): judged on C26 alone, VS20 and VS21 only guard."""
    gate = TierTests.gate

    def gate3(self, sc, origin="as_published"):
        return self.gate(sc, sources=("vs20", "vs21", "c26"), rules=[rule(origin=origin)])

    def test_not_applicable_without_c26(self):
        res, g = self.gate(scenario())
        self.assertEqual((g["status_exploratory"], g["label"]), ("not_applicable", "pass_strict"))

    def test_rule_passing_tier_3_only(self):
        # 40 TP + 6 FP clusters pooled (precision 0.87) fails tier 2; C26 holds 25 TP clusters and no FP.
        for origin in ("as_published", "tuned_pooled"):
            res, g = self.gate3(c26_split(scenario(tp_clusters=40, fp_clusters=6, members=1), 25), origin)
            self.assertEqual((g["status"], g["status_provisional"]), ("not_estimable", "fail"), origin)
            self.assertEqual((g["status_exploratory"], g["exploratory"]["failed"]), ("pass", []), origin)
            self.assertEqual(g["exploratory"]["criteria"]["n_clusters"]["value"], 25)
            self.assertIn("C26 alone", g["exploratory"]["criteria"]["heldout"]["note"])
            self.assertEqual(g["label"], "pass_exploratory_only", origin)
        with tempfile.TemporaryDirectory() as d:
            m.write_aggregates(d, res, "p3cs-20261127-abcdef1", [rule()], "")
            prop = json.loads((Path(d) / "review_proposal-p3cs-20261127-abcdef1.json").read_text(encoding="utf-8"))
            listed = prop["passed_exploratory_only"][0]
            self.assertEqual((listed["tier"], listed["criteria"]["n_clusters"]["threshold"]), ("exploratory", 20))
            row = next(r for r in m.read_rows(Path(d) / "rules_precision.csv") if r["gate_row"] == "True")
            self.assertEqual((row["status_exploratory"], row["label"]), ("pass", "pass_exploratory_only"))

    def test_tier_3_needs_20_c26_clusters(self):
        res, g = self.gate3(c26_split(scenario(tp_clusters=40, fp_clusters=6, members=1), 15))
        self.assertEqual((g["status_exploratory"], g["label"]), ("not_estimable", "not_estimable"))

    def test_tier_3_guard_on_seen_sets(self):
        # VS20 keeps 35 TP and 10 FP clusters: 45 fired, point precision 0.78 < 0.85.
        res, g = self.gate3(c26_split(scenario(tp_clusters=60, fp_clusters=10, members=1), 25))
        self.assertEqual(g["status_exploratory"], "fail")
        self.assertEqual(g["exploratory"]["failed"], ["per_source_precision"])
        self.assertEqual(set(g["exploratory"]["criteria"]["per_source_precision"]["value"]), {"vs20", "vs21"})

    def test_small_c26_gate_cell_is_blanked(self):
        res, g = self.gate3(c26_split(scenario(tp_clusters=40, fp_clusters=6, members=1), 3))
        self.assertEqual(g["exploratory"]["criteria"]["n_clusters"]["value"], 3)
        with tempfile.TemporaryDirectory() as d:
            m.write_aggregates(d, res, "p3cs-20261127-abcdef1", [rule()], "")
            doc = json.loads((Path(d) / "rules_precision.json").read_text(encoding="utf-8"))
            crit = doc["rules"][0]["exploratory"]["criteria"]
            self.assertEqual([crit[k]["value"] for k in ("n_clusters", "cluster_wilson_lo", "prior_stress", "decoy_fdr")],
                             [None] * 4)
            self.assertEqual(doc["rules"][0]["criteria"]["n_clusters"]["value"], 46)  # pooled kept

    def test_combined_label_order(self):
        cl = m.combined_label
        self.assertEqual(cl("pass", "pass", "pass"), "pass_strict")
        self.assertEqual(cl("fail", "pass", "pass"), "pass_provisional_only")
        self.assertEqual(cl("fail", "incomplete", "pass"), "pass_exploratory_only")
        self.assertEqual(cl("fail", "incomplete", "incomplete"), "incomplete_provisional")
        self.assertEqual(cl("not_estimable", "fail", "incomplete"), "incomplete_exploratory")
        self.assertEqual(cl("fail", "fail", "not_applicable"), "fail")

    def test_combined_label(self):
        self.assertEqual(m.combined_label("pass", "pass"), "pass_strict")
        self.assertEqual(m.combined_label("fail", "pass"), "pass_provisional_only")
        self.assertEqual(m.combined_label("not_estimable", "incomplete"), "incomplete_provisional")
        self.assertEqual(m.combined_label("incomplete", "incomplete"), "incomplete")
        self.assertEqual(m.combined_label("fail", "fail"), "fail")


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
        sfx = m.suffix("vs20", "vs21") if local else ""
        doc = m.write_aggregates(d, res, "p3cs-20261123-abcdef1", rules, sfx)
        m.write_report(d, doc, dict(prepared=dict(vs20=dict(provenance_unavailable=True)), p3_exclusions={}), sfx)
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
            self.assertIn("In-source recasts", report)
            self.assertIn("Tier 2 (provisional, registered 2026-11-24", report)
            prop = json.loads((Path(a) / "review_proposal-p3cs-20261123-abcdef1.json").read_text(encoding="utf-8"))
            self.assertEqual((prop["passed_strict"], prop["passed_provisional_only"]), ([], []))
            self.assertIn("author chooses which tier", prop["note"])

    def test_audit_catches_forbidden_columns_and_keys(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "x.csv").write_text("rule_id,first_block\n", encoding="utf-8")
            with self.assertRaises(AssertionError):
                m.audit_aggregates(d)
            (Path(d) / "x.csv").write_text("rule_id\nBSYNRYMUTXBXSQ-UHFFFAOYSA-N\n", encoding="utf-8")
            with self.assertRaises(AssertionError):
                m.audit_aggregates(d)

    def held_out_outputs(self, d):
        """Three sources: VS21 holds 3 fired clusters (6 fired records), C26 12 fired clusters with 2 FP clusters, VS20
        the rest, with one fired VS20 record in its own (Orbitrap) instrument stratum."""
        truth, by_rule, chim = scenario(fp_clusters=2)
        for rec in truth.values():
            c = rec["cluster"]
            if c in ("tp0", "tp1", "tp2") or c.startswith("n") and int(c[1:]) % 3 == 1:
                rec["source"] = "vs21"
            elif c in ("fp0", "fp1") or c.startswith("tp") and 3 <= int(c[2:]) < 13 or c.startswith("n") and int(c[1:]) % 3 == 2:
                rec["source"] = "c26"
        truth[next(u for u in sorted(truth) if truth[u]["cluster"] == "tp20")]["instrument"] = "Orbitrap/FT"
        rules = [rule(), rule("r2", measure=False)]
        res = m.measure_core(rules, truth, by_rule, {}, ["vs20", "vs21", "c26"], replicates=20, chimeric=chim)
        doc = m.write_aggregates(d, res, "p3cs-20261123-abcdef1", rules, m.suffix("vs20", "vs21", "c26"),
                                 extra=dict(sources=["vs20", "vs21", "c26"]))
        return res, doc

    def test_held_out_outputs_named_local_and_suppressed_by_denominator(self):
        with tempfile.TemporaryDirectory() as d:
            self.held_out_outputs(d)
            names = sorted(p.name for p in Path(d).iterdir())
            self.assertTrue(all(".local-vs21.local-c26." in n for n in names), names)
            rows = m.read_rows(Path(d) / "rules_precision.local-vs21.local-c26.csv")
            cell = {(r["source"], r["instrument"]): r for r in rows if (r["pass_id"], r["energy"]) == ("base", "all")}
            vs21, c26 = cell[("vs21", "all")], cell[("c26", "all")]
            # VS21: 3 fired clusters, a denominator below 5, so the cluster set goes; 6 fired records stay.
            self.assertEqual((vs21["c_fires"], vs21["c_tp"], vs21["c_fp"], vs21["c_precision"], vs21["c_wilson_lo"]),
                             ("", "", "", "", ""))
            self.assertEqual((vs21["fires"], vs21["tp"], vs21["fp"], vs21["precision"]), ("6", "6", "0", "1.0"))
            self.assertEqual(vs21["lr_plus"], "")
            self.assertNotEqual(vs21["n_universe"], "")
            # C26: 12 fired clusters; the 2 FP clusters are a numerator, not a denominator, and stay.
            self.assertEqual((c26["c_fires"], c26["c_fp"]), ("12", "2"))
            self.assertNotEqual(c26["c_precision"], "")
            # Development and pooled cells are never blanked, however small.
            for key in (("vs20", "all"), ("pooled", "all"), ("vs20", "Orbitrap/FT"), ("pooled", "Orbitrap/FT")):
                self.assertNotIn("", (cell[key]["fires"], cell[key]["c_fires"], cell[key]["c_precision"]), key)
            self.assertEqual(cell[("vs20", "Orbitrap/FT")]["fires"], "1")
            doc = json.loads((Path(d) / "rules_precision.local-vs21.local-c26.json").read_text(encoding="utf-8"))
            crit = doc["rules"][0]["criteria"]
            self.assertIsNotNone(crit["cluster_wilson_lo"]["value"])
            ps = crit["per_source_precision"]["value"]
            self.assertEqual(ps["vs21"], dict(c_fires=None, c_precision=None))
            self.assertEqual(ps["c26"]["c_fires"], 12)
            self.assertGreater(ps["vs20"]["c_fires"], 100)

    def test_suppress(self):
        row = dict(rule_id="r", version=1, source="vs21", fires=3, tp=3, fp=0, precision=1.0, n_positive=40,
                   n_universe=900, coverage=0.5, lr_plus=4.0)
        s = m.suppress(row)
        self.assertEqual((s["fires"], s["tp"], s["fp"], s["precision"], s["coverage"], s["lr_plus"]),
                         (None, None, None, None, 0.5, None))
        self.assertEqual((s["version"], s["n_universe"], s["n_positive"]), (1, 900, 40))
        for source in ("vs20", "pooled", None):
            self.assertEqual(m.suppress(dict(row, source=source)), dict(row, source=source), source)
        self.assertIsNone(m.suppress(dict(row, source="c26"))["precision"])
        self.assertIsNone(m.suppress({k: v for k, v in row.items() if k != "source"}, held_out=True)["fires"])

    def test_suppress_by_denominator_only(self):
        """Regression: a held-out cell is blanked only where its own denominator is below 5; a small numerator
        (fp = 3 of 200 fires) or a small count in another set does not blank the row."""
        row = dict(rule_id="r", version=1, source="vs21", fires=200, tp=197, fp=3, precision=0.985, c_fires=162, c_tp=160,
                   c_fp=2, c_precision=0.99, c_wilson_lo=0.96, c_coverage=0.8, n_universe=900, n_positive=4,
                   coverage=0.5, recall_universe=0.2, lr_plus=9.0,
                   decoy_fires_plus_5=2, decoy_fires_plus_9=0, decoy_fires_minus_7=0, decoy_fdr=0.003,
                   c_decoy_fires_plus_5=0, c_decoy_fires_plus_9=0, c_decoy_fires_minus_7=0, c_decoy_fdr=0.0)
        s = m.suppress(row)
        for k in ("n_positive", "coverage", "recall_universe", "lr_plus"):
            self.assertIsNone(s[k], k)
        for k in ("fires", "tp", "fp", "precision", "c_fires", "c_tp", "c_fp", "c_precision", "c_wilson_lo",
                  "c_coverage", "n_universe", "decoy_fires_plus_5", "decoy_fdr", "c_decoy_fdr"):
            self.assertEqual(s[k], row[k], k)
        # A small conditional population blanks the fired records it bounds, and their rates.
        s = m.suppress(dict(source="c26", n_cond=4, n_cond_positive=4, fires=4, tp=4, fp=0, precision=1.0, recall_cond=1.0))
        self.assertEqual(set(s.values()) - {"c26"}, {None})
        occ = m.suppress(dict(rule_id="r", source="vs21", n_positive=3, n_negative=897, n_negative_cond=40,
                              record_occupancy=0.04, c_positive=60, c_negative=4, c_negative_cond=1, c_record_occupancy=0.25))
        self.assertEqual((occ["n_positive"], occ["n_negative"], occ["record_occupancy"]), (None, 897, 0.04))
        self.assertEqual((occ["c_negative"], occ["c_negative_cond"], occ["c_record_occupancy"]), (None, None, None))

    def test_suppress_report_counts(self):
        d = {"vs21:pos:input": 1000, "vs21:pos:nominal_mass": 3, "vs21:pos:no_peaks": 20, "vs21:x": 0,
             "vs20:pos:nominal_mass": 2, "c26:bucket:x": 4, "c26/representatives/base": 1}
        s = m.suppress_counts(d)
        self.assertEqual(s, {"vs21:pos:input": 1000, "vs21:pos:nominal_mass": None, "vs21:pos:no_peaks": 20, "vs21:x": 0,
                             "vs20:pos:nominal_mass": 2, "c26:bucket:x": None, "c26/representatives/base": None})
        self.assertEqual(m.suppress_counts({"a": 2, "b": 3, "c": 50}), {"a": 2, "b": 3, "c": 50})
        self.assertEqual(m.suppress_counts({"a": 2, "b": 3, "c": 50}, held_out=True), {"a": None, "b": None, "c": 50})

    def test_vs21_report_suppresses_run_counts(self):
        with tempfile.TemporaryDirectory() as d:
            doc = self.outputs(d, local=True)
            run = dict(prepared=dict(vs21=dict(provenance_unavailable=False)),
                       exclusions={"vs21:pos:input": 5000, "vs21:pos:nominal_mass": 2},
                       p3_exclusions={"vs21:bucket:gnps": 4000, "vs21:bucket:riken": 3, "vs21:in_silico_records": 0},
                       instrument_classes={"QTOF-type -> QTOF-type": 900, "other -> unknown": 1},
                       chimeric={"vs21": dict(no_partner=4, no_base=0, records=3000)},
                       fire={"vs21/representatives/r-00": dict(counts={"base": dict(skipped=2, errors=0, skipped_by_reason={"skipped:profile": 2})})})
            m.write_report(d, doc, run, ".local-vs21")
            text = (Path(d) / "report.local-vs21.md").read_text(encoding="utf-8")
            self.assertIn('"vs21:pos:nominal_mass": null', text)
            self.assertIn('"riken": null', text)
            self.assertIn('"other -> unknown": 1', text)  # pooled over sources: never blanked
            self.assertIn('"no_partner": null', text)
            self.assertIn('"skipped": null', text)
            self.assertIn("in-silico detector counts suppressed", text)

    def test_report_lists_skip_counts(self):
        fire = {"vs20/representatives/r-00": dict(counts={"base": dict(skipped=2, errors=0, skipped_by_reason={"skipped:profile": 2}),
                                                          "mz_only_curated": dict(skipped=5, errors=0, skipped_by_reason={"skipped:adduct_unmapped": 5})}),
                "vs20/representatives/r-01": dict(counts={"base": dict(skipped=1, errors=0, skipped_by_reason={"skipped:profile": 1})})}
        sk = m.skip_counts(fire)
        self.assertEqual(sk["vs20/representatives/base"], {"errors": 0, "skipped": 3, "skipped:skipped:profile": 3})
        self.assertEqual(sk["vs20/representatives/mz_only_curated"]["skipped"], 5)
        with tempfile.TemporaryDirectory() as d:
            doc = self.outputs(d)
            m.write_report(d, doc, dict(fire=fire), "")
            self.assertIn("vs20/representatives/mz_only_curated", (Path(d) / "report.md").read_text(encoding="utf-8"))


def prepared_row(rid, source, cluster, motifs, mode="pos"):
    return dict(record_id=rid, source=source, key=f"K{rid}", formula="C1", motifs=sorted(motifs), cluster_id=cluster,
                mode=mode, adduct="[M+H]+", adduct_mapped=True, precursor=200.0, instrument="QTOF-type", energy="eV:<20",
                instrument_type="", collision_energy="", bucket="x", spectrum_count=1, msp_file="a.msp", msp_index=0,
                mz=[100.0], intensity=[1.0])


class InventoryTests(unittest.TestCase):
    def test_power_counts_by_truth_only(self):
        truth, _, _ = scenario(tp_clusters=10, members=2, neg_clusters=20)
        for uid in sorted(truth)[:6]:  # tp0..tp2, 6 motif records in 3 clusters
            truth[uid]["source"] = "c26"
        rules = [rule(), rule("r2", measure=False), rule("r3", source_checked=False), rule("r4", ion_mode="neg")]
        rows = m.power_counts(rules, truth, ["vs20", "vs21", "c26"])
        self.assertEqual({r["rule_id"] for r in rows}, {"r1", "r4"})
        by = {(r["rule_id"], r["source"]): r for r in rows}
        self.assertEqual({k: by[("r1", "c26")][k] for k in ("n_universe", "n_positive", "c_universe", "c_positive")},
                         dict(n_universe=6, n_positive=6, c_universe=3, c_positive=3))
        self.assertEqual((by[("r1", "pooled")]["n_positive"], by[("r1", "pooled")]["c_positive"]), (20 + 30, 10 + 30))
        self.assertEqual((by[("r1", "vs21")]["n_universe"], by[("r4", "pooled")]["n_universe"]), (0, 0))  # no neg records
        self.assertEqual(list(rows[0]), m.POWER_FIELDS)

    @unittest.skipUnless(HAS_ARROW, "pyarrow is required")
    def test_inventory_stage_needs_no_console_or_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "p3cs-20261127-abcdef1"
            out.mkdir()
            schema = m._schemas()[1]
            data = {"vs20": [prepared_row(0, "vs20", "a", {"M"}), prepared_row(2, "vs20", "b", set())],
                    "vs21": [prepared_row(1, "vs21", "c", {"M"})],
                    "c26": [prepared_row(4, "c26", "d", {"M"}), prepared_row(6, "c26", "d", {"M"}),
                            prepared_row(3, "c26", "e", set(), mode="neg")]}
            for source, rows in data.items():
                sink = p2.ParquetRows(m.row_path(out, "prepared", source), schema)
                for r in rows:
                    sink.add(r)
                sink.close()
            m.write_instrument_types(out, "c26", Counter({"LC-ESI-QTOF": 2, "Orbitrap Exploris": 1, "": 4}))
            run = dict(prepared={s: dict(representatives=len(r)) for s, r in data.items()},
                       exclusions={"vs20:pos:input": 9, "c26:pos:input": 7, "c26:vs21_overlap_records": 1},
                       p3_exclusions={"c26:excluded_ingest_lib": 5}, p2_parity={"c26": dict(result="not_checked")})
            args = SimpleNamespace(sources=["vs20", "vs21", "c26"], out=out, ruleset_dir=Path(d), console=None)
            saved = m.load_rules, m.vs21_guard
            m.load_rules = lambda d: ([rule(), rule("r2", measure=False)], None, {})
            m.vs21_guard = lambda a: True
            try:
                quiet(m.inventory, args, run)
            finally:
                m.load_rules, m.vs21_guard = saved
            names = sorted(p.name for p in out.iterdir() if p.name.startswith("inventory"))
            self.assertEqual(names, ["inventory.local-vs21.local-c26.json", "inventory.local-vs21.local-c26.tsv"])
            self.assertFalse((out / "fire").exists())
            doc = json.loads((out / names[0]).read_text(encoding="utf-8"))
            self.assertEqual(doc["prepare"]["c26"]["exclusions"], {"c26:pos:input": 7, "c26:vs21_overlap_records": 1})
            self.assertEqual(doc["prepare"]["c26"]["p3_exclusions"], {"c26:excluded_ingest_lib": 5})
            self.assertEqual(doc["prepare"]["vs20"]["stats"], dict(representatives=2))
            self.assertEqual([(r["instrument_type"], r["records"], r["instrument_class"]) for r in doc["instrument_types"]["c26"]],
                             [("", 4, "other"), ("LC-ESI-QTOF", 2, "QTOF-type"), ("Orbitrap Exploris", 1, "Orbitrap/FT")])
            power = {r["source"]: r for r in m.read_rows(out / names[1])}
            self.assertEqual([(power[s]["n_positive"], power[s]["c_positive"]) for s in ("pooled", "vs20", "vs21", "c26")],
                             [("4", "3"), ("1", "1"), ("1", "1"), ("2", "1")])
            self.assertEqual(power["c26"]["n_universe"], "2")  # the negative record is not in a positive-mode universe
            self.assertFalse(any("precision" in k or "fires" in k for r in doc["power"] for k in r))
            self.assertEqual((run["inventory"]["rules"], run["inventory"]["audit"]), (1, "pass"))

    def test_cli_accepts_the_inventory_stage(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(SystemExit):  # needs --ruleset-dir, like every other stage
                quiet(m.main, ["inventory", "--out", str(Path(d) / "p3cs-20261127-abcdef1")])
            with self.assertRaises(SystemExit):
                quiet(m.main, ["inventry", "--ruleset-dir", d])


def shard_dir(d, name="r-00", build="b1", engine="e1", errors=0, records=1, msp=None, local=False):
    out = Path(d) / name
    out.mkdir(parents=True, exist_ok=True)
    passes = {p: dict(evidence=f"Diagrules evidence-{p}-1.jsonl", refusals=[], counts=dict(errors=errors, skipped=0))
              for p in m.PASSES}
    for p in m.PASSES:
        (out / passes[p]["evidence"]).write_text("", encoding="utf-8")
    manifest = dict(rule_set=dict(sha256=m.EXPECTED_RULESET_SHA256, pinned=True, version=m.EXPECTED_RULESET_VERSION),
                    local_vs21=local, inputs=[dict(sha256=p2.sha256(msp), records=records)], passes=passes,
                    build_identity=build, engine_semantics=engine)
    (out / "Diagrules run-1.json").write_text(json.dumps(manifest), encoding="utf-8")
    return out


class ConsoleIdentityTests(unittest.TestCase):
    def test_shards_must_agree_and_match_the_freeze(self):
        with tempfile.TemporaryDirectory() as d:
            msp = Path(d) / "a.msp"
            msp.write_text("x", encoding="utf-8")
            r = {f"vs20/representatives/{i}": m.validate_shard(shard_dir(d, f"s{i}", msp=msp), msp, 1, False) for i in range(2)}
            self.assertEqual(m.check_shards(r), dict(build_identity="b1", engine_semantics="e1"))
            self.assertTrue(m.check_shards(r, dict(build_identity="b1", engine_semantics="e1")))
            with self.assertRaisesRegex(RuntimeError, "freeze"):
                m.check_shards(r, dict(build_identity="b2", engine_semantics="e1"))
            r["vs20/chimeric/x"] = m.validate_shard(shard_dir(d, "s9", build="b2", msp=msp), msp, 1, False)
            with self.assertRaisesRegex(RuntimeError, "disagree"):
                m.check_shards(r)

    def test_fire_fails_on_record_errors(self):
        with tempfile.TemporaryDirectory() as d:
            msp = Path(d) / "a.msp"
            msp.write_text("x", encoding="utf-8")
            r = {"vs20/representatives/0": m.validate_shard(shard_dir(d, errors=1, msp=msp), msp, 1, False)}
            with self.assertRaisesRegex(RuntimeError, "record errors"):
                m.check_shards(r)

    def test_reuse_refused_when_the_manifest_disagrees(self):
        with tempfile.TemporaryDirectory() as d:
            msp = Path(d) / "a.msp"
            msp.write_text("x", encoding="utf-8")
            out = shard_dir(d, msp=msp, records=1)
            self.assertEqual(m.run_console(Path(d) / "missing.exe", msp, Path(d) / "x.ini", out, 1, False), "reused")
            with self.assertRaisesRegex(RuntimeError, "Refusing to reuse"):
                m.run_console(Path(d) / "missing.exe", msp, Path(d) / "x.ini", out, 2, False)
            msp.write_text("y", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "Refusing to reuse"):
                m.run_console(Path(d) / "missing.exe", msp, Path(d) / "x.ini", out, 1, False)

    def test_console_files_cover_the_bin_folder(self):
        with tempfile.TemporaryDirectory() as d:
            b = Path(d) / "bin"
            (b / "lib" / ".svn").mkdir(parents=True)
            for name in ("app.exe", "app.exe.config", "lib/x.dll", "lib/.svn/y.dll", "readme.txt"):
                (b / name).write_bytes(name.encode())
            self.assertEqual(sorted(m.console_files(b / "app.exe")), ["app.exe", "app.exe.config", "lib/x.dll"])


class GuardTests(unittest.TestCase):
    def files(self, d):
        rs = Path(d) / "RuleSet"
        (rs / "rules").mkdir(parents=True)
        (rs / "rules/a.json").write_text("{}", encoding="utf-8")
        (Path(d) / "bin").mkdir()
        console = Path(d) / "bin" / "console.exe"
        console.write_bytes(b"MZ")
        (Path(d) / "bin" / "lib.dll").write_bytes(b"dll")
        return rs, console

    def freeze(self, d, rs, console, seed=1, bootstrap=10, name="freeze-20261123-abcdef1.json"):
        path = Path(d) / name
        path.write_text(json.dumps(dict(m.freeze_hashes(rs, console, seed, bootstrap),
                                        console_identity=dict(build_identity="b", engine_semantics="e"))), encoding="utf-8")
        return path

    def test_vs21_refused_without_flag_freeze_or_on_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            rs, console = self.files(d)
            with self.assertRaises(PermissionError):
                m.vs21_guard(SimpleNamespace(vs21=False, freeze=None, ruleset_dir=rs, console=console), d)
            with self.assertRaises(PermissionError):
                m.vs21_guard(SimpleNamespace(vs21=True, freeze=None, ruleset_dir=rs, console=console), d)
            freeze = self.freeze(d, rs, console)
            args = SimpleNamespace(vs21=True, freeze=freeze, ruleset_dir=rs, console=console, seed=1, bootstrap=10)
            self.assertTrue(m.vs21_guard(args, d))
            (rs / "rules/a.json").write_text('{"x": 1}', encoding="utf-8")
            with self.assertRaisesRegex(PermissionError, "ruleset_files"):
                m.vs21_guard(args, d)
            (rs / "rules/a.json").write_text("{}", encoding="utf-8")
            console.write_bytes(b"MZ2")
            with self.assertRaisesRegex(PermissionError, "console"):
                m.vs21_guard(args, d)

    def test_vs21_freeze_refuses_a_changed_insource_table(self):
        with tempfile.TemporaryDirectory() as d:
            rs, console = self.files(d)
            table = Path(d) / "insource-p3.tsv"
            table.write_text(INS_HEADER + INS_ROWS, encoding="utf-8")
            saved = m.INSOURCE_TABLE
            m.INSOURCE_TABLE = table
            try:
                args = SimpleNamespace(vs21=True, freeze=self.freeze(d, rs, console), ruleset_dir=rs, console=console,
                                       seed=1, bootstrap=10)
                self.assertTrue(m.vs21_guard(args, d))
                table.write_text(INS_HEADER + INS_ROWS + "*\t[M-H]-\t-CO2\t[M-H-CO2]-\ts\tn\n", encoding="utf-8")
                with self.assertRaisesRegex(PermissionError, "insource_table"):
                    m.vs21_guard(args, d)
            finally:
                m.INSOURCE_TABLE = saved

    def test_vs21_freeze_refuses_changed_criteria_of_either_tier(self):
        with tempfile.TemporaryDirectory() as d:
            rs, console = self.files(d)
            args = SimpleNamespace(vs21=True, freeze=self.freeze(d, rs, console), ruleset_dir=rs, console=console,
                                   seed=1, bootstrap=10)
            self.assertTrue(m.vs21_guard(args, d))
            for name, field, value in (("CRITERIA", "wilson_target", 0.94), ("CRITERIA_PROVISIONAL", "n_min", 25),
                                       ("CRITERIA_PROVISIONAL", "prior_stress_gated", True),
                                       ("CRITERIA_EXPLORATORY", "wilson_target", 0.79)):
                saved = getattr(m, name)
                setattr(m, name, dict(saved, **{field: value}))
                try:
                    with self.assertRaisesRegex(PermissionError, "criteria", msg=name + "." + field):
                        m.vs21_guard(args, d)
                finally:
                    setattr(m, name, saved)
            self.assertTrue(m.vs21_guard(args, d))

    def test_vs21_freeze_covers_dlls_parameters_and_path(self):
        """Regression: a changed DLL beside the console, another seed or bootstrap, or a freeze file outside the
        ws8-p3 root is refused."""
        with tempfile.TemporaryDirectory() as d:
            rs, console = self.files(d)
            freeze = self.freeze(d, rs, console)
            args = SimpleNamespace(vs21=True, freeze=freeze, ruleset_dir=rs, console=console, seed=1, bootstrap=10)
            self.assertTrue(m.vs21_guard(args, d))
            for field, value in (("seed", 2), ("bootstrap", 11)):
                with self.assertRaisesRegex(PermissionError, "parameters"):
                    m.vs21_guard(SimpleNamespace(**dict(vars(args), **{field: value})), d)
            frozen = json.loads(freeze.read_text(encoding="utf-8"))
            self.assertEqual(frozen["parameters"]["console_ini"], m.CONSOLE_INI)
            frozen["parameters"]["console_ini"] = "Ms2Tolerance=0.02"
            freeze.write_text(json.dumps(frozen), encoding="utf-8")
            with self.assertRaisesRegex(PermissionError, "parameters"):
                m.vs21_guard(args, d)
            freeze = self.freeze(d, rs, console)
            (Path(d) / "bin" / "lib.dll").write_bytes(b"other")
            with self.assertRaisesRegex(PermissionError, "console_files"):
                m.vs21_guard(args, d)
            (Path(d) / "bin" / "lib.dll").write_bytes(b"dll")
            other = self.freeze(d, rs, console, name="frozen.json")
            with self.assertRaisesRegex(PermissionError, "freeze-"):
                m.vs21_guard(SimpleNamespace(**dict(vars(args), freeze=other)), d)
            frozen = json.loads(freeze.read_text(encoding="utf-8"))
            del frozen["console_identity"]
            freeze.write_text(json.dumps(frozen), encoding="utf-8")
            with self.assertRaisesRegex(PermissionError, "console_identity"):
                m.vs21_guard(args, d)


    def test_freeze_pins_the_c26_parameters(self):
        with tempfile.TemporaryDirectory() as d:
            rs, console = self.files(d)
            args = SimpleNamespace(vs21=True, freeze=self.freeze(d, rs, console), ruleset_dir=rs, console=console,
                                   seed=1, bootstrap=10)
            self.assertTrue(m.vs21_guard(args, d))
            frozen = json.loads(args.freeze.read_text(encoding="utf-8"))["c26"]
            self.assertEqual((frozen["train_sha256"], frozen["ppm_max"], frozen["excluded_ingest_libs"]),
                             (m.C26_SHA256, 10, ["enveda-180"]))
            for name, value in (("C26_PPM_MAX", 20), ("C26_EXCLUDED_LIBS", ()), ("C26_SHA256", "0" * 64)):
                saved = getattr(m, name)
                setattr(m, name, value)
                try:
                    with self.assertRaisesRegex(PermissionError, "c26", msg=name):
                        m.vs21_guard(args, d)
                finally:
                    setattr(m, name, saved)
            self.assertTrue(m.vs21_guard(args, d))

    def test_cli_source_lists(self):
        with tempfile.TemporaryDirectory() as d:
            rs, _ = self.files(d)
            out = Path(d) / "p3cs-20261123-abcdef1"
            for sources in ("vs20,c26", "vs21,c26", "c26", "vs20,c26,vs21"):
                with self.assertRaises(SystemExit, msg=sources):
                    quiet(m.main, ["prepare", "--sources", sources, "--vs21", "--ruleset-dir", str(rs), "--out", str(out)])
            with self.assertRaises(SystemExit):  # the full list still needs --vs21 and a matching freeze
                quiet(m.main, ["prepare", "--sources", "vs20,vs21,c26", "--ruleset-dir", str(rs), "--out", str(out)])
            self.assertFalse(out.exists())

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
            self.assertIn("MsfinderConsoleApp.exe", run["console"]["files"])
            self.assertTrue(any(k.endswith(".dll") for k in run["console"]["files"]))
            self.assertEqual(m.console_identity(CONSOLE), run["console"]["identity"])
            self.assertGreater(run["insource"]["vs20"]["records"], 0)
            self.assertEqual(quiet(m.fire, args, run), None)  # every complete shard validates and is reused
            self.assertEqual({v["state"] for v in run["fire"].values()}, {"reused"})
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
            ins = [r for r in m.read_rows(out / "rules_synthetic.csv") if r["set"] == "insource"]
            self.assertTrue(any(r["transform"] == "all" for r in ins))
            self.assertTrue(any(f["set"] == "insource" for f in fires))


if __name__ == "__main__":
    unittest.main()
