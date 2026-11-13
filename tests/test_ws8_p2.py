"""WS8 P2 acceptance, using synthetic spectra and hand-picked public molecules."""
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from tools.ws8_p2 import p2_diagnostics as p2

HAS_RDKIT = importlib.util.find_spec("rdkit") is not None
P1_SHA = "5b47df5595050e59eb3a5b43fba110ad7094f175fd57166394b003242b31f71c"


class ObservationTests(unittest.TestCase):
    def setUp(self):
        self.d = dict(kind="product_ion", mz_or_loss=100.1, tolerance_da=0.01,
                      min_rel_intensity=0.01, precursor_types=["[M+H]+"])

    def observe(self, peaks, **kw):
        return p2.observe_one(peaks, 300.2, "[M+H]+", self.d, **kw)

    def test_tolerance_includes_edge(self):
        self.assertTrue(self.observe([(100.11, 1), (150.1, 100)])[0])
        self.assertFalse(self.observe([(100.11001, 1), (150.1, 100)])[0])
        self.assertFalse(self.observe([(100.106, 1), (150.1, 100)], tolerance=0.005)[0])

    def test_intensity_includes_floor(self):
        self.assertTrue(self.observe([(100.1, 1), (150.1, 100)])[0])
        self.assertFalse(self.observe([(100.1, 0.999), (150.1, 100)])[0])
        self.assertFalse(self.observe([(100.1, 4.99), (150.1, 100)], floor=0.05)[0])

    def test_precursor_plus_minus_one_excluded_from_base(self):
        self.assertTrue(self.observe([(100.1, 1), (150.1, 100), (300.2, 1000000), (301.2, 1000000)])[0])
        self.assertFalse(self.observe([(300.2, 1000)])[0])

    def test_neutral_loss_from_precursormz(self):
        self.d.update(kind="neutral_loss", mz_or_loss=176.032)
        self.assertTrue(self.observe([(124.168, 10), (150.1, 100)])[0])
        self.assertFalse(self.observe([(176.032, 10), (150.1, 100)])[0])
        self.d["mz_or_loss"] = 0.5
        self.assertFalse(self.observe([(299.7, 10), (150.1, 100)])[0])

    def test_precursor_type_gate(self):
        self.assertEqual(p2.observe_one([(100.1, 100)], 300.2, "[M+Na]+", self.d), (False, 0.0))

    def test_resolution_both_criteria(self):
        self.assertTrue(p2.nominal_mass([(100.0001, 1), (200.0009, 1)], [4, 4]))
        self.assertTrue(p2.nominal_mass([(100.1, 1), (200.2, 1)], [1, 1]))
        self.assertFalse(p2.nominal_mass([(100.12, 1), (200.21, 1)], [2, 2]))
        self.assertFalse(p2.nominal_mass([(100.12, 1), (200, 1)], [2, 0]))

    def test_parser_precision_and_allowlist(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "synthetic.msp"
            path.write_text('NAME: Synthetic\nSMILES: CCCC\nPRECURSORMZ: 300.2\nNum Peaks: 2\n100.1 2 "comment"\n200.120 3\n', encoding="utf-8")
            meta, peaks, decimals = next(p2.iter_msp(path))
            self.assertNotIn("name", meta)
            self.assertEqual(decimals, [1, 3])
            self.assertFalse(p2.nominal_mass(peaks, decimals))

    def test_formula_gate_counts_residue_heavy_elements(self):
        self.assertTrue(p2.formula_gate("C20H30O4S", "HSO4"))
        self.assertFalse(p2.formula_gate("C20H30O4", "HSO4"))
        self.assertFalse(p2.formula_gate("C20H30O3S", "HSO4"))
        self.assertTrue(p2.formula_gate("C20H30NO4P+", "C5H15NO4P"))
        self.assertTrue(p2.formula_gate("C20H30O4S-2", "HSO4"))

    def test_energy_units_and_instrument_strata(self):
        self.assertEqual(p2.energy_bin("35 eV"), "eV:20-<40")
        self.assertEqual(p2.energy_bin("NCE 40%"), "NCE:>=40")
        self.assertEqual(p2.energy_bin("10,20,40 eV"), "unparseable")
        self.assertEqual(p2.instrument_bin("LC-ESI-QTOF"), "QTOF-type")
        self.assertEqual(p2.instrument_bin("LC-ESI-ITFT"), "Orbitrap/FT")


class StatisticsTests(unittest.TestCase):
    def test_statistics_independent_values(self):
        s = p2.statistics(90, 10, 30, 870)
        self.assertEqual(s["precision"], 0.9)
        self.assertEqual(s["recall"], 0.75)
        self.assertEqual(s["fpr"], 10 / 880)
        self.assertEqual(s["prevalence"], 0.12)
        self.assertAlmostEqual(s["lr_plus"], 66)
        se = math.sqrt(1/90 - 1/120 + 1/10 - 1/880)
        self.assertAlmostEqual(s["lr_lo"], 66 * math.exp(-1.959963984540054 * se))
        self.assertAlmostEqual(s["prior_stress_precision"], .06*s["lr_lo"]/(.94+.06*s["lr_lo"]))
        self.assertAlmostEqual(s["wilson_lo"], 0.8256343, places=6)

    def test_zero_cell_correction_and_unestimable_denominators(self):
        s = p2.statistics(100, 0, 0, 900)
        self.assertAlmostEqual(s["lr_plus"], (100.5/101)/(.5/901))
        self.assertEqual(s["recall"], 1)
        self.assertEqual(s["fpr"], 0)
        self.assertIsNone(p2.statistics(0, 0, 0, 0)["lr_plus"])
        self.assertIsNone(p2.statistics(100, 0, 0, 0)["prior_stress_precision"])

    def test_rule_requires_all_conditions(self):
        good = p2.statistics(1000, 0, 0, 9000)
        self.assertEqual(p2.status(good), "pass")
        self.assertEqual(p2.status(p2.statistics(99, 0, 0, 900)), "not_estimable")
        self.assertEqual(p2.status(p2.statistics(100, 10, 0, 900)), "fail")
        self.assertEqual(p2.status(good, [dict(n_obs=30, precision=.94)]), "fail")
        self.assertEqual(p2.status(good, [dict(n_obs=29, precision=.1)]), "pass")
        self.assertEqual(p2.status(dict(good, prior_stress_precision=.8999)), "fail")
        self.assertEqual(p2.status(dict(good, wilson_lo=.9499)), "fail")

    def test_weighted_spectra_total_one_per_unit(self):
        counts = {}
        for obs in [True, True, False, False]:
            p2.add_confusion(counts, "unit", True, obs, .25)
        self.assertEqual(counts["unit"], [.5, 0, .5, 0])

    def test_mining_selection_only_public_scoring_only_heldout(self):
        public = [dict(mode="pos", motifs={"m"}, bins={("product_ion", 10010)}) for _ in range(200)]
        public += [dict(mode="pos", motifs=set(), bins={("neutral_loss", 5000)}) for _ in range(200)]
        selected = p2.select_candidates(public, ["m"])
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0]["vs20_n"], 200)
        heldout = [dict(mode="pos", motifs=set(), bins={("product_ion", 10010), ("product_ion", 20020)}) for _ in range(50)]
        scored = p2.score_candidates(selected, heldout)
        self.assertEqual(scored[0]["vs21_n_obs"], 50)
        self.assertEqual(scored[0]["vs21_precision"], 0)
        self.assertEqual(scored[0]["status"], "fail")
        self.assertEqual(p2.select_candidates(public, ["m"]), selected)
        self.assertEqual(p2.score_candidates(selected, [])[0]["status"], "not_estimable")
        self.assertNotIn(20020, [c["mass_bin"] for c in selected])

    def test_mining_collapses_all_precursor_types_per_compound(self):
        rows = [dict(key="synthetic", mode="pos", adduct=adduct, representative=True,
                     motifs=["m"], mz=[mz, 150.12], intensity=[100, 10], precursor=300.2)
                for adduct, mz in [("[M+H]+", 100.12), ("[M+2H]2+", 110.12)]]
        compounds = p2.mining_compounds(rows)
        self.assertEqual(len(compounds), 1)
        self.assertIn(("product_ion", 10012), compounds[0]["bins"])
        self.assertIn(("product_ion", 11012), compounds[0]["bins"])

    def test_csv_json_recomputation_and_mined_never_allowed(self):
        # A4 and A6 use the actual precision/choose/audit path on synthetic counts.
        if not HAS_RDKIT:
            self.skipTest("RDKit required to validate diagnostic table in precision stage")
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            d = p2.diagnostics()[0]
            rows = [dict(diagnostic_id=d["diagnostic_id"], ion_mode=d["ion_mode"], source=s,
                         precursor_type="all", variant="default", stratum="all", tp=1000, fp=0, fn=0, tn=9000)
                    for s in ("vs20", "vs21")]
            p2.write_csv(out / "confusion_counts.csv", rows)
            args, run = SimpleNamespace(out=out, sources=["vs20", "vs21"]), dict(inputs=[])
            p2.precision(args, run)
            # A planted passing mined diagnostic cannot be consumed by choose.
            p2.write_csv(out / "mined_candidates.csv", [dict(motif_id="m", diagnostic_id="mined_only", status="pass")])
            from unittest.mock import patch
            with patch.object(p2, "report"):
                p2.choose(args, run)
            result = json.loads((out / "diagnostics_precision.json").read_text())
            self.assertEqual(result["si"]["diagnostics_allowed"], [dict(diagnostic_id=d["diagnostic_id"], ion_mode=d["ion_mode"])])
            self.assertEqual(p2.audit_outputs(out), dict(P2_A4="pass", P2_A5="pass"))
            next(r for r in result["by_pair"] if r["diagnostic_id"] == d["diagnostic_id"])["lr_plus"] += 1
            p2.write_json(out / "diagnostics_precision.json", result)
            with self.assertRaises(AssertionError):
                p2.audit_outputs(out)

    def test_privacy_scan_rejects_identifying_column(self):
        # Scan source itself is not the input; fixture never contains a real identity.
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            p2.write_json(out / "diagnostics_precision.json", dict(by_pair=[]))
            (out / "diagnostics_precision.csv").write_text("diagnostic_id\n")
            (out / "bad.csv").write_text("smiles\nsynthetic\n")
            with self.assertRaises(AssertionError):
                p2.audit_outputs(out)


@unittest.skipUnless(HAS_RDKIT, "RDKit available only in the D5 chemistry interpreter")
class ChemistryTests(unittest.TestCase):
    def test_table_formula_masses_and_sources(self):
        from rdkit import Chem
        table = Chem.GetPeriodicTable()
        for r in p2.diagnostics():
            atoms = p2.formula_counts(r["formula"])
            expected = math.fsum(table.GetMassForIsotope(table.GetAtomicNumber(e), table.GetMostCommonIsotope(e)) * n for e, n in atoms.items())
            if r["kind"] == "product_ion":
                expected += -.000549 if r["ion_mode"] == "pos" else .000549
            self.assertAlmostEqual(expected, r["mz_or_loss"], delta=1e-5)
            self.assertTrue(r["source"].strip())

    def test_new_motifs_public_positive_negative(self):
        from rdkit import Chem
        rows = p2.motif_rows()
        extras = p2.read_tsv(p2.HERE / "motifs-p2-extra.tsv")
        self.assertEqual(len(rows), len({r["motif_id"] for r in rows}))
        # Public chemical definitions, not sampled from any library. Positives are
        # cited in the extra TSV (PubChem CID 11768891, 10917, 159864 and 190).
        # Negatives, respectively (stereochemistry is ignored):
        # https://pubchem.ncbi.nlm.nih.gov/compound/Methyl-beta-D-glucopyranoside
        # https://pubchem.ncbi.nlm.nih.gov/compound/Choline
        # https://pubchem.ncbi.nlm.nih.gov/compound/Taurine
        # https://pubchem.ncbi.nlm.nih.gov/compound/Guanine
        fixtures = {
            "pentopyranoside": ("COC1OCC(O)C(O)C1O", "COC1OC(CO)C(O)C(O)C1O"),
            "carnitine": ("C[N+](C)(C)CC(O)CC(=O)[O-]", "C[N+](C)(C)CCO"),
            "taurine_conjugate": ("CC(=O)NCCS(=O)(=O)O", "NCCS(=O)(=O)O"),
            "adenine": ("Nc1ncnc2[nH]cnc12", "Nc1nc2[nH]cnc2c(=O)[nH]1"),
        }
        for r in extras:
            q = Chem.MolFromSmarts(r["smarts"])
            self.assertIsNotNone(q)
            yes, no = fixtures[r["motif_id"]]
            self.assertTrue(p2.standardise(yes).HasSubstructMatch(q), r["motif_id"])
            self.assertFalse(p2.standardise(no).HasSubstructMatch(q), r["motif_id"])
            self.assertIn("PubChem", r["source"])
        # Expose P1's broad pentoside overlap without modifying that file.
        q = Chem.MolFromSmarts(next(r["smarts"] for r in rows if r["motif_id"] == "pentoside"))
        self.assertTrue(p2.standardise(fixtures["pentopyranoside"][1]).HasSubstructMatch(q))

    def test_adenine_matches_adenosine_public_cid_60961(self):
        from rdkit import Chem
        q = Chem.MolFromSmarts(next(r["smarts"] for r in p2.motif_rows() if r["motif_id"] == "adenine"))
        self.assertTrue(p2.standardise("Nc1ncnc2n(cnc12)C1OC(CO)C(O)C1O").HasSubstructMatch(q))

    def test_p1_motifs_byte_unchanged(self):
        self.assertEqual(p2.sha256(p2.P1_MOTIFS), P1_SHA)

    def test_formula_assignment_and_ambiguity(self):
        rows = [dict(mz_or_loss=18.01, kind="neutral_loss", ion_mode="pos"),
                dict(mz_or_loss=1.00, kind="neutral_loss", ion_mode="pos"),
                dict(mz_or_loss=500.00, kind="product_ion", ion_mode="pos")]
        p2.assign_formulas(rows)
        self.assertEqual(rows[0]["formula"], "H2O")
        self.assertEqual(rows[1]["formula_status"], "unassigned")
        self.assertEqual(rows[2]["formula_status"], "ambiguous")

    def test_prepare_representatives_resolution_dedup_and_weights(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            def record(smiles, peaks, adduct="[M+H]+", mode="Positive"):
                return f"SMILES: {smiles}\nPRECURSORMZ: 300.2\nPRECURSORTYPE: {adduct}\nIONMODE: {mode}\nNum Peaks: {len(peaks)}\n" + "\n".join(f"{m} {i}" for m, i in peaks) + "\n\n"
            public = root / "public-pos.msp"
            public.write_text(record("CCCC", [(100.12, 20), (200.23, 10)]) +
                              record("CCCC", [(100.12, 20), (200.23, 10), (250.12, 5)]) +
                              record("CCCCC", [(100, 10), (200, 20)]), encoding="utf-8")
            neg = root / "public-neg.msp"
            neg.write_text("", encoding="utf-8")
            (root / "MSMS-Private-pos-VS21.msp").write_text(
                record("CCCC", [(100.12, 10)]) + record("CCCCC", [(100.12, 10)]) +
                record("CCCCCC", [(100.12, 10)]), encoding="utf-8")
            (root / "MSMS-Private-neg-VS21.msp").write_text("", encoding="utf-8")
            out = root / "output"
            out.mkdir()
            args = SimpleNamespace(out=out, sources=["vs20", "vs21"], workers=1, limit=None)
            run = {}
            with patch.object(p2, "LIB_ROOT", root), patch.object(p2, "PUBLIC", [(public.name, p2.sha256(public)), (neg.name, p2.sha256(neg))]):
                p2.prepare(args, run)
            pubrows = list(p2.parquet_rows(p2.row_path(out, "prepared", "vs20")))
            private = list(p2.parquet_rows(p2.row_path(out, "prepared", "vs21")))
            self.assertEqual(len(pubrows), 2)
            self.assertEqual([len(r["mz"]) for r in pubrows if r["representative"]], [3])
            self.assertTrue(all(r["spectrum_count"] == 2 for r in pubrows))
            # Even a nominal-only VS20 compound excludes the corresponding held-out compound.
            self.assertEqual(len(private), 1)
            self.assertEqual(run["exclusions"]["vs21:vs20_overlap_records"], 2)
            self.assertEqual(run["exclusions"]["vs20:pos:nominal_mass"], 1)


if __name__ == "__main__":
    unittest.main()
