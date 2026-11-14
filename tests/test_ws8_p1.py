"""P1 acceptance on synthetic numbers and hand-picked public chemistry examples."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

from tools.ws8_p1 import mcs_core as core
from tools.ws8_p1 import p1_si_bar as p1
from tools.ws8_p1 import audit_p1 as audit

HAS_RDKIT = importlib.util.find_spec("rdkit") is not None


class StatisticsTests(unittest.TestCase):
    def test_pilot_samples_five_percent_of_sparse_cells(self):
        self.assertEqual(p1.sample_count(10, 600, True), 1)
        self.assertEqual(p1.sample_count(800, 600, True), 30)
        self.assertEqual(p1.sample_count(0, 150, True), 0)

    def test_wilson_known_values(self):
        self.assertAlmostEqual(p1.wilson(100, 100)[0], 0.9630065, places=6)
        self.assertEqual(p1.wilson(0, 0), (0.0, 1.0))
        self.assertAlmostEqual(p1.wilson(50, 100)[0], 0.4038315, places=6)

    def test_cumulative_uses_all_sizes_above_threshold(self):
        row = p1.cumulative([(1, False), (2, True), (3, True), (3, False)], 2)
        self.assertEqual((row["n"], row["success"]), (3, 2))
        self.assertEqual(p1.cumulative([], 1)["n"], 0)

    def test_chance_linear_quantile(self):
        self.assertAlmostEqual(p1.quantile(list(range(1, 101)), 0.95), 95.05)
        self.assertIsNone(p1.quantile([], 0.95))
        self.assertEqual(p1.quantile([7], 0.95), 7)

    def test_bar_requires_both_families_and_strict_chance(self):
        curves = []
        for family in p1.BAR_FAMILIES:
            for x in (5, 6, 7, 8):
                curves.append(dict(family=family, threshold=x, n=100, wilson_lo=0.97 if x >= (6 if family == "spectral_nn" else 7) else 0.9))
        bar = p1.choose_axis(curves, 5)
        self.assertEqual(bar["value"], 7)
        self.assertEqual(bar["deciding_family"], "same_formula_db_random")
        self.assertEqual(p1.choose_axis(curves, 7)["value"], 8)
        self.assertIsNone(p1.choose_axis(curves, 8))
        self.assertIsNone(p1.choose_axis(curves[:4], 1))

    def test_bar_rechecks_nonmonotonic_wilson_and_minimum_n(self):
        curves = [dict(family=f, threshold=x, n=99 if x == 3 else 100,
                       wilson_lo=0.96 if (f == "spectral_nn" and x == 1) or (f != "spectral_nn" and x >= 2) else 0.8)
                  for f in p1.BAR_FAMILIES for x in (1, 2, 3)]
        self.assertIsNone(p1.choose_axis(curves, 0))

    def test_parameter_identity(self):
        import hashlib
        digest = hashlib.sha256(json.dumps(core.mcs_params, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        self.assertEqual(digest, core.mcs_params_sha256)
        self.assertEqual(core.mcs_params["timeout"], 10)
        self.assertFalse(core.mcs_params["matchChiralTag"])

    def test_written_curves_choose_schema_and_acceptance(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            correct, chance = [], []
            for k in (2, 3, 4, 5):
                for axis, thresholds, q in (("heavy_atoms", [5, 6, 7], 5), ("fraction_smallest", [0.5, 0.55, 0.6], 0.5)):
                    chance.append(dict(kind="matched", source="all", lipid="all", heavy_bin="all", k=k, axis=axis, p95=q))
                    for family in p1.BAR_FAMILIES:
                        for threshold in thresholds:
                            correct.append(dict(source="all", lipid="all", heavy_bin="all", family=family, k=k,
                                                axis=axis, curve="cumulative", threshold=threshold, **p1.cumulative([(threshold, True)] * 100, threshold)))
            p1.write_csv(out / "curves_correctness.csv", correct)
            p1.write_csv(out / "curves_chance.csv", chance)
            p1.write_csv(out / "motif_precision.csv", [dict(motif_id="synthetic", family=f, n=100, precision=1,
                                                            wilson_lo=p1.wilson(100, 100)[0], status="pass") for f in p1.BAR_FAMILIES])
            args = SimpleNamespace(out=out, target=0.95, n_min=100)
            p1.choose(args, {"inputs": []})
            result = json.loads((out / "si_bar.json").read_text())
            for field in ("label", "run_id", "target", "interval", "n_min", "chance_quantile", "families_used_for_bar",
                          "mcs_params", "mcs_params_sha256", "motifs_file_sha256", "by_k", "global", "by_stratum", "strata_above_pooled", "motifs", "inputs"):
                self.assertIn(field, result)
            for k, row in result["by_k"].items():
                self.assertEqual(row["heavy_atoms"], 6)
                self.assertEqual(row["fraction"], 0.55)
                for axis in ("heavy_atoms", "fraction"):
                    x = row[axis]
                    chance_field = "chance95_atoms" if axis == "heavy_atoms" else "chance95_fraction"
                    self.assertGreater(x, row[chance_field])
                    selected = [r for r in p1.read_csv(out / "curves_correctness.csv") if r["k"] == k and
                                r["axis"] == ("heavy_atoms" if axis == "heavy_atoms" else "fraction_smallest") and float(r["threshold"]) == x]
                    self.assertTrue(all(int(r["n"]) >= 100 and float(r["wilson_lo"]) >= 0.95 for r in selected))
            self.assertEqual(result["global"]["heavy_atoms"], 6)
            audit.audit_bars(out, p1.read_csv(out / "curves_correctness.csv"), p1.read_csv(out / "curves_chance.csv"))

    def test_audit_rejects_invented_bar(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            result = {k: None for k in ("label", "run_id", "target", "interval", "n_min", "chance_quantile",
                                       "families_used_for_bar", "mcs_params", "mcs_params_sha256", "motifs_file_sha256",
                                       "global", "strata_above_pooled", "motifs", "inputs")}
            result.update(target=0.95, n_min=100, by_stratum={}, by_k={"2": {
                "heavy_atoms": 6, "fraction": None, "chance95_atoms": 5, "chance95_fraction": 0.5,
                "axis_status": {"heavy_atoms": "estimable", "fraction": "not_estimable"}}})
            p1.write_json(out / "si_bar.json", result)
            with self.assertRaisesRegex(AssertionError, "independently searched"):
                audit.audit_bars(out, [], [])


@unittest.skipUnless(HAS_RDKIT, "RDKit is optional; use the D5 interpreter")
class ChemistryTests(unittest.TestCase):
    def test_standardisation_fragments_stereo_isotopes_and_formula(self):
        mol = core.standardise("[Na+].[O-]C(=O)[C@@H]([13CH3])O")
        self.assertEqual(core.hill_formula(mol), "C3H6O3")
        self.assertEqual(core.first_block(mol), core.first_block(core.standardise("CC(O)C(=O)O")))

    def test_rejection_reasons(self):
        for smiles, reason in (("C", "heavy_atom_range"), ("CCCC[Fe]", "element"), ("[CH2]CCC", "radical"), ("invalid", "parse")):
            with self.subTest(reason=reason), self.assertRaisesRegex(core.RejectedStructure, reason):
                core.standardise(smiles)

    def test_truth_in_set_is_contained(self):
        truth = core.standardise("Cc1ccccc1O")
        smarts, size, timeout = core.mcs_over_set([truth, core.standardise("COc1ccccc1")])
        self.assertFalse(timeout)
        self.assertGreater(size, 0)
        self.assertTrue(core.contained(truth, smarts))

    def test_aromatic_benzene_does_not_match_cyclohexane(self):
        smarts, size, timeout = core.mcs_over_set([core.standardise("Cc1ccccc1"), core.standardise("Oc1ccccc1")])
        self.assertFalse(timeout)
        self.assertEqual(size, 6)
        self.assertFalse(core.contained(core.standardise("C1CCCCC1"), smarts))

    def test_complete_rings_and_ring_only(self):
        smarts, size, timeout = core.mcs_over_set([core.standardise("C1CCCCC1"), core.standardise("CCCCCC")])
        self.assertFalse(timeout)
        self.assertEqual(size, 0)

    def test_leave_compound_out_sets_and_nested_members(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            molecules = [core.standardise(s) for s in ("CCCCCCC", "CC(C)CCCC", "CCC(C)CCC", "CC(C)C(C)CC", "CC(C)(C)CCC", "CCC(C)(C)CC")]
            db = {core.first_block(m): p1.chem_record(__import__("rdkit").Chem.MolToSmiles(m)) for m in molecules}
            truth = dict(next(iter(db.values())), source="vs20", stratum="0:4-15", lipid=0, heavy_bin="4-15", motif_oversample=False)
            p1.dump(out / "prepared.local-vs21.pkl", {"db": db, "truths": [truth]})
            p1.dump(out / "neighbours.local-vs21.pkl", [dict(source="vs20", key=truth["key"], spectral_nn=[], spectral_nn_any=[], same_formula_spectral=[])])
            args = SimpleNamespace(out=out, sources=["vs20"], seed=20261113, pilot=True)
            p1.sets(args, {})
            rows = p1.load(out / "pending.local-vs21.pkl")
            for row in rows:
                if row["truth_key"] and row["family"] != "truth_in_set":
                    self.assertNotIn(row["truth_key"], row["member_keys"])
            random = [r for r in rows if r["truth_key"] and r["family"] == "same_formula_db_random"]
            self.assertEqual(len(random), 4)
            for a, b in zip(random, random[1:]):
                self.assertEqual(a["member_keys"], b["member_keys"][:a["k"]])

    def test_motifs_valid_and_public_positive_examples(self):
        from rdkit import Chem
        queries = p1.motifs()
        self.assertEqual(len(queries), 15)
        examples = {"indole": "c1ccc2[nH]ccc2c1", "catechol": "Oc1ccccc1O", "sulfate_ester": "CCOS(=O)(=O)O",
                    "phosphate_ester": "CCOP(=O)(O)O", "flavone": "O=c1cc(-c2ccccc2)oc2ccccc12",
                    "steroid": "CC12CCC3C(C1CCC2O)CCC4=CC(=CC=C34)O",
                    "glucuronide": "COC1OC(C(=O)O)C(O)C(O)C1O", "hexoside": "COC1OC(CO)C(O)C(O)C1O",
                    "pentoside": "COC1OCC(O)C(O)C1O", "deoxyhexoside": "COC1OC(C)C(O)C(O)C1O",
                    "gpc": "OCC(O)COP(=O)([O-])OCC[N+](C)(C)C", "gpe": "OCC(O)COP(=O)(O)OCCN",
                    "gps": "OCC(O)COP(=O)(O)OCC(N)C(=O)O", "gpi": "OCC(O)COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O",
                    "gpg": "OCC(O)COP(=O)(O)OCC(O)CO"}
        for motif, smiles in examples.items():
            self.assertTrue(Chem.MolFromSmiles(smiles).HasSubstructMatch(queries[motif]), motif)

    def test_timeout_partial_result_is_never_used(self):
        from unittest.mock import patch
        partial = SimpleNamespace(canceled=True, smartsString="[#6]", numAtoms=1)
        with patch("rdkit.Chem.rdFMCS.FindMCS", return_value=partial):
            self.assertEqual(core.mcs_over_set([core.standardise("CCCC"), core.standardise("CCCO")]), ("", 0, True))

    def test_spectral_scaffold_and_formula_families_exclude_truth(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            db_rows = [p1.chem_record(s) for s in ("CCc1ccccc1", "Cc1ccccc1C", "Cc1cccc(C)c1", "Cc1ccc(C)cc1", "Oc1ccccc1", "COc1ccccc1", "c1ccccc1")]
            db = {r["key"]: r for r in db_rows}
            truths = [dict(r, source="vs20", stratum="0:4-15", lipid=0, heavy_bin="4-15", motif_oversample=False) for r in (db_rows[0], db_rows[-1])]
            neighbours = []
            for truth in truths:
                candidates = [(0.9, r["key"], r["formula"], r["smiles"]) for r in db_rows if r["key"] != truth["key"]][:5]
                neighbours.append(dict(source="vs20", key=truth["key"], spectral_nn=candidates, spectral_nn_any=candidates,
                                       same_formula_spectral=[r for r in candidates if r[2] == truth["formula"]]))
            p1.dump(out / "prepared.local-vs21.pkl", {"db": db, "truths": truths})
            p1.dump(out / "neighbours.local-vs21.pkl", neighbours)
            p1.sets(SimpleNamespace(out=out, sources=["vs20"], seed=20261113, pilot=True), {})
            rows = [r for r in p1.load(out / "pending.local-vs21.pkl") if r["truth_key"]]
            for row in rows:
                if row["family"] != "truth_in_set":
                    self.assertNotIn(row["truth_key"], row["member_keys"])
                else:
                    self.assertIn(row["truth_key"], row["member_keys"])
            self.assertEqual({r["family"] for r in rows}, set(p1.FAMILIES))
            self.assertFalse(any(r["family"] == "spectral_nn_scaffold" and r["truth_key"] == db_rows[-1]["key"] for r in rows))


if __name__ == "__main__":
    unittest.main()
